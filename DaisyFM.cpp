#include "daisy_patch.h"
#include "daisysp.h"
#include "extra_fonts.h"

using namespace daisy;
using namespace daisysp;

DaisyPatch hw;

Parameter frequancyCtrl;
Parameter gainCtrl;
Parameter noiseCtrl;

FatFSInterface fsi;
SdmmcHandler   sdcard;
DIR dir;
FILINFO fil;

#define MAX_BUF_SIZE 2 * 1048576 // 4MB
int16_t DSY_SDRAM_BSS buffer[MAX_BUF_SIZE];
int16_t DSY_SDRAM_BSS buffer2[MAX_BUF_SIZE];

void ProcessControls();
void initFileSystem();
int initRadioPlayer(int sr);
float mapToFrequency(float normFreqCtrl);

float normFreqCtrl = 0.0f;
float gainCtrldB = 0.0f;
float noiseVariance = 0.0f;
int prevRegion = 0;
float centerFrequency = 5000.0f;
float outputGaindB = -60.0f;
//float sudoFreq = 87.5f;
static uint64_t gSamplesElapsed = 0;

inline double wrap_phase(double x, double len)
{
    const double q = std::floor(x / len);
    return x - q * len; // now in [0, len)
}

inline float gauss_approx(uint32_t& s){
    // xorshift32 PRNG
    auto rndu = [&](){
        s ^= s << 13; s ^= s >> 17; s ^= s << 5;
        return (s >> 8) * (1.0f/16777216.0f); // [0,1)
    };
    float sum = 0.f;
    for(int k=0;k<12;++k) sum += rndu();
    return sum - 6.f; // ~N(0,1)
}

struct BiquadDF2T
{
    // current (active) coefficients
    float cb0 = 1.f, cb1 = 0.f, cb2 = 0.f, ca1 = 0.f, ca2 = 0.f;

    // target coefficients
    float tb0 = 1.f, tb1 = 0.f, tb2 = 0.f, ta1 = 0.f, ta2 = 0.f;

    // per-sample deltas during a ramp
    float db0 = 0.f, db1 = 0.f, db2 = 0.f, da1 = 0.f, da2 = 0.f;

    // DF2T states
    float z1 = 0.f, z2 = 0.f;

    int rampSamples = 0;

    inline void reset()
    {
        z1 = z2 = 0.f;
    }

    // Set coefficients immediately (no ramp)
    inline void setImmediate(float b0, float b1, float b2, float a1, float a2)
    {
        cb0 = tb0 = b0; cb1 = tb1 = b1; cb2 = tb2 = b2;
        ca1 = ta1 = a1; ca2 = ta2 = a2;
        db0 = db1 = db2 = da1 = da2 = 0.f;
        rampSamples = 0;
    }

    // Set new target coefficients with a ramp
    inline void setTarget(float b0, float b1, float b2, float a1, float a2, int ramp)
    {
        tb0 = b0; tb1 = b1; tb2 = b2; ta1 = a1; ta2 = a2;
        rampSamples = ramp > 0 ? ramp : 0;

        if(rampSamples <= 0)
        {
            cb0 = tb0; cb1 = tb1; cb2 = tb2; ca1 = ta1; ca2 = ta2;
            db0 = db1 = db2 = da1 = da2 = 0.f;
        }
        else
        {
            db0 = (tb0 - cb0) / (float)rampSamples;
            db1 = (tb1 - cb1) / (float)rampSamples;
            db2 = (tb2 - cb2) / (float)rampSamples;
            da1 = (ta1 - ca1) / (float)rampSamples;
            da2 = (ta2 - ca2) / (float)rampSamples;
        }
    }

    inline float process(float x)
    {
        // ramp coeffs per sample (if needed)
        if(rampSamples > 0)
        {
            cb0 += db0; cb1 += db1; cb2 += db2; ca1 += da1; ca2 += da2;
            --rampSamples;
        }

        // Transposed Direct Form II
        float y = cb0 * x + z1;
        z1 = cb1 * x - ca1 * y + z2;
        z2 = cb2 * x - ca2 * y;
        return y;
    }
};

class LowPassFilter
{
	public:
	void setup(int sr)
    {
		sampleRate = (float)sr;
	}

	inline float process(float x)
    {
        return lp.process(x);
    }

	void setFrequency(float freq)
    {
        const float d = 1.414427157f;               // ~sqrt(2)
        auto designLP = [&](float fc, float& b0, float& b1, float& b2, float& a1, float& a2){
            float theta = 2.f * 3.14159265358979f * fc / sampleRate;
            float beta  = 0.5f * (1.f - 0.5f * d * sinf(theta)) / (1.f + 0.5f * d * sinf(theta));
            float gamma = (0.5f + beta) * cosf(theta);
            float a0    = 0.5f + beta - gamma;
            b0 = 0.5f * a0;  b1 = a0;  b2 = 0.5f * a0;
            a1 = -2.f * gamma; a2 = 2.f * beta;
        };

        float b0,b1,b2,a1,a2;

        // Low-pass edge at fHi
        designLP(freq, b0,b1,b2,a1,a2);
        lp.setImmediate(b0,b1,b2,a1,a2);
    }

	private:
    	float        sampleRate = 48000.f;
    	BiquadDF2T   lp;
};

class BandFilter
{
  public:
    void setup(int sr)
    {
        sampleRate = (float)sr;
        for(int i=0;i<3;++i){ hp[i].setImmediate(1.f, -2.f, 1.f, -2.f, 1.f); hp[i].reset(); }
        for(int i=0;i<3;++i){ lp[i].setImmediate(1.f,  2.f, 1.f, -2.f, 1.f); lp[i].reset(); }
    }

    inline float process(float x)
    {
        float y = x;
        // 4 × high-pass
        y = hp[0].process(y);
        y = hp[1].process(y);
        y = hp[2].process(y);
        //y = hp[3].process(y);
        // 4 × low-pass
        y = lp[0].process(y);
        y = lp[1].process(y);
        y = lp[2].process(y);
        //y = lp[3].process(y);
        return y;
    }

    void setCenterFrequency(float freq, float bandwidth, int rampSamples)
    {
        float nyq = 0.5f * sampleRate;
        float fLo = freq - 0.5f * bandwidth;
        float fHi = freq + 0.5f * bandwidth;
        if(fLo < 20.f)  fLo = 20.f;
        if(fHi > nyq*0.95f) fHi = nyq*0.95f;
        if(fLo >= fHi) { fLo = 0.5f * freq; fHi = 1.5f * freq; } // fallback

        const float d = 1.414427157f;               // ~sqrt(2)
        auto designLP = [&](float fc, float& b0, float& b1, float& b2, float& a1, float& a2){
            float theta = 2.f * 3.14159265358979f * fc / sampleRate;
            float beta  = 0.5f * (1.f - 0.5f * d * sinf(theta)) / (1.f + 0.5f * d * sinf(theta));
            float gamma = (0.5f + beta) * cosf(theta);
            float a0    = 0.5f + beta - gamma;
            b0 = 0.5f * a0;  b1 = a0;  b2 = 0.5f * a0;
            a1 = -2.f * gamma; a2 = 2.f * beta;
        };
        auto designHP = [&](float fc, float& b0, float& b1, float& b2, float& a1, float& a2){
            float theta = 2.f * 3.14159265358979f * fc / sampleRate;
            float beta  = 0.5f * (1.f - 0.5f * d * sinf(theta)) / (1.f + 0.5f * d * sinf(theta));
            float gamma = -(0.5f + beta) * cosf(theta);
            float a0    = 0.5f + beta - gamma;
            b0 =  0.5f * a0;  b1 = -a0;  b2 =  0.5f * a0;
            a1 = 2.f * gamma; a2 =  2.f * beta;
        };

        float b0,b1,b2,a1,a2;

        // High-pass edge at fLo
        designHP(fLo, b0,b1,b2,a1,a2);
        for(int i=0;i<3;++i)
            hp[i].setTarget(b0,b1,b2,a1,a2, rampSamples);

        // Low-pass edge at fHi
        designLP(fHi, b0,b1,b2,a1,a2);
        for(int i=0;i<3;++i)
            lp[i].setTarget(b0,b1,b2,a1,a2, rampSamples);
    }

    void reset() {
        for(int i=0;i<3;++i){ hp[i].reset(); lp[i].reset(); }
    }

  private:
    float        sampleRate = 48000.f;
    BiquadDF2T   hp[3];
    BiquadDF2T   lp[3];
};

class RadioStation {
	private:

	uint32_t bufferLength;
	size_t length;
	size_t position;
	int16_t *buffer_;
	FIL file;
	const TCHAR* filename;
	bool playing;
	int readIndex;
	float frac;
	int upsamplingFactor;
	int sampleRate;
	int currentFileIndex;
	float gain = 1.0f;

	Phasor carrierPhase;
	LowPassFilter inputFilter;
	float modulationIndex;
	float history;
	double start = 0.0;
	double pitch = 1.0;

	public:
	void Init(int16_t *buffer, uint32_t bufferLength, int sr){
		this->buffer_ = buffer;
		this->bufferLength = bufferLength;
		length = 0;
		playing = false;

		//readIndex = 0;
		frac = 0.0f;
		upsamplingFactor = 2;
		sampleRate = sr;
		modulationIndex = 500.0f;
		history = 0.0f;

		inputFilter.setup(sr);
		carrierPhase.Init(sampleRate);
	}

	void SetCarrierFreq(float freq){
		carrierPhase.SetFreq(freq);
	}

	void SetGain(float g){
		gain = g;
	}
	
	void Play(){
		playing = true;
	}
	void Stop(){
		playing = false;
		readIndex = 0;
		frac = 0.0f;
	}

	void debugPrint(){
		hw.seed.PrintLine("file: %d, readIndex: %d, length: %lu, frac: %d.%03d",
			currentFileIndex, readIndex, (unsigned long)length, (int)((frac - (int)frac) * 1000));
		//printf("RadioStation: playing: %d, readIndex: %d, length: %zu, frac: %f\n", playing, readIndex, length, frac);
	}

	float Stream(){

			if (!playing || length < 2) return 0.0f;

			const double L = (double)length;
			const double rate = pitch; 

			double phase = wrap_phase(start + (double)gSamplesElapsed * rate, L);

			// Linear interpolation
			size_t i0 = (size_t)phase;
			double frac = phase - (double)i0;
			size_t i1 = (i0 + 1 < length) ? (i0 + 1) : 0;

			float s0 = s162f(buffer_[i0]);
			float s1 = s162f(buffer_[i1]);
			return s0 + (s1 - s0) * (float)frac;
	}

	void StreamFM(float& out_i, float& out_q){
		float input = Stream();
		input = inputFilter.process(input);
		float theta = history + TWOPI_F * modulationIndex / sampleRate * input * gain;
		history = theta;
		float phs = carrierPhase.Process();
		out_i = sinf(TWOPI_F*phs + theta);
		out_q = cosf(TWOPI_F*phs + theta);
	}

	int SetFile(int fileIndex) {

		if (fileIndex == currentFileIndex)
		{
			return 0; //do nothing
		}
		switch (fileIndex) {
			case 0:
				filename = "radioStation-1.wav";
				break;
			case 1:
				filename = "radioStation-2.wav";
				break;
			case 2:
				filename = "radioStation-3.wav";
				break;
			case 3:
				filename = "radioStation-4.wav";
				break;
			case 4:
				filename = "radioStation-5.wav";
				break;
			case 5:
				filename = "radioStation-6.wav";
				break;
			default:
				filename = "radioStation-1.wav";
		}

		if(f_open(&file, filename, (FA_OPEN_EXISTING | FA_READ)) != FR_OK){
			return 1;
		}
		
		uint32_t data_offset = 0;
		char chunk_id[4];
		uint32_t chunk_size;
		UINT bytesread;
	
		f_lseek(&file, 12); // skip RIFF header
	
		while (true) {
			f_read(&file, chunk_id, 4, &bytesread);
			f_read(&file, &chunk_size, 4, &bytesread);
	
			if (strncmp(chunk_id, "data", 4) == 0) {
				data_offset = f_tell(&file); // current position = start of audio data
				break;
			} else {
				f_lseek(&file, f_tell(&file) + chunk_size); // skip this chunk
			}
		}
		f_lseek(&file, data_offset);
	
		if (f_read(&file, buffer_, MAX_BUF_SIZE * sizeof(int16_t), &bytesread) != FR_OK){
			return 2;
		}
	
		length = size_t(bytesread / 2);
		f_close(&file);
		currentFileIndex = fileIndex;
		return 0;
	}
};

class FMDemodulator {
	private:
		Phasor carrierPhase;
		//ComplexOsc carrierOsc;
		float modulationIndex;
		int sampleRate;
		float history_i;
		float history_q;
		BandFilter bandFilter_i;
		BandFilter bandFilter_q;
		LowPassFilter outputFilter;
		DcBlock dcBlock;

	public:

	void Init(float sr){
		//carrierOsc.SetFreq(sr, 6000.0f); // Default carrier frequency
		carrierPhase.Init(sr);
		carrierPhase.SetFreq(5000.0f); // Default carrier frequency
		bandFilter_i.setup(sr);
		bandFilter_q.setup(sr);
		outputFilter.setup(sr);
		outputFilter.setFrequency(8000.0f); // Default output filter frequency
		dcBlock.Init(sr);
		sampleRate = sr;
		modulationIndex = 500.0f;
	}

	void SetCarrierFreq(float freq){
		carrierPhase.SetFreq(freq);
		//carrierOsc.SetFreq(sampleRate, freq);
		bandFilter_i.setCenterFrequency(freq, 5000.0f, 64);
		bandFilter_q.setCenterFrequency(freq, 5000.0f, 64);
	}

	float Demodulate(float rx_i, float rx_q){
		// carrierOsc.Step(); // Update the oscillator phase
		// float c_i = carrierOsc.s; // Get the current cos value
		// float c_q = carrierOsc.c; // Get the current sin value

		float phs = carrierPhase.Process();
		float c_i = sinf(TWOPI_F*phs);
		float c_q = cosf(TWOPI_F*phs);

		float fltrx_i = bandFilter_i.process(rx_i);
		float fltrx_q = bandFilter_q.process(rx_q);

		//return fltrx_i;

		float cmplxmult_i = fltrx_i * c_i + fltrx_q * c_q;
		float cmplxmult_q = fltrx_q * c_i - fltrx_i * c_q;


		// --- Approach 1: Using complex multiplication and phase derivative (atan is expensive)
		// float deriv_i = cmplxmult_i * history_i + cmplxmult_q * history_q;
		// float deriv_q = cmplxmult_q * history_i - cmplxmult_i * history_q;

		// history_i = cmplxmult_i;
		// history_q = cmplxmult_q;

		//return atan2(deriv_q, deriv_i) * sampleRate / (TWOPI_F * modulationIndex);
		// float demod = atan2(deriv_q, deriv_i) * sampleRate / (TWOPI_F * modulationIndex);
		// float dc_blocked = dcBlock.Process(demod);
		// return SoftLimit(2.0f * outputFilter.process(dc_blocked));


		// --- Approach 2: Using classic discriminator: (without atan)
		float zi = cmplxmult_i, zq = cmplxmult_q;
		//float yi = zi*history_i + zq*history_q;          // real part (unused)
		float yq = zq*history_i - zi*history_q;          // imag part ~ Δphase
		float invpow = 1.0f / (zi*zi + zq*zq + 1e-6f);
		float demod = yq * invpow * 1.0f;         // choose scale to taste
		history_i = zi; history_q = zq;
		float dc_blocked = dcBlock.Process(demod);
		return SoftLimit(1.5f * outputFilter.process(dc_blocked));
	}
};

RadioStation radioStation1;
RadioStation radioStation2;
FMDemodulator radioDemodulator;

uint32_t seed_i = 0xA1B2C3D4u;  // any non-zero 32-bit seed
uint32_t seed_q = 0x5EED1234u;  // different non-zero seed

void AudioCallback(AudioHandle::InputBuffer in, AudioHandle::OutputBuffer out, size_t size)
{
	ProcessControls();

    for(size_t i = 0; i < size; i += 1)
    {
		float out_i, out_q, out_i2, out_q2;
		radioStation1.StreamFM(out_i, out_q);

		radioStation2.StreamFM(out_i2, out_q2);

		float awgn_i = gauss_approx(seed_i) * noiseVariance;
		float awgn_q = gauss_approx(seed_q) * noiseVariance;

		float output = radioDemodulator.Demodulate(out_i + out_i2 + awgn_i ,out_q + out_q2 + awgn_q);

		outputGaindB = 10.0f * log10f(output * output);
		
		out[0][i] = output; // Output demodulated signal to right channel
		gSamplesElapsed++;
	}
}

int main(void)
{
	hw.Init();
	hw.SetAudioBlockSize(64); // number of samples handled per callback
	hw.SetAudioSampleRate(SaiHandle::Config::SampleRate::SAI_48KHZ);

	frequancyCtrl.Init(hw.controls[0], 0.0f, 1.0f, Parameter::LINEAR);
	gainCtrl.Init(hw.controls[1], -10.0f, 20.0f, Parameter::LINEAR);
	noiseCtrl.Init(hw.controls[2], -60.0f, -20.0f, Parameter::LINEAR);

	initFileSystem();

	int sampleRate = hw.AudioSampleRate();
	initRadioPlayer(sampleRate);

	//hw.seed.StartLog(false);
	hw.StartAdc();
	hw.StartAudio(AudioCallback);

	//int counter = 0;
	while(1) {

		int region = floor(normFreqCtrl * 5.0f);
		if (region != prevRegion)
		{
			if (region % 2 == 0)
			{
				radioStation1.SetFile(region);
				radioStation2.SetFile(region + 1);
			} else 
			{
				radioStation1.SetFile(region + 1);
				radioStation2.SetFile(region);
			}
			prevRegion = region;
		}


		// GUI
		hw.display.Fill(false);

		hw.display.SetCursor(1, 0);
		std::string str  = "FM";
		char*       cstr = &str[0];
		hw.display.WriteString(cstr, Font_6x8, true);

		float sudoFreq = 87.5f + normFreqCtrl * 20.5f; // 87.5 MHz to 200 MHz
		int fracPart = (int)((sudoFreq - (int)sudoFreq) * 10);
		hw.display.SetCursor(18, 20);

		if (sudoFreq < 100.0f){
			str = " " + std::to_string((int) sudoFreq) + "." + std::to_string(fracPart);
		} else {
			str = std::to_string((int) sudoFreq) + "." + std::to_string(fracPart);
		}
		
    	hw.display.WriteString(cstr, digitalFont_16x26, true);
		hw.display.SetCursor(101, 38);
		hw.display.WriteString("MHz", Font_6x8, true);


		float barPerc = (int)(outputGaindB + 60.0f) / 60.0f; // -60dB to 0dB
		if (barPerc < 0.0f) barPerc = 0.0f;
		if (barPerc > 1.0f) barPerc = 1.0f;
		
		int barX = 5;
		int barBottomY = 60;
		int barHeight = 48;

		int meterdB = (int)(barPerc * (float)barHeight);

		hw.display.DrawLine(barX - 2, barBottomY - barHeight, barX + 2, barBottomY - barHeight, true); // 0 dB line
		hw.display.DrawRect(barX - 1, barBottomY - meterdB, barX + 1, barBottomY, true, true); // Draw the bar
		hw.display.DrawLine(barX - 2, barBottomY, barX + 2, barBottomY, true); // -60 dB line


		// hw.display.SetCursor(0, 0);

		// //hw.display.WriteString("!", gBitmap, true);
		// hw.display.WriteString("104.7", ddFont16x26, true);
		// //hw.display.WriteString("!", gBitmap, true);
		
		hw.display.Update();		
		hw.DelayMs(10);
	}
}

void initFileSystem()
{
	SdmmcHandler::Config sd_config;
	sd_config.Defaults();
	sd_config.speed = SdmmcHandler::Speed::STANDARD;
	sdcard.Init(sd_config);
	fsi.Init(FatFSInterface::Config::MEDIA_SD);
}

int initRadioPlayer(int sr)
{
	FRESULT result;
	result = f_mount(&fsi.GetSDFileSystem(), "/", 1);  // opt = 1 forces mount now
	f_opendir(&dir, fsi.GetSDPath());
	if (result != FR_OK){
		return -1;
	}

	radioDemodulator.Init(sr);
	radioStation1.Init(buffer, MAX_BUF_SIZE, sr);
	radioStation2.Init(buffer2, MAX_BUF_SIZE, sr);
	int res = radioStation1.SetFile(0);
	res = radioStation2.SetFile(1);

	radioStation1.SetCarrierFreq(6000.0f);
	radioStation2.SetCarrierFreq(18000.0f);
	radioDemodulator.SetCarrierFreq(6000.0f);

	radioStation1.Play();
	radioStation2.Play();

	return res;
}

float mapToFrequency(float normFreqCtrl)
{
	float centerFrequency = 6000.0f; // Default center frequenc
	int region = floor(normFreqCtrl * 5.0f);
	normFreqCtrl = normFreqCtrl * 5.0f - region;

	if (region % 2 == 0)
	{
		centerFrequency = 6000.0f + normFreqCtrl * 12000.0f;
	} else 
	{
		centerFrequency = 18000.0f - normFreqCtrl * 12000.0f;
	}
	//float centerFrequency = 20.0f + normFreqCtrl * 23980.0f;
	return centerFrequency;
}

void ProcessControls()
{
	hw.ProcessAllControls();

	normFreqCtrl = frequancyCtrl.Process();
	centerFrequency = mapToFrequency(normFreqCtrl);
	radioDemodulator.SetCarrierFreq(centerFrequency);

	gainCtrldB = gainCtrl.Process();
	radioStation1.SetGain(powf(10.0f, gainCtrldB / 20.0f));
	radioStation2.SetGain(powf(10.0f, gainCtrldB / 20.0f));

	noiseVariance = powf(10.0f, noiseCtrl.Process() / 20.0f);


}
