#include "daisy_patch.h"
#include "daisysp.h"

using namespace daisy;
using namespace daisysp;

extern uint64_t gSamplesElapsed;

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

	Phasor carrierPhase;
	LowPassFilter inputFilterL, inputFilterR;
	uint32_t maxBufferLength;
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
	float modulationIndex;
	float gain = 1.0f;
	float historyL;
	float historyR;
	double start = 0.0;
	double pitch = 1.0;

	public:
	void Init(int16_t *buffer, uint32_t maxBufferLength, int sr){
		this->buffer_ = buffer;
		this->maxBufferLength = maxBufferLength;
		length = 0;
		playing = false;
		frac = 0.0f;
		upsamplingFactor = 2;
		sampleRate = sr;
		modulationIndex = 500.0f;
		historyL = 0.0f;
		historyR = 0.0f;
		currentFileIndex = -1;

		inputFilterL.setup(sr);
		inputFilterR.setup(sr);
		carrierPhase.Init(sr);
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

	void Stream(float& sample_l, float& sample_r){
			const double L = (double)length;
			const double rate = pitch; 
			double phase = wrap_phase(start + (double)gSamplesElapsed * rate, L);

			// Linear interpolation
			size_t i0 = (size_t)phase;
			double frac = phase - (double)i0;
			size_t i1 = (i0 + 1 < length) ? (i0 + 1) : 0;

			float s0 = s162f(buffer_[2*i0]);
			float s1 = s162f(buffer_[2*i1]);
			sample_l = s0 + (s1 - s0) * (float)frac;

			s0 = s162f(buffer_[2*i0+1]);
			s1 = s162f(buffer_[2*i1+1]);
			sample_r = s0 + (s1 - s0) * (float)frac;
	}

	void StreamFM(float& out_i_l, float& out_q_l, float& out_i_r, float& out_q_r){
		float input_l, input_r;

		if (playing && length > 2){
			Stream(input_l, input_r);
		} else {
			input_l = 0.0f;
			input_r = 0.0f;
		}
		
		input_l = inputFilterL.process(input_l);
		input_r = inputFilterR.process(input_r);

		float phs = carrierPhase.Process();

		float thetaL = historyL + TWOPI_F * modulationIndex / sampleRate * input_l * gain;
		historyL = thetaL;
		out_i_l = sinf(TWOPI_F*phs + thetaL);
		out_q_l = cosf(TWOPI_F*phs + thetaL);

		float thetaR = historyR + TWOPI_F * modulationIndex / sampleRate * input_r * gain;
		historyR = thetaR;
		out_i_r = sinf(TWOPI_F*phs + thetaR);
		out_q_r = cosf(TWOPI_F*phs + thetaR);
	}

	int SetFile(int fileIndex) {
		if (fileIndex == currentFileIndex) {
			return 0; // do nothing
		}

		switch (fileIndex) {
			case 0: filename = "radioStation-1.wav"; break;
			case 1: filename = "radioStation-2.wav"; break;
			case 2: filename = "radioStation-3.wav"; break;
			case 3: filename = "radioStation-4.wav"; break;
			case 4: filename = "radioStation-5.wav"; break;
			case 5: filename = "radioStation-6.wav"; break;
			default: filename = "radioStation-1.wav"; break;
		}

		if (f_open(&file, filename, (FA_OPEN_EXISTING | FA_READ)) != FR_OK) {
			return 1;
		}

		// --- Read RIFF header (12 bytes) ---
		char riff_id[4];
		uint32_t riff_size = 0;
		char wave_id[4];
		UINT br = 0;

		if (f_read(&file, riff_id, 4, &br) != FR_OK || br != 4 ||
			f_read(&file, &riff_size, 4, &br) != FR_OK || br != 4 ||
			f_read(&file, wave_id, 4, &br) != FR_OK || br != 4) {
			f_close(&file);
			return 2; // I/O error
		}

		if (strncmp(riff_id, "RIFF", 4) != 0 || strncmp(wave_id, "WAVE", 4) != 0) {
			f_close(&file);
			return 3; // not a WAVE file
		}

		bool have_fmt = false;
		uint16_t num_channels = 0;
		uint16_t bits_per_sample = 0;
		uint16_t audio_format = 0;

		uint32_t data_size = 0;
		DWORD data_pos = 0;

		for (;;) {
			char chunk_id[4];
			uint32_t chunk_size = 0;

			if (f_read(&file, chunk_id, 4, &br) != FR_OK || br != 4) break;
			if (f_read(&file, &chunk_size, 4, &br) != FR_OK || br != 4) break;

			if (strncmp(chunk_id, "fmt ", 4) == 0) {
				uint8_t hdr[32]; // enough for common PCM/Extensible front part
				UINT toread = (chunk_size < sizeof(hdr)) ? chunk_size : (UINT)sizeof(hdr);
				if (f_read(&file, hdr, toread, &br) != FR_OK || br != toread) { f_close(&file); return 2; }

				if (chunk_size > toread) {
					f_lseek(&file, f_tell(&file) + (chunk_size - toread));
				}

				// parse the first 16 bytes (PCM core)
				if (chunk_size >= 16) {
					audio_format    = *(uint16_t*)(hdr + 0);
					num_channels    = *(uint16_t*)(hdr + 2);
					bits_per_sample = *(uint16_t*)(hdr + 14);
					have_fmt = true;
				}

				// pad byte if odd
				if (chunk_size & 1) f_lseek(&file, f_tell(&file) + 1);

			} else if (strncmp(chunk_id, "data", 4) == 0) {
				data_size = chunk_size;
				data_pos  = f_tell(&file); // start of audio frames
				// we can break here; optionally continue to scan other chunks if needed
				break;

			} else {
				// skip unknown chunk + pad byte if odd
				f_lseek(&file, f_tell(&file) + chunk_size + (chunk_size & 1));
			}
		}

		if (!have_fmt || data_pos == 0) {
			f_close(&file);
			return 4; // missing fmt or data
		}

		// --- Validate expected format: stereo, 16-bit PCM ---
		if (!(audio_format == 1 /*PCM*/ && num_channels == 2 && bits_per_sample == 16)) {
			f_close(&file);
			return 5; // unsupported format for this loader
		}

		// --- Read audio data (interleaved int16 LRLR...) ---
		f_lseek(&file, data_pos);

		// Limit read to both the buffer and the chunk size
		uint32_t max_bytes = maxBufferLength * sizeof(int16_t); // capacity in bytes
		uint32_t want_bytes = (data_size < max_bytes) ? data_size : max_bytes;

		if (f_read(&file, buffer_, want_bytes, &br) != FR_OK || br != want_bytes) {
			f_close(&file);
			return 8; // read error
		}

		length = size_t(br / sizeof(int16_t) / 2);       // total int16 samples (L+R interleaved)
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
