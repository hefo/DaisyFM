#include "daisy_patch.h"
#include "daisysp.h"
#include "extra_fonts.h"
#include "DaisyFM.h"

using namespace daisy;
using namespace daisysp;

DaisyPatch hw;

RadioStation radioStation1;
RadioStation radioStation2;
FMDemodulator radioDemodulator;

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
uint64_t gSamplesElapsed = 0;
uint32_t seed_i = 0xA1B2C3D4u;  // any non-zero 32-bit seed
uint32_t seed_q = 0x5EED1234u;  // different non-zero seed

float normFreqCtrl = 0.0f;
float gainCtrldB = 0.0f;
float noiseVariance = 0.0f;
float centerFrequency = 5000.0f;
float outputGaindB = -60.0f;
int prevRegion = 0;

void ProcessControls();
void InitFileSystem();
int InitRadioPlayer(int sr);
//float mapToFrequency(float normFreqCtrl);

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
	hw.SetAudioBlockSize(128); // number of samples handled per callback
	hw.SetAudioSampleRate(SaiHandle::Config::SampleRate::SAI_48KHZ);

	frequancyCtrl.Init(hw.controls[0], 0.0f, 1.0f, Parameter::LINEAR);
	gainCtrl.Init(hw.controls[1], -10.0f, 20.0f, Parameter::LINEAR);
	noiseCtrl.Init(hw.controls[2], -60.0f, -20.0f, Parameter::LINEAR);

	InitFileSystem();
	InitRadioPlayer(hw.AudioSampleRate());

	//hw.seed.StartLog(false);
	hw.StartAdc();
	hw.StartAudio(AudioCallback);

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

		// Draw GUI
		hw.display.Fill(false);

		hw.display.SetCursor(1, 0);
		std::string str  = "F M  S A M P L E R";
		char*       cstr = &str[0];
		hw.display.WriteString(cstr, Font_6x8, true);


		hw.display.SetCursor(12, 20);
		float sudoFreq = 87.5f + normFreqCtrl * 20.5f; // 87.5 MHz to 200 MHz
		int fracPart = (int)((sudoFreq - (int)sudoFreq) * 10);

		if (sudoFreq < 100.0f){
			str = " " + std::to_string((int) sudoFreq) + "." + std::to_string(fracPart);
		} else {
			str = std::to_string((int) sudoFreq) + "." + std::to_string(fracPart);
		}
		
    	hw.display.WriteString(cstr, digitalFont_16x26, true);
		hw.display.SetCursor(95, 38);
		hw.display.WriteString("MHz", Font_6x8, true);


		float barPerc = (int)(outputGaindB + 60.0f) / 60.0f; // -60dB to 0dB
		if (barPerc < 0.0f) barPerc = 0.0f;
		if (barPerc > 1.0f) barPerc = 1.0f;
		
		int barX = 123;
		int barBottomY = 57;
		int barHeight = 50;

		int meterdB = (int)(barPerc * (float)barHeight);

		hw.display.DrawLine(barX - 2, barBottomY - barHeight, barX + 2, barBottomY - barHeight, true); // 0 dB line
		hw.display.DrawRect(barX - 1, barBottomY - meterdB, barX + 1, barBottomY, true, true); // Draw the bar
		hw.display.DrawLine(barX - 2, barBottomY, barX + 2, barBottomY, true); // -60 dB line

		hw.display.Update();		
		hw.DelayMs(30);
	}
}

void InitFileSystem()
{
	SdmmcHandler::Config sd_config;
	sd_config.Defaults();
	sd_config.speed = SdmmcHandler::Speed::STANDARD;
	sdcard.Init(sd_config);
	fsi.Init(FatFSInterface::Config::MEDIA_SD);
}

int InitRadioPlayer(int sr)
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

	radioStation1.SetFile(0);
	radioStation1.SetCarrierFreq(6000.0f);
	radioStation1.Play();

	radioStation2.SetFile(1);
	radioStation2.SetCarrierFreq(18000.0f);
	radioStation2.Play();

	radioDemodulator.SetCarrierFreq(6000.0f);
	return 0;
}

inline float FrequencyMapping(float normFreqCtrl)
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

	return centerFrequency;
}

void ProcessControls()
{
	hw.ProcessAllControls();

	normFreqCtrl = frequancyCtrl.Process();
	centerFrequency = FrequencyMapping(normFreqCtrl);
	radioDemodulator.SetCarrierFreq(centerFrequency);

	gainCtrldB = gainCtrl.Process();
	radioStation1.SetGain(powf(10.0f, gainCtrldB / 20.0f));
	radioStation2.SetGain(powf(10.0f, gainCtrldB / 20.0f));

	noiseVariance = powf(10.0f, noiseCtrl.Process() / 20.0f);
}
