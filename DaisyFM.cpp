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
FMDemodulator radioDemodulator2;

Parameter frequancyCtrl;
Parameter gainCtrl;
Parameter noiseCtrl;

FatFSInterface fsi;
SdmmcHandler   sdcard;
DIR dir;
FILINFO fil;

#define MAX_BUF_SIZE 4 * 1048576 // 2 x 4MB
int16_t DSY_SDRAM_BSS buffer_1[MAX_BUF_SIZE];
int16_t DSY_SDRAM_BSS buffer_2[MAX_BUF_SIZE];
uint32_t gSamplesElapsed = 0;
uint32_t seed_i = 0xA1B2C3D4u;  // any non-zero 32-bit seed
uint32_t seed_q = 0x5EED1234u;  // different non-zero seed

float normFreqCtrl = 0.0f;
float gainCtrldB = 0.0f;
float noiseVariance = 0.0f;
float outputGaindB_l = -60.0f;
float outputGaindB_r = -60.0f;
int prevRegion = 0;
int pendingRegion = -1;
uint32_t regionHoldMs = 0;
static const uint32_t REGION_DEBOUNCE_MS = 80;

void ProcessControls();
void InitFileSystem();
void DrawDisplay();
int InitRadioPlayer(int sr);

void AudioCallback(AudioHandle::InputBuffer in, AudioHandle::OutputBuffer out, size_t size)
{
	ProcessControls();

	float peak_l = 0.0f, peak_r = 0.0f;
    for(size_t i = 0; i < size; i += 1)
    {
		float out_i, out_q, out_i2, out_q2;
		float out_iR, out_qR, out_i2R, out_q2R;
		radioStation1.StreamFM(out_i, out_q, out_iR, out_qR);
		radioStation2.StreamFM(out_i2, out_q2, out_i2R, out_q2R);

		float awgn_i = gauss_approx(seed_i) * noiseVariance;
		float awgn_q = gauss_approx(seed_q) * noiseVariance;

		float output_l = radioDemodulator.Demodulate(out_i + out_i2 + awgn_i ,out_q + out_q2 + awgn_q);
		//outputGaindB_l = 10.0f * log10f(output * output);

		float output_r = radioDemodulator2.Demodulate(out_iR + out_i2R + awgn_i, out_qR + out_q2R + awgn_q);
		//outputGaindB_r = 10.0f * log10f(output * output);

		// track peak amplitude for display purposes
		float a = fabsf(output_l);
    	if (a > peak_l) peak_l = a;
    	a = fabsf(output_r);
    	if (a > peak_r) peak_r = a;

		out[0][i] = output_l;
		out[1][i] = output_r;
		gSamplesElapsed++;
	}
	outputGaindB_l = 20.0f * log10f(peak_l + 1e-7f);
	outputGaindB_r = 20.0f * log10f(peak_r + 1e-7f);
}

int main(void)
{
	hw.Init();
	hw.SetAudioBlockSize(128); // number of samples handled per callback
	hw.SetAudioSampleRate(SaiHandle::Config::SampleRate::SAI_48KHZ);

	hw.display.Fill(false);
	hw.display.SetCursor(30, 28);
	hw.display.WriteString("Loading...", Font_6x8, true);
	hw.display.Update();

	frequancyCtrl.Init(hw.controls[0], 0.0f, 1.0f, Parameter::LINEAR);
	gainCtrl.Init(hw.controls[1], -10.0f, 20.0f, Parameter::LINEAR);
	noiseCtrl.Init(hw.controls[2], -60.0f, -20.0f, Parameter::LINEAR);

	InitFileSystem();
	InitRadioPlayer(hw.AudioSampleRate());

	//hw.seed.StartLog(false);
	hw.StartAdc();
	hw.StartAudio(AudioCallback);

	while(1) {
		int region = (int)floorf(normFreqCtrl * 5.0f);

		if (region != pendingRegion) {
			pendingRegion = region;
			regionHoldMs  = System::GetNow();
		}

		if (pendingRegion != prevRegion &&
		    (System::GetNow() - regionHoldMs) >= REGION_DEBOUNCE_MS)
		{
			int res = 0;
			if (pendingRegion % 2 == 0)
			{
				res = radioStation1.SetFile(pendingRegion);
				hw.seed.PrintLine("Result read 1: %d", res);
				res = radioStation2.SetFile(pendingRegion + 1);
				hw.seed.PrintLine("Result read 2: %d", res);
			} else
			{
				res = radioStation1.SetFile(pendingRegion + 1);
				hw.seed.PrintLine("Result read 1: %d", res);
				res = radioStation2.SetFile(pendingRegion);
				hw.seed.PrintLine("Result read 2: %d", res);
			}
			prevRegion = pendingRegion;
		}

		DrawDisplay();

		hw.DelayMs(10);
	}
}

void DrawDisplay()
{
	hw.display.Fill(false);

	hw.display.SetCursor(1, 0);
	//hw.display.WriteString(" F M  S A M P L E R", Font_6x8, true);
	hw.display.WriteString("   D A I S Y  F M", Font_6x8, true);
	hw.display.SetCursor(8, 20);
	float sudoFreq = 87.5f + normFreqCtrl * 20.5f; // 87.5 MHz to 200 MHz
	int fracPart = (int)((sudoFreq - (int)sudoFreq) * 10);

	char freqBuf[8];
	snprintf(freqBuf, sizeof(freqBuf), "%s%d.%d",
				sudoFreq < 100.0f ? " " : "",
				(int)sudoFreq, fracPart);
	hw.display.WriteString(freqBuf, digitalFont_16x26, true);
	hw.display.SetCursor(91, 38);
	hw.display.WriteString("MHz", Font_6x8, true);

	int barX = 123;
	int barBottomY = 57;
	int barHeight = 50;

	float barPerc_l = (int)(outputGaindB_l + 60.0f) / 60.0f; // -60dB to 0dB
	if (barPerc_l < 0.0f) barPerc_l = 0.0f;
	if (barPerc_l > 1.0f) barPerc_l = 1.0f;

	float barPerc_r = (int)(outputGaindB_r + 60.0f) / 60.0f; // -60dB to 0dB
	if (barPerc_r < 0.0f) barPerc_r = 0.0f;
	if (barPerc_r > 1.0f) barPerc_r = 1.0f;
	
	int meterdB_l = (int)(barPerc_l * (float)barHeight);
	int meterdB_r = (int)(barPerc_r * (float)barHeight);

	hw.display.DrawLine(barX - 3, barBottomY - barHeight, barX + 3, barBottomY - barHeight, true); // 0 dB line
	hw.display.DrawRect(barX - 2, barBottomY - meterdB_l, barX - 1, barBottomY, true, true); // Draw the bar
	hw.display.DrawRect(barX + 1, barBottomY - meterdB_r, barX + 2, barBottomY, true, true); // Draw the bar	
	hw.display.DrawLine(barX - 3, barBottomY, barX + 3, barBottomY, true); // -60 dB line

	// Tuner indicator
	const int tunerCenterX = 64;
	const int tunerBaseY   = 57;
	const float tunerScale = 100.0f;

	hw.display.DrawLine(1, tunerBaseY, 118, tunerBaseY, true);
	hw.display.DrawLine(tunerCenterX,     tunerBaseY - 3, tunerCenterX,     tunerBaseY + 5, true);
	hw.display.DrawLine(tunerCenterX + 1, tunerBaseY - 3, tunerCenterX + 1, tunerBaseY + 5, true);

	for(int i = 0; i <= 5; i++) {
		float stationNorm = (float)i / 5.0f;
		int markerX = tunerCenterX + (int)roundf((stationNorm - normFreqCtrl) * tunerScale);
		if(markerX >= 2 && markerX <= 117 && markerX != tunerCenterX && markerX != tunerCenterX + 1) {
			hw.display.DrawLine(markerX, tunerBaseY + 1, markerX, tunerBaseY + 5, true);
		}
	}

	hw.display.Update();
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
	FRESULT result = FR_NOT_READY;
	for (int attempt = 0; attempt < 5 && result != FR_OK; attempt++) {
    	result = f_mount(&fsi.GetSDFileSystem(), "/", 1);
    	if (result != FR_OK) hw.DelayMs(100);
	}
	if (result != FR_OK) {
    	// Show error on display and halt or return
    	return -1;
	}

	radioStation1.Init(buffer_1, MAX_BUF_SIZE, sr);
	radioStation1.SetFile(0);
	radioStation1.SetCarrierFreq(6000.0f);
	radioStation1.Play();

	radioStation2.Init(buffer_2, MAX_BUF_SIZE, sr);
	radioStation2.SetFile(1);
	radioStation2.SetCarrierFreq(18000.0f);
	radioStation2.Play();

	radioDemodulator.Init(sr);
	radioDemodulator.SetCarrierFreq(6000.0f);

	radioDemodulator2.Init(sr);
	radioDemodulator2.SetCarrierFreq(18000.0f);

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
	float centerFrequency = FrequencyMapping(normFreqCtrl);
	radioDemodulator.SetCarrierFreq(centerFrequency);
	radioDemodulator2.SetCarrierFreq(centerFrequency);

	gainCtrldB = gainCtrl.Process();
	radioStation1.SetGain(powf(10.0f, gainCtrldB / 20.0f));
	radioStation2.SetGain(powf(10.0f, gainCtrldB / 20.0f));

	noiseVariance = powf(10.0f, noiseCtrl.Process() / 20.0f);
}
