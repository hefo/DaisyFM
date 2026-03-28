#include "daisy_patch.h"
#include "daisysp.h"
#include "extra_fonts.h"
#include "DaisyFM.h"
#include <string.h>

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

#define MAX_BUF_SIZE  (4 * 1048576) // 8 MB per station (4 M int16 elements, ~44 s stereo 48 kHz)
#define NUM_STATIONS  6
int16_t DSY_SDRAM_BSS station_buffers[NUM_STATIONS][MAX_BUF_SIZE];

// Staging buffer in on-chip AXI-SRAM for SDMMC IDMA.
// SDMMC IDMA on STM32H7 reaches on-chip SRAM via a direct AHB path; writing
// directly to FMC-connected SDRAM goes through a slower AXI→FMC route that can
// stall. We DMA here first, then memcpy to SDRAM.
static uint8_t __attribute__((aligned(32))) sdReadStaging[64 * 1024];
size_t  station_lengths[NUM_STATIONS];
uint32_t gSamplesElapsed = 0;
uint32_t seed_i = 0xA1B2C3D4u;  // any non-zero 32-bit seed
uint32_t seed_q = 0x5EED1234u;  // different non-zero seed

//Frequency knob calibration
static constexpr float freqPotMin = 0.05f;
static constexpr float freqPotMax = 0.95f;

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
static int LoadWavFile(int fileIndex, int16_t* buf, uint32_t maxElements, size_t* outLength);

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
			int fileA, fileB;
			if (pendingRegion % 2 == 0) {
				fileA = pendingRegion;
				fileB = pendingRegion + 1;
			} else {
				fileA = pendingRegion + 1;
				fileB = pendingRegion;
			}
			if (fileA >= NUM_STATIONS) fileA = NUM_STATIONS - 1;
			if (fileB >= NUM_STATIONS) fileB = NUM_STATIONS - 1;

			radioStation1.SetBuffer(station_buffers[fileA], station_lengths[fileA]);
			radioStation2.SetBuffer(station_buffers[fileB], station_lengths[fileB]);
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

	hw.display.DrawLine(12, tunerBaseY, 115, tunerBaseY, true);
	hw.display.DrawLine(tunerCenterX,     tunerBaseY - 3, tunerCenterX,     tunerBaseY + 5, true);
	hw.display.DrawLine(tunerCenterX + 1, tunerBaseY - 3, tunerCenterX + 1, tunerBaseY + 5, true);

	for(int i = 0; i <= 5; i++) {
		float stationNorm = (float)i / 5.0f;
		int markerX = tunerCenterX + (int)roundf((stationNorm - normFreqCtrl) * tunerScale);
		if(markerX >= 13 && markerX <= 114 && markerX != tunerCenterX && markerX != tunerCenterX + 1) {
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

// Load a WAV file by index into buf (up to maxElements int16 values).
// Sets *outLength to the number of stereo sample-pairs loaded.
// Returns 0 on success, non-zero on error.
static int LoadWavFile(int fileIndex, int16_t* buf, uint32_t maxElements, size_t* outLength)
{
	const TCHAR* filename;
	switch (fileIndex) {
		case 0: filename = "radioStation-1.wav"; break;
		case 1: filename = "radioStation-2.wav"; break;
		case 2: filename = "radioStation-3.wav"; break;
		case 3: filename = "radioStation-4.wav"; break;
		case 4: filename = "radioStation-5.wav"; break;
		case 5: filename = "radioStation-6.wav"; break;
		default: return 1;
	}

	static FIL file;
	if (f_open(&file, filename, FA_OPEN_EXISTING | FA_READ) != FR_OK) return 1;

	UINT br = 0;
	char riff_id[4], wave_id[4];
	uint32_t riff_size = 0;
	if (f_read(&file, riff_id,  4, &br) != FR_OK || br != 4 ||
	    f_read(&file, &riff_size, 4, &br) != FR_OK || br != 4 ||
	    f_read(&file, wave_id,  4, &br) != FR_OK || br != 4) { f_close(&file); return 2; }
	if (strncmp(riff_id, "RIFF", 4) != 0 || strncmp(wave_id, "WAVE", 4) != 0)
		{ f_close(&file); return 3; }

	bool have_fmt = false;
	uint16_t num_channels = 0, bits_per_sample = 0, audio_format = 0;
	uint32_t data_size = 0;
	DWORD data_pos = 0;

	for (;;) {
		char chunk_id[4];
		uint32_t chunk_size = 0;
		if (f_read(&file, chunk_id,   4, &br) != FR_OK || br != 4) break;
		if (f_read(&file, &chunk_size, 4, &br) != FR_OK || br != 4) break;

		if (strncmp(chunk_id, "fmt ", 4) == 0) {
			uint8_t hdr[32];
			UINT toread = (chunk_size < sizeof(hdr)) ? chunk_size : (UINT)sizeof(hdr);
			if (f_read(&file, hdr, toread, &br) != FR_OK || br != toread)
				{ f_close(&file); return 2; }
			if (chunk_size > toread)
				f_lseek(&file, f_tell(&file) + (chunk_size - toread));
			if (chunk_size >= 16) {
				audio_format    = *(uint16_t*)(hdr + 0);
				num_channels    = *(uint16_t*)(hdr + 2);
				bits_per_sample = *(uint16_t*)(hdr + 14);
				have_fmt = true;
			}
			if (chunk_size & 1) f_lseek(&file, f_tell(&file) + 1);
		} else if (strncmp(chunk_id, "data", 4) == 0) {
			data_size = chunk_size;
			data_pos  = f_tell(&file);
			break;
		} else {
			f_lseek(&file, f_tell(&file) + chunk_size + (chunk_size & 1));
		}
	}

	if (!have_fmt || data_pos == 0) { f_close(&file); return 4; }
	if (!(audio_format == 1 && num_channels == 2 && bits_per_sample == 16))
		{ f_close(&file); return 5; }

	uint32_t max_bytes  = maxElements * sizeof(int16_t);
	uint32_t want_bytes = (data_size < max_bytes) ? data_size : max_bytes;
	want_bytes &= ~1u;
	if (want_bytes == 0) { f_close(&file); return 6; }

	// Read through on-chip SRAM staging buffer to avoid DMA directly into SDRAM.
	uint32_t totalRead = 0;
	while (totalRead < want_bytes) {
		uint32_t chunk = want_bytes - totalRead;
		if (chunk > (uint32_t)sizeof(sdReadStaging))
			chunk = (uint32_t)sizeof(sdReadStaging);
		UINT chunk_br = 0;
		if (f_read(&file, sdReadStaging, chunk, &chunk_br) != FR_OK || chunk_br == 0)
			break;
		memcpy((uint8_t*)buf + totalRead, sdReadStaging, chunk_br);
		totalRead += chunk_br;
		if (chunk_br < chunk) break; // EOF reached early
	}
	if (totalRead != want_bytes) { f_close(&file); return 8; }

	*outLength = size_t(totalRead / sizeof(int16_t) / 2);
	f_close(&file);
	return 0;
}

int InitRadioPlayer(int sr)
{
	FRESULT result = FR_NOT_READY;
	for (int attempt = 0; attempt < 5 && result != FR_OK; attempt++) {
    	result = f_mount(&fsi.GetSDFileSystem(), "/", 1);
    	if (result != FR_OK) hw.DelayMs(100);
	}
	if (result != FR_OK) return -1;

	// Preload all station files into SDRAM with progress shown on OLED
	for (int i = 0; i < NUM_STATIONS; i++) {
		hw.display.Fill(false);
		hw.display.SetCursor(20, 24);
		char msg[20];
		snprintf(msg, sizeof(msg), "Loading %d/%d", i + 1, NUM_STATIONS);
		hw.display.WriteString(msg, Font_6x8, true);
		hw.display.Update();

		station_lengths[i] = 0;
		LoadWavFile(i, station_buffers[i], MAX_BUF_SIZE, &station_lengths[i]);
	}

	radioStation1.Init(station_buffers[0], sr);
	radioStation1.SetCarrierFreq(6000.0f);
	radioStation1.SetBuffer(station_buffers[0], station_lengths[0]);

	radioStation2.Init(station_buffers[1], sr);
	radioStation2.SetCarrierFreq(18000.0f);
	radioStation2.SetBuffer(station_buffers[1], station_lengths[1]);

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
	normFreqCtrl = fclamp((normFreqCtrl - freqPotMin) / (freqPotMax - freqPotMin), 0.0f, 1.0f);
	float centerFrequency = FrequencyMapping(normFreqCtrl);
	radioDemodulator.SetCarrierFreq(centerFrequency);
	radioDemodulator2.SetCarrierFreq(centerFrequency);

	gainCtrldB = gainCtrl.Process();
	radioStation1.SetGain(powf(10.0f, gainCtrldB / 20.0f));
	radioStation2.SetGain(powf(10.0f, gainCtrldB / 20.0f));

	noiseVariance = powf(10.0f, noiseCtrl.Process() / 20.0f);
}
