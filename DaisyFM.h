#include "daisy_patch.h"
#include "daisysp.h"

using namespace daisy;
using namespace daisysp;

extern uint32_t gSamplesElapsed;

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

struct PreEmphasis {
    float K     = 7.2f;  // 2 * Fs * tau; tau=75e-6, Fs=48000
    float x_prev = 0.f;
    float y_prev = 0.f;
    void init(float fs, float tau = 75e-6f) { K = 2.f * fs * tau; }
    inline float process(float x) {
        float y = (1.f + K) * x + (1.f - K) * x_prev - y_prev;
        x_prev = x; y_prev = y;
        return y;
    }
};

struct DeEmphasis {
    float K     = 7.2f;
    float x_prev = 0.f;
    float y_prev = 0.f;
    void init(float fs, float tau = 75e-6f) { K = 2.f * fs * tau; }
    inline float process(float x) {
        float y = (x + x_prev - (1.f - K) * y_prev) / (1.f + K);
        x_prev = x; y_prev = y;
        return y;
    }
};

/* | fn (Hz) | ωn (rad/sample) | a1     | a2      | Character                    |
|---------|-----------------|--------|---------|------------------------------|
| 300     | 0.0393          | 0.0556 | 0.00154 | Slow, smooth — vintage feel  |
| 500     | 0.0654          | 0.0927 | 0.00428 | Balanced — recommended start |
| 1000    | 0.1309          | 0.1852 | 0.01713 | Fast, responsive              | */

struct PLL {
    float a1           = 0.1852f;  // proportional gain (fn=500 Hz, zeta=0.707)
    float a2           = 0.01713f; // integral gain
    float output_scale = 25.465f;  // Fs / (2*pi*kf); default kf=300 Hz
    float phi_v        = 0.f;      // VCO phase accumulator
    float s_int        = 0.f;      // loop integrator state

    void init(float fs, float kf = 300.f, float fn = 500.f, float zeta = 0.707f) {
        float wn   = TWOPI_F * fn / fs;
        a1         = 2.f * zeta * wn;
        a2         = wn * wn;
        output_scale = fs / (TWOPI_F * kf);
        phi_v = 0.f; s_int = 0.f;
    }

    inline float process(float I, float Q) {
        float cos_phi = cosf(phi_v);
        float sin_phi = sinf(phi_v);
        float eps_raw = Q * cos_phi - I * sin_phi;

        float amplitude = sqrtf(I * I + Q * Q) + 1e-10f;
        float eps = eps_raw / amplitude;

        s_int = s_int + a2 * eps;
        s_int = fclamp(s_int, -0.5f, 0.5f);  // anti-windup
        float u = a1 * eps + s_int;

        phi_v = phi_v + u;
        return u * output_scale;
    }
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
	PreEmphasis preEmphL, preEmphR;
	size_t length;
	size_t position;
	int16_t *buffer_;
	bool playing;
	int sampleRate;
	float modulationIndex;
	float gain = 1.0f;
	float historyL;
	float historyR;

	public:
	void Init(int16_t *buffer, int sr){
		buffer_ = buffer;
		length = 0;
		playing = false;
		sampleRate = sr;
		modulationIndex = 300.0f;
		historyL = 0.0f;
		historyR = 0.0f;

		preEmphL.init((float)sr);
		preEmphR.init((float)sr);
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
	}

	// Instantly switch to a pre-loaded buffer. ISR-safe: playing=false prevents
	// the audio callback from reading buffer_/length during the update window.
	void SetBuffer(int16_t* buf, size_t len) {
		playing = false;
		__DMB();
		buffer_ = buf;
		length  = len;
		__DMB();
		playing = true;
	}

	void Stream(float& sample_l, float& sample_r){
			//const double L = (double)length;
			//const double rate = pitch; 

			//double phase = wrap_phase(start + (double)gSamplesElapsed * rate, L);
			float phase = fmodf((float)gSamplesElapsed, (float)length);
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
		
		input_l = preEmphL.process(input_l);
		input_r = preEmphR.process(input_r);

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

};

class FMDemodulator {
	private:
		Phasor carrierPhase;
		int sampleRate;
		BandFilter bandFilter_i;
		BandFilter bandFilter_q;
		LowPassFilter basebandLP_i;   // baseband LP on I after carrier removal
		LowPassFilter basebandLP_q;   // baseband LP on Q after carrier removal
		PLL pll;
		DeEmphasis deEmph;
		LowPassFilter outputFilter;

	public:

	void Init(float sr){
		carrierPhase.Init(sr);
		carrierPhase.SetFreq(5000.0f);
		bandFilter_i.setup(sr);
		bandFilter_q.setup(sr);
		basebandLP_i.setup(sr);
		basebandLP_i.setFrequency(10000.f);
		basebandLP_q.setup(sr);
		basebandLP_q.setFrequency(10000.f);
		pll.init(sr, 300.f);          // kf=300 Hz, fn=500 Hz, zeta=0.707
		deEmph.init(sr);
		outputFilter.setup(sr);
		outputFilter.setFrequency(10000.f);
		sampleRate = sr;
	}

	void SetCarrierFreq(float freq){
		carrierPhase.SetFreq(freq);
		bandFilter_i.setCenterFrequency(freq, 10400.0f, 64);
		bandFilter_q.setCenterFrequency(freq, 10400.0f, 64);
	}

	float Demodulate(float rx_i, float rx_q){
		float phs = carrierPhase.Process();
		float c_i = sinf(TWOPI_F*phs);
		float c_q = cosf(TWOPI_F*phs);

		// IF bandpass
		float fltrx_i = bandFilter_i.process(rx_i);
		float fltrx_q = bandFilter_q.process(rx_q);

		// Carrier removal (complex multiply)
		float zi = fltrx_i * c_i + fltrx_q * c_q;
		float zq = fltrx_q * c_i - fltrx_i * c_q;

		// Baseband lowpass on I and Q
		zi = basebandLP_i.process(zi);
		zq = basebandLP_q.process(zq);

		// PLL discriminator
		float demod = pll.process(zi, zq);

		// De-emphasis
		demod = deEmph.process(demod);

		// Output lowpass + soft limit (no DC block; PLL has no DC)
		return SoftLimit(1.5f * outputFilter.process(demod));
	}
};
