#include "userosc.h"
#include "biquad.h"
#include "delay_line.h"

#define MIN_Q            M_SQRT2
#define MAX_Q            60
#define DELAY_LINE_SIZE  2048
#define DELAY_TAPS       7

static biquad bp;
static biquad *bpp = &bp;

static biquad bp2;
static biquad *bp2p = &bp2;

static delay_line dl;
static delay_line *dlp = &dl;

float dl_ram[DELAY_LINE_SIZE] __sdram;

static struct parameters {
  float reso;
  float ovrm;
  float ovrq;
  int   ovrt;
  float nois;
  float q;
} p;

void OSC_INIT(uint32_t platform, uint32_t api) {
  (void)platform; (void)api;
  
  biquad_flush(bpp);
  biquad_flush(bp2p);
  delay_line_init(dlp, dl_ram, DELAY_LINE_SIZE);
  
  /*
   * Agitate the feedback loop with a click -- or we won't get any
   * sound at all if the noise is turned down in the mix.
   */
  delay_line_write(dlp, 0.001f);
}

void OSC_CYCLE(const user_osc_param_t * const params,
               int32_t *yn,
               const uint32_t frames) {
  const uint8_t note = (params->pitch) >> 8;
  const uint8_t mod  = (params->pitch) & 0x00FF;

  const float f0 = osc_notehzf(note);
  const float f1 = osc_notehzf(note + 1);
  const float f  = clipmaxf(linintf(mod * k_note_mod_fscale, f0, f1),
			    k_note_max_hz);
  const float fo = clipmaxf((float)(p.ovrt) + 2 * f, k_note_max_hz);
  
  const uint32_t vtap1 = (uint32_t)((float)k_samplerate * 2 / f);
  const uint32_t vtap2 = vtap1 * 2;
  
  biquad_so_bp(bpp,
	       fasttanfullf(M_PI * biquad_wc(f,
					     k_samplerate_recipf)),
	       p.q);
  biquad_so_bp(bp2p,
	       fasttanfullf(M_PI * biquad_wc(fo,
					     k_samplerate_recipf)),
	       p.ovrq);

  const float noise_gain     = p.nois;
  const float reso_tap_gain  = p.reso / 2.0f;
  const float delay_tap_gain = ((1.0f - p.nois - p.reso)
				/ (float)(DELAY_TAPS));
  
  q31_t * __restrict y = (q31_t *)yn;
  const q31_t * y_e = y + frames;
  
  for (; y != y_e; ) {
    float sig = 0;
   
    sig += clip1m1f(noise_gain * _osc_white()
		    + reso_tap_gain  * delay_line_read(dlp, vtap1)
		    + reso_tap_gain  * delay_line_read(dlp, vtap2)
		    + delay_tap_gain * delay_line_read(dlp,   283)
		    + delay_tap_gain * delay_line_read(dlp,   419)
		    + delay_tap_gain * delay_line_read(dlp,   811)
		    + delay_tap_gain * delay_line_read(dlp,  1087)
		    + delay_tap_gain * delay_line_read(dlp,  1229)
		    + delay_tap_gain * delay_line_read(dlp,  1523)
		    + delay_tap_gain * delay_line_read(dlp,  1823));
    
    sig = clip1m1f(osc_sat_schetzenf(sig));
    sig = clip1m1f(osc_sat_schetzenf(sig));
    const float sigo = clip1m1f(biquad_process_so(bp2p, sig));
    sig = clip1m1f((1.0f - p.ovrm) * biquad_process_so(bpp, sig)
		   + p.ovrm * sigo);
    sig = clip1m1f(osc_sat_schetzenf(sig));
    sig = clip1m1f(osc_sat_schetzenf(sig));

    delay_line_write(dlp, sig);
    
    *(y++) = f32_to_q31(sig);
  }
}

void OSC_NOTEON(const user_osc_param_t * const params) {
  (void)params;
}

void OSC_NOTEOFF(const user_osc_param_t * const params) {
  (void)params;
}

void OSC_PARAM(uint16_t idx, uint16_t val) { 
  switch (idx) {

  case k_user_osc_param_id1:
    p.reso = (float)val / 100.0f * 0.5f;
    break;

  case k_user_osc_param_id2:
    p.ovrm = (float)val / 100.0f;
    break;

  case k_user_osc_param_id3:
    p.ovrq = ((float)MIN_Q
	      + (((float)MAX_Q - (float)MIN_Q)
		 * ((float)val / 100.0f)));
    break;

  case k_user_osc_param_id4:
    p.ovrt = val + 2;
    break;

  case k_user_osc_param_id5:
    p.nois = (float)val / 100.0f * 0.5f;
    break;

  case k_user_osc_param_shape:
    p.q = ((float)MIN_Q
	   + (((float)MAX_Q - (float)MIN_Q)
	      * (1.0f - param_val_to_f32(val))));
    break;

  case k_user_osc_param_shiftshape:
    (void)param_val_to_f32(val);
    break;
  }
}
