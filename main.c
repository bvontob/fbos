#include "userosc.h"
#include "biquad.h"
#include "delay_line.h"

#define MAX_Q            60
#define DELAY_LINE_SIZE  2048
#define WHITE_NOISE_GAIN 0.1f
#define DELAY_TAPS       9
#define DELAY_TAP_GAIN   ((1.0f - WHITE_NOISE_GAIN) / (float)(DELAY_TAPS))

static biquad bp;
static biquad *bpp = &bp;

static biquad dcf;
static biquad *dcfp = &dcf;

static delay_line dl;
static delay_line *dlp = &dl;

float dl_ram[DELAY_LINE_SIZE] __sdram;

static float q = 20.0f;

void OSC_INIT(uint32_t platform, uint32_t api) {
  (void)platform; (void)api;
  biquad_flush(bpp);
  biquad_so_dc(dcfp, 0.0f);
  biquad_flush(dcfp);
  delay_line_init(dlp, dl_ram, DELAY_LINE_SIZE);
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

  const uint32_t vtap1 = (uint32_t)((float)k_samplerate * 2 / f);
  const uint32_t vtap2 = vtap1 * 2;
  
  biquad_so_bp(bpp,
	       fasttanfullf(M_PI * biquad_wc(f,
					     k_samplerate_recipf)),
	       q);
  
  q31_t * __restrict y = (q31_t *)yn;
  const q31_t * y_e = y + frames;
  
  for (; y != y_e; ) {
    float sig = 0;
    
    sig += clip1m1f(WHITE_NOISE_GAIN * _osc_white()
		    + DELAY_TAP_GAIN * delay_line_read(dlp, vtap1)
		    + DELAY_TAP_GAIN * delay_line_read(dlp, vtap2)
		    + DELAY_TAP_GAIN * delay_line_read(dlp,   283)
		    + DELAY_TAP_GAIN * delay_line_read(dlp,   419)
		    // + DELAY_TAP_GAIN * delay_line_read(dlp,   547)
		    // + DELAY_TAP_GAIN * delay_line_read(dlp,   661)
		    + DELAY_TAP_GAIN * delay_line_read(dlp,   811)
		    // + DELAY_TAP_GAIN * delay_line_read(dlp,   947)
		    + DELAY_TAP_GAIN * delay_line_read(dlp,  1087)
		    + DELAY_TAP_GAIN * delay_line_read(dlp,  1229)
		    // + DELAY_TAP_GAIN * delay_line_read(dlp,  1381)
		    + DELAY_TAP_GAIN * delay_line_read(dlp,  1523)
		    // + DELAY_TAP_GAIN * delay_line_read(dlp,  1663)
		    + DELAY_TAP_GAIN * delay_line_read(dlp,  1823));
    /*    sig += clip1m1f(0.2f * _osc_white()
		    + 0.2f * delay_line_read(dlp, 289)
		    + 0.2f * delay_line_read(dlp, 449)
		    + 0.2f * delay_line_read(dlp, 912)
		    + 0.1f * delay_line_read(dlp, 1333)
		    + 0.1f * delay_line_read(dlp, 2031)); */
    
    sig = clip1m1f(osc_sat_schetzenf(sig));
    sig = clip1m1f(osc_sat_schetzenf(sig));
    sig = clip1m1f(biquad_process_so(bpp, sig));
    // sig = clip1m1f(biquad_process_so(dcfp, sig));
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
    break;

  case k_user_osc_param_shape:
    (void)param_val_to_f32(val);
    break;

  case k_user_osc_param_shiftshape:
    q = M_SQRT2 + (((float)MAX_Q - M_SQRT2) * (1.0f - param_val_to_f32(val)));
    break;
  }
}
