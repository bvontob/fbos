#include "float_math.h"

typedef struct biquad_coeffs_s { /* Default constructor set all to 0.0f */
  float ff0;
  float ff1;
  float ff2;
  float fb1;
  float fb2;
} biquad_coeffs;

typedef struct biquad_s { /* Default constructor can use flush() to 0.0f out both */
  biquad_coeffs c;
  float mZ1;
  float mZ2;
} biquad;

/**
 * Convert Hz frequency to radians
 *
 * @param   fc Frequency in Hz
 * @param   fsrecip Reciprocal of sampling frequency (1/Fs)
 */
static inline __attribute__((optimize("Ofast"),always_inline))
float biquad_wc(const float fc, const float fsrecip) {
  return fc * fsrecip;
}

/**
 * Calculate coefficients for single pole low pass filter.
 *
 * @param   pole Pole position in radians
 */
inline __attribute__((optimize("Ofast"),always_inline))
void biquad_pole_lp(biquad *bq, const float pole) {
  bq->c.ff0 = 1.f - pole;
  bq->c.fb1 = -pole;
  bq->c.fb2 = bq->c.ff2 = bq->c.ff1 = 0.f;
}

/**
 * Calculate coefficients for single pole high pass filter.
 *
 * @param   pole Pole position in radians
 */
inline __attribute__((optimize("Ofast"),always_inline))
void biquad_pole_hp(biquad *bq, const float pole) {
  bq->c.ff0 = 1.f - pole;
  bq->c.fb1 = pole;
  bq->c.fb2 = bq->c.ff2 = bq->c.ff1 = 0.f;
}

/**
 * Calculate coefficients for single pole DC filter.
 *
 * @param   pole Pole position in radians
 */
inline __attribute__((optimize("Ofast"),always_inline))
void biquad_fo_dc(biquad *bq, const float pole) {
  bq->c.ff0 = 1.f;
  bq->c.ff1 = -1.f;
  bq->c.fb1 = -pole;
  bq->c.fb2 = bq->c.ff2 = 0.f;
}

/**
 * Calculate coefficients for first order low pass filter.
 *
 * @param   k Tangent of PI x cutoff frequency in radians: tan(pi*wc)
 */
inline __attribute__((optimize("Ofast"),always_inline))
void biquad_fo_lp(biquad *bq, const float k) {
  const float kp1 = k + 1.f;
  const float km1 = k - 1.f;
  bq->c.ff0 = bq->c.ff1 = k / kp1;
  bq->c.fb1 = km1 / kp1;
  bq->c.fb2 = bq->c.ff2 = 0.f;
}

/**
 * Calculate coefficients for first order high pass filter.
 *
 * @param   k Tangent of PI x cutoff frequency in radians: tan(pi*wc)
 */
inline __attribute__((optimize("Ofast"),always_inline))
void biquad_fo_hp(biquad *bq, const float k) {
  // k = tan(pi*wc)
  const float kp1 = k + 1.f;
  const float km1 = k - 1.f;
  bq->c.ff0 = 1.f / kp1;
  bq->c.ff1 = -bq->c.ff0;
  bq->c.fb1 = km1 / kp1;
  bq->c.fb2 = bq->c.ff2 = 0.f;
}

/**
 * Calculate coefficients for first order all pass filter.
 *
 * @param   k Tangent of PI x cutoff frequency in radians: tan(pi*wc)
 */
inline __attribute__((optimize("Ofast"),always_inline))
void biquad_fo_ap(biquad *bq, const float k) {
  // k = tan(pi*wc)
  const float kp1 = k + 1.f;
  const float km1 = k - 1.f;
  bq->c.ff0 = bq->c.fb1 = km1 / kp1;
  bq->c.ff1 = 1.f;
  bq->c.fb2 = bq->c.ff2 = 0.f;
}

/**
 * Calculate coefficients for first order all pass filter.
 *
 * @param   wc cutoff frequency in radians
 *
 * @note Alternative implementation with no tangeant lookup
 */
inline __attribute__((optimize("Ofast"),always_inline))
void biquad_fo_ap2(biquad *bq, const float wc) {
  // Note: alternative implementation for use in phasers
  const float g1 = 1.f - wc;
  bq->c.ff0 = g1;
  bq->c.ff1 = -1;
  bq->c.fb1 = -g1;
  bq->c.fb2 = bq->c.ff2 = 0.f;
}

/**
 * Calculate coefficients for second order DC filter.
 *
 * @param   pole Pole position in radians
 */
inline __attribute__((optimize("Ofast"),always_inline))
void biquad_so_dc(biquad *bq, const float pole) {
  bq->c.ff0 = bq->c.ff2 = 1.f;
  bq->c.ff1 = 2.f;
  bq->c.fb1 = -2.f * pole;
  bq->c.fb2 = pole * pole;
}

/**
 * Calculate coefficients for second order low pass filter.
 *
 * @param   k Tangent of PI x cutoff frequency in radians: tan(pi*wc)
 * @param   q Resonance with flat response at q = sqrt(2)
 */
inline __attribute__((optimize("Ofast"),always_inline))
void biquad_so_lp(biquad *bq, const float k, const float q) {
  // k = tan(pi*wc)
  // flat response at q = sqrt(2)
  const float qk2 = q * k * k;
  const float qk2_k_q_r = 1.f / (qk2 + k + q);
  bq->c.ff0 = bq->c.ff2 = qk2 * qk2_k_q_r;
  bq->c.ff1 = 2.f * bq->c.ff0;
  bq->c.fb1 = 2.f * (qk2 - q) * qk2_k_q_r;
  bq->c.fb2 = (qk2 - k + q) * qk2_k_q_r;
}

/**
 * Calculate coefficients for second order high pass filter.
 *
 * @param   k Tangent of PI x cutoff frequency in radians: tan(pi*wc)
 * @param   q Resonance with flat response at q = sqrt(2)
 */
inline __attribute__((optimize("Ofast"),always_inline))
void biquad_so_hp(biquad *bq, const float k, const float q) {
  // k = tan(pi*wc)
  // flat response at q = sqrt(2)
  const float qk2 = q * k * k;
  const float qk2_k_q_r = 1.f / (qk2 + k + q);
  bq->c.ff0 = bq->c.ff2 = q * qk2_k_q_r;
  bq->c.ff1 = -2.f * bq->c.ff0;
  bq->c.fb1 = 2.f * (qk2 - q) * qk2_k_q_r;
  bq->c.fb2 = (qk2 - k + q) * qk2_k_q_r;
}

/**
 * Calculate coefficients for second order band pass filter.
 *
 * @param   k Tangent of PI x cutoff frequency in radians: tan(pi*wc)
 * @param   q Resonance with flat response at q = sqrt(2)
 */
inline __attribute__((optimize("Ofast"),always_inline))
void biquad_so_bp(biquad *bq, const float k, const float q) {
  // k = tan(pi*wc)
  // q is inverse of relative bandwidth (Fc / Fb)
  const float qk2 = q * k * k;
  const float qk2_k_q_r = 1.f / (qk2 + k + q);
  bq->c.ff0 = k * qk2_k_q_r;
  bq->c.ff1 = 0.f;
  bq->c.ff2 = -bq->c.ff0;
  bq->c.fb1 = 2.f * (qk2 - q) * qk2_k_q_r;
  bq->c.fb2 = (qk2 - k + q) * qk2_k_q_r;
}

/**
 * Calculate coefficients for second order band reject filter.
 *
 * @param   k Tangent of PI x cutoff frequency in radians: tan(pi*wc)
 * @param   q Resonance with flat response at q = sqrt(2)
 */
inline __attribute__((optimize("Ofast"),always_inline))
void biquad_so_br(biquad *bq, const float k, const float q) {
  // k = tan(pi*wc)
  // q is inverse of relative bandwidth (Fc / Fb)
  const float qk2 = q * k * k;
  const float qk2_k_q_r = 1.f / (qk2 + k + q);
  bq->c.ff0 = bq->c.ff2 = (qk2 + q) * qk2_k_q_r;
  bq->c.ff1 = bq->c.fb1 = 2.f * (qk2 - q) * qk2_k_q_r;
  bq->c.fb2 = (qk2 - k + q) * qk2_k_q_r;
}

/**
 * Calculate coefficients for second order all pass filter.
 *
 * @param   k Tangent of PI x cutoff frequency in radians: tan(pi*wc)
 * @param   q Inverse of relative bandwidth (Fc / Fb)
 */
inline __attribute__((optimize("Ofast"),always_inline))
void biquad_so_ap(biquad *bq, const float k, const float q) {
  // k = tan(pi*wc)
  // q is inverse of relative bandwidth (Fc / Fb)
  const float qk2 = q * k * k;
  const float qk2_k_q_r = 1.f / (qk2 + k + q);
  bq->c.ff0 = bq->c.fb2 = (qk2 - k + q) * qk2_k_q_r;
  bq->c.ff1 = bq->c.fb1 = 2.f * (qk2 - q) * qk2_k_q_r;
  bq->c.ff2 = 1.f;
}

/**
 * Calculate coefficients for second order all pass filter.
 *
 * @param   delta cos(2pi*wc)
 * @param   gamma tan(pi * wb)
 *
 * @note q is inverse of relative bandwidth (wc / wb)
 * @note Alternative implementation, so called "tunable" in DAFX second edition.
 */
inline __attribute__((optimize("Ofast"),always_inline))
void biquad_so_ap2(biquad *bq, const float delta, const float gamma) {
  // Note: Alternative implementation .. so called "tunable" in DAFX.
  // delta = cos(2pi*wc)
  const float c = (gamma - 1.f) / (gamma + 1.f);
  const float d = -delta;
  bq->c.ff0 = bq->c.fb2 = -c;
  bq->c.ff1 = bq->c.fb1 = d * (1.f - c);
  bq->c.ff2 = 1.f;
}

/**
 * Calculate coefficients for second order all pass filter.
 *
 * @param   delta cos(2pi*wc)
 * @param   radius 
 *
 * @note Another alternative implementation.
 */
inline __attribute__((optimize("Ofast"),always_inline))
void biquad_so_ap3(biquad *bq, const float delta, const float radius) {
  // Note: alternative implementation for use in phasers
  // delta = cos(2pi * wc)
  const float a1 = -2.f * radius * delta;
  const float a2 = radius * radius;
  bq->c.ff0 = bq->c.fb2 = a2;
  bq->c.ff1 = bq->c.fb1 = a1;
  bq->c.ff2 = 1.f;
}

/**
 * Flush internal delays
 */
inline __attribute__((optimize("Ofast"),always_inline))
void biquad_flush(biquad *bq) {
  bq->mZ1 = bq->mZ2 = 0;
}

/**
 * Second order processing of one sample
 *
 * @param xn  Input sample
 *
 * @return Output sample
 */
inline __attribute__((optimize("Ofast"),always_inline))
float biquad_process_so(biquad *bq, const float xn) {
  float acc = bq->c.ff0 * xn + bq->mZ1;
  bq->mZ1   = bq->c.ff1 * xn + bq->mZ2;
  bq->mZ2   = bq->c.ff2 * xn;
  bq->mZ1  -= bq->c.fb1 * acc;
  bq->mZ2  -= bq->c.fb2 * acc;
  return acc;
}

/**
 * First order processing of one sample
 *
 * @param xn  Input sample
 *
 * @return Output sample
 */
inline __attribute__((optimize("Ofast"),always_inline))
float biquad_process_fo(biquad *bq, const float xn) {
  float acc = bq->c.ff0 * xn + bq->mZ1;
  bq->mZ1   = bq->c.ff1 * xn;
  bq->mZ1  -= bq->c.fb1 * acc;
  return acc;
}
