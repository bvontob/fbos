#include "float_math.h"
#include "int_math.h"
#include "buffer_ops.h"

#ifndef __sdram
#define __sdram __attribute__((section(".sdram")))
#endif

typedef struct delay_line_s {
  float   *line;
  float    fracz;
  size_t   size;
  size_t   mask;
  uint32_t wr_idx;
} delay_line;
  
static inline __attribute__((optimize("Ofast"), always_inline))
void delay_line_clear(delay_line *dl) {
  buf_clr_f32((float *)dl->line, dl->size);
}

/**
 * Constructor with explicit memory area to use as backing buffer for
 * delay line.
 *
 * @param ram Pointer to memory buffer
 * @param line_size Size in float of memory buffer (must be power of 2)
 */
static inline __attribute__((optimize("Ofast"), always_inline))
void delay_line_init(delay_line *dl, float *ram, size_t line_size) {
  dl->line   = ram;
  dl->fracz  = 0;
  dl->size   = nextpow2_u32(line_size);
  dl->mask   = dl->size - 1;
  dl->wr_idx = 0;
  delay_line_clear(dl);
}

/**
 * Write a single sample to the head of the delay line
 *
 * @param s Sample to write
 */
static inline __attribute__((optimize("Ofast"), always_inline))
void delay_line_write(delay_line *dl, const float s) {
  dl->line[(dl->wr_idx--) & dl->mask] = s;
}

/**
 * Read a single sample from the delay line at given position from
 * current write index.
 *
 * @param pos Offset from write index
 * @return Sample at given position from write index
 */
static inline __attribute__((optimize("Ofast"), always_inline))
float delay_line_read(const delay_line *dl, const uint32_t pos) {
  return dl->line[(dl->wr_idx + pos) & dl->mask];
}

/**
 * Read a sample from the delay line at a fractional position from
 * current write index.
 *
 * @param pos Offset from write index as floating point.
 * @return Interpolated sample at given fractional position from write index
 */
static inline __attribute__((optimize("Ofast"), always_inline))
float delay_line_read_frac(const delay_line *dl, const float pos) {
  const uint32_t base = (uint32_t)pos;
  const float frac = pos - base;
  const float s0 = delay_line_read(dl, base);
  const float s1 = delay_line_read(dl, base + 1);
  return linintf(frac, s0, s1);
}

/**
 * Read a sample from the delay line at a position from current write
 * index with interpolation from last read.
 *
 * @param pos Offset from write index
 * @param frac Interpolation from last read pair.
 * @return Interpolation of last read sample and sample at given position from write index.
 */
static inline __attribute__((optimize("Ofast"), always_inline))
float delay_line_read_fracz(delay_line *dl,
			    const uint32_t pos,
			    const float frac) {
  const float s0 = delay_line_read(dl, pos);
  const float y = linintf(frac, s0, dl->fracz);
  dl->fracz = s0;
  return y;
}
