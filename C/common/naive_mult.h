#ifndef NAIVE_MULT
#define NAIVE_MULT

#include <stddef.h>

#include "tools.h"

// ================================

// Multiplying size-len polynomials stored at src1 and src2 in in R[x] / (x^len - twiddle)
// where R = ring.
// The resulting polynomial is stored at des.
void naive_mulR(
    void *des,
    const void *src1, const void *src2,
    size_t len, const void *twiddle,
    struct ring ring
    );

// Multiplying size-len polynomials stored at src1 and src2 in R[x] where R = ring.
// The resulting polynomial is stored at des.
void naive_mul_long(
    void *des,
    const void *src1, const void *src2,
    size_t len,
    struct ring ring
    );

// Point-wise multiplication of src1[len * jump] by src2[len] over R where R = ring.
// In particular, for i in {0, ..., len - 1} and j in {0, ..., jump - 1},
// src1[i * jump + j] is multiplied by src2[i].
// The destination contains the same number of elements as src1.
void point_mul(
    void *des,
    const void *src1, const void *src2,
    size_t len, size_t jump,
    struct ring ring
    );

#endif

