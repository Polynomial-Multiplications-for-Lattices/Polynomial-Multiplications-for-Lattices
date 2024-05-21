#ifndef TOOLS_H
#define TOOLS_H

#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

// ================================
// Introduce this file...

struct compress_profile{
    size_t array_n;
    size_t ntt_n;
    size_t log_ntt_n;
    size_t compressed_layers;
    size_t merged_layers[16];
};

// This structure
struct commutative_ring{
    // sizeZ is refers to the size in bytes of an element in the ring.
    size_t sizeZ;
    // memberZ maps *src to its representative in the ring and store the result in des.
    void (*memberZ)(void *des, void *src);
    // addZ adds up *src1 and *src2, and stores the result in des.
    void (*addZ)(void *des, void *src1, void *src2);
    // subZ subtract *src2 from *src1, and stores the result in des.
    void (*subZ)(void *des, void *src1, void *src2);
    // mulZ multiplies *src1 and *src2, and stores the result in des.
    void (*mulZ)(void *des, void *src1, void *src2);
    // expZ computes (*src)^e, and stores the result in des.
    void (*expZ)(void *des, void *src, size_t e);
};

// ================================
// We also provide several commonly used functions.

// ================================
// Reducing elements to mod mod_v.

void cmod_int16(void *des, void *src, void *mod);
void cmod_int32(void *des, void *src, void *mod);
void cmod_int64(void *des, void *src, void *mod);

// ================================
// Addition in mod mov_v.

void addmod_int16(void *des, void *src1, void *src2, void *mod);
void addmod_int32(void *des, void *src1, void *src2, void *mod);

// ================================
// Subtraction in mod mov_v.

void submod_int16(void *des, void *src1, void *src2, void *mod);
void submod_int32(void *des, void *src1, void *src2, void *mod);

// ================================
// Multiplication in mod mov_v.

void mulmod_int16(void *des, void *src1, void *src2, void *mod);
void mulmod_int32(void *des, void *src1, void *src2, void *mod);

// ================================
// Exponentiation in mod mov_v.

void expmod_int16(void *des, void *src, size_t e, void *mod);
void expmod_int32(void *des, void *src, size_t e, void *mod);

// ================================
// In-place bit-reversing the array.
void bitreverse(void *src, size_t len, size_t size);

#endif

