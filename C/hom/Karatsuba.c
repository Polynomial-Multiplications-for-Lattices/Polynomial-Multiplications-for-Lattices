
#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>

#include "tools.h"
#include "naive_mult.h"

// ================
// This file demonstrates recursive Karatsuba with symmetric inputs.
// We compute the product of two size-96 polynomials in Z_{2^64}[x].

// ================
// Theory.

// ================
// A small example.

// ================
// Optimization guide.
/*

1. Instead of computing one layer at a time, try to compute multiple ones and save memory operations.

*/

// ================
// Applications to lattice-based cryptosystems.

#define ARRAY_N 96

// ================
// Z_{2^64}

void memberZ(void *des, void *src){
    *(uint64_t*)des = *(uint64_t*)src;
}

void addZ(void *des, void *src1, void *src2){
    *(uint64_t*)des = (*(uint64_t*)src1) + (*(uint64_t*)src2);
}

void subZ(void *des, void *src1, void *src2){
    *(uint64_t*)des = (*(uint64_t*)src1) - (*(uint64_t*)src2);
}

void mulZ(void *des, void *src1, void *src2){
    *(uint64_t*)des = (*(uint64_t*)src1) * (*(uint64_t*)src2);
}

void expZ(void *des, void *src, size_t e){

    uint64_t src_v = *(uint64_t*)src;
    uint64_t tmp_v;

    tmp_v = 1;
    for(; e; e >>= 1){
        if(e & 1){
            tmp_v = tmp_v * src_v;
        }
        src_v = src_v * src_v;
    }

    memmove(des, &tmp_v, sizeof(uint64_t));
}

struct ring coeff_ring = {
    .sizeZ = sizeof(uint64_t),
    .memberZ = memberZ,
    .addZ = addZ,
    .subZ = subZ,
    .mulZ = mulZ,
    .expZ = expZ
};

// ================

// len must be even.
static
void karatsuba_eval(void *des, void *src, size_t len, struct ring ring){

    for(size_t i = 0; i < (len / 2); i++){
        ring.addZ(des + i * ring.sizeZ, src + i * ring.sizeZ, src + ((len / 2) + i) * ring.sizeZ);
    }

}

// len must be even.
static
void karatsuba_interpol(void *des, void *src, size_t len, struct ring ring){

    // Interpolation.
    for(size_t i = 0; i < len - 1; i++){
        ring.subZ(src + i * ring.sizeZ, src + i * ring.sizeZ, des + i * ring.sizeZ);
        ring.subZ(src + i * ring.sizeZ, src + i * ring.sizeZ, des + (len + i) * ring.sizeZ);
    }

    // Sum up the overlapped parts.
    for(size_t i = 0; i < len - 1; i++){
        ring.addZ(des + ((len / 2) + i) * ring.sizeZ, des + ((len / 2) + i) * ring.sizeZ, src + i * ring.sizeZ);
    }

}

// threshold | len,
// len / threshold must be a power of two.
static
void karatsuba_recur(void *des, void *src1, void *src2, size_t len, size_t threshold, struct ring ring){

    // If len <= threshold, we apply the naive long multiplication.
    if(len <= threshold){
        naive_mul_long(des, src1, src2, len, ring);
        return;
    }

    // Declare buffers for the evaluation.
    char src1mid[(len / 2) * ring.sizeZ], src2mid[(len / 2) * ring.sizeZ];
    char desmid[(len - 1) * ring.sizeZ];

    // Evaluating half-size polynomials at 1.
    karatsuba_eval(src1mid, src1, len, ring);
    karatsuba_eval(src2mid, src2, len, ring);

    memset(des, 0, (2 * len - 1) * ring.sizeZ);

    // Karatsuba for the point 0.
    karatsuba_recur(des, src1, src2, len / 2, threshold, ring);
    // Karatsuba for the point \infty.
    karatsuba_recur(des + len * ring.sizeZ, src1 + (len / 2) * ring.sizeZ, src2 + (len / 2) * ring.sizeZ, len / 2, threshold, ring);
    // Karatsuba for the point 1.
    karatsuba_recur(desmid, src1mid, src2mid, len / 2, threshold, ring);

    // Apply Karatsuba interpolation and sum up the overlapped parts.
    karatsuba_interpol(des, desmid, len, ring);

}

int main(void){

    uint64_t src1[ARRAY_N], src2[ARRAY_N];
    uint64_t ref[2 * ARRAY_N], res[2 * ARRAY_N];

    for(size_t i = 0; i < ARRAY_N; i++){
        src1[i] = rand();
        src2[i] = rand();
    }

    // Compute the reference.
    naive_mul_long(ref, src1, src2, ARRAY_N, coeff_ring);

    // Apply recursive Karatsuba.
    karatsuba_recur(res, src1, src2, ARRAY_N, 6, coeff_ring);

    for(size_t i = 0; i < 2 * ARRAY_N - 1; i++){
        assert(ref[i] == res[i]);
    }

    printf("Test finished!\n");

}







