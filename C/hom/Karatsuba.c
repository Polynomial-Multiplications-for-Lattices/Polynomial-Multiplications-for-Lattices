
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
// We compute the product of two size-96 polynomials in Z_{2^32}[x].

// ================
// Theory.
// Given two size-n polynomials in R[x], we wish to compute their product in R[x].
// For simplicity, we illustrate the idea when n is even.
// Karatsuba converts the computing task into three polynomial multiplications with input size n/2.
// We illustrate below with the most simple case.

// ================
// The most simple case.
// Consider the case n = 2, we wish to compute (a0 + a1 x) (b0 + b1 x) in R[x].
// For a0 + a1 x, we form the following terms:
//   1. a0
//   2. a0 + a1
//   3. a1
// in R.
// Term 1. has x-degree 0, term 2. has x-degree 1, and term 3. has x-degree 2. We will later justify this.
// For b0 + b1 x, we also form the similar terms and compute the following products:
//   1. a0 b0
//   2. (a0 + a1) (b0 + b1)
//   3. a1 b1
// in R.
// Now, we subtract terms 1. and 3. from 2. and denote the results as follows:
//   1. c0 = a0 b0
//   2. c1 = (a0 + a1) (b0 + b1) - a0 b0
//   3. c2 = a1 b1
// and find c0 + c1 x + c2 x^2 = (a0 + a1 x) (b0 + b1 x) in R[x].
// This is the reason why we associate the x-degrees to each of the terms in the beginning of Karatsuba.
// The association of the x-degrees plays an important role when n is greater than 2 as illustrated in the next example.

// ================
// Another example.
// Goal: compute (a0 + a1 x + a2 x^2 + a3 x^3) (b0 + b1 x + b2 x^2 + b3 x^3) in R[x].
// For a0 + a1 x + a2 x^2 + a3 x^3), we form the following
//   1. a0 + a1 x
//   2. (a0 + a2) + (a1 + a3) x
//   3. a2 + a3 x
// in R[x]
// where term 1. has x-degree 0, term 2. has x-degree 2, and term 3. has x-degree 4.
// We also form the similar terms for b0 + b1 x + b2 x^2 + b3 x^3 and compute the following:
//   1. (a0 + a1 x) (b0 + b1 x)
//   2. ((a0 + a2) + (a1 + a3) x) ((b0 + b2) + (b1 + b3) x)
//   3. (a2 + a3 x) (b2 + b3 x)
// in R[x].
// Next, we subtract terms 1. and 3. from 2. as before and denote the results as follows:
//   1. c0 + c1 x + c0' x^2 = (a0 + a1 x) (b0 + b1 x)
//   2. c2 + c3 x + c2' x^2 = ((a0 + a2) + (a1 + a3) x) ((b0 + b2) + (b1 + b3) x)
//                            - (a0 + a1 x) (b0 + b1 x) - (a2 + a3 x) (b2 + b3 x)
//   3. c4 + c5 x + c4' x^2 = (a2 + a3 x) (b2 + b3 x)
// in R[x].
// To combine them into the target size-7 polynomial, we position the terms at the right positions according to
// the associated x-degree.
// Graphically, we sum up the following rows:
//    c0,  c1, c0',   0,   0,   0,   0
//     0,   0,  c2,  c3, c2',   0,   0
//     0,   0,   0,   0,  c4,  c5, c4'
// One can verify that
// c0 + c1 x + (c2 + c0') x^2 + c3 x^3 + (c4 + c2') x^4 + c5 x^5 + c4' x^6 =
// (a0 + a1 x + a2 x^2 + a3 x^3) (b0 + b1 x + b2 x^2 + b3 x^3) in R[x].

// ================
// Algebraic view.

// ================
// Estimating the cost.
// For multiplying two size-n polynomials with Karatsuba, we decompose the multiplication task into three
// polynomial multiplications of size-(n/2). If we continueously apply Karatsuba until n <= 1,
// eventually, we have 3^log_2 n = n^(log_2 3) multiplications in R.
// Furthermore, we also need Theta(n^(log_2 3)) number of additions/subtractions.
// In summary, we need Theta(n^(log_2 3)) number of multiplications and additions/subtractions and R.

// ================
// Optimization guide.
/*

1. For the recursive Karatsuba,
   instead of computing one layer at a time, try to compute multiple layers at once and save memory operations.

*/

// ================
// Applications to lattice-based cryptosystems.

// ARRAY_N must be even.
#define ARRAY_N 96

// ================
// Z_{2^32}

void memberZ(void *des, const void *src){
    *(int32_t*)des = *(int32_t*)src;
}

void addZ(void *des, const void *src1, const void *src2){
    *(int32_t*)des = (*(int32_t*)src1) + (*(int32_t*)src2);
}

void subZ(void *des, const void *src1, const void *src2){
    *(int32_t*)des = (*(int32_t*)src1) - (*(int32_t*)src2);
}

void mulZ(void *des, const void *src1, const void *src2){
    *(int32_t*)des = (*(int32_t*)src1) * (*(int32_t*)src2);
}

void expZ(void *des, const void *src, size_t e){

    int32_t src_v = *(int32_t*)src;
    int32_t tmp_v;

    tmp_v = 1;
    for(; e; e >>= 1){
        if(e & 1){
            tmp_v = tmp_v * src_v;
        }
        src_v = src_v * src_v;
    }

    memmove(des, &tmp_v, sizeof(int32_t));
}

struct ring coeff_ring = {
    .sizeZ = sizeof(int32_t),
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

    int32_t src1[ARRAY_N], src2[ARRAY_N];
    int32_t ref[2 * ARRAY_N], res[2 * ARRAY_N];

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







