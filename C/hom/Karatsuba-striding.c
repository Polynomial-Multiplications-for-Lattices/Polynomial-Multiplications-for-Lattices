
#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>

#include "tools.h"
#include "naive_mult.h"

// ================
// This file demonstrates 2-layer Karatsuba with symmetric inputs.
// We compute the product of two size-256 polynomials in Z_{2^32}[x].

// ================
// Theory.
// We recall that Karatsuba computes the product (a0 + a1 x) (b0 + b1 x) in R[x] from the products
// 1. a0 b0
// 2. a1 b1
// 3. (a0 + a1) (b0 + b1)
// in R.
// We can apply the same idea to products in R[x] / (x^(2^k) + 1) and generally, any polynomial ring of the form R[x] / (f(x^2)).
// Recall that the polynomial ring R[x] / (f(x^2)) contains the polynomial ring R[y] / (f(y)) as a subring via the embedding y -> x^2.
// Conceptually, if we have an efficient polynomial multiplication in R[y] / (f(y)), then we have an efficient polynomial multiplication
// in R[x] / (f(x^2)) with an extension of the inversion of the embedding y -> x^2.
// Formally, we introduce the relation x^2 - y and rewrite the polynomial ring R[x] / (f(x^2)) as
// (R[y] / (f(y)) )[x] / (x^2 - y).
// If we apply Karatsuba in x, we have three polynomial multiplications in R[y] / (f(y)).
// This file demonstrates the implementation for f(x) = x^ARRAY_N + 1 with an even ARRAY_N.

// ================
// A small example.
// Goal: compute (a0 + a1 x + a2 x^2 + a3 x^3) (b0 + b1 x + b2 x^2 + b3 x^3) in R[x] / (x^4 + 1).
// Map
// a0 + a1 x + a2 x^2 + a3 x^3 in R[x] / (x^4 + 1)
// to a0 + a1 x + y(a2 + a3 x) in (R[y] / (y^2 + 1) ) / [x] / (x^2 - y)
// to a0 + a2 y + (a1 + a3 y) x in (R[y] / (y^2 + 1) ) / [x] / (x^2 - y)
// Apply Karatsuba in x, we have three terms:
//     1. a0 + a2 y
//     2. a1 + a3 y
//     3. (a0 + a1) + (a2 + a3) y
// Do the same for b0 + b1 x + b2 x^2 + b3 x^3 and multiply the corresponding terms.
// We have three products:
//     1. (a0 + a2 y) (b0 + b2 y)
//     2. (a1 + a3 y) (b1 + b3 y)
//     3. ((a0 + a1) + (a2 + a3) y) ((b0 + b1) + (b2 + b3) y)
// in R[y] / (y^2 + 1).
// Inverting Karatsuba gives:
//     1. c0 + c1 y = (a0 + a2 y) (b0 + b2 y)
//     2. c2 + c3 y = (a1 + a3 y) (b1 + b3 y)
//     3. c4 + c5 y = ((a0 + a1) + (a2 + a3) y) ((b0 + b1) + (b2 + b3) y) - (c0 + c1 y) - (c2 + c3 y)
// We combine the terms by summing up the following rows:
// c0, 0, c1, 0
// 0, c4, 0, c5
// -c3, 0, c2, 0
// and get
// c0 - c3, c4, c1 + c2, c5.
// One can verify that this is the product
// (a0 + a1 x + a2 x^2 + a3 x^3) (b0 + b1 x + b2 x^2 + b3 x^3) in R[x] / (x^4 + 1).

// ================
// Algebraic view.

// ================
// Optimization guide.

// ================
// Applications to lattice-based cryptosystems.

// ARRAY_N must be even.
#define ARRAY_N 192

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

static
void negacyclic_Karatsuba_striding_recur(int32_t *des, const int32_t *src1, const int32_t *src2, size_t len, size_t threshold){

    const int32_t twiddle = -1;

    if(len <= threshold){
        naive_mulR(des, src1, src2, len, &twiddle, coeff_ring);
        return;
    }

    int32_t src1lo[len / 2], src1hi[len / 2], src1mid[len / 2];
    int32_t src2lo[len / 2], src2hi[len / 2], src2mid[len / 2];
    int32_t reslo[len / 2], reshi[len / 2], resmid[len / 2];

    for(size_t i = 0; i < len / 2; i++){
        src1lo[i] = src1[2 * i + 0];
        src1hi[i] = src1[2 * i + 1];
        src1mid[i] = src1lo[i] + src1hi[i];
    }
    for(size_t i = 0; i < len / 2; i++){
        src2lo[i] = src2[2 * i + 0];
        src2hi[i] = src2[2 * i + 1];
        src2mid[i] = src2lo[i] + src2hi[i];
    }

    negacyclic_Karatsuba_striding_recur(reslo, src1lo, src2lo, len / 2, threshold);
    negacyclic_Karatsuba_striding_recur(reshi, src1hi, src2hi, len / 2, threshold);
    negacyclic_Karatsuba_striding_recur(resmid, src1mid, src2mid, len / 2, threshold);

    for(size_t i = 0; i < len / 2; i++){
        resmid[i] = resmid[i] - reslo[i] - reshi[i];
    }

    for(size_t i = 0; i < len / 2; i++){
        des[2 * i + 0] = reslo[i];
        des[2 * i + 1] = resmid[i];
    }

    des[0] -= reshi[len / 2 - 1];
    for(size_t i = 1; i < len / 2; i++){
        des[2 * i] += reshi[i - 1];
    }

}


int main(void){

    int32_t poly1[ARRAY_N], poly2[ARRAY_N];
    int32_t ref[ARRAY_N], res[ARRAY_N];

    const int32_t twiddle = -1;

    for(size_t i = 0; i < ARRAY_N; i++){
        poly1[i] = rand();
        poly2[i] = rand();
    }

    // Compute the product in Z_{2^32}[x] / (x^ARRAY_N + 1).
    naive_mulR(ref, poly1, poly2, ARRAY_N, &twiddle, coeff_ring);

    // Compute the product in Z_{2^32} [x] / (x^ARRAY_N + 1) via striding followed by two layers of Karatsuba.
    negacyclic_Karatsuba_striding_recur(res, poly1, poly2, ARRAY_N, 4);

    // Test for correctness.
    for(size_t i = 0; i < ARRAY_N; i++){
        assert(ref[i] == res[i]);
    }

    printf("Test finished!\n");


}




