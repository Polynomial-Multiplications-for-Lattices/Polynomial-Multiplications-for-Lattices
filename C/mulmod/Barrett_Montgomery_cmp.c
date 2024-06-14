
#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>

#include "tools.h"

// ================
// This file demonstrates that Barrett multiplication and the accumulative variant of Montgomery
// multiplication compute the same thing with careful choices of integer approximations.

// ================
// Theory.
// Let a and b be the operands that we wish to multiply, Q be the modulus, R > Q be
// the size of the arithmetic, and approx_0 and approx_1 be integer approximations.
// We denote mod^approx_0 and mod^approx_1 the corresponding modular reduction defined by:
// for all z in Z, z mod^approx_0 Q = z - approx_0(z / Q) Q
//                 z mod^apporx_1 Q = z - approx_1(z / Q) Q
// Barrett--Montgomery correspondence states the following:
// a b - approx_1( a approx_0(b R / Q) / R ) Q =
// ( a (b R mod^approx_0 Q) + (a (b R mod^approx_0 Q) Qprime mod^approx_1 R ) Q ) / R.
// We have Barrett multiplication on the left-hand side and the accumulative variant of Montgomery multiplication
// on the right-hand side if we choose
// approx_0 = approx_1 = round.
// We test with this case, but the identity holds for arbitrary integer approximations.

// ================
// Proof.

// R = 2^32 below
#define Q 8380417
// RmodQ = R mod^+- Q
#define RmodQ (-4186625)
// Qprime = -Q^{-1} mod^+- R
#define Qprime (-58728449)

#define NTESTS 1000

// ================
// Definition of Z_Q with signed arithmetic.
// See "tools.h" for explanations.

int32_t mod = Q;

void memberZ(void *des, void *src){
    cmod_int32(des, src, &mod);
}

void addZ(void *des, void *src1, void *src2){
    addmod_int32(des, src1, src2, &mod);
}

void subZ(void *des, void *src1, void *src2){
    submod_int32(des, src1, src2, &mod);
}

void mulZ(void *des, void *src1, void *src2){
    mulmod_int32(des, src1, src2, &mod);
}

void expZ(void *des, void *src, size_t e){
    expmod_int32(des, src, e, &mod);
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
// Multiplication instructions

// mulhi computes the long product of a and b, rounds the result to the high part,
// and returns the high part.
int32_t mulhir(int32_t a, int32_t b){
    return (int32_t)(((int64_t)a * b + 0x80000000) >> 32);
}

// mullo computes the low part of the long product of a and b.
int32_t mullo(int32_t a, int32_t b){
    return a * b;
}

// mullong computes the long product of a and b.
int64_t mullong(int32_t a, int32_t b){
    return (int64_t)a * b;
}

// Return the high part of a.
int32_t gethi(int64_t a){
    return (int32_t)(a >> 32);
}

// round( a R / Q )
// a R mod^+- Q = a R - round(a R / Q) Q
// => round(a R / Q) Q = a R - (a R mod^+- Q)
// => round(a R / Q) = (a R mod^+- Q) * (-Q^(-1) mod^+- R) mod^+- R
int32_t get_barrett_hi(int32_t a, int32_t rmodq, int32_t qprime){

    int32_t t;

    // a * RmodQ mod^+- Q
    coeff_ring.mulZ(&t, &a, &rmodq);

    // (a * RmodQ mod^+- Q) * Qprime mod^+- R
    return t * qprime;

}

// Barrett multiplication.
int32_t barrett_mul(int32_t a, int32_t b, int32_t q, int32_t rmodq, int32_t qprime){

    int32_t lo, hi;
    int32_t bhi;

    // lo = a * b mod^+- R
    lo = a * b;
    // bhi = round(b R / Q)
    bhi = get_barrett_hi(b, rmodq, qprime);
    // hi = round(a * round(b R / Q) / R)
    hi = mulhir(a, bhi);

    // lo = (a * b mod^+- R) - round(a * round(b R / Q) / R) * Q
    return lo - hi * q;

}

// The accumulative variant of Montgomery multiplication.
int32_t montgomery_acc_mul(int32_t a, int32_t b, int32_t q, int32_t qprime){

    int64_t prod;
    int32_t lo;

    // prod = a * b
    prod = mullong(a, b);
    // lo = a * b * Qprime mod^+- R
    lo = mullo(prod, qprime);
    // prod = a * b + (a * b * Qprime mod^+- R) * Q
    prod += mullong(lo, q);

    // prod = (a * b + (a * b * Qprime mod^+- R) Q) / R
    return gethi(prod);

}

int main(void){

    int32_t a, b, t, q, qprime, rmodq;
    int32_t res_montgomery, res_barrett;

    q = Q;
    qprime = Qprime;
    rmodq = RmodQ;

    for(size_t i = 0; i < NTESTS; i++){

        // Generate random elements in Z_Q.
        t = rand() % Q;
        coeff_ring.memberZ(&a, &t);
        t = rand() % Q;
        coeff_ring.memberZ(&b, &t);

        // Precompute b R mod^+- Q.
        coeff_ring.mulZ(&t, &b, &rmodq);

        // Call Barrett multiplication.
        res_barrett = barrett_mul(a, b, q, rmodq, qprime);
        // Call the accumulative variant of Montgomery multiplication.
        res_montgomery = montgomery_acc_mul(a, t, q, qprime);

        // Compare if the results of Barrett and Montgomery multiplications are the same.
        assert(res_montgomery == res_barrett);

    }

    printf("Test finished!\n");

}


