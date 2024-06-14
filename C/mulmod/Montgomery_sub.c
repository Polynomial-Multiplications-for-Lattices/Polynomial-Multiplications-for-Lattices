
#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>

#include "tools.h"

// ================
// This file demonstrates the subtractive variant of the signed Montgomery multiplication.
// Let a and b be the operands that we wish to multiply, Q be the modulus, and R > Q be
// the size of the arithmetic.
// Montgomery multiplication computes a value that is equivalent to a b R^(-1) mod^+- Q.
// If b is known, we replace it with b R mod^+- Q and then Montgomery multiplication computes
// a value that is equivalent to a b mod^+- Q.

// ================
// Theory.
// We compute a b / R - ( a b Q^(-1) mod^+- R) Q / R in the subtractive variant of Montgomery multiplication.
// Observe that a b - ( a b Q^(-1) mod^+- R) Q is equivalent to 0 modulo R and equivalent to a b modulo Q,
// (a b - ( a b Q^(-1) mod^+- R ) Q ) / R is an integer that is equivalent to a b R^(-1) mod^+- Q.
// Moreover, we also know that a b mod^+- R = (a b Q(-1) mod^+- R) Q mod^+- R.
// This implies (a b - ( a b Q^(-1) mod^+- R ) Q ) / R = a b / R - ( a b Q^(-1) mod^+- R) Q / R.
// Finally, a b / R - ( a b Q^(-1) mod^+- R) Q / R is reduced as we have already seen in the
// accumulative variant.

// R = 2^32 below
#define Q 8380417
// RmodQ = R mod^+- Q
#define RmodQ (-4186625)
// Qprime = Q^{-1} mod^+- R
#define Qprime (58728449)

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

// mulhi computes the high part of the long product of a and b.
int32_t mulhi(int32_t a, int32_t b){
    return (int32_t)(((int64_t)a * b) >> 32);
}

// mulhi computes the low part of the long product of a and b.
int32_t mullo(int32_t a, int32_t b){
    return a * b;
}

// The subtractive variant of Montgomery multiplication.
int32_t montgomery_sub_mul(int32_t a, int32_t b, int32_t q, int32_t qprime){

    int32_t lo, hi;

    // hi = a * b / R
    hi = mulhi(a, b);
    // lo = b * Qprime mod^+- R
    lo = mullo(b, qprime);
    // lo = a * b * Qprime mod^+- R
    lo = mullo(a, lo);
    // hi = a * b / R - (a * b * Qprime mod^+- R) Q / R
    hi -= mulhi(lo, q);

    return hi;

}

// The subtractive variant of Montgomery multiplication with precomputation.
int32_t montgomery_sub_mul_pre(int32_t a, int32_t b, int32_t bqprime, int32_t q){

    int32_t lo, hi;

    // hi = a * b / R
    hi = mulhi(a, b);
    // lo = a * (b Qprime mod^+- R) mod^+- R
    lo = mullo(a, bqprime);
    // hi = a * b / R - (a * b Qprime mod^+- R) Q / R
    hi -= mulhi(lo, q);

    return hi;

}

int main(void){

    int32_t a, b, t, q, qprime, rmodq;
    int32_t ref, res;

    q = Q;
    qprime = Qprime;
    rmodq = RmodQ;

    for(size_t i = 0; i < NTESTS; i++){

        // Generate random elements in Z_Q.
        t = rand() % Q;
        coeff_ring.memberZ(&a, &t);
        t = rand() % Q;
        coeff_ring.memberZ(&b, &t);

        // Compute the product of a and b modulo Q.
        coeff_ring.mulZ(&ref, &a, &b);

        // Compute a value equivalent to the product of a and b with the subtractive variant of
        // Montgomery multiplication.
        res = montgomery_sub_mul(a, b, q, qprime);

        // Map the value to Z_Q.
        // Notice that this step is needed only when we want the canonical representations of the
        // values.
        coeff_ring.mulZ(&res, &res, &rmodq);

        // Compare the resulting values.
        assert(ref == res);

    }

    for(size_t i = 0; i < NTESTS; i++){

        // Generate random elements in Z_Q.
        t = rand() % Q;
        coeff_ring.memberZ(&a, &t);
        t = rand() % Q;
        coeff_ring.memberZ(&b, &t);

        // Compute the product of a and b modulo Q.
        coeff_ring.mulZ(&ref, &a, &b);

        // Compute a value equivalent to the product of a and b with the subtractive variant of
        // Montgomery multiplication with precomputation.
        res = montgomery_sub_mul_pre(a, b, mullo(b, qprime), q);

        // Map the value to Z_Q.
        // Notice that this step is needed only when we want the canonical representations of the
        // values.
        coeff_ring.mulZ(&res, &res, &rmodq);

        // Compare the resulting values.
        assert(ref == res);

    }

    printf("Test finished!\n");

}


