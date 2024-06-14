
#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>

#include "tools.h"

// ================
// This file demonstrates the accumulative variant of the signed Montgomery multiplication.
// Let a and b be the operands that we wish to multiply, Q be the modulus, and R > Q be
// the size of the arithmetic.
// Montgomery multiplication computes a value that is equivalent to a b R^(-1) mod^+- Q.
// If b is known, we replace it with b R mod^+- Q and then Montgomery multiplication computes
// a value that is equivalent to a b mod^+- Q.

// ================
// Theory.
// Observe that a b + (- a b Q^(-1) mod^+- R) Q is
// equivalent to 0 modulo R and
// equivalent to a b modulo Q,
// (a b + ( - a b Q^(-1) mod^+- R ) Q ) / R is an integer that is
// equivalent to a b R^(-1) mod^+- Q.
// It remains to show that (a b + ( - a b Q^(-1) mod^+- R ) Q ) / R is reduced.
// Taking the absolute value yields the upper bound Q / 2 + |a b| / R.
// If |b| > Q /2, we have Q / 2 (1 + |a| / R).

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

// mullong computes the long product of a and b.
int64_t mullong(int32_t a, int32_t b){
    return (int64_t)a * b;
}

// mullo computes the low part of the long product of a and b.
int32_t mullo(int32_t a, int32_t b){
    return a * b;
}

// Return the high part of a.
int32_t gethi(int64_t a){
    return (int32_t)(a >> 32);
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

        // Compute a value equivalent to the product of a and b with the accumulative variant of
        // Montgomery multiplication.
        res = montgomery_acc_mul(a, b, q, qprime);

        // Map the value to Z_Q.
        // Notice that this step is needed only when we want the canonical representations of the
        // values.
        coeff_ring.mulZ(&res, &res, &rmodq);

        // Compare the resulting values.
        assert(ref == res);

    }

    printf("Test finished!\n");

}





