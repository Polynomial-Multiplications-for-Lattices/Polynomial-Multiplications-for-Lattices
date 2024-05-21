
#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>

#include "tools.h"

// ================
// This file demonstrates signed Barrett multiplication.
// Let a and b be the operands that we wish to multiply, Q be the modulus, and R > Q be
// the arithmetic precision.
// Barrett multiplication computes a value that is reasonably close to a b mod^+- Q
// by approximating the quotient a b / Q and subtract approx(a b / Q) Q from a b.
// As long as approx(a b / Q) is close to a b / Q, |a b - approx(a b / Q) Q| < R / 2.


// ================
// Theory.
// Observe that a b mod^+- Q = a b - round(a b / Q) Q, if we replace round(a b / Q) with
// a function admit efficient computation with a reasonably tolerable error delta,
// then the result is off by delta Q. Therefore, we just need to ensure
// (delta + 1 / 2) Q < R / 2.
// See Barrett_Montgomery_cmp for a formal proof.

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

struct commutative_ring coeff_ring = {
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

    lo = a * b;
    bhi = get_barrett_hi(b, rmodq, qprime);
    hi = mulhir(a, bhi);

    return lo - hi * q;

}

int main(void){

    int32_t a, b, t, q, qprime, rmodq;
    int32_t ref, res;

    q = Q;
    qprime = Qprime;
    rmodq = RmodQ;

    for(size_t i = 0; i < NTESTS; i++){

        t = rand() % Q;
        coeff_ring.memberZ(&a, &t);
        t = rand() % Q;
        coeff_ring.memberZ(&b, &t);

        coeff_ring.mulZ(&ref, &a, &b);

        res = barrett_mul(a, b, q, rmodq, qprime);

        coeff_ring.memberZ(&res, &res);

        assert(ref == res);

    }

    printf("Test finished!\n");

}
