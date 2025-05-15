
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
// the size of the arithmetic.
// Barrett multiplication computes a value that is reasonably close to a b mod^+- Q
// by approximating the quotient a b / Q and subtract it from a b.
// As long as the approximation is close to a b / Q, the absolute value of the result is smaller than R / 2.

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

void memberZ(void *des, const void *src){
    cmod_int32(des, src, &mod);
}

void addZ(void *des, const void *src1, const void *src2){
    addmod_int32(des, src1, src2, &mod);
}

void subZ(void *des, const void *src1, const void *src2){
    submod_int32(des, src1, src2, &mod);
}

void mulZ(void *des, const void *src1, const void *src2){
    mulmod_int32(des, src1, src2, &mod);
}

void expZ(void *des, const void *src, size_t e){
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
// Let bhi = round(b R  / Q). This function computes a b - round(a bhi / R) Q with
// (a b - round(a bhi / R) Q ) mod^+- R
// = ( (a b mod^+- R) - round(a bhi / R) Q ) mod^+- R.
// As long as | a b - round(a bhi / R) Q | < R / 2,
// reducting modulo R yields the same result as an integer.
int32_t barrett_mul(int32_t a, int32_t b, int32_t q, int32_t rmodq, int32_t qprime){

    int32_t lo, hi;
    int32_t bhi;

    // bhi = round(b R / Q)
    bhi = get_barrett_hi(b, rmodq, qprime);
    // lo = a b mod^+- R
    lo = a * b;
    // hi = round(a bhi / R)
    hi = mulhir(a, bhi);

    // ( (a b mod^+- R) - round(a bhi / R ) Q ) mod^+- R
    // = a b - round(a bhi / R) Q (since |a b - round(a bhi / R) Q| < R / 2)
    return lo - hi * q;

}

// Barrett multiplication with precomputed bhi = round(b R / Q).
int32_t barrett_mul_pre(int32_t a, int32_t b, int32_t bhi, int32_t q){

    int32_t lo, hi;

    // lo = a b mod^+- R
    lo = a * b;
    // hi = round(a bhi / R)
    hi = mulhir(a, bhi);

     // ( (a b mod^+- R) - round(a bhi / R ) Q ) mod^+- R
    // = a b - round(a bhi / R) Q (since |a b - round(a bhi / R) Q| < R / 2)
    return lo - hi * q;

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

        // Compute a value equivalent to the product of a and b with Barrett multiplication.
        res = barrett_mul(a, b, q, rmodq, qprime);

        // Map the value to Z_Q.
        // Notice that this step is needed only when we want the canonical representations of the
        // values.
        coeff_ring.memberZ(&res, &res);

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

        // Assuming b is known, we precompute round(b R / Q) and compute a value equivalent to
        // the product of a and b with Barrett multiplication.
        res = barrett_mul_pre(a, b, get_barrett_hi(b, rmodq, qprime), q);

        // Map the value to Z_Q.
        // Notice that this step is needed only when we want the canonical representations of the
        // values.
        coeff_ring.memberZ(&res, &res);

        // Compare the resulting values.
        assert(ref == res);

    }

    printf("Test finished!\n");

}
