
#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>

#include "tools.h"

// ================
// This file demonstrates signed Barrett multiplication with a non-conventional integer approximation suitable
// for multi-limb arithmetic.
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
// Let b be a constant. We define the following integer approximation:
// For all r, approx(r) = a_{r, h} b_h + floor( a_{r, l} b_h / sqrt(R) ) + floor( a_{r, h} b_l / sqrt(R) )
// for a_{r, l} + a_{r, h} sqrt(R) = r R / round(b R / Q) and b_l + b_h sqrt(R) = round(b R / Q).
// For -Q / 2 <= b < Q / 2 and -R/2 <= a_{r, l} + a_{r, h} sqrt(R) < R/2,
// we have |r - approx(r)| <= 3, implying that the absolute value of result of the product
// is smaller than or equal to 7 R / 2.

// ================
// Proof.
// |r - approx(r)|
// = | (a_{r, l} + a_{r, h} sqrt(R)) (b_l + b_h sqrt(R)) / R
// - ( a_{r, h} b_h + floor( a_{r, l} b_h / sqrt(R) ) + floor( a_{r, h} b_l / sqrt(R) ) ) |
// = | a_{r, l} b_l / R + ( a_{r, h} b_l / sqrt(R) - floor( a_{r, h} b_l / sqrt(R) ) ) +
//                        ( a_{r, l} b_h / sqrt(R) - floor( a_{r, l} b_r / sqrt(R) ) ) |
// <= 3.

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

// mulhi_approx computes an approximation of the high part of the long product of a and b.
// Let's write a = alo + ahi * sqrt(R), b = blo + bhi * sqrt(R).
// mulhi_approx computes ahi * bhi + floor(ahi * blo / sqrt(R)) + floor(alo * bhi / sqrt(R)).
// One can show that the results is within +- of floor(a * b / R).
// We denote the output as approx(a * b / R).
int32_t mulhi_approx(int32_t a, int32_t b){

    int32_t alo, ahi;
    int32_t blo, bhi;
    int32_t res;

    alo = (int32_t)((uint16_t)a);
    ahi = a >> 16;
    blo = (int32_t)((uint16_t)b);
    bhi = b >> 16;

    res = ahi * bhi;
    res += (ahi * blo) >> 16;
    res += (alo * bhi) >> 16;

    return res;

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

    // lo = a * b mod^+- R
    lo = a * b;
    // bhi = round( b R / Q)
    bhi = get_barrett_hi(b, rmodq, qprime);
    // hi = approx(a * round(b R / Q) / R)
    hi = mulhi_approx(a, bhi);

    // lo = (a * b mod^+- R) - approx(a * round(b R / Q) / R) * Q
    return lo - hi * q;

}

// Barrett multiplication with precomputation.
int32_t barrett_mul_pre(int32_t a, int32_t b, int32_t bhi, int32_t q){

    int32_t lo, hi;

    // lo = a * b mod^+- R
    lo = a * b;
    // hi = apporx(a * bhi / R)
    hi = mulhi_approx(a, bhi);

    // lo = (a * b mod^+- R) - approx(a * bhi / R)
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
