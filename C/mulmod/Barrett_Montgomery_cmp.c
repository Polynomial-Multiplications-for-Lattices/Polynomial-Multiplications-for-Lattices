
#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>

#include "tools.h"

// R = 2^32 below
#define Q 8380417
// RmodQ = R mod^+- Q
#define RmodQ (-4186625)
// Qprime = -Q^{-1} mod^+- R
#define Qprime (-58728449)

#define NTESTS 1000

// ================
// Z_Q

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

    lo = a * b;
    bhi = get_barrett_hi(b, rmodq, qprime);
    hi = mulhir(a, bhi);

    return lo - hi * q;

}

// The accumulative variant of Montgomery multiplication.
int32_t montgomery_acc_mul(int32_t a, int32_t b, int32_t q, int32_t qprime){

    int64_t prod;
    int32_t lo;

    prod = mullong(a, b);
    lo = mullo(prod, qprime);
    prod += mullong(lo, q);

    return gethi(prod);

}

int main(void){

    int32_t a, b, t, q, qprime, rmodq;
    int32_t res_montgomery, res_barrett;

    q = Q;
    qprime = Qprime;
    rmodq = RmodQ;

    for(size_t i = 0; i < NTESTS; i++){

        t = rand() % Q;
        coeff_ring.memberZ(&a, &t);
        t = rand() % Q;
        coeff_ring.memberZ(&b, &t);

        coeff_ring.mulZ(&t, &b, &rmodq);

        res_barrett = barrett_mul(a, b, q, rmodq, qprime);
        res_montgomery = montgomery_acc_mul(a, t, q, qprime);

        // Compare if the results of Barrett and Montgomery multiplications are the same.
        assert(res_montgomery == res_barrett);

    }

    printf("Test finished!\n");

}


