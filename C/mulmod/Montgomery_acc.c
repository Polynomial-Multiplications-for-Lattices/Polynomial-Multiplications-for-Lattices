
#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>

#include "tools.h"
#include "naive_mult.h"

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

int64_t mullong(int32_t a, int32_t b){
    return (int64_t)a * b;
}

int32_t mullo(int32_t a, int32_t b){
    return a * b;
}

int32_t gethi(int64_t a){
    return (int32_t)(a >> 32);
}

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

        res = montgomery_acc_mul(a, b, q, qprime);

        coeff_ring.mulZ(&res, &res, &rmodq);

        assert(ref == res);

    }

    printf("Test finished!\n");

}




