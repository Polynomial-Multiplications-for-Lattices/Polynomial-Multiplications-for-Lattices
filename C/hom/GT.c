
#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>

#include "tools.h"
#include "naive_mult.h"
#include "gen_table.h"
#include "ntt_c.h"

/*

TODO:

Apply Cooley--Tukey FFT to the radix-2 part.

*/

// ================
// This file demonstrates the isomorphism Z_Q[x] / (x^1536 - 1) \cong
// Z_Q[z] / (z^3 - 1) \otimes Z_Q[y] / (y^512 - 1).

// ================
// Optimization guide.

/*

TBA

*/

#define ARRAY_N 1536

#define Q (7681)

// ================
// Z_Q

int16_t mod = Q;

void memberZ(void *des, void *src){
    cmod_int16(des, src, &mod);
}

void addZ(void *des, void *src1, void *src2){
    addmod_int16(des, src1, src2, &mod);
}

void subZ(void *des, void *src1, void *src2){
    submod_int16(des, src1, src2, &mod);
}

void mulZ(void *des, void *src1, void *src2){
    mulmod_int16(des, src1, src2, &mod);
}

void expZ(void *des, void *src, size_t e){
    expmod_int16(des, src, e, &mod);
}

struct commutative_ring coeff_ring = {
    .sizeZ = sizeof(int16_t),
    .memberZ = memberZ,
    .addZ = addZ,
    .subZ = subZ,
    .mulZ = mulZ,
    .expZ = expZ
};

// ================
// Z_Q[z] / (z^3 - 1)

void memberZ_convol(void *des, void *src){

    size_t size = sizeof(int16_t);

    for(size_t i = 0; i < 3; i++){
        cmod_int16(des + i * size, src + i * size, &mod);
    }

}

void addZ_convol(void *des, void *src1, void *src2){

    size_t size = sizeof(int16_t);

    for(size_t i = 0; i < 3; i++){
        addmod_int16(des + i * size, src1 + i * size, src2 + i * size, &mod);
    }

}

void subZ_convol(void *des, void *src1, void *src2){

    size_t size = sizeof(int16_t);

    for(size_t i = 0; i < 3; i++){
        submod_int16(des + i * size, src1 + i * size, src2 + i * size, &mod);
    }

}

void mulZ_convol(void *des, void *src1, void *src2){

    size_t twiddle = 1;

    naive_mulR(des,
        src1, src2, 3, &twiddle, coeff_ring);

}

void expZ_convol(void *des, void *src, size_t e){
    return;
}

struct commutative_ring convol_ring = {
    .sizeZ = sizeof(int16_t) * 3,
    .memberZ = memberZ_convol,
    .addZ = addZ_convol,
    .subZ = subZ_convol,
    .mulZ = mulZ_convol,
    .expZ = expZ_convol
};

int main(void){

    int16_t poly1[ARRAY_N], poly2[ARRAY_N];
    int16_t ref[ARRAY_N], res[ARRAY_N];

    int16_t poly1_NTT[ARRAY_N], poly2_NTT[ARRAY_N];
    int16_t res_NTT[ARRAY_N];

    int16_t twiddle_convol[3];

    int16_t twiddle, t;

    for(size_t i = 0; i < ARRAY_N; i++){
        t = rand();
        coeff_ring.memberZ(poly1 + i, &t);
        t = rand();
        coeff_ring.memberZ(poly2 + i, &t);
    }

// ================

    // Compute the product in Z_Q[x] / (x^1536 - 1).
    twiddle = 1;
    naive_mulR(ref,
        poly1, poly2, ARRAY_N, &twiddle, coeff_ring);

// ================

    // Permute so we have Z_Q[x] / (x^1536 - 1) \cong
    // Z_Q[y] / (y^512 - 1) \otimes Z_Q[z] / (z^3 - 1).
    for(size_t i = 0; i < ARRAY_N; i++){
        poly1_NTT[(i % 512) * 3 + (i % 3)] = poly1[i];
        poly2_NTT[(i % 512) * 3 + (i % 3)] = poly2[i];
    }

// ================

    // Compute the product in Z_Q[y] / (y^512 - 1) \otimes Z_Q[z] / (z^3 - 1).
    twiddle_convol[0] = 1;
    twiddle_convol[1] = 0;
    twiddle_convol[2] = 0;
    naive_mulR(res_NTT,
        poly1_NTT, poly2_NTT, 512, &twiddle_convol, convol_ring);

// ================

    // Permute so we have
    // Z_Q[y] / (y^512 - 1) \otimes Z_Q[z] / (z^3 - 1)
    // \cong Z_Q[x] / (x^1536 - 1).
    for(size_t i = 0; i < ARRAY_N; i++){
        res[i] = res_NTT[(i % 512) * 3 + (i % 3)];
    }

    for(size_t i = 0; i < ARRAY_N; i++){
        assert(ref[i] == res[i]);
    }

    printf("Test finished!\n");

}











