
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

// ================
// This file demonstrate the correctness of Nussbaumer + Cooley--Tukey
// for Z_Q[x] / (x^256 + 1) for Q = 1, 2, 4, ..., 2^27.

// ================
// Optimization guide.
/*

1. Notice that all the twiddle factors are negacyclic shifts.
   Currently, they are polynomial multiplications.
   Change them into negacyclic shifts.

2. After applying Nussbaumer and Cooley--Tukey, the remaining computing tasks are
   32 polynomial multiplications in Z_{32 Q}[y] / (y^16 + 1).
   Design fast computations for them.

*/

#define ARRAY_N 256
#define INNER_N 16

// Q = 1, 2, 4, ..., 2^27
#define Q (1 << 27)
#define COEFF_SIZE (sizeof(int32_t))

int32_t mod = Q;

// ================
// Z_{2^32}

void memberZ(void *des, void *src){
    *(int32_t*)des = *(int32_t*)src;
}

void addZ(void *des, void *src1, void *src2){
    *(int32_t*)des = (*(int32_t*)src1) + (*(int32_t*)src2);
}

void subZ(void *des, void *src1, void *src2){
    *(int32_t*)des = (*(int32_t*)src1) - (*(int32_t*)src2);
}

void mulZ(void *des, void *src1, void *src2){
    *(int32_t*)des = (*(int32_t*)src1) * (*(int32_t*)src2);
}

void expZ(void *des, void *src, size_t e){

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

struct commutative_ring coeff_ring = {
    .sizeZ = sizeof(int32_t),
    .memberZ = memberZ,
    .addZ = addZ,
    .subZ = subZ,
    .mulZ = mulZ,
    .expZ = expZ
};

// ================
// Z_{2^32}[y] / (y^16 + 1)

void memberZ_negacyclic(void *des, void *src){
    for(size_t i = 0; i < INNER_N; i++){
        coeff_ring.memberZ(des + i * COEFF_SIZE, src + i * COEFF_SIZE);
    }
}

void addZ_negacyclic(void *des, void *src1, void *src2){
    for(size_t i = 0; i < INNER_N; i++){
        coeff_ring.addZ(des + i * COEFF_SIZE, src1 + i * COEFF_SIZE, src2 + i * COEFF_SIZE);
    }
}

void subZ_negacyclic(void *des, void *src1, void *src2){
    for(size_t i = 0; i < INNER_N; i++){
        coeff_ring.subZ(des + i * COEFF_SIZE, src1 + i * COEFF_SIZE, src2 + i * COEFF_SIZE);
    }
}

void mulZ_negacyclic(void *des, void *src1, void *src2){
    int32_t twiddle = -1;
    naive_mulR(des, src1, src2, INNER_N, &twiddle, coeff_ring);
}

void expZ_negacyclic(void *des, void *src, size_t e){

    int32_t src_v[INNER_N];
    int32_t tmp_v[INNER_N];
    int32_t twiddle = -1;

    memmove(src_v, src, INNER_N * coeff_ring.sizeZ);

    memset(tmp_v, 0, INNER_N * coeff_ring.sizeZ);
    tmp_v[0] = 1;

    for(; e; e >>= 1){
        if(e & 1){
            naive_mulR(tmp_v, tmp_v, src_v, INNER_N, &twiddle, coeff_ring);
        }
        naive_mulR(src_v, src_v, src_v, INNER_N, &twiddle, coeff_ring);
    }

    memmove(des, tmp_v, INNER_N * coeff_ring.sizeZ);
}

struct commutative_ring negacyclic_ring = {
    .sizeZ = COEFF_SIZE * INNER_N,
    .memberZ = memberZ_negacyclic,
    .addZ = addZ_negacyclic,
    .subZ = subZ_negacyclic,
    .mulZ = mulZ_negacyclic,
    .expZ = expZ_negacyclic
};

// ================
// buffers for twiddle factors in (Z_{2^32}[y] / (y^16 + 1))[x] / (x^32 - 1)

int32_t streamlined_twiddle_symbolic_CT_table[32][INNER_N];
int32_t streamlined_twiddle_symbolic_CT_itable[32][INNER_N];

int main(void){

    int32_t poly1[ARRAY_N], poly2[ARRAY_N];
    int32_t poly1_NTT[32 * INNER_N], poly2_NTT[32 * INNER_N];
    int32_t res_NTT[32 * INNER_N];
    int32_t ref[ARRAY_N], res[ARRAY_N];

    int32_t twiddle, t;

    int32_t twiddle_negacyclic[INNER_N];
    int32_t scale_negacyclic[INNER_N];
    int32_t zeta_negacyclic[INNER_N];

// ================

    for(size_t i = 0; i < ARRAY_N; i++){
        t = rand();
        coeff_ring.memberZ(poly1 + i, &t);
        t = rand();
        coeff_ring.memberZ(poly2 + i, &t);
    }

// ================

    // Compute poly1 * poly2 in Z_{2^32}[x] / (x^256 + 1).
    twiddle = -1;
    naive_mulR(ref,
        poly1, poly2, ARRAY_N, &twiddle, coeff_ring);
    // Reduce the coefficient ring from Z_{2^32} to Z_Q.
    for(size_t i = 0; i < ARRAY_N; i++){
        cmod_int32(ref + i, ref + i, &mod);
    }

// Nussbaumer for Z_Q[x]/ (x^256 + 1)
// ================
// Notice that starting from this step, we compute the 2^5-multiple of the
// product. Therefore, Q must not larger than 2^32 / 2^5 = 2^27.
// This is why Q = 1, 2, 4, ..., 2^27 are the only options.

    // Initialization
    for(size_t i = 0; i < 32; i++){
        for(size_t j = 0; j < 16; j++){
            poly1_NTT[i * 16 + j] = 0;
            poly2_NTT[i * 16 + j] = 0;
            res_NTT[i * 16 + j] = 0;
        }
    }

    // Z_{2^32}[x] / (x^256 + 1)
    // to
    // Z_{2^32}[x, y] / (x^16 - y, y^16 + 1)
    // to
    // ( Z_{2^32}[y] / (y^16 + 1) ) [x] / (x^16 - y)
    // to
    // ( Z_{2^32}[y] / (y^16 + 1) ) [x] / (x^32 - 1)
    for(size_t i = 0; i < 16; i++){
        for(size_t j = 0; j < 16; j++){
            poly1_NTT[j * 16 + i] = poly1[i * 16 + j];
            poly2_NTT[j * 16 + i] = poly2[i * 16 + j];
        }
    }

// ========

    // Specify the layer-merging strategy.
    // For simplicity, we compute one layer at a time.
    struct compress_profile profile = {
    32, 32, 5, 5, {1, 1, 1, 1, 1}
    };

    // Initialize constants for generating twiddle factors.
    memset(twiddle_negacyclic, 0, negacyclic_ring.sizeZ);
    twiddle_negacyclic[1] = 1;
    memset(scale_negacyclic, 0, negacyclic_ring.sizeZ);
    scale_negacyclic[0] = 1;
    memset(zeta_negacyclic, 0, negacyclic_ring.sizeZ);
    zeta_negacyclic[0] = 1;

    // Generate the twiddle factors for computing FFTs in x.
    gen_streamlined_DWT_table(streamlined_twiddle_symbolic_CT_table,
        scale_negacyclic, twiddle_negacyclic, zeta_negacyclic, profile, 0, negacyclic_ring);

    // Apply symbolic FFTs.
    // Now we have prod_i ( Z_{2^32}[y] / (y^16 + 1) ) [x] / (x - y^i)
    CT_NTT(poly1_NTT, streamlined_twiddle_symbolic_CT_table, profile, negacyclic_ring);
    CT_NTT(poly2_NTT, streamlined_twiddle_symbolic_CT_table, profile, negacyclic_ring);

    // Compute the products in prod_i ( Z_{2^32}[y] / (y^16 + 1) ) [x] / (x - y^i)
    twiddle = -1;
    for(size_t i = 0; i < 32; i++){
        naive_mulR(res_NTT + i * INNER_N, poly1_NTT + i * INNER_N, poly2_NTT + i * INNER_N, 16, &twiddle, coeff_ring);
    }

    // Initialize constants for generating twiddle factors.
    memset(twiddle_negacyclic, 0, negacyclic_ring.sizeZ);
    twiddle_negacyclic[15] = -1;
    memset(scale_negacyclic, 0, negacyclic_ring.sizeZ);
    scale_negacyclic[0] = 1;

    // Generate twiddle factors.
    gen_streamlined_inv_CT_table(streamlined_twiddle_symbolic_CT_itable,
        scale_negacyclic, twiddle_negacyclic, profile, 0, negacyclic_ring);

    // Apply the inverse of symbolic FFT.
    CT_iNTT(res_NTT, streamlined_twiddle_symbolic_CT_itable, profile, negacyclic_ring);

    // At this point, we compute the 2^5-multiple of the desired result.
    // Notice that this step is commutative with the follow-up steps.
    for(size_t i = 0; i < 32 * INNER_N; i++){
        res_NTT[i] >>= 5;
    }

// ========

    // ( Z_{2^32}[y] / (y^16 + 1) ) [x] / (x^32 - 1)
    // to
    // ( Z_{2^32}[y] / (y^16 + 1) ) [x] / (x^16 - y)
    for(size_t j = 1; j < 16; j++){
        for(size_t i = 0; i < 16; i++){
            coeff_ring.addZ(res_NTT + i * 16 + j, res_NTT + i * 16 + j, res_NTT + (i + 16) * 16 + (j - 1));
        }
    }

    for(size_t i = 0; i < 16; i++){
        coeff_ring.subZ(res_NTT + i * 16, res_NTT + i * 16, res_NTT + (i + 16) * 16 + 15);
    }

    // ( Z_{2^32}[y] / (y^16 + 1) ) [x] / (x^16 - y)
    // to
    // Z_{2^32}[x] / (x^256 + 1)
    for(size_t i = 0; i < 16; i++){
        for(size_t j = 0; j < 16; j++){
            res[i * 16 + j] = res_NTT[j * 16 + i];
        }
    }

// ========

    // Z_{2^32}[x] / (x^256 + 1)
    // to
    // Z_Q[x] / (x^256 + 1)
    for(size_t i = 0; i < ARRAY_N; i++){
        cmod_int32(res + i, res + i, &mod);
    }

// ================

    for(size_t i = 0; i < ARRAY_N; i++){
        assert(ref[i] == res[i]);
    }

    printf("Test finished!\n");

}









