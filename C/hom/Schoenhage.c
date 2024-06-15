
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
// This file demonstrate the correctness of Schoenhage + Cooley--Tukey
// for Z_Q[x] / (x^256 - 1) with Q = 1, 2, 4, ..., 2^27.

// ================
// Optimization guide.
/*

1. Notice that all the twiddle factors are negacyclic shifts.
   Currently, they are polynomial multiplications.
   Change them into negacyclic shifts.

2. After applying Schoenhage and Cooley--Tukey, the remaining computing tasks are
   32 polynomial multiplications in Z_{32 Q}[x] / (x^16 + 1).
   Design fast computations for them.

*/

// ================
// Applications to lattice-based cryptosystems.

#define ARRAY_N 256
#define INNER_N 16
#define TWIDDLE_POS 1

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

struct ring coeff_ring = {
    .sizeZ = sizeof(int32_t),
    .memberZ = memberZ,
    .addZ = addZ,
    .subZ = subZ,
    .mulZ = mulZ,
    .expZ = expZ
};

// ================
// Z_{2^32}[x] / (x^INNER_N + 1)

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

struct ring negacyclic_ring = {
    .sizeZ = COEFF_SIZE * INNER_N,
    .memberZ = memberZ_negacyclic,
    .addZ = addZ_negacyclic,
    .subZ = subZ_negacyclic,
    .mulZ = mulZ_negacyclic,
    .expZ = expZ_negacyclic
};

// ================
// buffers for twiddle factors in (Z_{2^32}[x] / (x^INNER_N + 1))[y] / (y^32 - 1)

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

    // Compute poly1 * poly2 in Z_{2^32}[x] / (x^ARRAY_N - 1).
    twiddle = 1;
    naive_mulR(ref,
        poly1, poly2, ARRAY_N, &twiddle, coeff_ring);
    // Reduce the coefficient ring from Z_{2^32} to Z_Q.
    for(size_t i = 0; i < ARRAY_N; i++){
        cmod_int32(ref + i, ref + i, &mod);
    }

// Schoenhage for Z_Q[x]/ (x^ARRAY_N - 1)
// ================
// Notice that starting from this step, we compute the 2^5-multiple of the
// product. Therefore, Q must not larger than 2^32 / 2^5 = 2^27.
// This is why Q = 1, 2, 4, ..., 2^27 are the only options.

    // Initialization
    memset(poly1_NTT, 0, sizeof(poly1_NTT));
    memset(poly2_NTT, 0, sizeof(poly2_NTT));
    memset(res_NTT, 0, sizeof(res_NTT));

    // Z_Q[x] / (x^ARRAY_N - 1)
    // to
    // Z_Q[x, y] / (x^(INNER_N / 2) - y, y^32 - 1)
    // to
    // ( Z_Q[x] / (x^(INNER_N / 2) - y) ) [y] / (y^32 - 1)
    // to
    // ( Z_Q[x] / (x^INNER_N + 1) ) [y] / (y^32 - 1)
    for(size_t i = 0; i < 32; i++){
        for(size_t j = 0; j < (INNER_N / 2); j++){
            poly1_NTT[i * INNER_N + j] = poly1[i * (INNER_N / 2) + j];
            poly2_NTT[i * INNER_N + j] = poly2[i * (INNER_N / 2) + j];
        }
    }

// ================

    // Specify the layer-merging strategy.
    // For simplicity, we compute one layer at a time.
    struct compress_profile profile = {
        32, 32, 5, 5, {1, 1, 1, 1, 1}
    };

    // Initialize constants for generating twiddle factors.
    memset(twiddle_negacyclic, 0, negacyclic_ring.sizeZ);
    twiddle_negacyclic[TWIDDLE_POS] = 1;
    memset(scale_negacyclic, 0, negacyclic_ring.sizeZ);
    scale_negacyclic[0] = 1;
    memset(zeta_negacyclic, 0, negacyclic_ring.sizeZ);
    zeta_negacyclic[0] = 1;

    // Generate the twiddle factors for computing FFTs in y.
    gen_streamlined_DWT_table(streamlined_twiddle_symbolic_CT_table,
        scale_negacyclic, twiddle_negacyclic, zeta_negacyclic, profile, 0, negacyclic_ring);

    // Apply symbolic FFTs.
    // Now we have prod_i ( Z_{2^32}[x] / (x^INNER_N + 1) ) [y] / (y - x^i).
    CT_NTT(poly1_NTT, streamlined_twiddle_symbolic_CT_table, profile, negacyclic_ring);
    CT_NTT(poly2_NTT, streamlined_twiddle_symbolic_CT_table, profile, negacyclic_ring);

    // Compute the products in prod_i ( Z_{2^32}[x] / (x^INNER_N + 1) ) [y] / (y - x^i)
    twiddle = -1;
    for(size_t i = 0; i < 32; i++){
        naive_mulR(res_NTT + i * INNER_N, poly1_NTT + i * INNER_N, poly2_NTT + i * INNER_N, INNER_N, &twiddle, coeff_ring);
    }

    // Initialize constants for generating twiddle factors.
    memset(twiddle_negacyclic, 0, negacyclic_ring.sizeZ);
    twiddle_negacyclic[INNER_N - TWIDDLE_POS] = -1;
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

// ================

    // ( Z_Q[x] / (x^INNER_N + 1) ) [y] / (y^32 - 1)
    // to
    // ( Z_Q[x] / (x^(INNER_N / 2) - y) ) [y] / (y^32 - 1)
    for(size_t i = 1; i < 32; i++){
        for(size_t j = 0; j < (INNER_N / 2); j++){
            coeff_ring.addZ(res_NTT + i * INNER_N + j, res_NTT + i * INNER_N + j, res_NTT + (i - 1) * INNER_N + j + (INNER_N / 2));
        }
    }

    for(size_t j = 0; j < (INNER_N / 2); j++){
        coeff_ring.addZ(res_NTT + 0 * INNER_N + j, res_NTT + 0 * INNER_N + j, res_NTT + 31 * INNER_N + j + (INNER_N / 2));
    }

    // ( Z_Q[x] / (x^(INNER_N / 2) - y) ) [y] / (y^32 - 1)
    // to
    // Z_{2^32}[x] / (x^ARRAY_N - 1)
    for(size_t i = 0; i < 32; i++){
        for(size_t j = 0; j < (INNER_N / 2); j++){
            res[i * (INNER_N / 2) + j] = res_NTT[i * INNER_N + j];
        }
    }

// ================

    // Z_{2^32}[x] / (x^ARRAY_N - 1)
    // to
    // Z_Q[x] / (x^ARRAY_N - 1)
    for(size_t i = 0; i < ARRAY_N; i++){
        cmod_int32(res + i, res + i, &mod);
    }

// ================

    for(size_t i = 0; i < ARRAY_N; i++){
        assert(ref[i] == res[i]);
    }

    printf("Test finished!\n");

}









