
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
// This file demonstrates the uses of Fermat number transform multiplying
// polynomials in Z_{65537} [x] / (x^64 - 1).

// ================
// Optimization guide.
/*

1. Notice that twiddle factors in the initial layers are powers of two. Implement the twiddle factor multiplications
   with shifts and see if they are faster than generic modular multiplications.

*/

// ================
// Applications to lattice-based cryptosystems.

#define ARRAY_N 64
#define NTT_N 64
#define LOGNTT_N 6

#define Q (65537)

#define OMEGA (-4080)
#define OMEGA_INV (-2040)

// ================
// Z_{65537}

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

int32_t streamlined_NTT_table[NTT_N - 1];

int32_t streamlined_iNTT_table[NTT_N - 1];

int32_t streamlined_twiddle_table[(NTT_N - 1)];

struct compress_profile profile;

// ================
// Tables of squares and square roots.

int32_t __sq[1u << 17];
int32_t __sqrt[1u << 17];

int main(void){

    int32_t poly1[ARRAY_N], poly2[ARRAY_N];
    int32_t ref[ARRAY_N], res[ARRAY_N];

    int32_t omega, zeta, twiddle, scale, t;

    for(size_t i = 0; i < ARRAY_N; i++){
        t = rand();
        coeff_ring.memberZ(poly1 + i, &t);
        t = rand();
        coeff_ring.memberZ(poly2 + i, &t);
    }

// ================
// Compute the product in Z_{65537}[x] / (x^256 - 1).

    twiddle = 1;
    naive_mulR(ref,
        poly1, poly2, ARRAY_N, &twiddle, coeff_ring);

// ================
// Compute sqrt(2).

    // Initialization.
    for(size_t i = 0; i < Q; i++){
        __sq[i] = __sqrt[i] = -1;
    }

    for(size_t i = 0; i < Q; i++){
        t = i;
        mulmod_int32(__sq + i, &t, &t, &mod);
    }

    for(size_t i = 0; i < Q; i++){
        t = __sq[i];
        if(t < 0){
            t += Q;
        }
        __sqrt[t] = i;
        cmod_int32(__sqrt + t, __sqrt + t, &mod);
    }

// ================
// Specify the layer-merging strategy.

    profile.array_n = ARRAY_N;
    profile.ntt_n = NTT_N;
    profile.log_ntt_n = LOGNTT_N;

    profile.compressed_layers = LOGNTT_N;
    for(size_t i = 0; i < profile.compressed_layers; i++){
        profile.merged_layers[i] = 1;
    }

// ================
// Generate twiddle factors for FFT.

    zeta = 1;
    omega = OMEGA;
    scale = 1;
    gen_streamlined_DWT_table(streamlined_twiddle_table,
        &scale, &omega, &zeta, profile, 0, coeff_ring);

// ================
// Apply Cooley--Tukey FFT.

    compressed_CT_NTT(poly1,
        0, LOGNTT_N - 1, streamlined_twiddle_table, profile, coeff_ring);
    compressed_CT_NTT(poly2,
        0, LOGNTT_N - 1, streamlined_twiddle_table, profile, coeff_ring);

// ================

    point_mul(res, poly1, poly2, ARRAY_N, 1, coeff_ring);

// ================
// Generate twiddle factors for the inverse of Cooley--Tukey FFT.

    zeta = 1;
    omega = OMEGA_INV;
    scale = 1;
    gen_streamlined_DWT_table(streamlined_twiddle_table,
        &scale, &omega, &zeta, profile, 0, coeff_ring);

// ================
// Apply the inverse of Cooley--Tukey FFT.

    compressed_GS_iNTT(res,
        0, LOGNTT_N - 1, streamlined_twiddle_table, profile, coeff_ring);

// ================
// Multiply the scale to the reference.

    scale = NTT_N;
    for(size_t i = 0; i < ARRAY_N; i++){
        coeff_ring.mulZ(ref + i, ref + i, &scale);
    }

// ================

    for(size_t i = 0; i < ARRAY_N; i++){
        assert(ref[i] == res[i]);
    }

    printf("Test finished!\n");

}








