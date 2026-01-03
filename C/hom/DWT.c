
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
// This file computes the discrete weighted transformation (DWT) and its inversion for Z_Q[x] / (x^512 + 1)
// via Cooley--Tukey and Gentlemand--Sande FFT.

// ================
// Theory.

// ================
// A small example.

// ================
// Optimization guide.
/*

 See the file DWT_merged_layers.c

*/

// ================
// Applications to lattice-based cryptosystems.
// Generally speaking, DWT is definable for polynomial rings of the form R[x] / (x^n - zeta^n) as long as the following
// hold:
// 1. The positive integer n, encoded as the sum of n copies of the identity of R, is invertible.
// 2. There is a principal n-th root of unity in R. A root of unity w is called a principal n-th root of unity if
//    Phi_n(w) = 0 in R where Phi_n(x) is the n-th cyclotomic polynomial.
// For the DWT to be invertible when R is commutative, we also require zeta to be invertible.
// If R is non-commutative, we ask zeta to commute with all the elements in R (so zeta belongs
// to the center of R by definition).
// When R takes the form Z_Q, the definability of an invertible DWT reduces to
// 1. n | gcd(q_1 - 1, ..., q_d - 1) where Q = prod_i q_i (see [Pol71]).
// 2. zeta must be invertible in R.
// We summarize below real-world examples for the power-of-two-size DWTs.
// 1. Kyber
//    - Polynomial ring: Z_3329[x] / (x^256 + 1).
//    - DWT: Z_3329[x] / (x^256 + 1) with size 128. This transformation is written into the specification of Kyber.
// 2. Dilithium:
//    - Polynomial ring: Z_8380417[x] / (x^256 + 1).
//    - DWT: Z_8380417[x] / (x^256 + 1) with size 256. This transformation is written into the specification of Dilithium.
// 3. Saber:
//    - Polynomial ring: Z_8192[x] / (x^256 + 1) where one of the input polynomials has small coefficients.
//    - DWT:
//      (a) Z_25166081[x] / (x^256 + 1) ([CHK+21]).
//      (b) Z_20972417[x] / (x^256 + 1) ([CHK+21]).
//      (c) Z_{3329 x 7681}[x] / (x^256 + 1) ([ACC+22]).
//      (d) (Z_3329 x Z_7681)[x] / (x^256 + 1) ([ACC+22]).

// ================
// Below are the parameters for this file.
// We demonstrate how to compute products of two polynomials in Z_12289[x] / (x^512 + 1) with size-512 DWT.
// The DWT is implemented with Cooley--Tukey FFT and its inverse is implemented with Gentleman--Sande FFT.

#define ARRAY_N 512
#define NTT_N 512
#define LOGNTT_N 9

#define Q (12289)

// OMEGA is a principal (2 NTT_N)-th root of unity in Z_Q (OMEGA^NTT_N = -1 in Z_Q since NTT_N is a power of two).
#define OMEGA (49)
#define OMEGA_INV (1254)

// ================
// Z_Q

int16_t mod = Q;

void memberZ(void *des, const void *src){
    cmod_int16(des, src, &mod);
}

void addZ(void *des, const void *src1, const void *src2){
    addmod_int16(des, src1, src2, &mod);
}

void subZ(void *des, const void *src1, const void *src2){
    submod_int16(des, src1, src2, &mod);
}

void mulZ(void *des, const void *src1, const void *src2){
    mulmod_int16(des, src1, src2, &mod);
}

void expZ(void *des, const void *src, size_t e){
    expmod_int16(des, src, e, &mod);
}

struct ring coeff_ring = {
    .sizeZ = sizeof(int16_t),
    .memberZ = memberZ,
    .addZ = addZ,
    .subZ = subZ,
    .mulZ = mulZ,
    .expZ = expZ
};

// ================

int16_t streamlined_NTT_table[NTT_N - 1] = {
-1479, -5146, 4043, -1305, 722, 5736, -4134, 3542, -3504, -2545, 3621, -1646, 1212, 3195, 5860, -4821, 2639, -2625, -949, -563, -2975, -3006, -2744, 5728, -4591, 5023, 5828, -3328, -5777, -4978, 1351, 2319, -1170, -955, -790, -3201, 3014, 5086, -1326, 4846, -2747, -3135, 3712, 4805, -3553, -1062, -2294, 3091, -81, -4320, -1000, -2963, -4896, -3051, 2366, -1177, -4255, -1635, -2768, -140, -1853, -4611, -726, 1260, 4388, 4632, -5755, 2426, 334, 1428, 1696, 2013, -3289, 729, 3241, 2881, 3284, -5092, -2089, -3694, -5179, -1759, -3707, 3382, -355, -2548, -4231, 3637, 3459, 145, -5542, -2731, -3932, -4890, -5911, -2842, 480, 1022, 9, -2468, 339, 5791, 544, -1673, 4278, -5331, -4989, -4177, -3584, 1381, -2525, -953, -3748, 827, 5767, 2476, 118, 2197, -5067, 3949, -3296, 4452, 2396, -4354, 130, 2837, -5374, 2401, 442, -5101, -1067, 390, 773, -3833, 3778, 354, 4861, -2912, 5698, 5012, -2481, 2859, -1045, 1017, -4885, 1632, -5084, 27, -3066, -3763, -1440, 1537, 242, 4714, -4143, -2678, 3704, 5019, -545, 1002, 5011, 5088, -4284, -4976, -1607, -3780, -875, -2437, 3646, 6022, 2987, -2566, -2187, -6039, -2422, -1065, 2143, -404, -4645, 1168, 5277, -1207, 3248, 493, -4096, -5444, 2381, -4337, -435, 1378, 1912, 2166, 3915, -113, -4919, -160, 3149, -3, 4437, 3636, 4938, 5291, 2704, -1426, -4654, 1663, -1777, 3364, 1689, 4057, -3271, -2847, -4414, 2174, 4372, -5042, -2305, 4053, 2645, 5195, -2780, -4895, 1484, -3247, -2686, -3978, -2969, -2370, 2865, 5332, 3510, 1630, -2126, 5407, 3186, -1153, -2884, -2249, -4048, -2399, -3400, -5191, -3136, -3000, 671, 3016, 243, -5559, 420, -2178, 1544, 3985, 4905, 3531, 476, 49, 1263, 5915, 1483, -2500, -1489, -1583, -5942, 1512, 350, -1815, 5383, 5369, -2057, -3202, 4493, -2738, -5868, -5735, 2655, -3009, 1693, 174, 723, -1975, -3757, 347, 2925, -3315, -426, 1858, 4754, 3030, 4115, 2361, -1843, 2908, 218, 3434, -3529, 3963, 576, 6142, -2447, 1954, -2051, -2882, -1805, 3991, -3969, -2767, 156, 2281, 5876, -2031, 5333, 3772, 418, 5908, -453, 5429, -4774, -4737, 1293, 295, 6099, 5766, 652, -4016, 4077, -3762, -2919, 325, -1404, -1146, -948, 5990, 1159, -3728, -4049, 3329, 4298, -168, 2692, 5961, -5106, -1962, 1594, -6122, -2555, -5184, -1200, 1360, 3956, -6119, 5297, -4079, -1058, 922, 441, 1958, 4322, 1112, 2078, 4046, 709, -3150, 1319, 4240, -3570, -6065, -835, 2459, 683, 3656, -64, -1566, 5782, -2948, -2503, -3123, -1747, -3054, -5486, -4433, -5919, 3834, -5257, -5241, -2920, -4169, -3127, -5468, 1010, -3482, 787, 5057, 4698, 4780, -3445, -192, 1321, 4912, -2049, 677, -5874, -6055, -3336, 1323, -2766, -52, 3174, 1579, -431, -2505, 5906, 3957, -2839, 151, -2127, -58, -241, 3532, -1003, 1956, -5009, -885, -6008, 3477, -5681, 142, -1105, -2844, 3438, -975, 4212, -3029, -5594, 4782, 5886, -4213, 504, 2302, -605, -421, -4080, 3602, 6068, -3600, 3263, 6077, -4624, -4467, -4789, -5537, 4749, 4449, -5456, -147, -3789, 6118, -3818, 1190, -2683, 3860, 5445, -4536, -1050, 5079, -3262, 2169, -522, -4324, 4916, -4075, 5315, -1278, -2344, 1973, -5574, -3514, -1041, 5925, -1018, 654, 3565, 1702, 1987, -5529, 5206, 3199, -56, 6136, -5862, -5415, -3643, 4948, -6137, 400, -1728, 5339, 5446, 3710, 6093, 468, -3988, 316, -382, -2033, -3998, 3879, 1922, -1359, -5435, 973, -1254
};

int16_t streamlined_iNTT_table[NTT_N - 1] = {
1479, -4043, 5146, 4134, -5736, -722, 1305, -5860, -3195, -1212, 1646, -3621, 2545, 3504, -3542, -1351, 4978, 5777, 3328, -5828, -5023, 4591, -5728, 2744, 3006, 2975, 563, 949, 2625, -2639, 4821, 726, 4611, 1853, 140, 2768, 1635, 4255, 1177, -2366, 3051, 4896, 2963, 1000, 4320, 81, -3091, 2294, 1062, 3553, -4805, -3712, 3135, 2747, -4846, 1326, -5086, -3014, 3201, 790, 955, 1170, -2319, 5374, -2837, -130, 4354, -2396, -4452, 3296, -3949, 5067, -2197, -118, -2476, -5767, -827, 3748, 953, 2525, -1381, 3584, 4177, 4989, 5331, -4278, 1673, -544, -5791, -339, 2468, -9, -1022, -480, 2842, 5911, 4890, 3932, 2731, 5542, -145, -3459, -3637, 4231, 2548, 355, -3382, 3707, 1759, 5179, 3694, 2089, 5092, -3284, -2881, -3241, -729, 3289, -2013, -1696, -1428, -334, -2426, 5755, -4632, -4388, -1260, -476, -3531, -4905, -3985, -1544, 2178, -420, 5559, -243, -3016, -671, 3000, 3136, 5191, 3400, 2399, 4048, 2249, 2884, 1153, -3186, -5407, 2126, -1630, -3510, -5332, -2865, 2370, 2969, 3978, 2686, 3247, -1484, 4895, 2780, -5195, -2645, -4053, 2305, 5042, -4372, -2174, 4414, 2847, 3271, -4057, -1689, -3364, 1777, -1663, 4654, 1426, -2704, -5291, -4938, -3636, -4437, 3, -3149, 160, 4919, 113, -3915, -2166, -1912, -1378, 435, 4337, -2381, 5444, 4096, -493, -3248, 1207, -5277, -1168, 4645, 404, -2143, 1065, 2422, 6039, 2187, 2566, -2987, -6022, -3646, 2437, 875, 3780, 1607, 4976, 4284, -5088, -5011, -1002, 545, -5019, -3704, 2678, 4143, -4714, -242, -1537, 1440, 3763, 3066, -27, 5084, -1632, 4885, -1017, 1045, -2859, 2481, -5012, -5698, 2912, -4861, -354, -3778, 3833, -773, -390, 1067, 5101, -442, -2401, 1254, -973, 5435, 1359, -1922, -3879, 3998, 2033, 382, -316, 3988, -468, -6093, -3710, -5446, -5339, 1728, -400, 6137, -4948, 3643, 5415, 5862, -6136, 56, -3199, -5206, 5529, -1987, -1702, -3565, -654, 1018, -5925, 1041, 3514, 5574, -1973, 2344, 1278, -5315, 4075, -4916, 4324, 522, -2169, 3262, -5079, 1050, 4536, -5445, -3860, 2683, -1190, 3818, -6118, 3789, 147, 5456, -4449, -4749, 5537, 4789, 4467, 4624, -6077, -3263, 3600, -6068, -3602, 4080, 421, 605, -2302, -504, 4213, -5886, -4782, 5594, 3029, -4212, 975, -3438, 2844, 1105, -142, 5681, -3477, 6008, 885, 5009, -1956, 1003, -3532, 241, 58, 2127, -151, 2839, -3957, -5906, 2505, 431, -1579, -3174, 52, 2766, -1323, 3336, 6055, 5874, -677, 2049, -4912, -1321, 192, 3445, -4780, -4698, -5057, -787, 3482, -1010, 5468, 3127, 4169, 2920, 5241, 5257, -3834, 5919, 4433, 5486, 3054, 1747, 3123, 2503, 2948, -5782, 1566, 64, -3656, -683, -2459, 835, 6065, 3570, -4240, -1319, 3150, -709, -4046, -2078, -1112, -4322, -1958, -441, -922, 1058, 4079, -5297, 6119, -3956, -1360, 1200, 5184, 2555, 6122, -1594, 1962, 5106, -5961, -2692, 168, -4298, -3329, 4049, 3728, -1159, -5990, 948, 1146, 1404, -325, 2919, 3762, -4077, 4016, -652, -5766, -6099, -295, -1293, 4737, 4774, -5429, 453, -5908, -418, -3772, -5333, 2031, -5876, -2281, -156, 2767, 3969, -3991, 1805, 2882, 2051, -1954, 2447, -6142, -576, -3963, 3529, -3434, -218, -2908, 1843, -2361, -4115, -3030, -4754, -1858, 426, 3315, -2925, -347, 3757, 1975, -723, -174, -1693, 3009, -2655, 5735, 5868, 2738, -4493, 3202, 2057, -5369, -5383, 1815, -350, -1512, 5942, 1583, 1489, 2500, -1483, -5915, -1263, -49
};

int16_t streamlined_twiddle_table[(NTT_N - 1)];

int main(void){

    int16_t poly1[ARRAY_N], poly2[ARRAY_N];
    int16_t ref[ARRAY_N], res[ARRAY_N];

    int16_t omega, zeta, twiddle, scale, t;

    for(size_t i = 0; i < ARRAY_N; i++){
        t = rand();
        coeff_ring.memberZ(poly1 + i, &t);
        t = rand();
        coeff_ring.memberZ(poly2 + i, &t);
    }

// ================
// Compute the product in Z_Q[x] / (x^512 + 1).

    twiddle = -1;
    naive_mulR(ref,
        poly1, poly2, ARRAY_N, &twiddle, coeff_ring);

// ================
// Specify the layer-merging strategy.

    struct compress_profile profile = {
        ARRAY_N, NTT_N, LOGNTT_N, LOGNTT_N
    };

    for(size_t i = 0; i < profile.compressed_layers; i++){
        profile.merged_layers[i] = 1;
    }

// ================
// Generate twiddle factors for Cooley--Tukey FFT.

    zeta = OMEGA;
    coeff_ring.expZ(&omega, &zeta, 2);
    scale = 1;
    gen_streamlined_DWT_table(streamlined_twiddle_table,
        &scale, &omega, &zeta, profile, 0, coeff_ring);

    //

    assert(memcmp(streamlined_NTT_table, streamlined_twiddle_table, (NTT_N - 1) * sizeof(int16_t)) == 0);

// ================
// Apply Cooley--Tukey FFT.

    CT_NTT(poly1, streamlined_NTT_table, profile, coeff_ring);
    CT_NTT(poly2, streamlined_NTT_table, profile, coeff_ring);

// ================

    point_mul(res, poly1, poly2, ARRAY_N, 1, coeff_ring);

// ================
// Generate twiddle factors for the inverse via Gentlemans--Sande FFT.

    zeta = OMEGA_INV;
    coeff_ring.expZ(&omega, &zeta, 2);
    scale = 1;
    gen_streamlined_DWT_table(streamlined_twiddle_table,
        &scale, &omega, &zeta, profile, 0, coeff_ring);

    assert(memcmp(streamlined_iNTT_table, streamlined_twiddle_table, (NTT_N - 1) * sizeof(int16_t)) == 0);

// ================
// Apply Gentleman--Sande FFT.

    GS_iNTT(res, streamlined_iNTT_table, profile, coeff_ring);

// ================
// Multiply the scale to reference.

    scale = NTT_N;
    for(size_t i = 0; i < ARRAY_N; i++){
        coeff_ring.mulZ(ref + i, ref + i, &scale);
    }

    for(size_t i = 0; i < ARRAY_N; i++){
        assert(ref[i] == res[i]);
    }

    printf("Test finished!\n");

}








