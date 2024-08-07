
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
// Optimization guide.
/*

  There are several things that frequently come to experts' mind for a fixed transformation.
  This includes the choices of modular multiplications, the precision of arithmetic, the number of memory
  operations, and the availibity of vectorization. We demonstrate how memory operations are optimized in the literature.

  Recall that a radix-2 Cooley--Tukey FFT is constructed from butterfly operations, each mapping two input coefficients
  to two output coefficients. The most straightforward way is to load two coefficients from memory,
  apply a butterfly, store the resulting two coefficients to memory, and iterate through the input coefficeints
  until all the coefficients are updated. Such memory access pattern frequently arose in the reference implementations
  where readibility is the main purpose.

  For saving memory operations, we can instead load four coefficients, compute two layers of butterflies,
  and store the resulting four coefficients to memory. In this way, we save a load-store pair for
  each coefficients. We call the idea layer-merging. One can certainly extrapolate layer-merging to higher powers of two.
  The limiting factor is the register pressure.

  Register pressure. Below we summarize the recommanded layer-merging strategies for each of the ISAs/extensions.
  - Armv7-M:
    - 14 general purpose registers (32-bit).
    - 2-layer-merging.
  - Armv7-M with FPU:
    - 14 general purpose registers (32-bit).
    - 32 single-precision floating-point registers.
    - 3-layer-merging, and 4-layer-merging if unavoidable (instead of two 2-layer-merging).
  - Armv8-A:
    - 32 SIMD registers (128-bit).
    - 4-layer-merging.
  - x86-64 AVX2:
    - 16 ymm registers (256-bit).
    - 3-layer-merging.

*/

// ================
// Applications to lattice-based cryptosystems.
//

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
-1479, -5146, 4043, -1305, 722, 5736, -4134, 3542, -4821, 2639, -3504, -2625, -949, -2545, -563, -2975, 3621, -3006, -2744, -1646, 5728, -4591, 1212, 5023, 5828, 3195, -3328, -5777, 5860, -4978, 1351, 2319, 1260, 4388, -1170, 4632, -5755, -955, 2426, 334, -790, 1428, 1696, -3201, 2013, -3289, 3014, 729, 3241, 5086, 2881, 3284, -1326, -5092, -2089, 4846, -3694, -5179, -2747, -1759, -3707, -3135, 3382, -355, 3712, -2548, -4231, 4805, 3637, 3459, -3553, 145, -5542, -1062, -2731, -3932, -2294, -4890, -5911, 3091, -2842, 480, -81, 1022, 9, -4320, -2468, 339, -1000, 5791, 544, -2963, -1673, 4278, -4896, -5331, -4989, -3051, -4177, -3584, 2366, 1381, -2525, -1177, -953, -3748, -4255, 827, 5767, -1635, 2476, 118, -2768, 2197, -5067, -140, 3949, -3296, -1853, 4452, 2396, -4611, -4354, 130, -726, 2837, -5374, 2401, 49, 1263, 442, 5915, 1483, -5101, -2500, -1489, -1067, -1583, -5942, 390, 1512, 350, 773, -1815, 5383, -3833, 5369, -2057, 3778, -3202, 4493, 354, -2738, -5868, 4861, -5735, 2655, -2912, -3009, 1693, 5698, 174, 723, 5012, -1975, -3757, -2481, 347, 2925, 2859, -3315, -426, -1045, 1858, 4754, 1017, 3030, 4115, -4885, 2361, -1843, 1632, 2908, 218, -5084, 3434, -3529, 27, 3963, 576, -3066, 6142, -2447, -3763, 1954, -2051, -1440, -2882, -1805, 1537, 3991, -3969, 242, -2767, 156, 4714, 2281, 5876, -4143, -2031, 5333, -2678, 3772, 418, 3704, 5908, -453, 5019, 5429, -4774, -545, -4737, 1293, 1002, 295, 6099, 5011, 5766, 652, 5088, -4016, 4077, -4284, -3762, -2919, -4976, 325, -1404, -1607, -1146, -948, -3780, 5990, 1159, -875, -3728, -4049, -2437, 3329, 4298, 3646, -168, 2692, 6022, 5961, -5106, 2987, -1962, 1594, -2566, -6122, -2555, -2187, -5184, -1200, -6039, 1360, 3956, -2422, -6119, 5297, -1065, -4079, -1058, 2143, 922, 441, -404, 1958, 4322, -4645, 1112, 2078, 1168, 4046, 709, 5277, -3150, 1319, -1207, 4240, -3570, 3248, -6065, -835, 493, 2459, 683, -4096, 3656, -64, -5444, -1566, 5782, 2381, -2948, -2503, -4337, -3123, -1747, -435, -3054, -5486, 1378, -4433, -5919, 1912, 3834, -5257, 2166, -5241, -2920, 3915, -4169, -3127, -113, -5468, 1010, -4919, -3482, 787, -160, 5057, 4698, 3149, 4780, -3445, -3, -192, 1321, 4437, 4912, -2049, 3636, 677, -5874, 4938, -6055, -3336, 5291, 1323, -2766, 2704, -52, 3174, -1426, 1579, -431, -4654, -2505, 5906, 1663, 3957, -2839, -1777, 151, -2127, 3364, -58, -241, 1689, 3532, -1003, 4057, 1956, -5009, -3271, -885, -6008, -2847, 3477, -5681, -4414, 142, -1105, 2174, -2844, 3438, 4372, -975, 4212, -5042, -3029, -5594, -2305, 4782, 5886, 4053, -4213, 504, 2645, 2302, -605, 5195, -421, -4080, -2780, 3602, 6068, -4895, -3600, 3263, 1484, 6077, -4624, -3247, -4467, -4789, -2686, -5537, 4749, -3978, 4449, -5456, -2969, -147, -3789, -2370, 6118, -3818, 2865, 1190, -2683, 5332, 3860, 5445, 3510, -4536, -1050, 1630, 5079, -3262, -2126, 2169, -522, 5407, -4324, 4916, 3186, -4075, 5315, -1153, -1278, -2344, -2884, 1973, -5574, -2249, -3514, -1041, -4048, 5925, -1018, -2399, 654, 3565, -3400, 1702, 1987, -5191, -5529, 5206, -3136, 3199, -56, -3000, 6136, -5862, 671, -5415, -3643, 3016, 4948, -6137, 243, 400, -1728, -5559, 5339, 5446, 420, 3710, 6093, -2178, 468, -3988, 1544, 316, -382, 3985, -2033, -3998, 4905, 3879, 1922, 3531, -1359, -5435, 476, 973, -1254
};

int16_t streamlined_iNTT_table[NTT_N - 1] = {
1479, -4043, 5146, 4134, -5736, -722, 1305, -5860, -1351, 4978, -3195, 5777, 3328, -1212, -5828, -5023, 1646, 4591, -5728, -3621, 2744, 3006, 2545, 2975, 563, 3504, 949, 2625, -3542, -2639, 4821, 726, 5374, -2837, 4611, -130, 4354, 1853, -2396, -4452, 140, 3296, -3949, 2768, 5067, -2197, 1635, -118, -2476, 4255, -5767, -827, 1177, 3748, 953, -2366, 2525, -1381, 3051, 3584, 4177, 4896, 4989, 5331, 2963, -4278, 1673, 1000, -544, -5791, 4320, -339, 2468, 81, -9, -1022, -3091, -480, 2842, 2294, 5911, 4890, 1062, 3932, 2731, 3553, 5542, -145, -4805, -3459, -3637, -3712, 4231, 2548, 3135, 355, -3382, 2747, 3707, 1759, -4846, 5179, 3694, 1326, 2089, 5092, -5086, -3284, -2881, -3014, -3241, -729, 3201, 3289, -2013, 790, -1696, -1428, 955, -334, -2426, 1170, 5755, -4632, -2319, -4388, -1260, -476, 1254, -973, -3531, 5435, 1359, -4905, -1922, -3879, -3985, 3998, 2033, -1544, 382, -316, 2178, 3988, -468, -420, -6093, -3710, 5559, -5446, -5339, -243, 1728, -400, -3016, 6137, -4948, -671, 3643, 5415, 3000, 5862, -6136, 3136, 56, -3199, 5191, -5206, 5529, 3400, -1987, -1702, 2399, -3565, -654, 4048, 1018, -5925, 2249, 1041, 3514, 2884, 5574, -1973, 1153, 2344, 1278, -3186, -5315, 4075, -5407, -4916, 4324, 2126, 522, -2169, -1630, 3262, -5079, -3510, 1050, 4536, -5332, -5445, -3860, -2865, 2683, -1190, 2370, 3818, -6118, 2969, 3789, 147, 3978, 5456, -4449, 2686, -4749, 5537, 3247, 4789, 4467, -1484, 4624, -6077, 4895, -3263, 3600, 2780, -6068, -3602, -5195, 4080, 421, -2645, 605, -2302, -4053, -504, 4213, 2305, -5886, -4782, 5042, 5594, 3029, -4372, -4212, 975, -2174, -3438, 2844, 4414, 1105, -142, 2847, 5681, -3477, 3271, 6008, 885, -4057, 5009, -1956, -1689, 1003, -3532, -3364, 241, 58, 1777, 2127, -151, -1663, 2839, -3957, 4654, -5906, 2505, 1426, 431, -1579, -2704, -3174, 52, -5291, 2766, -1323, -4938, 3336, 6055, -3636, 5874, -677, -4437, 2049, -4912, 3, -1321, 192, -3149, 3445, -4780, 160, -4698, -5057, 4919, -787, 3482, 113, -1010, 5468, -3915, 3127, 4169, -2166, 2920, 5241, -1912, 5257, -3834, -1378, 5919, 4433, 435, 5486, 3054, 4337, 1747, 3123, -2381, 2503, 2948, 5444, -5782, 1566, 4096, 64, -3656, -493, -683, -2459, -3248, 835, 6065, 1207, 3570, -4240, -5277, -1319, 3150, -1168, -709, -4046, 4645, -2078, -1112, 404, -4322, -1958, -2143, -441, -922, 1065, 1058, 4079, 2422, -5297, 6119, 6039, -3956, -1360, 2187, 1200, 5184, 2566, 2555, 6122, -2987, -1594, 1962, -6022, 5106, -5961, -3646, -2692, 168, 2437, -4298, -3329, 875, 4049, 3728, 3780, -1159, -5990, 1607, 948, 1146, 4976, 1404, -325, 4284, 2919, 3762, -5088, -4077, 4016, -5011, -652, -5766, -1002, -6099, -295, 545, -1293, 4737, -5019, 4774, -5429, -3704, 453, -5908, 2678, -418, -3772, 4143, -5333, 2031, -4714, -5876, -2281, -242, -156, 2767, -1537, 3969, -3991, 1440, 1805, 2882, 3763, 2051, -1954, 3066, 2447, -6142, -27, -576, -3963, 5084, 3529, -3434, -1632, -218, -2908, 4885, 1843, -2361, -1017, -4115, -3030, 1045, -4754, -1858, -2859, 426, 3315, 2481, -2925, -347, -5012, 3757, 1975, -5698, -723, -174, 2912, -1693, 3009, -4861, -2655, 5735, -354, 5868, 2738, -3778, -4493, 3202, 3833, 2057, -5369, -773, -5383, 1815, -390, -350, -1512, 1067, 5942, 1583, 5101, 1489, 2500, -442, -1483, -5915, -2401, -1263, -49
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
        ARRAY_N, NTT_N, LOGNTT_N, 4, {3, 2, 2, 2}
    };

// ================
// Generate twiddle factors for Cooley--Tukey FFT.

    zeta = OMEGA;
    coeff_ring.expZ(&omega, &zeta, 2);
    scale = 1;
    gen_streamlined_DWT_table(streamlined_twiddle_table,
        &scale, &omega, &zeta, profile, 0, coeff_ring);

    assert(memcmp(streamlined_twiddle_table, streamlined_NTT_table, (NTT_N - 1) * sizeof(int16_t)) == 0);

// ================
// Apply Cooley--Tukey FFT.

    compressed_CT_NTT(poly1,
        0, profile.compressed_layers - 1, streamlined_twiddle_table, profile, coeff_ring);
    compressed_CT_NTT(poly2,
        0, profile.compressed_layers - 1, streamlined_twiddle_table, profile, coeff_ring);

// ================

    point_mul(res, poly1, poly2, ARRAY_N, 1, coeff_ring);

// ================
// Generate twiddle factors for the inverse via Gentlemans--Sande FFT.

    zeta = OMEGA_INV;
    coeff_ring.expZ(&omega, &zeta, 2);
    scale = 1;
    gen_streamlined_DWT_table(streamlined_twiddle_table,
        &scale, &omega, &zeta, profile, 0, coeff_ring);

    assert(memcmp(streamlined_twiddle_table, streamlined_iNTT_table, (NTT_N - 1) * sizeof(int16_t)) == 0);

// ================
// Apply Gentleman--Sande FFT.

    compressed_GS_iNTT(res,
        0, profile.compressed_layers - 1, streamlined_twiddle_table, profile, coeff_ring);

// ================
// Multiply the scale to reference.

    scale = 512;
    for(size_t i = 0; i < ARRAY_N; i++){
        coeff_ring.mulZ(ref + i, ref + i, &scale);
    }

    for(size_t i = 0; i < ARRAY_N; i++){
        assert(ref[i] == res[i]);
    }

    printf("Test finished!\n");

}








