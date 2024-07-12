
#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>

#include "tools.h"
#include "naive_mult.h"

// ================
// This file demonstrates polynomial multiplication of size-4m polynomials in Z_Q[x]
// via Toom-4 with the point set {0, 1, -1, 2, -2, 1/2, \infty}.

// ================
// Theory.

// ================
// A small example.

// ================
// Optimization guide.
/*

1. While applying the matrices, make use of barrel shifter to multiply
   by the correct constants. For example, 3 * a can be implemented as
   a + (a << 1) and 9 * a can beimplemented as a + (a << 3).
   With such optimization, we save the memory operations loading the matrices.

2. Notice that point set is carefully chosen. In principle, when an integer z
   is chosen, we also choose -z. So evaluating a polynomial at {z, -z} will be faster
   by first evaluating for the odds and evens individually and applying an add-sub pair.

*/

// ================
// Applications to lattice-based cryptosystems.

#define ARRAY_N 256

// Q = 1, 2, 4, ..., 2^29
#define Q (1 << 29)

int32_t mod = Q;

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
// Toom-4 matrices
/*

// ========

T
=

 1,  0,  0,  0,  0,   0,  0
 1,  1,  1,  1,  1,   1,  1
 1, -1,  1, -1,  1,  -1,  1
 1,  2,  4,  8, 16,  32, 64
 1, -2,  4, -8, 16, -32, 64
64, 32, 16,  8,  4,   2,  1
 0,  0,  0,  0,  0,   0,  1

// ========

T^{-1}
=

   1,     0,     0,     0,     0,     0,    0
  -2,  -2/3,  -2/9,  1/36,  1/60,  2/45,   -2
-5/4,   2/3,   2/3, -1/24, -1/24,     0,    4
 5/2,   3/2, -7/18, -1/18,     0, -1/18,  5/2
 1/4,  -1/6,  -1/6,  1/24,  1/24,     0,   -5
-1/2,  -1/3,   1/9,  1/36, -1/60,  1/90, -1/2
   0,     0,     0,     0,     0,     0,    1

=

diag(1, 1/4, 1/8, 1/2, 1/8, 1/4, 1)

   1,    0,    0,    0,     0,    0,   0
  -8, -8/3, -8/9,  1/9,  1/15, 8/45,  -8
 -10, 16/3, 16/3, -1/3,  -1/3,    0,  32
   5,    3, -7/9, -1/9,     0, -1/9,   5
   2, -4/3, -4/3,  1/3,   1/3,    0, -40
  -2, -4/3,  4/9,  1/9, -1/15, 2/45,  -2
   0,    0,    0,    0,     0,    0,   1

=

 4,     0,    0,    0,     0,     0,    0
-8,  -4/3, -4/9,  2/9,  2/15,  4/45,   -4
-5,   4/3,  4/3, -1/3,  -1/3,     0,    8
10,     3, -7/9, -4/9,     0,  -1/9,    5
 1,  -1/3, -1/3,  1/3,   1/3,     0,  -10
-2,  -2/3,  2/9,  2/9, -2/15,  1/45,   -1
 0,     0,    0,    0,     0,     0,    2

diag(1/4, 1/2, 1/2, 1/8, 1/8, 1/2, 1/2)
*/

// ================
// inverses modulo 2^32

// 3^(-1) = -1431655765
// 5^(-1) = -858993459
// 9^(-1) = 954437177
// 15^(-1) = -286331153
// 45^(-1) = -1527099483

// ================

// Toom-4 evaluation full matrix
int32_t TC4[7][7] = {
{ 1,  0,  0,  0,  0,   0,  0},
{ 1,  1,  1,  1,  1,   1,  1},
{ 1, -1,  1, -1,  1,  -1,  1},
{ 1,  2,  4,  8, 16,  32, 64},
{ 1, -2,  4, -8, 16, -32, 64},
{64, 32, 16,  8,  4,   2,  1},
{ 0,  0,  0,  0,  0,   0,  1}
};

// Toom-4 evaluation matrix
// This is also referred as Hom-V in TMVP.
int32_t TC4_trunc[7][7] = {
{ 1,  0,  0,  0, 0, 0, 0},
{ 1,  1,  1,  1, 0, 0, 0},
{ 1, -1,  1, -1, 0, 0, 0},
{ 1,  2,  4,  8, 0, 0, 0},
{ 1, -2,  4, -8, 0, 0, 0},
{ 8,  4,  2,  1, 0, 0, 0},
{ 0,  0,  0,  1, 0, 0, 0}
};

// Toom-4 inversion matrix
int32_t iTC4[7][7] = {
{   1,           0,           0,           0,           0,          0,   0},
{  -8, -1431655768,   954437176,   954437177,  -286331153,  668106024,  -8},
{ -10, -1431655760, -1431655760,  1431655765,  1431655765,          0,  32},
{   5,           3,  1908874353,  -954437177,           0, -954437177,   5},
{   2,  1431655764,  1431655764, -1431655765, -1431655765,          0, -40},
{  -2,  1431655764,  -477218588,   954437177,   286331153, 1240768330,  -2},
{   0,           0,           0,           0,           0,          0,   1}
};

// matrix-vector multiplication
static
void matrix_vector_mul(int32_t *des, int32_t *srcM, int32_t *srcV, size_t len){

    int32_t buff[len];

    memset(buff, 0, sizeof(buff));

    for(size_t i = 0; i < len; i++){
        for(size_t k = 0; k < len; k++){
            buff[i] += srcM[i * len + k] * srcV[k];
        }
    }

    memmove(des, buff, len * sizeof(int32_t));

}

// len must be a 4-multiple.
// This function computes the product of two size-len polynomials in Z_{2^32}[x]
// using Toom-4 with the point set {0, 1, -1, 2, -2, 1/2, \infty}.
// Notice that matrices are modifed to ensure the well-defineness over Z_{2^32}.
static
void TC_striding_mul(int32_t *des, int32_t *src1, int32_t *src2, size_t len){

    int32_t src1_extended[7][len / 4], src2_extended[7][len / 4];
    int32_t res[7][len / 4];
    int32_t TC4_buff[7];
    int32_t twiddle;

    memset(src1_extended, 0, sizeof(src1_extended));
    memset(src2_extended, 0, sizeof(src2_extended));

    for(size_t i = 0; i < len / 4; i++){
        for(size_t j = 0; j < 4; j++){
            src1_extended[j][i] = src1[i * 4 + j];
        }
    }

    for(size_t i = 0; i < len / 4; i++){
        for(size_t j = 0; j < 4; j++){
            src2_extended[j][i] = src2[i * 4 + j];
        }
    }

    // Apply Toom-4 evaluation matrix.
    for(size_t i = 0; i < len / 4; i++){
        for(size_t j = 0; j < 7; j++){
            TC4_buff[j] = src1_extended[j][i];
        }
        matrix_vector_mul(TC4_buff, (int32_t*)&TC4_trunc[0][0], TC4_buff, 7);
        for(size_t j = 0; j < 7; j++){
            src1_extended[j][i] = TC4_buff[j];
        }
    }

    for(size_t i = 0; i < len / 4; i++){
        for(size_t j = 0; j < 7; j++){
            TC4_buff[j] = src2_extended[j][i];
        }
        matrix_vector_mul(TC4_buff, (int32_t*)&TC4_trunc[0][0], TC4_buff, 7);
        for(size_t j = 0; j < 7; j++){
            src2_extended[j][i] = TC4_buff[j];
        }
    }

    // Compute small-dimensional products.
    twiddle = -1;
    for(size_t i = 0; i < 7; i++){
        naive_mulR((int32_t*)&res[i][0], (int32_t*)&src1_extended[i][0], (int32_t*)&src2_extended[i][0], len / 4, &twiddle, coeff_ring);
    }

    // Apply Toom-4 inversion matrix.
    for(size_t i = 0; i < len / 4; i++){
        for(size_t j = 0; j < 7; j++){
            TC4_buff[j] = res[j][i];
        }
        matrix_vector_mul(TC4_buff, (int32_t*)&iTC4[0][0], TC4_buff, 7);
        // Multiply by powers of two.
        res[0][i] = TC4_buff[0];
        res[1][i] = TC4_buff[1] >> 2;
        res[2][i] = TC4_buff[2] >> 3;
        res[3][i] = TC4_buff[3] >> 1;
        res[4][i] = TC4_buff[4] >> 3;
        res[5][i] = TC4_buff[5] >> 2;
        res[6][i] = TC4_buff[6];
    }

    // Export the result.
    for(size_t i = 0; i < len / 4; i++){
        des[i * 4 + 3] = res[3][i];
    }

    for(size_t i = 0; i < len / 4 - 1; i++){
        for(size_t j = 4; j < 7; j++){
            des[(i + 1) * 4 + j - 4] = res[j - 4][i + 1] + res[j][i];
        }
    }

    for(size_t j = 4; j < 7; j++){
        des[j - 4] = res[j - 4][0] - res[j][len / 4 - 1];
    }

}

int main(void){

    int32_t poly1[ARRAY_N], poly2[ARRAY_N];
    int32_t ref[ARRAY_N], res[ARRAY_N];

    int32_t t;
    int32_t twiddle;

    for(size_t i = 0; i < ARRAY_N; i++){
        t = rand();
        coeff_ring.memberZ(poly1 + i, &t);
        t = rand();
        coeff_ring.memberZ(poly2 + i, &t);
    }

    twiddle = -1;
    // Compute the product in Z_{2^32}[x].
    naive_mulR(ref, poly1, poly2, ARRAY_N, &twiddle, coeff_ring);
    // Reduce from Z_{2^32} to Z_Q.
    for(size_t i = 0; i < ARRAY_N; i++){
        cmod_int32(ref + i, ref + i, &mod);
    }

    // Compute the product in Z_{2^32}[x] / (x^ARRAY_N + 1) via
    // striding followed by Toom-4 with the point set
    // {0, 1, -1, 2, -2, 1/2, \infty}.
    TC_striding_mul(res, poly1, poly2, ARRAY_N);
    // Reduce from Z_{2^32} to Z_Q.
    for(size_t i = 0; i < ARRAY_N; i++){
        cmod_int32(res + i, res + i, &mod);
    }

    for(size_t i = 0; i < ARRAY_N; i++){
        assert(ref[i] == res[i]);
    }

    printf("Test finished!\n");



}
























