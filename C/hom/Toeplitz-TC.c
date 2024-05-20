
#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>

#include "tools.h"
#include "naive_mult.h"

/*

TODO:

Add Toeplitz matrix-vector product derivation from an algebra homomorphism.

*/

// ================
// This file demonstrates polynomial multiplication in Z_Q[x] / (x^4m + 1)
// via Toeplitz matrix-vector product built upon Toom-4 with the point set
// {0, 1, -1, 2, -2, 1/2, \infty}.

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

struct commutative_ring coeff_ring = {
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
// TMVP matrices built upon the above Toom-4
/*

(T^(-1))^*

=

diag(1/4, 1/2, 1/2, 1/8, 1/8, 1/2, 1/2)

diag(1, 1/3, 1/9, 1/9, 1/15, 1/45)

4,   -8,   -5,   10,    1,    -2, 0
0,   -4,    4,    9,   -1,    -2, 0
0,   -4,   12,   -7,   -3,     2, 0
0,    2,   -3,   -4,    3,     2, 0
0,    2,   -5,    0,    5,    -2, 0
0,    4,    0,   -5,    0,     1, 0
0,   -4,    8,    5,  -10,    -1, 2

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

// transpose of Toom-4 evaluation matrix
int32_t TC4_trunc_T[7][7] = {
{ 1,  1,  1,  1,  1,  8,  0},
{ 0,  1, -1,  2, -2,  4,  0},
{ 0,  1,  1,  4,  4,  2,  0},
{ 0,  1, -1,  8, -8,  1,  1},
{ 0,  0,  0,  0,  0,  0,  0},
{ 0,  0,  0,  0,  0,  0,  0},
{ 0,  0,  0,  0,  0,  0,  0}
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

// scaling of Hom-M from Toom-4 inversion matrix
int32_t iTC4_T_modified_scale[7] = {
1, -1431655765, 954437177, 954437177, -286331153, -1527099483, 1
};

// Hom-M from Toom-4 inversion matrix
int32_t iTC4_T_modified[7][7] = {
{4,   -8,   -5,   10,    1,    -2, 0},
{0,   -4,    4,    9,   -1,    -2, 0},
{0,   -4,   12,   -7,   -3,     2, 0},
{0,    2,   -3,   -4,    3,     2, 0},
{0,    2,   -5,    0,    5,    -2, 0},
{0,    4,    0,   -5,    0,     1, 0},
{0,   -4,    8,    5,  -10,    -1, 2}
};

// Hom-I from Toom-4 evaluation matrix
// We need to multiply the scales 1/8, 1/4, 1/2, 1/2 at the end.
int32_t TC4_trunc_T_modified[7][7] = {
{ 2,  4,  4,  1,  1, 32,  0},
{ 0,  2, -2,  1, -1,  8,  0},
{ 0,  1,  1,  1,  1,  2,  0},
{ 0,  1, -1,  2, -2,  1,  1},
{ 0,  0,  0,  0,  0,  0,  0},
{ 0,  0,  0,  0,  0,  0,  0},
{ 0,  0,  0,  0,  0,  0,  0}
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

// TMVP where the matrix is stored in the compressed format.
static
void TMVP(int32_t *des, int32_t *srcM, int32_t *srcV, size_t len){

    for(size_t i = 0; i < len; i++){
        des[i] = 0;
        for(size_t j = 0; j < len; j++){
            des[i] += srcM[len - 1 - i + j] * srcV[j];
        }
    }

}

// This function illustrate how to compute Z_Q[x] / (x^len + 1)
// with Toeplitz transformation built upon Toom-4 with the point set
// {0, 1, -1, 2, -2, 1/2, \infty} where
// len should be a multiple of 4.
static
void TMVP_TC4_negacyclic_mul(int32_t *des, int32_t *src1, int32_t *src2, size_t len){

    int32_t src1_V[4][len / 4];
    int32_t src1_V_full[7][len / 4];

    int32_t src2_Toeplitz[2 * len];
    int32_t src2_Toeplitz_full[7][len / 2];

    int32_t res_V[4][len / 4];
    int32_t res_V_full[7][len / 4];
    int32_t buff1[7], buff2[7], buff3[7];

    // Copy.
    memmove(&src1_V[0][0], src1, len * sizeof(int32_t));

    // Construct the compressed format of the Toeplitz matrix from
    // the multiplication map b -> a b mod x^len + 1.
    for(size_t i = 0; i < len; i++){
        src2_Toeplitz[i] = src2[len - 1 - i];
    }
    for(size_t i = len; i < 2 * len - 1; i++){
        src2_Toeplitz[i] = -src2[2 * len - 1 - i];
    }

    // Copying.
    for(size_t i = 0; i < 7; i++){
        memmove(&src2_Toeplitz_full[i][0], src2_Toeplitz + i * 4, ((len / 2) - 1) * sizeof(int32_t));
    }

    // Apply Hom-V.
    for(size_t i = 0; i < len / 4; i++){
        memset(buff1, 0, 7 * sizeof(int32_t));
        for(size_t j = 0; j < 4; j++){
            buff1[j] = src1_V[j][i];
        }
        matrix_vector_mul(buff3, (int32_t*)&TC4_trunc[0][0], buff1, 7);
        for(size_t j = 0; j < 7; j++){
            src1_V_full[j][i] = buff3[j];
        }
    }

    // Apply Hom-M.
    for(size_t i = 0; i < len / 2 - 1; i++){
        for(size_t j = 0; j < 7; j++){
            buff2[j] = src2_Toeplitz_full[j][i];
        }
        matrix_vector_mul(buff3, (int32_t*)&iTC4_T_modified[0][0], buff2, 7);
        for(size_t j = 0; j < 7; j++){
            buff3[j] *= iTC4_T_modified_scale[j];
        }
        for(size_t j = 0; j < 7; j++){
            src2_Toeplitz_full[j][i] = buff3[j];
        }
    }

    // Apply small-dimensional TMVP.
    for(size_t i = 0; i < 7; i++){
        TMVP((int32_t*)&res_V_full[i][0], (int32_t*)&src2_Toeplitz_full[i][0], (int32_t*)&src1_V_full[i][0], 4);
    }

    // Apply Hom-I.
    for(size_t i = 0; i < len / 4; i++){
        for(size_t j = 0; j < 7; j++){
            buff2[j] = res_V_full[j][i];
        }
        matrix_vector_mul(buff3, (int32_t*)&TC4_trunc_T_modified[0][0], buff2, 7);
        res_V[3][i] = buff3[0] >> 3;
        res_V[2][i] = buff3[1] >> 2;
        res_V[1][i] = buff3[2] >> 1;
        res_V[0][i] = buff3[3] >> 1;
    }

    memmove(des, (int32_t*)&res_V[0][0], len * sizeof(int32_t));

}

int main(void){

    int32_t poly1[16], poly2[16];
    int32_t ref[16], res[16];

    int32_t twiddle, t;

    for(size_t i = 0; i < 4; i++){
        t = rand();
        coeff_ring.memberZ(poly1 + i, &t);
        t = rand();
        coeff_ring.memberZ(poly2 + i, &t);
    }

    // Compute the product in Z_{2^32}[x] / (x^16 + 1).
    twiddle = -1;
    naive_mulR(ref, poly1, poly2, 16, &twiddle, coeff_ring);
    // Reduce from Z_{2^32} to Z_Q.
    for(size_t i = 0; i < 16; i++){
        cmod_int32(ref + i, ref + i, &mod);
    }

    // Compute the product in Z_{2^32}[x] / (x^16 + 1) with TMVP built upon the
    // Toom-4 with the point set {0, 1, -1, 2, -2, 1/2, \infty}.
    TMVP_TC4_negacyclic_mul(res, poly1, poly2, 16);
    // Reduce from Z_{2^32} to Z_Q.
    for(size_t i = 0; i < 16; i++){
        cmod_int32(res + i, res + i, &mod);
    }

    for(size_t i = 0; i < 16; i++){
        assert(ref[i] == res[i]);
    }

    printf("Test finished!\n");



}
























