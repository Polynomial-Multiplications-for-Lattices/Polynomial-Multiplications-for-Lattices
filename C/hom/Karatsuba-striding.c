
#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>

#include "tools.h"
#include "naive_mult.h"

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

struct commutative_ring coeff_ring = {
    .sizeZ = sizeof(int32_t),
    .memberZ = memberZ,
    .addZ = addZ,
    .subZ = subZ,
    .mulZ = mulZ,
    .expZ = expZ
};

// Multiply two size-64 polynomials in R[x] / (x^64 + 1) via
// R[x] / (x^64 + 1)
// to
// (R[y] / (y^16 + 1)) / (x^4 - y)
// to
// (R[y] / (y^16 + 1)) / (x^7) (2 layers of Karatsuba)
static
void negacyclic_Karatsuba_striding(int32_t *des, const int32_t *src1, const int32_t *src2, size_t len){

    int32_t res_p00[len / 4], res_p01[len / 4], res_p02[len / 4],
            res_p10[len / 4], res_p11[len / 4], res_p12[len / 4],
            res_p20[len / 4], res_p21[len / 4], res_p22[len / 4];

    int32_t p00, p01, p02,
            p10, p11, p12,
            p20, p21, p22;

    int32_t _p01, _p11, _p21;

    int32_t buff[3];

    memset(res_p00, 0, sizeof(res_p00));
    memset(res_p01, 0, sizeof(res_p01));
    memset(res_p02, 0, sizeof(res_p02));
    memset(res_p10, 0, sizeof(res_p10));
    memset(res_p11, 0, sizeof(res_p11));
    memset(res_p12, 0, sizeof(res_p12));
    memset(res_p20, 0, sizeof(res_p20));
    memset(res_p21, 0, sizeof(res_p21));
    memset(res_p22, 0, sizeof(res_p22));

    for(size_t i = 0; i < len / 4; i++){

        // Load 4, cache 9.
        p00 = src1[4 * i + 0];
        p02 = src1[4 * i + 1];
        p20 = src1[4 * i + 2];
        p22 = src1[4 * i + 3];

        p01 = p00 + p02;
        p21 = p20 + p22;

        p10 = p00 + p20;
        p12 = p02 + p22;

        p11 = p10 + p12;

        for(size_t j = 0; j < len / 4 - i; j++){

            res_p00[i + j]           += p00 * src2[4 * j + 0];
            res_p02[i + j]           += p02 * src2[4 * j + 1];
            _p01                      = src2[4 * j + 0] + src2[4 * j + 1];
            res_p01[i + j]           += p01 * _p01;

            res_p20[i + j]           += p20 * src2[4 * j + 2];
            res_p22[i + j]           += p22 * src2[4 * j + 3];
            _p21                      = src2[4 * j + 2] + src2[4 * j + 3];
            res_p21[i + j]           += p21 * _p21;

            res_p10[i + j]           += p10 * (src2[4 * j + 0] + src2[4 * j + 2]);
            res_p12[i + j]           += p12 * (src2[4 * j + 1] + src2[4 * j + 3]);
            _p11                      = _p01 + _p21;
            res_p11[i + j]           += p11 * _p11;

        }

        for(size_t j = len / 4 - i; j < len / 4; j++){

            res_p00[i + j - len / 4] -= p00 * src2[4 * j + 0];
            res_p02[i + j - len / 4] -= p02 * src2[4 * j + 1];
            _p01                      = src2[4 * j + 0] + src2[4 * j + 1];
            res_p01[i + j - len / 4] -= p01 * _p01;

            res_p20[i + j - len / 4] -= p20 * src2[4 * j + 2];
            res_p22[i + j - len / 4] -= p22 * src2[4 * j + 3];
            _p21                      = src2[4 * j + 2] + src2[4 * j + 3];
            res_p21[i + j - len / 4] -= p21 * _p21;

            res_p10[i + j - len / 4] -= p10 * (src2[4 * j + 0] + src2[4 * j + 2]);
            res_p12[i + j - len / 4] -= p12 * (src2[4 * j + 1] + src2[4 * j + 3]);
            _p11                      = _p01 + _p21;
            res_p11[i + j - len / 4] -= p11 * _p11;

        }

    }

    // Load 9, store 6.
    for(size_t i = 0; i < len / 4; i++){
        res_p01[i] = res_p01[i] - res_p00[i] - res_p02[i];
        res_p11[i] = res_p11[i] - res_p10[i] - res_p12[i];
        res_p21[i] = res_p21[i] - res_p20[i] - res_p22[i];

        res_p10[i] = res_p10[i] - res_p00[i] - res_p20[i];
        res_p11[i] = res_p11[i] - res_p01[i] - res_p21[i];
        res_p12[i] = res_p12[i] - res_p02[i] - res_p22[i];

        res_p10[i] += res_p02[i];
        res_p20[i] += res_p12[i];
    }

    // Load 11, store 4, Cache 3.
    des[0]  = res_p00[0] - res_p20[len / 4 - 1];
    des[1]  = res_p01[0] - res_p21[len / 4 - 1];
    des[2]  = res_p10[0] - res_p22[len / 4 - 1];
    des[3]  = res_p11[0];
    buff[0] = res_p20[0];
    buff[1] = res_p21[0];
    buff[2] = res_p22[0];
    // Load 7, store 4, cache 3.
    for(size_t i = 1; i < len / 4 - 1; i++){
        des[4 * i + 0] = buff[0] + res_p00[i];
        des[4 * i + 1] = buff[1] + res_p01[i];
        des[4 * i + 2] = buff[2] + res_p10[i];
        des[4 * i + 3] =           res_p11[i];
        buff[0]        =           res_p20[i];
        buff[1]        =           res_p21[i];
        buff[2]        =           res_p22[i];
    }
    // Load 4, store 4.
    des[len - 4] = buff[0] + res_p00[len / 4 - 1];
    des[len - 3] = buff[1] + res_p01[len / 4 - 1];
    des[len - 2] = buff[2] + res_p10[len / 4 - 1];
    des[len - 1] =           res_p11[len / 4 - 1];

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

    // Compute the product in Z_{2^32} [x] / (x^ARRAY_N + 1)
    // via striding followed by two layers of Karatsuba.
    negacyclic_Karatsuba_striding(res, poly1, poly2, ARRAY_N);
    // Reduce from Z_{2^32} to Z_Q.
    for(size_t i = 0; i < ARRAY_N; i++){
        cmod_int32(res + i, res + i, &mod);
    }

    for(size_t i = 0; i < ARRAY_N; i++){
        assert(ref[i] == res[i]);
    }

    printf("Test finished!\n");


}




