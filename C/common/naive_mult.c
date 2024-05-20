
#include <stdlib.h>
#include <memory.h>

#include "naive_mult.h"

// ================================
// Multiplying polynomials in in R[x] / (x^len - twiddle)
// where R is the ring defined by mod, addmod, mulmod.
void naive_mulR(
    void *des,
    void *src1, void *src2,
    size_t len, void *twiddle,
    struct commutative_ring ring
    ){

    char buff[(len << 1) * ring.sizeZ];
    char tmp[ring.sizeZ];

    memset(buff, 0, (len << 1) * ring.sizeZ);

    for(size_t i = 0; i < len; i++){
        for(size_t j = 0; j < len; j++){
            ring.mulZ(tmp, src1 + i * ring.sizeZ, src2 + j * ring.sizeZ);
            ring.addZ(buff + (i + j) * ring.sizeZ, buff + (i + j) * ring.sizeZ, tmp);
        }
    }

    for(size_t i = ((len - 1) << 1); i >= len; i--){
        ring.mulZ(tmp, buff + i * ring.sizeZ, twiddle);
        ring.addZ(des + (i - len) * ring.sizeZ, buff + (i - len) * ring.sizeZ, tmp);
    }
    memcpy(des + (len - 1) * ring.sizeZ, buff + (len - 1) * ring.sizeZ, ring.sizeZ);

}

// ================================
// Multiplying size-len polynomials in R[x]
// where R is the ring defined by mod, addmod, mulmod.
void naive_mul_long(
    void *des,
    void *src1, void *src2,
    size_t len,
    struct commutative_ring ring
    ){

    char buff[(len << 1) * ring.sizeZ];
    char tmp[ring.sizeZ];

    memset(buff, 0, (len << 1) * ring.sizeZ);

    for(size_t i = 0; i < len; i++){
        for(size_t j = 0; j < len; j++){
            ring.mulZ(tmp, src1 + i * ring.sizeZ, src2 + j * ring.sizeZ);
            ring.addZ(buff + (i + j) * ring.sizeZ, buff + (i + j) * ring.sizeZ, tmp);
        }
    }
    memcpy(des, buff, (2 * len - 1) * ring.sizeZ);

}

// ================================
// Point-wise multiplication of src1[len * jump] by src2[len].
// In particular, for i in {0, ..., len - 1} and j in {0, ..., jump - 1},
// src1[i * jump + j] is multiplied by src2[i].
void point_mul(
    void *des,
    void *src1, void *src2,
    size_t len, size_t jump,
    struct commutative_ring ring
    ){

    for(size_t i = 0; i < len; i++){
        for(size_t j = 0; j < jump; j++){
            ring.mulZ(des + (i * jump + j) * ring.sizeZ, src1 + (i * jump + j) * ring.sizeZ, src2 + i * ring.sizeZ);
        }
    }

}







