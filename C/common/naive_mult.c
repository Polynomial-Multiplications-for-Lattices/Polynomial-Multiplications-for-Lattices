
#include <memory.h>

#include "naive_mult.h"

// ================================

// Multiplying size-len polynomials stored at src1 and src2 in in R[x] / (x^len - twiddle)
// where R = ring.
void naive_mulR(
    void *des,
    const void *src1, const void *src2,
    size_t len, const void *twiddle,
    struct ring ring
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

// Multiplying size-len polynomials stored at src1 and src2 in R[x] where R = ring.
void naive_mul_long(
    void *des,
    const void *src1, const void *src2,
    size_t len,
    struct ring ring
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

// Point-wise multiplication of src1[len * jump] by src2[len] over R where R = ring.
void point_mul(
    void *des,
    const void *src1, const void *src2,
    size_t len, size_t jump,
    struct ring ring
    ){

    for(size_t i = 0; i < len; i++){
        for(size_t j = 0; j < jump; j++){
            ring.mulZ(des + (i * jump + j) * ring.sizeZ, src1 + (i * jump + j) * ring.sizeZ, src2 + i * ring.sizeZ);
        }
    }

}







