
#include <memory.h>

#include "tools.h"

// ================================
// Candidates for memberZ.

void cmod_int16(void *des, const void *src, const void *mod){
    int16_t mod_v = *(int16_t*)mod;
    int16_t t = (*(int16_t*)src) % mod_v;
    if(t < -(mod_v >> 1)){
        t += mod_v;
    }
    if(t > (mod_v >> 1)){
        t -= mod_v;
    }
    *(int16_t*)des = t;
}

void cmod_int32(void *des, const void *src, const void *mod){
    int32_t mod_v = *(int32_t*)mod;
    int32_t t = (*(int32_t*)src) % mod_v;
    if(t < -(mod_v >> 1)){
        t += mod_v;
    }
    if(t > (mod_v >> 1)){
        t -= mod_v;
    }
    *(int32_t*)des = t;
}

void cmod_int64(void *des, const void *src, const void *mod){
    int64_t mod_v = *(int64_t*)mod;
    int64_t t = (*(int64_t*)src) % mod_v;
    if(t < -(mod_v >> 1)){
        t += mod_v;
    }
    if(t > (mod_v >> 1)){
        t -= mod_v;
    }
    *(int64_t*)des = t;
}

// ================================
// Candidates for addZ.

void addmod_int16(void *des, const void *src1, const void *src2, const void *mod){

    int32_t tmp_v, mod_v, des_v;

    tmp_v = (int32_t)(*(int16_t*)src1) + (int32_t)(*(int16_t*)src2);
    mod_v = (int32_t)(*(int16_t*)mod);

    cmod_int32(&des_v, &tmp_v, &mod_v);

    *(int16_t*)des = (int16_t)des_v;

}

void addmod_int32(void *des, const void *src1, const void *src2, const void *mod){

    int64_t tmp_v, mod_v, des_v;

    tmp_v = (int64_t)(*(int32_t*)src1) + (int64_t)(*(int32_t*)src2);
    mod_v = (int64_t)(*(int32_t*)mod);

    cmod_int64(&des_v, &tmp_v, &mod_v);

    *(int32_t*)des = (int32_t)des_v;

}

// ================================
// Candidates for subZ.

void submod_int16(void *des, const void *src1, const void *src2, const void *mod){

    int32_t tmp_v, mod_v, des_v;

    tmp_v = (int32_t)(*(int16_t*)src1) - (int32_t)(*(int16_t*)src2);
    mod_v = (int32_t)(*(int16_t*)mod);

    cmod_int32(&des_v, &tmp_v, &mod_v);

    *(int16_t*)des = (int16_t)des_v;

}

void submod_int32(void *des, const void *src1, const void *src2, const void *mod){

    int64_t tmp_v, mod_v, des_v;

    tmp_v = (int64_t)(*(int32_t*)src1) - (int64_t)(*(int32_t*)src2);
    mod_v = (int64_t)(*(int32_t*)mod);

    cmod_int64(&des_v, &tmp_v, &mod_v);

    *(int32_t*)des = (int32_t)des_v;

}

// ================================
// Candidates for mulZ.

void mulmod_int16(void *des, const void *src1, const void *src2, const void *mod){

    int32_t tmp_v, mod_v, des_v;

    tmp_v = (int32_t)(*(int16_t*)src1) * (int32_t)(*(int16_t*)src2);
    mod_v = (int32_t)(*(int16_t*)mod);

    cmod_int32(&des_v, &tmp_v, &mod_v);

    *(int16_t*)des = (int16_t)des_v;

}

void mulmod_int32(void *des, const void *src1, const void *src2, const void *mod){

    int64_t tmp_v, mod_v, des_v;

    tmp_v = (int64_t)(*(int32_t*)src1) * (int64_t)(*(int32_t*)src2);
    mod_v = (int64_t)(*(int32_t*)mod);

    cmod_int64(&des_v, &tmp_v, &mod_v);

    *(int32_t*)des = (int32_t)des_v;

}

// ================================
// Candidates for expZ.

void expmod_int16(void *des, const void *src, size_t e, const void *mod){

    int16_t src_v = *(int16_t*)src;
    int16_t tmp_v;

    tmp_v = 1;
    for(; e; e >>= 1){
        if(e & 1){
            mulmod_int16(&tmp_v, &tmp_v, &src_v, mod);
        }
        mulmod_int16(&src_v, &src_v, &src_v, mod);
    }

    memcpy(des, &tmp_v, sizeof(int16_t));

}

void expmod_int32(void *des, const void *src, size_t e, const void *mod){

    int32_t src_v = *(int32_t*)src;
    int32_t tmp_v;

    tmp_v = 1;
    for(; e; e >>= 1){
        if(e & 1){
            mulmod_int32(&tmp_v, &tmp_v, &src_v, mod);
        }
        mulmod_int32(&src_v, &src_v, &src_v, mod);
    }

    memcpy(des, &tmp_v, sizeof(int32_t));

}

// ================================
// In-place bit-reversal.

void bitreverse(void *src, size_t len, size_t size){

    char tmp[size];

    for(size_t i = 0, j = 0; i < len; i++){
        if(i < j){
            memcpy(tmp, src + i * size, size);
            memcpy(src + i * size, src + j * size, size);
            memcpy(src + j * size, tmp, size);
        }
        for(size_t k = len >> 1; (j ^= k) < k; k >>= 1);
    }

}

#if defined(__x86_64__) || defined(__aarch64__)

// void cmod_int128(void *des, void *src, void *mod){
//     __int128 mod_v = *(__int128*)mod;
//     __int128 t = (*(__int128*)src) % mod_v;
//     if(t >= (mod_v >> 1)){
//         t -= mod_v;
//     }
//     if(t < -(mod_v >> 1)){
//         t += mod_v;
//     }
//     *(__int128*)des = t;
// }

// void addmod_int128(void *des, void *src1, void *src2, void *mod){

//     __int128 tmp_v, mod_v, des_v;

//     tmp_v = (__int128)(*(int64_t*)src1) + (__int128)(*(int64_t*)src2);
//     mod_v = (__int128)(*(int64_t*)mod);

//     cmod_int128(&des_v, &tmp_v, &mod_v);

//     *(int64_t*)des = (int64_t)des_v;

// }

#endif













