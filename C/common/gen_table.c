#include <stdlib.h>
#include <memory.h>

#include "tools.h"
#include "gen_table.h"

// ================================
// Generate twiddle factors for cyclic NTT with Cooley-Tukey butterflies.
void gen_CT_table(
    void *des,
    void *scale, void *omega,
    struct compress_profile _profile,
    struct commutative_ring ring
    ){

    char zeta[ring.sizeZ];
    char twiddle[ring.sizeZ];

    memcpy(zeta, omega, ring.sizeZ);

    memcpy(twiddle, scale, ring.sizeZ);
    for(size_t i = 0; i < (_profile.ntt_n >> 1); i++){
        memcpy(des, twiddle, ring.sizeZ);
        des += ring.sizeZ;
        ring.mulZ(twiddle, twiddle, zeta);
    }

    des -= ring.sizeZ * (_profile.ntt_n >> 1);
    bitreverse(des, _profile.ntt_n >> 1, ring.sizeZ);

}

// ================================
// Generate twiddle factors for DWT with Cooley-Tukey butterflies.
void gen_DWT_table(
    void *des,
    void *scale, void *omega, void *zeta,
    struct compress_profile _profile,
    struct commutative_ring ring
    ){

    char buff[_profile.ntt_n * ring.sizeZ];
    char zeta_buff[_profile.log_ntt_n * ring.sizeZ];


    gen_CT_table(buff, scale, omega, _profile, ring);

    memcpy(zeta_buff + (_profile.log_ntt_n - 1) * ring.sizeZ, zeta, ring.sizeZ);

    for(ssize_t i = _profile.log_ntt_n - 2; i >= 0; i--){
        ring.expZ(zeta_buff + i * ring.sizeZ, zeta_buff + (i + 1) * ring.sizeZ, 2);
    }

    for(size_t i = 0; i < _profile.log_ntt_n; i++){
        for(size_t j = 0; j < (1u << i); j++){
            ring.mulZ(des + j * ring.sizeZ, buff + j * ring.sizeZ, zeta_buff + i * ring.sizeZ);
        }
        des += (1u << i) * ring.sizeZ;
    }


}

// ================================
// Generate twiddle factors for cyclic iNTT with Cooley-Tukey butterflies.
void gen_inv_CT_table(
    void *des,
    void *scale, void *omega,
    struct compress_profile _profile,
    struct commutative_ring ring
    ){

    char zeta[ring.sizeZ];
    char twiddle[ring.sizeZ];

    for(size_t level = 0; level < _profile.log_ntt_n; level++){
        ring.expZ(zeta, omega, (1u << _profile.log_ntt_n) >> (level + 1));
        memcpy(twiddle, scale, ring.sizeZ);
        for(size_t i = 0; i < (1u << level); i++){
            memcpy(des, twiddle, ring.sizeZ);
            des += ring.sizeZ;
            ring.mulZ(twiddle, twiddle, zeta);
        }
    }

}

// ================================
// Generate twiddle factors for DWT with Cooley-Tukey butterflies.
// The table is re-ordered according to _profile.
void gen_streamlined_DWT_table(
    void *des,
    void *scale, void *omega, void *zeta,
    struct compress_profile _profile, bool pad,
    struct commutative_ring ring
    ){

    size_t start_level;

    char tmp[_profile.ntt_n * ring.sizeZ];
    void *level_ptr[_profile.log_ntt_n];

    gen_DWT_table(
        tmp, scale, omega, zeta,
        _profile,
        ring
    );

    for(size_t i = 0; i < _profile.log_ntt_n; i++){
        level_ptr[i] = tmp + ring.sizeZ * ((1 << i) - 1);
    }

    start_level = 0;
    for(size_t i = 0; i < _profile.compressed_layers; i++){
        for(size_t j = 0; j < (1u << start_level); j++){
            if(pad){
                memset(des, 0, ring.sizeZ);
                des += ring.sizeZ;
            }
            for(size_t k = 0; k < (_profile.merged_layers[i]); k++){
                for(size_t h = 0; h < (1u << k); h++){
                    memcpy(des,
                        level_ptr[start_level + k] + (j * (1 << k) + h) * ring.sizeZ,
                        ring.sizeZ);
                    des += ring.sizeZ;
                }
            }
        }
    start_level += (_profile.merged_layers)[i];
    }

}

// ================================
// Generate twiddle factors for cyclic iNTT with Cooley-Tukey butterflies.
// The table is re-ordered according to _profile.
void gen_streamlined_inv_CT_table(
    void *des,
    void *scale, void *omega,
    struct compress_profile _profile, bool pad,
    struct commutative_ring ring
    ){

    char zeta[ring.sizeZ];
    size_t start_level;

    char tmp[_profile.ntt_n * ring.sizeZ];

    void *level_ptr[_profile.log_ntt_n];

    memcpy(zeta, omega, ring.sizeZ);

    gen_inv_CT_table(
        tmp, scale, zeta,
        _profile,
        ring
    );

    for(size_t i = 0; i < _profile.log_ntt_n; i++){
        level_ptr[i] = tmp + ring.sizeZ * ((1 << i) - 1);
    }

    start_level = 0;
    for(size_t i = 0; i < _profile.compressed_layers; i++){
        for(size_t j = 0; j < (1u << start_level); j++){
            if(pad){
                memset(des, 0, ring.sizeZ);
                des += ring.sizeZ;
            }
            for(size_t k = 0; k < (_profile.merged_layers[i]); k++){
                for(size_t h = 0; h < (1u << k); h++){
                    memcpy(
                        des,
                        level_ptr[start_level + k] + (j + (h << start_level)) * ring.sizeZ,
                        ring.sizeZ);
                    des += ring.sizeZ;
                }
            }
        }
        start_level += (_profile.merged_layers)[i];
    }

}

// ================================
// Generate twiddle factors for twisting (x^NTT_N - omega^NTT_N) to (x^NTT_N - 1).
void gen_twist_table(
    void *des,
    void *scale, void *omega,
    struct compress_profile _profile,
    struct commutative_ring ring
    ){

    char zeta[ring.sizeZ];
    char twiddle[ring.sizeZ];

    memcpy(zeta, omega, ring.sizeZ);

    memcpy(twiddle, scale, ring.sizeZ);
    for(size_t i = 0; i < _profile.ntt_n; i++){
        memcpy(des, twiddle, ring.sizeZ);
        des += ring.sizeZ;
        ring.mulZ(twiddle, twiddle, zeta);
    }

}

// ================================
// Generate twiddle factors for multiplication in x^(ARRAY_N / NTT_N) +- omega^i.
void gen_mul_table(
    void *des,
    void *scale, void *omega,
    struct compress_profile _profile,
    struct commutative_ring ring
    ){

    char zeta[ring.sizeZ];
    char twiddle[ring.sizeZ];

    memcpy(zeta, omega, ring.sizeZ);

    memcpy(twiddle, scale, ring.sizeZ);
    for(size_t i = 0; i < (_profile.ntt_n >> 1); i++){
        memcpy(des, twiddle, ring.sizeZ);
        des += ring.sizeZ;
        ring.mulZ(twiddle, twiddle, zeta);
    }

    des -= (_profile.ntt_n >> 1) * ring.sizeZ;

    bitreverse(des, _profile.ntt_n >> 1, ring.sizeZ);

}








