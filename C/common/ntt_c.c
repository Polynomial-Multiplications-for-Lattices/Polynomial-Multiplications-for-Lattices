
#include <memory.h>
#include <sys/types.h>

#include "tools.h"
#include "ntt_c.h"

// ================================
// Cooley-Tukey butterfly.
void CT_butterfly(
    void *src,
    size_t indx_a, size_t indx_b,
    void *twiddle,
    struct ring ring
    ){

    char tmp[ring.sizeZ];

    ring.mulZ(tmp, src + indx_b * ring.sizeZ, twiddle);
    ring.subZ(src + indx_b * ring.sizeZ, src + indx_a * ring.sizeZ, tmp);
    ring.addZ(src + indx_a * ring.sizeZ, src + indx_a * ring.sizeZ, tmp);

}

// ================================
// Gentleman-Sande butterfly.
void GS_butterfly(
    void *src,
    size_t indx_a, size_t indx_b,
    void *twiddle,
    struct ring ring
    ){

    char tmp[ring.sizeZ];

    ring.subZ(tmp, src + indx_a * ring.sizeZ, src + indx_b * ring.sizeZ);
    ring.addZ(src + indx_a * ring.sizeZ, src + indx_a * ring.sizeZ, src + indx_b * ring.sizeZ);
    ring.mulZ(src + indx_b * ring.sizeZ, tmp, twiddle);

}

// ================================
void CT_NTT_core(
    void *src,
    size_t level,
    void *_root_table,
    struct compress_profile _profile,
    struct ring ring
    ){

    size_t step;
    void *real_root_table;

    step = (_profile.array_n) >> (level + 1);
    real_root_table = _root_table + ((1u << level) - 1) * ring.sizeZ;

    for(size_t i = 0; i < _profile.array_n; i += 2 * step){
        for(size_t j = 0; j < step; j++){
            CT_butterfly(src + (i + j) * ring.sizeZ, 0, step, real_root_table, ring);
        }
        real_root_table += ring.sizeZ;
    }

}

// ================================
void CT_iNTT_core(
    void *src,
    size_t level,
    void *_root_table,
    struct compress_profile _profile,
    struct ring ring
    ){

    size_t step;
    void *real_root_table;

    step = 1u << level;
    real_root_table = _root_table + ((1u << level) - 1) * ring.sizeZ;

    for(size_t i = 0; i < step; i++){
        for(size_t j = 0; j < _profile.array_n; j += 2 * step){
            CT_butterfly(src + (i + j) * ring.sizeZ, 0, step, real_root_table, ring);
        }
        real_root_table += ring.sizeZ;
    }

}

// ================================
void GS_iNTT_core(
    void *src,
    size_t level,
    void *_root_table,
    struct compress_profile _profile,
    struct ring ring
    ){

    size_t step;
    void *real_root_table;

    step = (_profile.array_n) >> (level + 1);
    real_root_table = _root_table + ((1U << level) - 1) * ring.sizeZ;

    for(size_t i = 0; i < _profile.array_n; i += 2 * step){
        for(size_t j = 0; j < step; j++){
            GS_butterfly(src + (i + j) * ring.sizeZ, 0, step, real_root_table, ring);
        }
        real_root_table += ring.sizeZ;
    }

}

// ================================
void CT_NTT(
    void *src,
    void *_root_table,
    struct compress_profile _profile,
    struct ring ring
    ){

    for(size_t i = 0; i < _profile.log_ntt_n; i++){
        CT_NTT_core(src, i, _root_table, _profile, ring);
    }

}

// ================================
void CT_iNTT(
    void *src,
    void *_root_table,
    struct compress_profile _profile,
    struct ring ring
    ){

    for(size_t i = 0; i < _profile.log_ntt_n; i++){
        CT_iNTT_core(src, i, _root_table, _profile, ring);
    }

}

// ================================
void GS_iNTT(
    void *src,
    void *_root_table,
    struct compress_profile _profile,
    struct ring ring
    ){

    for(ssize_t i = _profile.log_ntt_n - 1; i >= 0; i--){
        GS_iNTT_core(src, i, _root_table, _profile, ring);
    }

}

// ================================

// ================================
// Multi-layer Cooley-Tukey butterfly for the forward transformation.
void m_layer_CT_butterfly(
    void *src,
    size_t layers, size_t step,
    void *_root_table,
    struct ring ring
    ){

    size_t real_count, real_step, twiddle_count, jump;
    void *real_root_table;

    for(size_t i = 0; i < layers; i++){

        twiddle_count = 1u << i;

        real_count = 1u << (layers - 1 - i);

        jump = step << (layers - i);

        real_root_table = _root_table + ((1u << i) - 1) * ring.sizeZ;

        real_step = step << (layers - 1 - i);

        for(size_t k = 0; k < real_count; k++){
            for(size_t j = 0; j < twiddle_count; j++){
                CT_butterfly(
                    src + (j * jump + k * step) * ring.sizeZ,
                    0, real_step,
                    real_root_table + j * ring.sizeZ,
                    ring
                    );
            }
        }

    }

}

// ================================
// Multi-layer Cooley-Tukey butterfly for the inverse transformation.
void m_layer_CT_ibutterfly(
    void *src,
    size_t layers, size_t step,
    void *_root_table,
    struct ring ring
    ){

    size_t real_count, real_step, twiddle_count, jump;
    void *real_root_table;

    for(size_t i = 0; i < layers; i++){

        twiddle_count = 1u << i;

        real_count = 1u << (layers - 1 - i);

        jump = step << (i + 1);

        real_root_table = _root_table + ((1u << i) - 1) * ring.sizeZ;

        real_step = step << i;

        for(size_t k = 0; k < real_count; k++){
            for(size_t j = 0; j < twiddle_count; j++){
                CT_butterfly(
                    src + (j * step + k * jump) * ring.sizeZ,
                    0, real_step,
                    real_root_table + j * ring.sizeZ,
                    ring
                    );
            }
        }

    }

}

// ================================
// Multi-layer Gentleman-Sande butterfly for the inverse transformation.
void m_layer_GS_ibutterfly(
    void *src,
    size_t layers, size_t step,
    void *_root_table,
    struct ring ring
    ){

    size_t real_count, real_step, twiddle_count, jump;
    void *real_root_table;

    for(ssize_t i = layers - 1; i >= 0; i--){

        twiddle_count = 1u << i;

        real_count = 1u << (layers - 1 - i);

        jump = step << (layers - i);

        real_root_table = _root_table + ((1u << i) - 1) * ring.sizeZ;

        real_step = step << (layers - 1 - i);

        for(size_t k = 0; k < real_count; k++){
            for(size_t j = 0; j < twiddle_count; j++){
                GS_butterfly(
                    src + (j * jump + k * step) * ring.sizeZ,
                    0, real_step,
                    real_root_table + j * ring.sizeZ,
                    ring
                    );
            }
        }

    }

}

// ================================
// NTT with Cooley-Tukey butterfly.
// We must use m_layer_CT_butterfly here.
void compressed_CT_NTT(
    void *src,
    size_t start_level, size_t end_level,
    void *_root_table,
    struct compress_profile _profile,
    struct ring ring
    ){

    size_t step, offset;
    size_t real_start_level, real_end_level;
    void *real_root_table;
    size_t *level_indx;

    real_start_level = 0;
    for(size_t i = 0; i < start_level; i++){
        real_start_level += (_profile.merged_layers)[i];
    }

    real_end_level = real_start_level;
    for(size_t i = start_level; i < end_level; i++){
        real_end_level += (_profile.merged_layers)[i];
    }

    level_indx = (_profile.merged_layers) + start_level;

    for(size_t level = real_start_level; level <= real_end_level; level += *(level_indx++)){

        step = _profile.array_n >> (level + (*level_indx));

        offset = 0;

        real_root_table = _root_table + ((1u << level) - 1) * ring.sizeZ;

        for(size_t count = 0; count < (1u << level); count++){

            for(size_t i = 0; i < step; i++){
                m_layer_CT_butterfly(
                    src + (offset + i) * ring.sizeZ,
                    *level_indx, step,
                    real_root_table,
                    ring
                    );
            }

            offset += _profile.array_n >> level;

            real_root_table += ((1u << (*level_indx)) - 1) * ring.sizeZ;

        }

    }

}

// ================================
// iNTT with Cooley-Tukey butterfly.
// We must use m_layer_CT_ibutterfly here.
void compressed_CT_iNTT(
    void *src,
    size_t start_level, size_t end_level,
    void *_root_table,
    struct compress_profile _profile,
    struct ring ring
    ){

    size_t step, offset;
    size_t real_start_level, real_end_level;
    void *real_root_table;
    size_t *level_indx;

    real_start_level = 0;
    for(size_t i = 0; i < start_level; i++){
        real_start_level += (_profile.merged_layers)[i];
    }

    real_end_level = real_start_level;
    for(size_t i = start_level; i < end_level; i++){
        real_end_level += (_profile.merged_layers)[i];
    }

    real_root_table = _root_table + ((1u << real_start_level) - 1) * ring.sizeZ;

    level_indx = (_profile.merged_layers) + start_level;

    for(size_t level = real_start_level; level <= real_end_level; level += *(level_indx++)){

        step = (_profile.array_n >> _profile.log_ntt_n) << level;

        for(size_t count = 0; count < (1u << level); count++){

            offset = count * (_profile.array_n >> _profile.log_ntt_n);

            for(size_t i = 0; i < (_profile.ntt_n >> ((*level_indx) + level)); i++){

                for(size_t j = 0; j < (_profile.array_n >> _profile.log_ntt_n); j++){
                    m_layer_CT_ibutterfly(
                        src + (offset + j) * ring.sizeZ,
                        *level_indx, step,
                        real_root_table,
                        ring
                        );
                }

                offset += (_profile.array_n >> _profile.log_ntt_n) << ((*level_indx) + level);

            }

            real_root_table += ((1u << (*level_indx)) - 1) * ring.sizeZ;

        }

    }

}

// ================================
// iNTT with Gentleman-Sande butterfly.
// We must use m_layer_GS_ibutterfly here.
void compressed_GS_iNTT(
    void *src,
    size_t start_level, size_t end_level,
    void *_root_table,
    struct compress_profile _profile,
    struct ring ring
    ){

    size_t step, offset;
    ssize_t real_start_level, real_end_level;
    void *real_root_table;
    size_t *level_indx;

    real_start_level = 0;
    for(size_t i = 0; i < start_level; i++){
        real_start_level += (_profile.merged_layers)[i];
    }

    real_end_level = real_start_level;
    for(size_t i = start_level; i < end_level; i++){
        real_end_level += (_profile.merged_layers)[i];
    }

    level_indx = (_profile.merged_layers) + end_level;

    for(ssize_t level = real_end_level; level >= real_start_level; level -= *(level_indx-- -1) ){

        step = _profile.array_n >> (level + (*(level_indx)));

        offset = 0;

        real_root_table = _root_table + ((1u << level) - 1) * ring.sizeZ;

        for(size_t count = 0; count < (1u << level); count++){

            for(size_t i = 0; i < step; i++){
                m_layer_GS_ibutterfly(
                    src + (offset + i) * ring.sizeZ,
                    *level_indx, step,
                    real_root_table,
                    ring
                    );
            }

            offset += _profile.array_n >> level;

            real_root_table += ((1u << (*level_indx)) - 1) * ring.sizeZ;

        }

    }

}












