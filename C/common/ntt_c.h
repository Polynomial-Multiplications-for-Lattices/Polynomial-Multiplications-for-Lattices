#ifndef NTT_C_H
#define NTT_C_H

#include "tools.h"

// ================================
// Butterfly operations.

// Cooley-Tukey butterfly.
// This function computes (src[indx_a] + (*twiddle) src[indx_b], src[indx_a] - (*twiddle) src[indx_b])
// and stores the result to src + indx_a and src + indx_b.
void CT_butterfly(
    void *src,
    size_t indx_a, size_t indx_b,
    void *twiddle,
    struct ring ring
    );

// Gentleman-Sande butterfly.
// This function computes (src[indx_a] + src[indx_b], (src[indx_a] - src[indx_b]) (*twiddle) )
// and stores the result to src + indx_a and src + indx_b.
void GS_butterfly(
    void *src,
    size_t indx_a, size_t indx_b,
    void *twiddle,
    struct ring ring
    );

// ================================
// Core operations computing one layer of butterflies.

// This function computes the level-th layer of Cooley--Tukey butterflies in the NTT.
void CT_NTT_core(
    void *src,
    size_t level,
    void *_root_table,
    struct compress_profile _profile,
    struct ring ring
    );

// This function computes the level-th layer of Cooley--Tukey butterflies in the iNTT.
void CT_iNTT_core(
    void *src,
    size_t level,
    void *_root_table,
    struct compress_profile _profile,
    struct ring ring
    );

// This function computes the level-th layer of Gentleman--Sande butterflies in the iNTT.
void GS_iNTT_core(
    void *src,
    size_t level,
    void *_root_table,
    struct compress_profile _profile,
    struct ring ring
    );

// ================================
// NTT computations without layer-merging.

void CT_NTT(
    void *src,
    void *_root_table,
    struct compress_profile _profile,
    struct ring ring
    );

void CT_iNTT(
    void *src,
    void *_root_table,
    struct compress_profile _profile,
    struct ring ring
    );

void GS_iNTT(
    void *src,
    void *_root_table,
    struct compress_profile _profile,
    struct ring ring
    );

// ================================
// Multi-layer butterly.

// Multi-layer Cooley-Tukey butterfly for the forward transformation.
void m_layer_CT_butterfly(
    void *src,
    size_t layers, size_t step,
    void *_root_table,
    struct ring ring
    );

// Multi-layer Cooley-Tukey butterfly for the inverse transformation.
void m_layer_CT_ibutterfly(
    void *src,
    size_t layers, size_t step,
    void *_root_table,
    struct ring ring
    );

// Multi-layer Gentleman-Sande butterfly for the inverse transformation.
void m_layer_GS_ibutterfly(
    void *src,
    size_t layers, size_t step,
    void *_root_table,
    struct ring ring
    );

// ================================
// NTT with custom layer-merging.

// NTT with Cooley-Tukey butterfly.
// We must use m_layer_CT_butterfly here.
void compressed_CT_NTT(
    void *src,
    size_t start_level, size_t end_level,
    void *_root_table,
    struct compress_profile _profile,
    struct ring ring
    );

// iNTT with Cooley-Tukey butterfly.
// We must use m_layer_CT_ibutterfly here.
void compressed_CT_iNTT(
    void *src,
    size_t start_level, size_t end_level,
    void *_root_table,
    struct compress_profile _profile,
    struct ring ring
    );

// iNTT with Gentleman-Sande butterfly.
// We must use m_layer_GS_ibutterfly here.
void compressed_GS_iNTT(
    void *src,
    size_t start_level, size_t end_level,
    void *_root_table,
    struct compress_profile _profile,
    struct ring ring
    );

#endif

