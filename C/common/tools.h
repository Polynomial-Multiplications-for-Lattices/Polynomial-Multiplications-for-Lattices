#ifndef TOOLS_H
#define TOOLS_H

#include <stdint.h>
#include <stddef.h>

// ================================
// Structure compress_profile

// Explain what this structure is about...
// This structure contains the information of a layer-merging strategy. The structure controls
// the computation flow of functions in gen_table.c and ntt_c.c.
// - array_n
//      - This refers to the length of the array.
// - ntt_n
//      - This refers to the length of the NTT. ntt_n must be a power-of-two factor of array_n.
// - log_ntt_n
//      - log_ntt_n must be the base-2 logarithm of ntt_n. It determines the total number of
//      - layers prior to the compression.
// - compressed_layers
//      -
// - merged_layers[16]
//      -
struct compress_profile {
    size_t array_n;
    size_t ntt_n;
    size_t log_ntt_n;
    size_t compressed_layers;
    size_t merged_layers[16];
};

// ================================
// Structure ring

// This structure defines a ring with each element occupying sizeZ bytes in memory.
// The rest are all functions defining a ring along with several auxiliary helper functions.
// Notice that all elements of the structure are user-defined. We give recommendations on how
// each functions should be implemented, but you can do whatever you want as long as you know
// what you are doing.
// - sizeZ
//      - sizeZ defines the number of bytes occupied by an element.
// - memberZ
//      - This pointer points to the membership function mapping an element to its representative
//        in the ring.
//      - Ideally, the image should be a canonical representation of the ring.
// - addZ
//      - This pointer points to a function implementing an addition function for the ring.
//      - In the simplest form, the result should be representative of the ring ring.
//      - However, you can also implement it in an unreduced fashion as long as the result is
//        still in sizeZ bytes.
// - subZ
//      - This pointer points to a function implementing a subtraction function for the ring.
//      - See addZ for further information.
// - mulZ
//      - This pointer points to a function implementing a multiplication function for the ring.
//      - See addZ for further information.
// - expZ
//      - This pointer points to a function implementing an exponentiation function as repeated
//        multiplication for the ring.
//      - This function is only used for generating the tables of twiddle factors. You can skip
//        it if there is no need to generate the tables of twiddle factors.
//      - See addZ for further information.
struct ring {
    // sizeZ is refers to the size in bytes of an element in the ring.
    size_t sizeZ;
    // memberZ maps *src to its representative in the ring and store the result in des.
    void (*memberZ)(void *des, const void *src);
    // addZ adds up *src1 and *src2, and stores the result in des.
    void (*addZ)(void *des, const void *src1, const void *src2);
    // subZ subtract *src2 from *src1, and stores the result in des.
    void (*subZ)(void *des, const void *src1, const void *src2);
    // mulZ multiplies *src1 and *src2, and stores the result in des.
    void (*mulZ)(void *des, const void *src1, const void *src2);
    // expZ computes (*src)^e, and stores the result in des.
    void (*expZ)(void *des, const void *src, size_t e);
};

// ================================
// We also provide several commonly used functions.
// By default, Z_Q is defined as the set of integers in [-Q / 2, Q / 2).

// ================================
// Candidates for memberZ.

// This function assumes sizeZ = 2 and maps the element to the representative
// in the ring Z_{*mod} with signed representation.
void cmod_int16(void *des, const void *src, const void *mod);
// This function assumes sizeZ = 4 and maps the element to the representative
// in the ring Z_{*mod} with signed representation.
void cmod_int32(void *des, const void *src, const void *mod);
// This function assumes sizeZ = 8 and maps the element to the representative
// in the ring Z_{*mod} with signed representation.
void cmod_int64(void *des, const void *src, const void *mod);

// ================================
// Candidates for addZ.

// This function assumes sizeZ = 2. It sums *src1 and *src2 and maps the result to the
// representative in the ring Z_{*mod} with signed representation.
void addmod_int16(void *des, const void *src1, const void *src2, const void *mod);
// This function assumes sizeZ = 4. It sums *src1 and *src2 and maps the result to the
// representative in the ring Z_{*mod} with signed representation.
void addmod_int32(void *des, const void *src1, const void *src2, const void *mod);

// ================================
// Candidates for subZ.

// This function assumes sizeZ = 2. It subtracts *src2 from *src1 and maps the result to the
// representative in the ring Z_{*mod} with signed representation.
void submod_int16(void *des, const void *src1, const void *src2, const void *mod);
// This function assumes sizeZ = 4. It subtracts *src2 from *src1 and maps the result to the
// representative in the ring Z_{*mod} with signed representation.
void submod_int32(void *des, const void *src1, const void *src2, const void *mod);

// ================================
// Candidates for mulZ.

// This function assumes sizeZ = 2. It multiplies *src1 and *src2 and maps the result to the
// representative in the ring Z_{*mod} with signed representation.
void mulmod_int16(void *des, const void *src1, const void *src2, const void *mod);
// This function assumes sizeZ = 4. It multiplies *src1 and *src2 and maps the result to the
// representative in the ring Z_{*mod} with signed representation.
void mulmod_int32(void *des, const void *src1, const void *src2, const void *mod);

// ================================
// Candidates for expZ.

// This function assumes sizeZ = 2. It computes (*src)^e and maps the result to the
// representative in the ring Z_{*mod} with signed representation.
void expmod_int16(void *des, const void *src, size_t e, const void *mod);
// This function assumes sizeZ = 4. It computes (*src)^e and maps the result to the
// representative in the ring Z_{*mod} with signed representation.
void expmod_int32(void *des, const void *src, size_t e, const void *mod);

// ================================
// In-place bit-reversal.

// This function is used for generating tables of twiddle factors.
// See gen_table.c for example usages.
void bitreverse(void *src, size_t len, size_t size);

#endif

