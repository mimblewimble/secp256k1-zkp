/**********************************************************************
 * Copyright (c) 2013, 2014, 2017 Pieter Wuille, Andrew Poelstra,     *
 *   Peter Dettmann                                                   *
 * Distributed under the MIT software license, see the accompanying   *
 * file COPYING or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#ifndef _SECP256K1_ECMULT_IMPL_H_
#define _SECP256K1_ECMULT_IMPL_H_

#include <string.h>

#include "group.h"
#include "scalar.h"
#include "ecmult.h"

#if defined(EXHAUSTIVE_TEST_ORDER)
/* We need to lower these values for exhaustive tests because
 * the tables cannot have infinities in them (this breaks the
 * affine-isomorphism stuff which tracks z-ratios) */
#  if EXHAUSTIVE_TEST_ORDER > 128
#    define WINDOW_A 5
#    define WINDOW_G 8
#  elif EXHAUSTIVE_TEST_ORDER > 8
#    define WINDOW_A 4
#    define WINDOW_G 4
#  else
#    define WINDOW_A 2
#    define WINDOW_G 2
#  endif
#else
/* optimal for 128-bit and 256-bit exponents. */
#define WINDOW_A 5
/** larger numbers may result in slightly better performance, at the cost of
    exponentially larger precomputed tables. */
#ifdef USE_ENDOMORPHISM
/** Two tables for window size 15: 1.375 MiB. */
#define WINDOW_G 15
#else
/** One table for window size 16: 1.375 MiB. */
#define WINDOW_G 16
#endif
#endif

/** The number of entries a table with precomputed multiples needs to have. */
#define ECMULT_TABLE_SIZE(w) (1 << ((w)-2))

/** Fill a table 'prej' with precomputed odd multiples of a. Prej will contain
 *  the values [1*a,3*a,...,(2*n-1)*a], so it space for n values. zr[0] will
 *  contain prej[0].z / a.z. The other zr[i] values = prej[i].z / prej[i-1].z.
 *  Prej's Z values are undefined, except for the last value.
 */
static void secp256k1_ecmult_odd_multiples_table(int n, secp256k1_gej *prej, secp256k1_fe *zr, const secp256k1_gej *a) {
    secp256k1_gej d;
    secp256k1_ge a_ge, d_ge;
    int i;

    VERIFY_CHECK(!a->infinity);

    secp256k1_gej_double_var(&d, a, NULL);

    /*
     * Perform the additions on an isomorphism where 'd' is affine: drop the z coordinate
     * of 'd', and scale the 1P starting value's x/y coordinates without changing its z.
     */
    d_ge.x = d.x;
    d_ge.y = d.y;
    d_ge.infinity = 0;

    secp256k1_ge_set_gej_zinv(&a_ge, a, &d.z);
    prej[0].x = a_ge.x;
    prej[0].y = a_ge.y;
    prej[0].z = a->z;
    prej[0].infinity = 0;

    zr[0] = d.z;
    for (i = 1; i < n; i++) {
        secp256k1_gej_add_ge_var(&prej[i], &prej[i-1], &d_ge, &zr[i]);
    }

    /*
     * Each point in 'prej' has a z coordinate too small by a factor of 'd.z'. Only
     * the final point's z coordinate is actually used though, so just update that.
     */
    secp256k1_fe_mul(&prej[n-1].z, &prej[n-1].z, &d.z);
}

/** Fill a table 'pre' with precomputed odd multiples of a.
 *
 *  There are two versions of this function:
 *  - secp256k1_ecmult_odd_multiples_table_globalz_windowa which brings its
 *    resulting point set to a single constant Z denominator, stores the X and Y
 *    coordinates as ge_storage points in pre, and stores the global Z in rz.
 *    It only operates on tables sized for WINDOW_A wnaf multiples.
 *  - secp256k1_ecmult_odd_multiples_table_storage_var, which converts its
 *    resulting point set to actually affine points, and stores those in pre.
 *    It operates on tables of any size, but uses heap-allocated temporaries.
 *
 *  To compute a*P + b*G, we compute a table for P using the first function,
 *  and for G using the second (which requires an inverse, but it only needs to
 *  happen once).
 */
static void secp256k1_ecmult_odd_multiples_table_globalz_windowa(secp256k1_ge *pre, secp256k1_fe *globalz, const secp256k1_gej *a) {
    secp256k1_gej prej[ECMULT_TABLE_SIZE(WINDOW_A)];
    secp256k1_fe zr[ECMULT_TABLE_SIZE(WINDOW_A)];

    /* Compute the odd multiples in Jacobian form. */
    secp256k1_ecmult_odd_multiples_table(ECMULT_TABLE_SIZE(WINDOW_A), prej, zr, a);
    /* Bring them to the same Z denominator. */
    secp256k1_ge_globalz_set_table_gej(ECMULT_TABLE_SIZE(WINDOW_A), pre, globalz, prej, zr);
}

static void secp256k1_ecmult_odd_multiples_table_storage_var(int n, secp256k1_ge_storage *pre, const secp256k1_gej *a, const secp256k1_callback *cb) {
    secp256k1_gej *prej = (secp256k1_gej*)checked_malloc(cb, sizeof(secp256k1_gej) * n);
    secp256k1_ge *prea = (secp256k1_ge*)checked_malloc(cb, sizeof(secp256k1_ge) * n);
    secp256k1_fe *zr = (secp256k1_fe*)checked_malloc(cb, sizeof(secp256k1_fe) * n);
    int i;

    /* Compute the odd multiples in Jacobian form. */
    secp256k1_ecmult_odd_multiples_table(n, prej, zr, a);
    /* Convert them in batch to affine coordinates. */
    secp256k1_ge_set_table_gej_var(prea, prej, zr, n);
    /* Convert them to compact storage form. */
    for (i = 0; i < n; i++) {
        secp256k1_ge_to_storage(&pre[i], &prea[i]);
    }

    free(prea);
    free(prej);
    free(zr);
}

/** The following two macro retrieves a particular odd multiple from a table
 *  of precomputed multiples. */
#define ECMULT_TABLE_GET_GE(r,pre,n,w) do { \
    VERIFY_CHECK(((n) & 1) == 1); \
    VERIFY_CHECK((n) >= -((1 << ((w)-1)) - 1)); \
    VERIFY_CHECK((n) <=  ((1 << ((w)-1)) - 1)); \
    if ((n) > 0) { \
        *(r) = (pre)[((n)-1)/2]; \
    } else { \
        secp256k1_ge_neg((r), &(pre)[(-(n)-1)/2]); \
    } \
} while(0)

#define ECMULT_TABLE_GET_GE_STORAGE(r,pre,n,w) do { \
    VERIFY_CHECK(((n) & 1) == 1); \
    VERIFY_CHECK((n) >= -((1 << ((w)-1)) - 1)); \
    VERIFY_CHECK((n) <=  ((1 << ((w)-1)) - 1)); \
    if ((n) > 0) { \
        secp256k1_ge_from_storage((r), &(pre)[((n)-1)/2]); \
    } else { \
        secp256k1_ge_from_storage((r), &(pre)[(-(n)-1)/2]); \
        secp256k1_ge_neg((r), (r)); \
    } \
} while(0)

static void secp256k1_ecmult_context_init(secp256k1_ecmult_context *ctx) {
    ctx->pre_g = NULL;
#ifdef USE_ENDOMORPHISM
    ctx->pre_g_128 = NULL;
#endif
}

static void secp256k1_ecmult_context_build(secp256k1_ecmult_context *ctx, const secp256k1_callback *cb) {
    secp256k1_gej gj;

    if (ctx->pre_g != NULL) {
        return;
    }

    /* get the generator */
    secp256k1_gej_set_ge(&gj, &secp256k1_ge_const_g);

    ctx->pre_g = (secp256k1_ge_storage (*)[])checked_malloc(cb, sizeof((*ctx->pre_g)[0]) * ECMULT_TABLE_SIZE(WINDOW_G));

    /* precompute the tables with odd multiples */
    secp256k1_ecmult_odd_multiples_table_storage_var(ECMULT_TABLE_SIZE(WINDOW_G), *ctx->pre_g, &gj, cb);

#ifdef USE_ENDOMORPHISM
    {
        secp256k1_gej g_128j;
        int i;

        ctx->pre_g_128 = (secp256k1_ge_storage (*)[])checked_malloc(cb, sizeof((*ctx->pre_g_128)[0]) * ECMULT_TABLE_SIZE(WINDOW_G));

        /* calculate 2^128*generator */
        g_128j = gj;
        for (i = 0; i < 128; i++) {
            secp256k1_gej_double_var(&g_128j, &g_128j, NULL);
        }
        secp256k1_ecmult_odd_multiples_table_storage_var(ECMULT_TABLE_SIZE(WINDOW_G), *ctx->pre_g_128, &g_128j, cb);
    }
#endif
}

static void secp256k1_ecmult_context_clone(secp256k1_ecmult_context *dst,
                                           const secp256k1_ecmult_context *src, const secp256k1_callback *cb) {
    if (src->pre_g == NULL) {
        dst->pre_g = NULL;
    } else {
        size_t size = sizeof((*dst->pre_g)[0]) * ECMULT_TABLE_SIZE(WINDOW_G);
        dst->pre_g = (secp256k1_ge_storage (*)[])checked_malloc(cb, size);
        memcpy(dst->pre_g, src->pre_g, size);
    }
#ifdef USE_ENDOMORPHISM
    if (src->pre_g_128 == NULL) {
        dst->pre_g_128 = NULL;
    } else {
        size_t size = sizeof((*dst->pre_g_128)[0]) * ECMULT_TABLE_SIZE(WINDOW_G);
        dst->pre_g_128 = (secp256k1_ge_storage (*)[])checked_malloc(cb, size);
        memcpy(dst->pre_g_128, src->pre_g_128, size);
    }
#endif
}

static int secp256k1_ecmult_context_is_built(const secp256k1_ecmult_context *ctx) {
    return ctx->pre_g != NULL;
}

static void secp256k1_ecmult_context_clear(secp256k1_ecmult_context *ctx) {
    free(ctx->pre_g);
#ifdef USE_ENDOMORPHISM
    free(ctx->pre_g_128);
#endif
    secp256k1_ecmult_context_init(ctx);
}

/** Convert a number to WNAF notation. The number becomes represented by sum(2^i * wnaf[i], i=0..bits),
 *  with the following guarantees:
 *  - each wnaf[i] is either 0, or an odd integer between -(1<<(w-1) - 1) and (1<<(w-1) - 1)
 *  - two non-zero entries in wnaf are separated by at least w-1 zeroes.
 *  - the number of set values in wnaf is returned. This number is at most 256, and at most one more
 *    than the number of bits in the (absolute value) of the input.
 */
static int secp256k1_ecmult_wnaf(int *wnaf, int len, const secp256k1_scalar *a, int w) {
    secp256k1_scalar s = *a;
    int last_set_bit = -1;
    int bit = 0;
    int sign = 1;
    int carry = 0;

    VERIFY_CHECK(wnaf != NULL);
    VERIFY_CHECK(0 <= len && len <= 256);
    VERIFY_CHECK(a != NULL);
    VERIFY_CHECK(2 <= w && w <= 31);

    memset(wnaf, 0, len * sizeof(wnaf[0]));

    if (secp256k1_scalar_get_bits(&s, 255, 1)) {
        secp256k1_scalar_negate(&s, &s);
        sign = -1;
    }

    while (bit < len) {
        int now;
        int word;
        if (secp256k1_scalar_get_bits(&s, bit, 1) == (unsigned int)carry) {
            bit++;
            continue;
        }

        now = w;
        if (now > len - bit) {
            now = len - bit;
        }

        word = secp256k1_scalar_get_bits_var(&s, bit, now) + carry;

        carry = (word >> (w-1)) & 1;
        word -= carry << w;

        wnaf[bit] = sign * word;
        last_set_bit = bit;

        bit += now;
    }
#ifdef VERIFY
    CHECK(carry == 0);
    while (bit < 256) {
        CHECK(secp256k1_scalar_get_bits(&s, bit++, 1) == 0);
    } 
#endif
    return last_set_bit + 1;
}

static void secp256k1_ecmult(const secp256k1_ecmult_context *ctx, secp256k1_gej *r, const secp256k1_gej *a, const secp256k1_scalar *na, const secp256k1_scalar *ng) {
    secp256k1_ge pre_a[ECMULT_TABLE_SIZE(WINDOW_A)];
    secp256k1_ge tmpa;
    secp256k1_fe Z;
#ifdef USE_ENDOMORPHISM
    secp256k1_ge pre_a_lam[ECMULT_TABLE_SIZE(WINDOW_A)];
    secp256k1_scalar na_1, na_lam;
    /* Splitted G factors. */
    secp256k1_scalar ng_1, ng_128;
    int wnaf_na_1[130];
    int wnaf_na_lam[130];
    int bits_na_1;
    int bits_na_lam;
    int wnaf_ng_1[129];
    int bits_ng_1;
    int wnaf_ng_128[129];
    int bits_ng_128;
#else
    int wnaf_na[256];
    int bits_na;
    int wnaf_ng[256];
    int bits_ng;
#endif
    int i;
    int bits;

#ifdef USE_ENDOMORPHISM
    /* split na into na_1 and na_lam (where na = na_1 + na_lam*lambda, and na_1 and na_lam are ~128 bit) */
    secp256k1_scalar_split_lambda(&na_1, &na_lam, na);

    /* build wnaf representation for na_1 and na_lam. */
    bits_na_1   = secp256k1_ecmult_wnaf(wnaf_na_1,   130, &na_1,   WINDOW_A);
    bits_na_lam = secp256k1_ecmult_wnaf(wnaf_na_lam, 130, &na_lam, WINDOW_A);
    VERIFY_CHECK(bits_na_1 <= 130);
    VERIFY_CHECK(bits_na_lam <= 130);
    bits = bits_na_1;
    if (bits_na_lam > bits) {
        bits = bits_na_lam;
    }
#else
    /* build wnaf representation for na. */
    bits_na     = secp256k1_ecmult_wnaf(wnaf_na,     256, na,      WINDOW_A);
    bits = bits_na;
#endif

    /* Calculate odd multiples of a.
     * All multiples are brought to the same Z 'denominator', which is stored
     * in Z. Due to secp256k1' isomorphism we can do all operations pretending
     * that the Z coordinate was 1, use affine addition formulae, and correct
     * the Z coordinate of the result once at the end.
     * The exception is the precomputed G table points, which are actually
     * affine. Compared to the base used for other points, they have a Z ratio
     * of 1/Z, so we can use secp256k1_gej_add_zinv_var, which uses the same
     * isomorphism to efficiently add with a known Z inverse.
     */
    secp256k1_ecmult_odd_multiples_table_globalz_windowa(pre_a, &Z, a);

#ifdef USE_ENDOMORPHISM
    for (i = 0; i < ECMULT_TABLE_SIZE(WINDOW_A); i++) {
        secp256k1_ge_mul_lambda(&pre_a_lam[i], &pre_a[i]);
    }

    /* split ng into ng_1 and ng_128 (where gn = gn_1 + gn_128*2^128, and gn_1 and gn_128 are ~128 bit) */
    secp256k1_scalar_split_128(&ng_1, &ng_128, ng);

    /* Build wnaf representation for ng_1 and ng_128 */
    bits_ng_1   = secp256k1_ecmult_wnaf(wnaf_ng_1,   129, &ng_1,   WINDOW_G);
    bits_ng_128 = secp256k1_ecmult_wnaf(wnaf_ng_128, 129, &ng_128, WINDOW_G);
    if (bits_ng_1 > bits) {
        bits = bits_ng_1;
    }
    if (bits_ng_128 > bits) {
        bits = bits_ng_128;
    }
#else
    bits_ng     = secp256k1_ecmult_wnaf(wnaf_ng,     256, ng,      WINDOW_G);
    if (bits_ng > bits) {
        bits = bits_ng;
    }
#endif

    secp256k1_gej_set_infinity(r);

    for (i = bits - 1; i >= 0; i--) {
        int n;
        secp256k1_gej_double_var(r, r, NULL);
#ifdef USE_ENDOMORPHISM
        if (i < bits_na_1 && (n = wnaf_na_1[i])) {
            ECMULT_TABLE_GET_GE(&tmpa, pre_a, n, WINDOW_A);
            secp256k1_gej_add_ge_var(r, r, &tmpa, NULL);
        }
        if (i < bits_na_lam && (n = wnaf_na_lam[i])) {
            ECMULT_TABLE_GET_GE(&tmpa, pre_a_lam, n, WINDOW_A);
            secp256k1_gej_add_ge_var(r, r, &tmpa, NULL);
        }
        if (i < bits_ng_1 && (n = wnaf_ng_1[i])) {
            ECMULT_TABLE_GET_GE_STORAGE(&tmpa, *ctx->pre_g, n, WINDOW_G);
            secp256k1_gej_add_zinv_var(r, r, &tmpa, &Z);
        }
        if (i < bits_ng_128 && (n = wnaf_ng_128[i])) {
            ECMULT_TABLE_GET_GE_STORAGE(&tmpa, *ctx->pre_g_128, n, WINDOW_G);
            secp256k1_gej_add_zinv_var(r, r, &tmpa, &Z);
        }
#else
        if (i < bits_na && (n = wnaf_na[i])) {
            ECMULT_TABLE_GET_GE(&tmpa, pre_a, n, WINDOW_A);
            secp256k1_gej_add_ge_var(r, r, &tmpa, NULL);
        }
        if (i < bits_ng && (n = wnaf_ng[i])) {
            ECMULT_TABLE_GET_GE_STORAGE(&tmpa, *ctx->pre_g, n, WINDOW_G);
            secp256k1_gej_add_zinv_var(r, r, &tmpa, &Z);
        }
#endif
    }

    if (!r->infinity) {
        secp256k1_fe_mul(&r->z, &r->z, &Z);
    }
}

/* begin ecmult_multi */

typedef struct {
    uint32_t *tree;
    const secp256k1_scalar *scalars;
    size_t size;
} secp256k1_scalar_heap;

static void secp256k1_sift_down(secp256k1_scalar_heap *heap, size_t node, uint32_t index) {
    uint32_t child_index, other_index;
    size_t child, other, half_size = heap->size >> 1;
    const secp256k1_scalar *sc = heap->scalars;

    while (node < half_size) {
        /* Initially assume the left child is the larger child */
        child = (node << 1) + 1;
        child_index = heap->tree[child];

        /* If there is a right child, check whether it's larger than the left */
        other = child + 1;
        if (other < heap->size) {
            other_index = heap->tree[other];
            if (secp256k1_scalar_cmp_var(&sc[other_index], &sc[child_index]) > 0) {
                child = other;
                child_index = other_index;
            }
        }

        /* If the current node is larger than its largest child, stop at this level */
        if (secp256k1_scalar_cmp_var(&sc[index], &sc[child_index]) > 0) {
            break;
        }

        /* Move the larger child up, and recurse from its previous position */
        heap->tree[node] = child_index;
        node = child;
    }

    heap->tree[node] = index;
}

static void secp256k1_sift_up(secp256k1_scalar_heap *heap, size_t node, uint32_t index) {
    size_t parent;
    uint32_t parent_index;
    const secp256k1_scalar *sc = heap->scalars;

    while (node > 0) {
        parent = (node - 1) >> 1;
        parent_index = heap->tree[parent];

        /* If the current node is not larger than its parent, stop at this level */
        if (secp256k1_scalar_cmp_var(&sc[index], &sc[parent_index]) <= 0) {
            break;
        }

        /* Move the parent down, and recurse from its previous position */
        heap->tree[node] = parent_index;
        node = parent;
    }

    heap->tree[node] = index;
}

static void secp256k1_sift_floyd(secp256k1_scalar_heap *heap, size_t node, uint32_t index) {
    uint32_t child_index, other_index;
    size_t child, other, half_size = heap->size >> 1;
    const secp256k1_scalar *sc = heap->scalars;

    while (node < half_size) {
        /* Initially assume the left child is the larger child */
        child = (node << 1) + 1;
        child_index = heap->tree[child];

        /* If there is a right child, check whether it's larger than the left */
        other = child + 1;
        if (other < heap->size) {
            other_index = heap->tree[other];
            if (secp256k1_scalar_cmp_var(&sc[other_index], &sc[child_index]) > 0) {
                child = other;
                child_index = other_index;
            }
        }

        /* Move the larger child up, and recurse from its previous position */
        heap->tree[node] = child_index;
        node = child;
    }

    secp256k1_sift_up(heap, node, index);
}

SECP256K1_INLINE static void secp256k1_heapify(secp256k1_scalar_heap *heap) {
    size_t root = heap->size >> 1;;
    while (root-- > 0) {
        secp256k1_sift_down(heap, root, heap->tree[root]);
    }
}

static void secp256k1_heap_initialize(secp256k1_scalar_heap *heap, uint32_t *tree, const secp256k1_scalar *scalars, const secp256k1_gej *pt, size_t n) {
    size_t i, size = 0;

    heap->tree = tree;
    for (i = 0; i < n; ++i) {
        if (!secp256k1_scalar_is_zero(&scalars[i]) && !secp256k1_gej_is_infinity(&pt[i])) {
            heap->tree[size++] = i;
        }
    }

    heap->scalars = scalars;
    heap->size = size;

    secp256k1_heapify(heap);
}

SECP256K1_INLINE static uint32_t secp256k1_replace(secp256k1_scalar_heap *heap, uint32_t index) {
    uint32_t result = heap->tree[0];
    VERIFY_CHECK(heap->size > 0);
    secp256k1_sift_floyd(heap, 0, index);
    return result;
}

SECP256K1_INLINE static uint32_t secp256k1_heap_remove(secp256k1_scalar_heap *heap) {
    uint32_t result = heap->tree[0];
    VERIFY_CHECK(heap->size > 0);
    if (--heap->size > 0) {
        secp256k1_sift_down(heap, 0, heap->tree[heap->size]);
    }
    return result;
}

/** Multi-multiply: R = sum_i ni * Ai */
static void secp256k1_ecmult_multi_bos_coster(uint32_t *tree_space, secp256k1_gej *r, secp256k1_scalar *sc, secp256k1_gej *pt, size_t n) {
    secp256k1_scalar_heap heap;
    uint32_t first, second;

    secp256k1_gej_set_infinity(r);
    secp256k1_heap_initialize(&heap, tree_space, sc, pt, n);

    if (heap.size == 0) {
        return;
    }

    first = secp256k1_heap_remove(&heap);

    while (heap.size > 0) {
        second = heap.tree[0];        

        do {
            /* Observe that nX + mY = (n-m)X + m(X + Y), and if n > m this transformation
             * reduces the magnitude of the larger scalar, significantly because X and Y are
             * chosen to be the two largest values, and therefore will be similar in magnitude.
             * So by repeating this we will quickly zero out all but one exponent, which will
             * be small. */
            secp256k1_gej_add_var(&pt[second], &pt[first], &pt[second], NULL);  /* Y -> X + Y */
            secp256k1_scalar_numsub(&sc[first], &sc[first], &sc[second]);  /* n -> n - m */

            if (secp256k1_scalar_cmp_var(&sc[first], &sc[second]) < 0) {
                break;
            }

            /* For pathological inputs, n and m may not be similar in magnitude (e.g. if
             * n ~ 2^256 and m ~ 1. In this case the above step will not reduce the magnitude
             * of the larger scalar, which we detect with the above condition. In this case
             * we simply halve the scalar and double its point, ensuring we make progress. */
            if (secp256k1_scalar_shr_int(&sc[first], 1) == 1) {
                secp256k1_gej_add_var(r, r, &pt[first], NULL);
            }
            secp256k1_gej_double_var(&pt[first], &pt[first], NULL);
        }
        while (secp256k1_scalar_cmp_var(&sc[first], &sc[second]) >= 0);

        if (secp256k1_scalar_is_zero(&sc[first])) {
            first = secp256k1_heap_remove(&heap);
        } else {
            first = secp256k1_replace(&heap, first);
        }
    }

    VERIFY_CHECK(!secp256k1_scalar_is_zero(&sc[first]));

    /* Now the desired result is heap_sc[0] * heap_pt[0], and for random scalars it is
     * very likely that heap_sc[0] = 1, and extremely likely heap_sc[0] < 5. (After
     * about 100k trials I saw around 200 2's and one 3.) So use a binary ladder rather
     * than any heavy machinery to finish it off. */
    for (;;) {
        if (secp256k1_scalar_shr_int(&sc[first], 1) == 1) {
            secp256k1_gej_add_var(r, r, &pt[first], NULL);
            if (secp256k1_scalar_is_zero(&sc[first])) {
                break;
            }
        }
        secp256k1_gej_double_var(&pt[first], &pt[first], NULL);
    }
}

#ifdef USE_ENDOMORPHISM
SECP256K1_INLINE static void secp256k1_ecmult_endo_split(secp256k1_scalar *s1, secp256k1_scalar *s2, secp256k1_gej *p1, secp256k1_gej *p2) {
    secp256k1_scalar tmp = *s1;
    secp256k1_scalar_split_lambda(s1, s2, &tmp);
    secp256k1_gej_mul_lambda(p2, p1);

    if (secp256k1_scalar_is_high(s1)) {
        secp256k1_scalar_negate(s1, s1);
        secp256k1_gej_neg(p1, p1);
    }
    if (secp256k1_scalar_is_high(s2)) {
        secp256k1_scalar_negate(s2, s2);
        secp256k1_gej_neg(p2, p2);
    }
}
#endif

static int secp256k1_ecmult_multi(const secp256k1_ecmult_context *ctx, secp256k1_scratch *scratch, const secp256k1_callback* error_callback, secp256k1_gej *r, const secp256k1_scalar *inp_g_sc, secp256k1_ecmult_multi_callback cb, void *cbdata, size_t n) {
    const size_t entry_size = sizeof(secp256k1_gej) + sizeof(secp256k1_scalar) + sizeof(size_t);
    const size_t max_entries = secp256k1_scratch_max_allocation(scratch) / entry_size;
    /* Use 2(n+1) with the endomorphism, n+1 without, when calculating batch sizes.
     * The reason for +1 is that Bos-Coster requires we add the G scalar to the list of
     * other scalars. */
#ifdef USE_ENDOMORPHISM
    const size_t n_batches = (2*n + max_entries + 1) / max_entries;
    size_t entries_per_batch = (2*n + n_batches + 1) / n_batches;
#else
    const size_t n_batches = (n + max_entries) / max_entries;
    size_t entries_per_batch = (n + n_batches) / n_batches;
#endif

    secp256k1_gej tmp;
    secp256k1_gej *pt;
    secp256k1_scalar *sc;
    uint32_t *tree_space;
    size_t idx = 0;
    size_t point_idx = 0;

    /* Attempt to allocate sufficient space for Bos-Coster */
    while (!secp256k1_scratch_resize(scratch, error_callback, entries_per_batch * entry_size)) {
        entries_per_batch /= 2;
        if (entries_per_batch < 2) {
            return 0;
        }
    }
    secp256k1_scratch_reset(scratch);
    pt = (secp256k1_gej *) secp256k1_scratch_alloc(scratch, entries_per_batch * sizeof(*pt));
    sc = (secp256k1_scalar *) secp256k1_scratch_alloc(scratch, entries_per_batch * sizeof(*sc));
    tree_space = (uint32_t *) secp256k1_scratch_alloc(scratch, entries_per_batch * sizeof(*tree_space));
    VERIFY_CHECK(pt != NULL);
    VERIFY_CHECK(sc != NULL);
    VERIFY_CHECK(ctx != NULL);
    VERIFY_CHECK(tree_space != NULL);

    sc[0] = *inp_g_sc;
    secp256k1_gej_set_ge(&pt[0], &secp256k1_ge_const_g);
    idx++;
#ifdef USE_ENDOMORPHISM
    secp256k1_ecmult_endo_split(&sc[0], &sc[1], &pt[0], &pt[1]);
    idx++;
#endif

    secp256k1_gej_set_infinity(r);
    while (point_idx < n) {
        if (!cb(&sc[idx], &pt[idx], point_idx, cbdata)) {
            return 0;
        }
        idx++;
#ifdef USE_ENDOMORPHISM
        secp256k1_ecmult_endo_split(&sc[idx - 1], &sc[idx], &pt[idx - 1], &pt[idx]);
        idx++;
        if (idx >= entries_per_batch - 1) {
#else
        if (idx >= entries_per_batch) {
#endif
            secp256k1_ecmult_multi_bos_coster(tree_space, &tmp, sc, pt, idx);
            secp256k1_gej_add_var(r, r, &tmp, NULL);
            idx = 0;
        }
        point_idx++;
    }
    secp256k1_ecmult_multi_bos_coster(tree_space, &tmp, sc, pt, idx);
    secp256k1_gej_add_var(r, r, &tmp, NULL);
    return 1;
}

#endif
