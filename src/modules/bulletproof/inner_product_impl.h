/**********************************************************************
 * Copyright (c) 2017 Andrew Poelstra                                 *
 * Distributed under the MIT software license, see the accompanying   *
 * file COPYING or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#ifndef SECP256K1_MODULE_BULLETPROOF_INNER_PRODUCT_IMPL
#define SECP256K1_MODULE_BULLETPROOF_INNER_PRODUCT_IMPL

#include "group.h"
#include "scalar.h"

#include "modules/bulletproof/main_impl.h"
#include "modules/bulletproof/generators.h"
#include "modules/bulletproof/util.h"

#define POPCOUNT(x)	(__builtin_popcountl((unsigned long)(x)))  /* TODO make these portable */
#define CTZ(x)	(__builtin_ctzl((unsigned long)(x)))

typedef int (secp256k1_bulletproof_vfy_callback)(secp256k1_scalar *sc, secp256k1_ge *pt, secp256k1_scalar *randomizer, size_t idx, void *data);

typedef struct {
    const unsigned char *serialized_points;
    size_t n_extra_ser_points;
    secp256k1_scalar a;
    secp256k1_scalar b;
    secp256k1_scalar dot;
    secp256k1_scalar p_offs;
    unsigned char commit[32];
    secp256k1_bulletproof_vfy_callback *rangeproof_cb;
    void *rangeproof_cb_data;
    size_t n_extra_rangeproof_points;
} secp256k1_bulletproof_innerproduct_context;

typedef struct {
    const secp256k1_bulletproof_innerproduct_context *proof;
    secp256k1_scalar xsq[SECP256K1_BULLETPROOF_MAX_DEPTH + 1];
    secp256k1_scalar xsqinv[SECP256K1_BULLETPROOF_MAX_DEPTH + 1];
    secp256k1_scalar xcache[SECP256K1_BULLETPROOF_MAX_DEPTH + 1];
    secp256k1_scalar xsqinv_mask;
} secp256k1_bulletproof_innerproduct_vfy_data;

typedef struct {
    size_t n_proofs;
    const secp256k1_ge *geng;
    const secp256k1_ge *genh;
    size_t n_gen_pairs;
    size_t ceil_lg_gens;
    secp256k1_scalar randomizer[MAX_BATCH_QTY];
    secp256k1_bulletproof_innerproduct_vfy_data proof[MAX_BATCH_QTY];
} secp256k1_bulletproof_innerproduct_vfy_ecmult_context;

/* Bulletproof rangeproof verification comes down to a single multiexponentiation of the form
 *
 *   P + (c-a*b)*x*G - sum_{i=1}^n [a*s'_i*G_i + b*s_i*H_i] + sum_{i=1}^log2(n) [x_i^-2 L_i + x_i^2 R_i
 *
 * which will equal infinity if the rangeproof is correct. Here
 *   - `G_i` and `H_i` are standard NUMS generators. `G` is the standard secp256k1 generator.
 *   - `P` and `c` are inputs to the proof, which claims that there exist `a_i` and `b_i`, `i` ranging
 *     from 0 to `n-1`, such that `P = sum_i [a_i G_i + b_i H_i]` and that `<{a_i}, {b_i}> = c`.
 *   - `a`, `b`, `L_i` and `R_i`are auxillary components of the proof, where `i` ranges from 0 to `log2(n)-1`.
 *   - `x_i = H(x_{i-1} || L_i || R_i)`, where `x_{-1}` is passed through the `commit` variable and
 *     must be a commitment to `P` and `c`.
 *   - `x` is a hash of `commit` and is used to rerandomize `c`. See Protocol 2 vs Protocol 1 in the paper.
 *   - `s_i` and `s'_i` are computed as follows.
 *
 * For each `i` between 0 and `n-1` inclusive, let `b_{ij}` be -1 (1) if the `j`th bit of `i` is zero (one).
 * Here `j` ranges from 0 to `log2(n)-1`. Then for each such `i` we define
 *   - `s_i = prod_j x_j^{b_{ij}}`
 *   - `s'_i = 1/s_i`
 *
 * Alternately we can define `s_i` and `s'_i` recursively as follows:
 *   - `s_0 = s`_{n - 1} = 1 / prod_j x_j`
 *   - `s_i = s'_{n - 1 - i} = s_{i - 2^j} * x_j^2` where `j = i & (i - 1)` is `i` with its least significant 1 set to 0.
 *
 * Our ecmult_multi function takes `(c - a*b)*x` directly and multiplies this by `G`. For every other
 * (scalar, point) pair it calls the following callback function, which takes an index and outputs a
 * pair. The function therefore has three regimes:
 *
 * For the first `2n` invocations, it alternately returns `(s'_{n - i}, G_{n - i})` and `(s_i, H_i)`,
 * where `i` is `floor(idx / 2)`. The reason for the funny indexing is that we use the above recursive
 * definition of `s_i` and `s'_i` which produces each element with only a single scalar multiplication,
 * but in this mixed order. (We start with an array of `x_j^2` for each `x_j`.)
 *
 * As a side-effect, whenever `n - i = 2^j` for some `j`, `s_i = x_j^{-1} * prod_{j' != j} x_{j'}`,
 * so `x_j^{-2} = s_i*s_0`. Therefore we compute an array of inverse squares during this computation,
 * using only one multiplication per. We will need it in the following step.
 *
 * For the next `2*log2(n)` invocations it alternately returns `(x_i^-2, L_i)` and `(x_i^2, R_i)`
 * where `i` is `idx - 2*n`.
 *
 * For the remaining invocations it passes through to another callback, `rangeproof_cb_data` which
 * computes `P`. The reason for this is that in practice `P` is usually defined by another multiexp
 * rather than being a known point, and it is more efficient to compute one exponentiation.
 *
 */

/* For the G and H generators, we choose the ith generator with a scalar computed from the
 * L/R hashes as follows: prod_{j=1}^m x_j^{e_j}, where each exponent e_j is either -1 or 1.
 * The choice directly maps to the bits of i: for the G generators, a 0 bit means e_j is 1
 * and a 1 bit means e_j is -1. For the H generators it is the opposite. Finally, each of the
 * G scalars is further multiplied by -a, while each of the H scalars is further multiplied
 * by -b.
 *
 * These scalars are computed starting from I, the inverse of the product of every x_j, which
 * is then selectively multiplied by x_j^2 for whichever j's are needed. As it turns out, by
 * caching logarithmically many scalars, this can always be done by multiplying one of the
 * cached values by a single x_j, rather than starting from I and doing multiple multiplications.
 */

static int secp256k1_bulletproof_innerproduct_vfy_ecmult_callback(secp256k1_scalar *sc, secp256k1_ge *pt, size_t idx, void *data) {
    secp256k1_bulletproof_innerproduct_vfy_ecmult_context *ctx = (secp256k1_bulletproof_innerproduct_vfy_ecmult_context *) data;

    /* First 2N points use the standard Gi, Hi generators, and the scalars can be aggregated across proofs  */
    if (idx < 2 * ctx->n_gen_pairs) {
        size_t i;
        /* TODO zero this point when appropriate for non-2^n numbers of pairs */
        if (idx < ctx->n_gen_pairs) {
            *pt = ctx->geng[idx];
        } else {
            *pt = ctx->genh[idx - ctx->n_gen_pairs];
        }

        secp256k1_scalar_clear(sc);
        for (i = 0; i < ctx->n_proofs; i++) {
            const size_t cache_idx = POPCOUNT(idx);
            secp256k1_scalar term;
            VERIFY_CHECK(cache_idx < SECP256K1_BULLETPROOF_MAX_DEPTH);
            /* Compute the normal inner-product scalar... */
            if (cache_idx > 0) {
                if (idx < ctx->n_gen_pairs) {
                    const size_t xsq_idx = ctx->ceil_lg_gens - 1 - CTZ(idx);
                    secp256k1_scalar_mul(&ctx->proof[i].xcache[cache_idx], &ctx->proof[i].xcache[cache_idx - 1], &ctx->proof[i].xsq[xsq_idx]);
                    if (idx == ctx->n_gen_pairs - 1) {
                        ctx->proof[i].xcache[0] = ctx->proof[i].xcache[cache_idx];
                    }
                } else if (idx == ctx->n_gen_pairs) {
                    const size_t prev_cache_idx = POPCOUNT(idx - 1);
                    secp256k1_scalar_mul(
                        &ctx->proof[i].xcache[cache_idx],
                        &ctx->proof[i].xcache[prev_cache_idx],
                        &ctx->proof[i].xsqinv[prev_cache_idx]
                    );
                } else {
                    const size_t xsq_idx = ctx->ceil_lg_gens - 1 - CTZ(idx);
                    secp256k1_scalar_mul(&ctx->proof[i].xcache[cache_idx], &ctx->proof[i].xcache[cache_idx - 1], &ctx->proof[i].xsqinv[xsq_idx]);
                }
            }
            term = ctx->proof[i].xcache[cache_idx];

            /* When going through the G generators, compute the x-inverses as side effects */
            if (idx < ctx->n_gen_pairs && POPCOUNT(idx) == ctx->ceil_lg_gens - 1) {  /* if the scalar has only one 0, i.e. only one inverse... */
                const size_t xsqinv_idx = ctx->ceil_lg_gens - 1 - CTZ(~idx);
                /* ...multiply it by the total inverse, to get x_j^-2 */
                secp256k1_scalar_mul(&ctx->proof[i].xsqinv[xsqinv_idx], &ctx->proof[i].xcache[cache_idx], &ctx->proof[i].xsqinv_mask);
            }

            /* ...add whatever offset the rangeproof wants... */
            if (ctx->proof[i].proof->rangeproof_cb != NULL) {
                secp256k1_scalar rangeproof_offset;
                if ((ctx->proof[i].proof->rangeproof_cb)(&rangeproof_offset, NULL, &ctx->randomizer[i], idx, ctx->proof[i].proof->rangeproof_cb_data) == 0) {
                    return 0;
                }
                /* evil hack to change H_i to H_i' */
                if (idx >= ctx->n_gen_pairs) {
                    secp256k1_scalar_mul(&term, &term, (secp256k1_scalar *) ctx->proof[i].proof->rangeproof_cb_data);
                }
                secp256k1_scalar_add(&term, &term, &rangeproof_offset);
            }

            secp256k1_scalar_add(sc, sc, &term);
        }
    /* Next 2lgN points are the L and R vectors */
    } else if (idx < 2 * (ctx->n_gen_pairs + ctx->ceil_lg_gens * ctx->n_proofs)) {
        size_t real_idx = idx / 2 - ctx->n_gen_pairs;
        const size_t proof_idx = real_idx / ctx->ceil_lg_gens;
        real_idx = real_idx % ctx->ceil_lg_gens;
        if (idx % 2 == 0) {
            secp256k1_bulletproof_deserialize_point(
                pt,
                ctx->proof[proof_idx].proof->serialized_points,
                ctx->proof[proof_idx].proof->n_extra_ser_points + real_idx,
                ctx->proof[proof_idx].proof->n_extra_ser_points + 2 * ctx->ceil_lg_gens
            );
            *sc = ctx->proof[proof_idx].xsq[real_idx];
        } else {
            secp256k1_bulletproof_deserialize_point(
                pt,
                ctx->proof[proof_idx].proof->serialized_points,
                ctx->proof[proof_idx].proof->n_extra_ser_points + real_idx + ctx->ceil_lg_gens,
                ctx->proof[proof_idx].proof->n_extra_ser_points + 2 * ctx->ceil_lg_gens
            );
            *sc = ctx->proof[proof_idx].xsqinv[real_idx];
        }
        secp256k1_scalar_mul(sc, sc, &ctx->randomizer[proof_idx]);
    /* Remaining points are whatever the rangeproof wants */
    } else if (idx == 2 * (ctx->n_gen_pairs + ctx->ceil_lg_gens * ctx->n_proofs)) {
        /* Special case: the first extra point is independent of the proof, for both rangeproof and circuit */
        size_t i;
        secp256k1_scalar_clear(sc);
        for (i = 0; i < ctx->n_proofs; i++) {
            secp256k1_scalar term;
            if ((ctx->proof[i].proof->rangeproof_cb)(&term, pt, &ctx->randomizer[i], 2 * (ctx->n_gen_pairs + ctx->ceil_lg_gens), ctx->proof[i].proof->rangeproof_cb_data) == 0) {
                return 0;
            }
            secp256k1_scalar_add(sc, sc, &term);
        }
    } else {
        size_t proof_idx = 0;
        size_t real_idx = idx - 2 * (ctx->n_gen_pairs + ctx->ceil_lg_gens * ctx->n_proofs) - 1;
        while (real_idx >= ctx->proof[proof_idx].proof->n_extra_rangeproof_points - 1) {
            real_idx -= ctx->proof[proof_idx].proof->n_extra_rangeproof_points - 1;
            proof_idx++;
            VERIFY_CHECK(proof_idx < ctx->n_proofs);
        }
        if ((ctx->proof[proof_idx].proof->rangeproof_cb)(sc, pt, &ctx->randomizer[proof_idx], 2 * (ctx->n_gen_pairs + ctx->ceil_lg_gens), ctx->proof[proof_idx].proof->rangeproof_cb_data) == 0) {
            return 0;
        }
    }

    return 1;
}

/* nb For security it is essential that `commit_inp` already commit to all data
 *    needed to compute `P`. We do not hash it in during verification since `P`
 *    may be specified indirectly as a bunch of scalar offsets.
 */
static int secp256k1_bulletproof_inner_product_verify_impl(const secp256k1_context* ctx, const secp256k1_ecmult_context *ecmult_ctx, secp256k1_scratch *scratch, const secp256k1_ge *geng, const secp256k1_ge *genh, size_t n_gens, const secp256k1_bulletproof_innerproduct_context *proof, size_t n_proofs) {
    secp256k1_sha256 sha256;
    secp256k1_bulletproof_innerproduct_vfy_ecmult_context ecmult_data;
    unsigned char commit[32];
    size_t total_n_points = 2 * n_gens + 1; /* +1 for shared G */
    secp256k1_gej r;
    secp256k1_scalar p_offs;
    size_t i;

    ecmult_data.n_proofs = n_proofs;
    ecmult_data.geng = geng;
    ecmult_data.genh = genh;
    ecmult_data.n_gen_pairs = n_gens;
    ecmult_data.ceil_lg_gens = CTZ(n_gens);
    /* Seed RNG for per-proof randomizers */
    secp256k1_sha256_initialize(&sha256);
    for (i = 0; i < n_proofs; i++) {
        const size_t n_points = 2 * ecmult_data.ceil_lg_gens + proof[i].n_extra_ser_points;
        const size_t sp_len = 32 * n_points + (n_points + 7) / 8;
        secp256k1_sha256_write(&sha256, proof[i].serialized_points, sp_len);
        secp256k1_scalar_get_b32(commit, &proof[i].a);
        secp256k1_sha256_write(&sha256, commit, 32);
        secp256k1_scalar_get_b32(commit, &proof[i].b);
        secp256k1_sha256_write(&sha256, commit, 32);
        secp256k1_scalar_get_b32(commit, &proof[i].dot);
        secp256k1_sha256_write(&sha256, commit, 32);
        secp256k1_scalar_get_b32(commit, &proof[i].p_offs);
        secp256k1_sha256_write(&sha256, commit, 32);
    }
    secp256k1_sha256_finalize(&sha256, commit);

    secp256k1_scalar_clear(&p_offs);
    for (i = 0; i < n_proofs; i++) {
        secp256k1_scalar negprod;
        secp256k1_scalar x;
        int overflow;
        size_t j;

        total_n_points += 2 * ecmult_data.ceil_lg_gens + proof[i].n_extra_rangeproof_points - 1; /* -1 for shared G */

        ecmult_data.proof[i].proof = &proof[i];
        /* set per-proof randomizer */
        secp256k1_sha256_initialize(&sha256);
        secp256k1_sha256_write(&sha256, commit, 32);
        secp256k1_sha256_finalize(&sha256, commit);
        secp256k1_scalar_set_b32(&ecmult_data.randomizer[i], commit, &overflow);
        if (overflow || secp256k1_scalar_is_zero(&ecmult_data.randomizer[i])) {
            /* cryptographically unreachable */
            return 0;
        }

        /* Compute x*(dot - a*b) for each proof */
        secp256k1_scalar_set_b32(&x, proof[i].commit, &overflow);
        if (overflow || secp256k1_scalar_is_zero(&x)) {
            return 0;
        }
        secp256k1_scalar_mul(&negprod, &proof[i].a, &proof[i].b);
        secp256k1_scalar_negate(&negprod, &negprod);
        secp256k1_scalar_add(&negprod, &negprod, &proof[i].dot);
        secp256k1_scalar_mul(&x, &x, &negprod);
        secp256k1_scalar_add(&x, &x, &proof[i].p_offs);

        secp256k1_scalar_mul(&x, &x, &ecmult_data.randomizer[i]);
        secp256k1_scalar_add(&p_offs, &p_offs, &x);

        /* Compute the inverse product and the array of squares; the rest will be filled
         * in by the callback during the multiexp. */
        memcpy(commit, proof[i].commit, 32);
        ecmult_data.proof[i].xcache[0] = ecmult_data.randomizer[i];
        for (j = 0; j < ecmult_data.ceil_lg_gens; j++) {
            secp256k1_scalar xi;
            const size_t lidx = proof[i].n_extra_ser_points + j;
            const size_t ridx = proof[i].n_extra_ser_points + j + ecmult_data.ceil_lg_gens;
            const size_t bitveclen = (proof[i].n_extra_ser_points + 2 * ecmult_data.ceil_lg_gens + 7) / 8;
            const unsigned char lrparity = 2 * !!(proof[i].serialized_points[lidx / 8] & (1 << (lidx % 8))) + !!(proof[i].serialized_points[ridx / 8] & (1 << (ridx % 8)));
            /* Map commit -> H(commit || LR parity || Lx || Rx), compute xi from it */
            secp256k1_sha256_initialize(&sha256);
            secp256k1_sha256_write(&sha256, commit, 32);
            secp256k1_sha256_write(&sha256, &lrparity, 1);
            secp256k1_sha256_write(&sha256, &proof[i].serialized_points[32 * lidx + bitveclen], 32);
            secp256k1_sha256_write(&sha256, &proof[i].serialized_points[32 * ridx + bitveclen], 32);
            secp256k1_sha256_finalize(&sha256, commit);

            secp256k1_scalar_set_b32(&xi, commit, &overflow);
            if (overflow || secp256k1_scalar_is_zero(&xi)) {
                return 0;
            }
            secp256k1_scalar_mul(&ecmult_data.proof[i].xcache[0], &ecmult_data.proof[i].xcache[0], &xi);
            secp256k1_scalar_sqr(&ecmult_data.proof[i].xsq[j], &xi);
        }
        /* Compute (-a * r * x1 * ... * xn)^-1 which will be used to mask out individual x_i^-2's */
        secp256k1_scalar_negate(&ecmult_data.proof[i].xsqinv_mask, &proof[i].a);
        secp256k1_scalar_mul(&ecmult_data.proof[i].xsqinv_mask, &ecmult_data.proof[i].xsqinv_mask, &ecmult_data.proof[i].xcache[0]);
        secp256k1_scalar_inverse_var(&ecmult_data.proof[i].xsqinv_mask, &ecmult_data.proof[i].xsqinv_mask);

        /* Extract a*b^-1 which is used to switch from multiplication by b to multiplication by a. Use negprod as dummy. */
        secp256k1_scalar_mul(&negprod, &ecmult_data.proof[i].xsqinv_mask, &ecmult_data.proof[i].xcache[0]);  /* -a^-1 */
        secp256k1_scalar_mul(&ecmult_data.proof[i].xsqinv[j], &negprod, &proof[i].b);  /* -a^-1*b */
        secp256k1_scalar_negate(&ecmult_data.proof[i].xsqinv[j], &ecmult_data.proof[i].xsqinv[j]); /* a^-1*b */

        /* Extract -a * r * (x1 * ... * xn)^-1 which is our first coefficient */
        secp256k1_scalar_mul(&negprod, &ecmult_data.randomizer[i], &proof[i].a); /* r*a */
        secp256k1_scalar_sqr(&negprod, &negprod); /* (r*a)^2 */
        secp256k1_scalar_mul(&ecmult_data.proof[i].xcache[0], &ecmult_data.proof[i].xsqinv_mask, &negprod);  /* -a * r * (x1 * x2 * ... * xn)^-1 */
    }

    /* Do the multiexp */
    if (secp256k1_ecmult_multi_var(ecmult_ctx, scratch, &ctx->error_callback, &r, &p_offs, secp256k1_bulletproof_innerproduct_vfy_ecmult_callback, (void *) &ecmult_data, total_n_points) != 1) {
        return 0;
    }
    return secp256k1_gej_is_infinity(&r);
}

static void secp256k1_scalar_dot_product(secp256k1_scalar *r, const secp256k1_scalar *a, const secp256k1_scalar *b, size_t n) {
    secp256k1_scalar_clear(r);
    while(n--) {
        secp256k1_scalar term;
        secp256k1_scalar_mul(&term, &a[n], &b[n]);
        secp256k1_scalar_add(r, r, &term);
    }
}

typedef struct {
    secp256k1_gej *pt1;
    const secp256k1_scalar *sc1;
    secp256k1_gej *pt2;
    const secp256k1_scalar *sc2;
    size_t n;
} secp256k1_bulletproof_innerproduct_pf_ecmult_context;

static int secp256k1_bulletproof_innerproduct_pf_ecmult_callback(secp256k1_scalar *sc, secp256k1_ge *pt, size_t idx, void *data) {
    secp256k1_bulletproof_innerproduct_pf_ecmult_context *ctx = (secp256k1_bulletproof_innerproduct_pf_ecmult_context *) data;
    if (idx < ctx->n) {
        secp256k1_ge_set_gej(pt, &ctx->pt1[idx]);
        *sc = ctx->sc1[idx];
    } else {
        secp256k1_ge_set_gej(pt, &ctx->pt2[idx - ctx->n]);
        *sc = ctx->sc2[idx - ctx->n];
    }
    return 1;
}
/* These proofs are not zero-knowledge. There is no need to worry about constant timeness.
 * Having said that, we avoid using ecmult_multi because it'd clobber our scratch space.
 * `commit_inp` must contain 256 bits of randomness, it is used immediately as a randomizer.
 */
static int secp256k1_bulletproof_inner_product_prove_impl(const secp256k1_context* ctx, const secp256k1_ecmult_context *ecmult_ctx, secp256k1_scratch *scratch, secp256k1_scalar *final_a, secp256k1_scalar *final_b, secp256k1_ge *lpt_arr, secp256k1_ge *rpt_arr, const size_t n, secp256k1_ecmult_multi_callback *cb, void *cb_data, const unsigned char *commit_inp) {
    const size_t depth = secp256k1_ceil_lg(n);
    const size_t next_pow2 = 1ull << depth;
    secp256k1_scalar zero;
    size_t i;
    unsigned char commit[32];
    secp256k1_scalar *a_arr;
    secp256k1_scalar *b_arr;
    secp256k1_gej *geng;
    secp256k1_gej *genh;
    secp256k1_scalar ux;
    int overflow;
    secp256k1_scratch *myscr;

    /* Special-case the depth-0 proof */
    if (depth == 0) {
        secp256k1_ge tmpge;
        cb(final_a, &tmpge, 0, cb_data);
        cb(final_b, &tmpge, 1, cb_data);
        return 1;
    }

    /* setup */
    VERIFY_CHECK(depth < SECP256K1_BULLETPROOF_MAX_DEPTH);

    secp256k1_scalar_clear(&zero);

    if (!secp256k1_scratch_resize(scratch, &ctx->error_callback, 2 * next_pow2 * (sizeof(secp256k1_scalar) + sizeof(secp256k1_gej)), 4)) {
        return 0;
    }

    VERIFY_CHECK(n == next_pow2);
    secp256k1_scratch_reset(scratch);
    a_arr = (secp256k1_scalar*)secp256k1_scratch_alloc(scratch, next_pow2 * sizeof(secp256k1_scalar));
    b_arr = (secp256k1_scalar*)secp256k1_scratch_alloc(scratch, next_pow2 * sizeof(secp256k1_scalar));
    geng = (secp256k1_gej*)secp256k1_scratch_alloc(scratch, next_pow2 * sizeof(secp256k1_gej));
    genh = (secp256k1_gej*)secp256k1_scratch_alloc(scratch, next_pow2 * sizeof(secp256k1_gej));
    VERIFY_CHECK(a_arr != NULL);
    VERIFY_CHECK(b_arr != NULL);
    VERIFY_CHECK(geng != NULL);
    VERIFY_CHECK(genh != NULL);

    myscr = secp256k1_scratch_create(NULL, 100000000, 100000000);  /* 100M should be waay overkill */

    memcpy(commit, commit_inp, 32);
    for (i = 0; i < n; i++) {
        secp256k1_ge tmpge;
        cb(&a_arr[i], &tmpge, 2*i, cb_data);
        secp256k1_gej_set_ge(&geng[i], &tmpge);
        cb(&b_arr[i], &tmpge, 2*i+1, cb_data);
        secp256k1_gej_set_ge(&genh[i], &tmpge);
    }

    /* Protocol 2: obtain randomizer */
    secp256k1_scalar_set_b32(&ux, commit, &overflow);
    if (overflow || secp256k1_scalar_is_zero(&ux)) {
        return 0;
    }

    /* Protocol 1: Iterate, halving vector size until it is 1 */
    for (i = 0; i < depth; i++) {
        const size_t halfwidth = n >> (i + 1);
        secp256k1_bulletproof_innerproduct_pf_ecmult_context pfdata;
        secp256k1_gej tmplj, tmprj;
        secp256k1_scalar x, xinv;
        secp256k1_scalar tmps;
        size_t j;

        /* L */
        secp256k1_scalar_dot_product(&tmps, &a_arr[0], &b_arr[halfwidth], halfwidth);
        secp256k1_scalar_mul(&tmps, &tmps, &ux);
        pfdata.pt1 = &geng[halfwidth];
        pfdata.sc1 = &a_arr[0];
        pfdata.pt2 = &genh[0];
        pfdata.sc2 = &b_arr[halfwidth];
        pfdata.n = halfwidth;
        secp256k1_ecmult_multi_var(ecmult_ctx, myscr, &ctx->error_callback, &tmplj, &tmps, &secp256k1_bulletproof_innerproduct_pf_ecmult_callback, (void *) &pfdata, halfwidth * 2);
        secp256k1_ge_set_gej(&lpt_arr[i], &tmplj);
        /* R */
        secp256k1_scalar_dot_product(&tmps, &b_arr[0], &a_arr[halfwidth], halfwidth);
        secp256k1_scalar_mul(&tmps, &tmps, &ux);
        pfdata.pt1 = &geng[0];
        pfdata.sc1 = &a_arr[halfwidth];
        pfdata.pt2 = &genh[halfwidth];
        pfdata.sc2 = &b_arr[0];
        pfdata.n = halfwidth;
        secp256k1_ecmult_multi_var(ecmult_ctx, myscr, &ctx->error_callback, &tmprj, &tmps, &secp256k1_bulletproof_innerproduct_pf_ecmult_callback, (void *) &pfdata, halfwidth * 2);
        secp256k1_ge_set_gej(&rpt_arr[i], &tmprj);

        /* x, x^2, x^-1, x^-2 */
        secp256k1_bulletproof_update_commit(commit, &lpt_arr[i], &rpt_arr[i]);
        secp256k1_scalar_set_b32(&x, commit, &overflow);
        if (overflow || secp256k1_scalar_is_zero(&x)) {
            return 0;
        }
        secp256k1_scalar_inverse_var(&xinv, &x);

        /* update generators and scalar array */
        for (j = 0; j < halfwidth; j++) {
            secp256k1_gej tmp1j, tmp2j;

            secp256k1_ecmult(ecmult_ctx, &tmp1j, &geng[j], &xinv, &zero);
            secp256k1_ecmult(ecmult_ctx, &tmp2j, &geng[j + halfwidth], &x, &zero);
            secp256k1_gej_add_var(&geng[j], &tmp1j, &tmp2j, NULL);

            secp256k1_ecmult(ecmult_ctx, &tmp1j, &genh[j], &x, &zero);
            secp256k1_ecmult(ecmult_ctx, &tmp2j, &genh[j + halfwidth], &xinv, &zero);
            secp256k1_gej_add_var(&genh[j], &tmp1j, &tmp2j, NULL);

            secp256k1_scalar_mul(&a_arr[j], &a_arr[j], &x);
            secp256k1_scalar_mul(&tmps, &a_arr[j + halfwidth], &xinv);
            secp256k1_scalar_add(&a_arr[j], &a_arr[j], &tmps);

            secp256k1_scalar_mul(&b_arr[j], &b_arr[j], &xinv);
            secp256k1_scalar_mul(&tmps, &b_arr[j + halfwidth], &x);
            secp256k1_scalar_add(&b_arr[j], &b_arr[j], &tmps);
        }
    }

    *final_a = a_arr[0];
    *final_b = b_arr[0];
    secp256k1_scratch_space_destroy(myscr);
    return 1;
}

#endif
