/**********************************************************************
 * Copyright (c) 2017 Andrew Poelstra                                 *
 * Distributed under the MIT software license, see the accompanying   *
 * file COPYING or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#ifndef SECP256K1_MODULE_BULLETPROOF_RANGEPROOF_IMPL
#define SECP256K1_MODULE_BULLETPROOF_RANGEPROOF_IMPL

#include "modules/bulletproof/inner_product_impl.h"
#include "modules/bulletproof/util.h"
#include "group.h"

#define MAX_NBITS	64

typedef struct {
    secp256k1_scalar mul_by;
    secp256k1_scalar yinv;
    secp256k1_scalar z;
    secp256k1_scalar z_randomized;
    secp256k1_scalar zsq;
    secp256k1_scalar g_exponent;
    secp256k1_scalar negz;
    secp256k1_scalar x;
    secp256k1_ge a;
    secp256k1_ge s;
    size_t n;
    /* eq (61) stuff */
    size_t count;
    secp256k1_scalar randomizer61;
    secp256k1_scalar y;
    secp256k1_scalar t;
    const secp256k1_ge *asset;
    const secp256k1_ge *commit;
    size_t n_commits;
    secp256k1_ge t1;
    secp256k1_ge t2;
} secp256k1_bulletproof_vfy_ecmult_context;

static int secp256k1_bulletproof_rangeproof_vfy_callback(secp256k1_scalar *sc, secp256k1_ge *pt, secp256k1_scalar *randomizer, size_t idx, void *data) {
    secp256k1_bulletproof_vfy_ecmult_context *ctx = (secp256k1_bulletproof_vfy_ecmult_context *) data;

    if (idx == 0) {
        secp256k1_scalar_mul(&ctx->g_exponent, &ctx->negz, randomizer);
        secp256k1_scalar_mul(&ctx->z_randomized, &ctx->z, randomizer);
    }

    if (idx < ctx->n) {
        *sc = ctx->g_exponent;
    } else if (idx < 2 * ctx->n) {
        const size_t nbits = ctx->n / ctx->n_commits;
        const size_t commit_idx = (idx - ctx->n) / nbits;
        const size_t bit_idx = (idx - ctx->n) % nbits;

        if (bit_idx == 0) {
            size_t i;
            secp256k1_scalar_sqr(&ctx->zsq, &ctx->z);
            for (i = 0; i < commit_idx; i++) {
                secp256k1_scalar_mul(&ctx->zsq, &ctx->zsq, &ctx->z);
            }
            secp256k1_scalar_mul(&ctx->zsq, &ctx->zsq, &ctx->mul_by);
            secp256k1_scalar_mul(&ctx->zsq, &ctx->zsq, &ctx->yinv);
            secp256k1_scalar_mul(&ctx->zsq, &ctx->zsq, randomizer);
        }
        secp256k1_scalar_add(sc, &ctx->zsq, &ctx->z_randomized);

        secp256k1_scalar_mul(&ctx->zsq, &ctx->zsq, &ctx->yinv);
        secp256k1_scalar_add(&ctx->zsq, &ctx->zsq, &ctx->zsq);
        secp256k1_scalar_mul(&ctx->mul_by, &ctx->mul_by, &ctx->yinv);
    } else {
        switch(ctx->count) {
        /* S^x in eq (62) */
        case 0:
            *sc = ctx->x;
            *pt = ctx->s;
            break;
        /* A in eq (62) */
        case 1:
            *pt = ctx->a;
            secp256k1_scalar_set_int(sc, 1);
            break;
        /* G^[k(y, z) + sum_i y^i - t] from eq (61) */
        case 2: {
            size_t i;
            secp256k1_scalar yn;
            secp256k1_scalar twosum;
            secp256k1_scalar tmp;

            secp256k1_scalar_clear(&twosum);
            secp256k1_scalar_clear(&yn);
            secp256k1_scalar_set_int(&tmp, 1);

            secp256k1_scalar_sqr(&ctx->zsq, &ctx->z);  /* need to re-set this */
            secp256k1_scalar_negate(sc, &ctx->zsq);  /* -z^2 */
            secp256k1_scalar_add(sc, sc, &ctx->z);   /* z - z^2 */

            for (i = 0; i < ctx->n_commits; i++) {
                const size_t nbits = ctx->n / ctx->n_commits;
                secp256k1_scalar negzn;
                secp256k1_scalar twon;
                size_t j;

                secp256k1_scalar_clear(&twon);
                for (j = 0; j < nbits; j++) {

                    secp256k1_scalar_mul(&yn, &yn, &ctx->y);
                    secp256k1_scalar_add(&twon, &twon, &twon);

                    secp256k1_scalar_add(&yn, &yn, &tmp);
                    secp256k1_scalar_add(&twon, &twon, &tmp);
                }

                secp256k1_scalar_mul(&negzn, &ctx->zsq, &ctx->negz);
                for (j = 0; j < i; j++) {
                    secp256k1_scalar_mul(&negzn, &negzn, &ctx->z);
                }
                secp256k1_scalar_mul(&twon, &twon, &negzn);
                secp256k1_scalar_add(&twosum, &twosum, &twon);
            }  /* yn = 1 + y + ... + y^(n-1); twosum = (z^3 + ... + z^{2 + n_commits})(1 + 2 + ... + 2^(n-1)) */


            secp256k1_scalar_mul(sc, sc, &yn);    /* (z - z^2)(1 + ... + y^(n-1)) */
            secp256k1_scalar_add(sc, sc, &twosum);  /* (z - z^2)(1 + ... + y^(n-1)) - z^3(1 + ... + 2^(n-1)) */
            secp256k1_scalar_negate(&tmp, &ctx->t);
            secp256k1_scalar_add(sc, sc, &tmp);    /* (z - z^2)(1 + ... + y^n) - z^3(1 + ... + 2^n) - t */
            secp256k1_scalar_mul(sc, sc, &ctx->randomizer61);
            *pt = *ctx->asset;
            break;
        }
        /* T1^x in eq (61) */
        case 3:
            secp256k1_scalar_mul(sc, &ctx->x, &ctx->randomizer61);
            *pt = ctx->t1;
            break;
        /* T2^x^2 in eq (61) */
        case 4:
            secp256k1_scalar_sqr(sc, &ctx->x);
            secp256k1_scalar_mul(sc, sc, &ctx->randomizer61);
            *pt = ctx->t2;
            break;
        /* V^z^2 in eq (61) */
        default:
            VERIFY_CHECK(ctx->count < 5 + ctx->n_commits);

            secp256k1_scalar_mul(sc, &ctx->zsq, &ctx->randomizer61);
            secp256k1_scalar_mul(&ctx->zsq, &ctx->zsq, &ctx->z);
            *pt = ctx->commit[ctx->count - 5];
            break;
        }
        secp256k1_scalar_mul(sc, sc, randomizer);
        ctx->count++;
    }
    return 1;
}

static int secp256k1_bulletproof_rangeproof_verify_impl(const secp256k1_context *ctx, const secp256k1_ecmult_context *ecmult_ctx, secp256k1_scratch *scratch, const unsigned char **proof, size_t *plen, size_t n_proofs, size_t nbits, const secp256k1_ge **commitp, size_t n_commits, const secp256k1_ge *genp, const secp256k1_ge *geng, const secp256k1_ge *genh, const unsigned char *extra_commit, size_t extra_commit_len) {
    const size_t depth = secp256k1_ceil_lg(nbits * n_commits);
    secp256k1_bulletproof_vfy_ecmult_context ecmult_data[MAX_BATCH_QTY];
    secp256k1_bulletproof_innerproduct_context innp_ctx[MAX_BATCH_QTY];
    size_t i;

    /* sanity-check input */
    if (POPCOUNT(nbits) != 1 || nbits > MAX_NBITS) {
        return 0;
    }
    for (i = 0; i < n_proofs; i++) {
        if (plen[i] != (9 + 2*depth) * 32 + (4 + 2*depth + 7) / 8) {
            return 0;
        }
        if (depth > SECP256K1_BULLETPROOF_MAX_DEPTH || plen[i] > SECP256K1_BULLETPROOF_MAX_PROOF) {
            return 0;
        }
    }

    for (i = 0; i < n_proofs; i++) {
        secp256k1_sha256 sha256;
        unsigned char commit[32] = {0};
        unsigned char randomizer61[32] = {0};  /* randomizer for eq (61) so we can add it to eq (62) to save a separate multiexp */
        secp256k1_scalar taux, mu, a, b;
        secp256k1_ge age, sge;
        int overflow;
        size_t j;

        /* Commit to all input data: pedersen commit, asset generator, extra_commit */
        for (j = 0; j < n_commits; j++) {
            secp256k1_bulletproof_update_commit(commit, &commitp[i][j], genp);
        }
        if (extra_commit != NULL) {
            secp256k1_sha256_initialize(&sha256);
            secp256k1_sha256_write(&sha256, commit, 32);
            secp256k1_sha256_write(&sha256, extra_commit, extra_commit_len);
            secp256k1_sha256_finalize(&sha256, commit);
        }

        /* Compute y, z, x */
        secp256k1_bulletproof_deserialize_point(&age, &proof[i][160], 0, 4 + 2*depth);
        secp256k1_bulletproof_deserialize_point(&sge, &proof[i][160], 1, 4 + 2*depth);

        secp256k1_bulletproof_update_commit(commit, &age, &sge);
        secp256k1_scalar_set_b32(&ecmult_data[i].y, commit, &overflow);
        if (overflow || secp256k1_scalar_is_zero(&ecmult_data[i].y)) {
            return 0;
        }
        secp256k1_bulletproof_update_commit(commit, &age, &sge);
        secp256k1_scalar_set_b32(&ecmult_data[i].z, commit, &overflow);
        if (overflow || secp256k1_scalar_is_zero(&ecmult_data[i].z)) {
            return 0;
        }

        secp256k1_bulletproof_deserialize_point(&ecmult_data[i].t1, &proof[i][160], 2, 4 + 2*depth);
        secp256k1_bulletproof_deserialize_point(&ecmult_data[i].t2, &proof[i][160], 3, 4 + 2*depth);

        secp256k1_bulletproof_update_commit(commit, &ecmult_data[i].t1, &ecmult_data[i].t2);
        secp256k1_scalar_set_b32(&ecmult_data[i].x, commit, &overflow);
        if (overflow || secp256k1_scalar_is_zero(&ecmult_data[i].x)) {
            return 0;
        }

        /* compute exponent offsets */
        ecmult_data[i].mul_by = ecmult_data[i].y;
        secp256k1_scalar_inverse_var(&ecmult_data[i].yinv, &ecmult_data[i].y);  /* TODO somehow batch this w the inner-product argument inverse */
        secp256k1_scalar_sqr(&ecmult_data[i].zsq, &ecmult_data[i].z);
        secp256k1_scalar_negate(&ecmult_data[i].negz, &ecmult_data[i].z);

        /* Update commit with remaining data for the inner product prof */
        secp256k1_sha256_initialize(&sha256);
        secp256k1_sha256_write(&sha256, commit, 32);
        secp256k1_sha256_write(&sha256, &proof[i][0], 96);
        secp256k1_sha256_finalize(&sha256, commit);

        secp256k1_sha256_initialize(&sha256);
        secp256k1_sha256_write(&sha256, commit, 32);
        secp256k1_sha256_finalize(&sha256, randomizer61);
        secp256k1_scalar_set_b32(&ecmult_data[i].randomizer61, randomizer61, &overflow);
        if (overflow || secp256k1_scalar_is_zero(&ecmult_data[i].randomizer61)) {
            return 0;
        }

        /* Deserialize everything else */
        secp256k1_scalar_set_b32(&ecmult_data[i].t, &proof[i][0], &overflow);
        if (overflow || secp256k1_scalar_is_zero(&ecmult_data[i].t)) {
            return 0;
        }
        secp256k1_scalar_set_b32(&taux, &proof[i][32], &overflow);
        if (overflow || secp256k1_scalar_is_zero(&taux)) {
            return 0;
        }
        secp256k1_scalar_set_b32(&mu, &proof[i][64], &overflow);
        if (overflow || secp256k1_scalar_is_zero(&mu)) {
            return 0;
        }
        secp256k1_scalar_set_b32(&a, &proof[i][96], &overflow);
        if (overflow || secp256k1_scalar_is_zero(&a)) {
    	    return 0;
        }
        secp256k1_scalar_set_b32(&b, &proof[i][128], &overflow);
        if (overflow || secp256k1_scalar_is_zero(&b)) {
            return 0;
        }

        /* Verify inner product proof */
        ecmult_data[i].a = age;
        ecmult_data[i].s = sge;
        ecmult_data[i].n = nbits * n_commits;
        ecmult_data[i].count = 0;
        ecmult_data[i].asset = genp;
        ecmult_data[i].commit = commitp[i];
        ecmult_data[i].n_commits = n_commits;
        secp256k1_scalar_mul(&taux, &taux, &ecmult_data[i].randomizer61);
        secp256k1_scalar_add(&mu, &mu, &taux);

        innp_ctx[i].serialized_points = &proof[i][160];
        innp_ctx[i].n_extra_ser_points = 4;
        innp_ctx[i].a = a;
        innp_ctx[i].b = b;
        innp_ctx[i].dot = ecmult_data[i].t;
        innp_ctx[i].p_offs = mu;
        memcpy(innp_ctx[i].commit, commit, 32);
        innp_ctx[i].rangeproof_cb = secp256k1_bulletproof_rangeproof_vfy_callback;
        innp_ctx[i].rangeproof_cb_data = (void *) &ecmult_data[i];
        innp_ctx[i].n_extra_rangeproof_points = 5 + n_commits;
    }

    return secp256k1_bulletproof_inner_product_verify_impl(ctx, ecmult_ctx, scratch, geng, genh, 1u << depth, innp_ctx, n_proofs);
}

typedef struct {
    secp256k1_rfc6979_hmac_sha256 rng;
    secp256k1_scalar y;
    secp256k1_scalar z;
    secp256k1_scalar yn;
    secp256k1_scalar z22n;
    const uint64_t *val;
    size_t n_vals;
    size_t nbits;
    size_t count;
} secp256k1_bulletproof_lr_generator;

static void secp256k1_lr_generator_init(secp256k1_bulletproof_lr_generator *generator, const secp256k1_rfc6979_hmac_sha256 *rng, const secp256k1_scalar *y, const secp256k1_scalar *z, size_t nbits, const uint64_t *val, size_t n_vals) {
    memcpy(&generator->rng, rng, sizeof(*rng));
    generator->y = *y;
    generator->z = *z;
    secp256k1_scalar_set_int(&generator->yn, 1);
    generator->nbits = nbits;
    generator->val = val;
    generator->n_vals = n_vals;
    generator->count = 0;
}

static void secp256k1_lr_generator_finalize(secp256k1_bulletproof_lr_generator *generator) {
    secp256k1_rfc6979_hmac_sha256_finalize(&generator->rng);
}

static void secp256k1_lr_generate(secp256k1_bulletproof_lr_generator *generator, secp256k1_scalar *lout, secp256k1_scalar *rout, const secp256k1_scalar *x) {
    const size_t commit_idx = generator->count / generator->nbits;
    const size_t bit_idx = generator->count % generator->nbits;
    const int bit = (generator->val[commit_idx] >> bit_idx) & 1;
    secp256k1_scalar sl, sr;
    secp256k1_scalar negz;

    if (bit_idx == 0) {
        size_t i;
        secp256k1_scalar_sqr(&generator->z22n, &generator->z);
        for (i = 0; i < commit_idx; i++) {
            secp256k1_scalar_mul(&generator->z22n, &generator->z22n, &generator->z);
        }
    }

    secp256k1_bulletproof_genrand_pair(&generator->rng, &sl, &sr);
    secp256k1_scalar_mul(&sl, &sl, x);
    secp256k1_scalar_mul(&sr, &sr, x);

    secp256k1_scalar_set_int(lout, bit);
    secp256k1_scalar_negate(&negz, &generator->z);
    secp256k1_scalar_add(lout, lout, &negz);
    secp256k1_scalar_add(lout, lout, &sl);

    secp256k1_scalar_set_int(rout, 1 - bit);
    secp256k1_scalar_negate(rout, rout);
    secp256k1_scalar_add(rout, rout, &generator->z);
    secp256k1_scalar_add(rout, rout, &sr);
    secp256k1_scalar_mul(rout, rout, &generator->yn);
    secp256k1_scalar_add(rout, rout, &generator->z22n);

    generator->count++;
    secp256k1_scalar_mul(&generator->yn, &generator->yn, &generator->y);
    secp256k1_scalar_add(&generator->z22n, &generator->z22n, &generator->z22n);
}

typedef struct {
    secp256k1_scalar yinv;
    secp256k1_scalar x;
    secp256k1_scalar mul;
    secp256k1_scalar cache;
    secp256k1_bulletproof_lr_generator lr_gen;
    const secp256k1_ge *geng;
    const secp256k1_ge *genh;
} secp256k1_bulletproof_abgh_data;

static int secp256k1_bulletproof_abgh_callback(secp256k1_scalar *sc, secp256k1_ge *pt, size_t idx, void *data) {
    secp256k1_bulletproof_abgh_data *ctx = (secp256k1_bulletproof_abgh_data *) data;
    const int is_g = idx % 2 == 0;

    if (is_g) {
        secp256k1_lr_generate(&ctx->lr_gen, sc, &ctx->cache, &ctx->x);
        *pt = ctx->geng[idx / 2];
    } else {
        /* Map h -> h' (eqn 59) */
        secp256k1_gej genhj;
        secp256k1_ecmult_const(&genhj, &ctx->genh[idx / 2], &ctx->mul, 256);
        *sc = ctx->cache;
        secp256k1_ge_set_gej(pt, &genhj);

        secp256k1_scalar_mul(&ctx->mul, &ctx->mul, &ctx->yinv);
    }

    return 1;
}

/* Proof format: t, tau_x, mu, a, b, A, S, T_1, T_2, {L_i}, {R_i}
 *               5 scalar + [4 + 2log(n)] ge
 *
 * The non-bold `h` in the Bulletproofs paper corresponds to our secp256k1_ge_const_g
 * while the non-bold `g` corresponds to the asset type `genp`.
 */
static int secp256k1_bulletproof_rangeproof_prove_impl(const secp256k1_context *ctx, const secp256k1_ecmult_gen_context *ecmult_gen_ctx, const secp256k1_ecmult_context *ecmult_ctx, secp256k1_scratch *scratch, unsigned char *proof, size_t *plen, const size_t nbits, const uint64_t *value, const secp256k1_scalar *blind, const secp256k1_ge *commitp, size_t n_commits, const secp256k1_ge *genp, const secp256k1_ge *geng, const secp256k1_ge *genh, const unsigned char *nonce, const unsigned char *extra_commit, size_t extra_commit_len) {
    secp256k1_bulletproof_lr_generator lr_gen;
    secp256k1_bulletproof_abgh_data abgh_data;
    secp256k1_scalar zero;
    secp256k1_sha256 sha256;
    const size_t depth = secp256k1_ceil_lg(nbits * n_commits);
    unsigned char commit[32] = {0};
    secp256k1_scalar alpha, rho;
    secp256k1_scalar t, t0, t1, t2;
    secp256k1_scalar tau1, tau2, taux, mu;
    secp256k1_gej tj[2];      /* T_1, T_2 */
    secp256k1_scalar y;
    secp256k1_scalar z, zsq;
    secp256k1_scalar x, xsq;
    secp256k1_scalar tmps;
    secp256k1_gej aj, sj;
    secp256k1_gej tmpj;
    size_t i, j;
    int overflow;
    unsigned char rngseed[64];
    secp256k1_rfc6979_hmac_sha256 rng;
    /* inner product proof variables */
    secp256k1_scalar a, b;
    secp256k1_ge out_pt[4 + 2*SECP256K1_BULLETPROOF_MAX_DEPTH];

    if (POPCOUNT(nbits) != 1 || nbits > MAX_NBITS) {
        return 0;
    }
    for (i = 0; i < n_commits; i++) {
        if (nbits < 64 && value[i] >= (1ull << nbits)) {
            return 0;
        }
    }
    if (*plen < (9 + 2*depth) * 32 + (4 + 2*depth + 7) / 8) {
        return 0;
    } else {
        *plen = (9 + 2*depth) * 32 + (4 + 2*depth + 7) / 8;
    }

    secp256k1_scalar_clear(&zero);

    /* Commit to all input data: pedersen commit, asset generator, extra_commit */
    for (i = 0; i < n_commits; i++) {
        secp256k1_bulletproof_update_commit(commit, &commitp[i], genp); /* TODO be less stupid about this */
    }
    if (extra_commit != NULL) {
        secp256k1_sha256_initialize(&sha256);
        secp256k1_sha256_write(&sha256, commit, 32);
        secp256k1_sha256_write(&sha256, extra_commit, extra_commit_len);
        secp256k1_sha256_finalize(&sha256, commit);
    }

    memcpy(rngseed, nonce, 32);
    memcpy(rngseed+32, commit, 32);
    secp256k1_rfc6979_hmac_sha256_initialize(&rng, rngseed, 64);
    memset(rngseed, 0, 32);
    secp256k1_bulletproof_genrand_pair(&rng, &alpha, &rho);
    secp256k1_bulletproof_genrand_pair(&rng, &tau1, &tau2);

    /* Compute A and S */
    lr_gen.rng = rng;
    secp256k1_ecmult_gen(ecmult_gen_ctx, &aj, &alpha);
    secp256k1_ecmult_gen(ecmult_gen_ctx, &sj, &rho);
    for (i = 0; i < n_commits; i++) {
        for (j = 0; j < nbits; j++) {
            secp256k1_scalar sl, sr;
            size_t al = !!(value[i] & (1ull << j));
            secp256k1_ge aterm = genh[i * nbits + j];
            secp256k1_ge sterm;
            secp256k1_gej stermj;

            secp256k1_bulletproof_genrand_pair(&lr_gen.rng, &sl, &sr);

            secp256k1_ge_neg(&aterm, &aterm);
            secp256k1_fe_cmov(&aterm.x, &geng[i * nbits + j].x, al);
            secp256k1_fe_cmov(&aterm.y, &geng[i * nbits + j].y, al);

            secp256k1_gej_add_ge(&aj, &aj, &aterm);

            secp256k1_ecmult_const(&stermj, &geng[i * nbits + j], &sl, 256);
            secp256k1_ge_set_gej(&sterm, &stermj);
            secp256k1_gej_add_ge(&sj, &sj, &sterm);
            secp256k1_ecmult_const(&stermj, &genh[i * nbits + j], &sr, 256);
            secp256k1_ge_set_gej(&sterm, &stermj);
            secp256k1_gej_add_ge(&sj, &sj, &sterm);
        }
    }
    secp256k1_rfc6979_hmac_sha256_finalize(&lr_gen.rng);

    /* get challenges y and z */
    secp256k1_ge_set_gej(&out_pt[0], &aj);
    secp256k1_ge_set_gej(&out_pt[1], &sj);

    secp256k1_bulletproof_update_commit(commit, &out_pt[0], &out_pt[1]);
    secp256k1_scalar_set_b32(&y, commit, &overflow);
    if (overflow || secp256k1_scalar_is_zero(&y)) {
        return 0;
    }
    secp256k1_bulletproof_update_commit(commit, &out_pt[0], &out_pt[1]); /* TODO rehashing A and S to get a second challenge is overkill */
    secp256k1_scalar_set_b32(&z, commit, &overflow);
    if (overflow || secp256k1_scalar_is_zero(&z)) {
        return 0;
    }
    secp256k1_scalar_sqr(&zsq, &z);

    /* Compute coefficients t0, t1, t2 of the <l, r> polynomial */
    /* t0 = l(0) dot r(0) */
    secp256k1_lr_generator_init(&lr_gen, &rng, &y, &z, nbits, value, n_commits);
    secp256k1_scalar_clear(&t0);
    for (i = 0; i < nbits * n_commits; i++) {
        secp256k1_scalar l, r;
        secp256k1_lr_generate(&lr_gen, &l, &r, &zero);
        secp256k1_scalar_mul(&l, &l, &r);
        secp256k1_scalar_add(&t0, &t0, &l);
    }
    secp256k1_lr_generator_finalize(&lr_gen);

    /* A = t0 + t1 + t2 = l(1) dot r(1) */
    secp256k1_lr_generator_init(&lr_gen, &rng, &y, &z, nbits, value, n_commits);
    secp256k1_scalar_clear(&t1);
    for (i = 0; i < nbits * n_commits; i++) {
        secp256k1_scalar one;
        secp256k1_scalar l, r;
        secp256k1_scalar_set_int(&one, 1);
        secp256k1_lr_generate(&lr_gen, &l, &r, &one);
        secp256k1_scalar_mul(&l, &l, &r);
        secp256k1_scalar_add(&t1, &t1, &l);
    }
    secp256k1_lr_generator_finalize(&lr_gen);

    /* B = t0 - t1 + t2 = l(-1) dot r(-1) */
    secp256k1_lr_generator_init(&lr_gen, &rng, &y, &z, nbits, value, n_commits);
    secp256k1_scalar_clear(&t2);
    for (i = 0; i < nbits * n_commits; i++) {
        secp256k1_scalar negone;
        secp256k1_scalar l, r;
        secp256k1_scalar_set_int(&negone, 1);
        secp256k1_scalar_negate(&negone, &negone);
        secp256k1_lr_generate(&lr_gen, &l, &r, &negone);
        secp256k1_scalar_mul(&l, &l, &r);
        secp256k1_scalar_add(&t2, &t2, &l);
    }
    secp256k1_lr_generator_finalize(&lr_gen);

    /* t1 = (A - B)/2 */
    secp256k1_scalar_set_int(&tmps, 2);
    secp256k1_scalar_inverse_var(&tmps, &tmps);
    secp256k1_scalar_negate(&t2, &t2);
    secp256k1_scalar_add(&t1, &t1, &t2);
    secp256k1_scalar_mul(&t1, &t1, &tmps);

    /* t2 = -(-B + t0) + t1 */
    secp256k1_scalar_add(&t2, &t2, &t0);
    secp256k1_scalar_negate(&t2, &t2);
    secp256k1_scalar_add(&t2, &t2, &t1);

    /* Compute Ti = t_i*A + tau_i*G for i = 1,2 */
    secp256k1_gej_set_ge(&tmpj, genp);
    secp256k1_ecmult(ecmult_ctx, &tj[0], &tmpj, &t1, &tau1);
    secp256k1_ecmult(ecmult_ctx, &tj[1], &tmpj, &t2, &tau2);
    secp256k1_ge_set_gej(&out_pt[2], &tj[0]);
    secp256k1_ge_set_gej(&out_pt[3], &tj[1]);

    /* get challenge x */
    secp256k1_bulletproof_update_commit(commit, &out_pt[2], &out_pt[3]);
    secp256k1_scalar_set_b32(&x, commit, &overflow);
    if (overflow || secp256k1_scalar_is_zero(&x)) {
        return 0;
    }
    secp256k1_scalar_sqr(&xsq, &x);

    /* compute tau_x, mu and t */
    secp256k1_scalar_mul(&taux, &tau1, &x);
    secp256k1_scalar_mul(&tmps, &tau2, &xsq);
    secp256k1_scalar_add(&taux, &taux, &tmps);
    for (i = 0; i < n_commits; i++) {
        secp256k1_scalar_mul(&tmps, &zsq, &blind[i]);
        secp256k1_scalar_add(&taux, &taux, &tmps);
        secp256k1_scalar_mul(&zsq, &zsq, &z);
    }

    secp256k1_scalar_mul(&mu, &rho, &x);
    secp256k1_scalar_add(&mu, &mu, &alpha);

    secp256k1_scalar_mul(&tmps, &t2, &xsq);
    secp256k1_scalar_mul(&t, &t1, &x);
    secp256k1_scalar_add(&t, &t, &tmps);
    secp256k1_scalar_add(&t, &t, &t0);

    /* Negate taux and mu so the verifier doesn't have to */
    secp256k1_scalar_negate(&taux, &taux);
    secp256k1_scalar_negate(&mu, &mu);

    /* Mix these scalars into the hash so the input to the inner product proof is fixed */
    secp256k1_sha256_initialize(&sha256);
    secp256k1_sha256_write(&sha256, commit, 32);
    secp256k1_scalar_get_b32(commit, &t);
    secp256k1_sha256_write(&sha256, commit, 32);
    secp256k1_scalar_get_b32(commit, &taux);
    secp256k1_sha256_write(&sha256, commit, 32);
    secp256k1_scalar_get_b32(commit, &mu);
    secp256k1_sha256_write(&sha256, commit, 32);
    secp256k1_sha256_finalize(&sha256, commit);

    /* Compute l and r, do inner product proof */
    secp256k1_scalar_inverse_var(&abgh_data.yinv, &y);
    abgh_data.x = x;
    secp256k1_scalar_set_int(&abgh_data.mul, 1);
    abgh_data.geng = geng;
    abgh_data.genh = genh;
    secp256k1_lr_generator_init(&abgh_data.lr_gen, &rng, &y, &z, nbits, value, n_commits);
    secp256k1_bulletproof_inner_product_prove_impl(ctx, ecmult_ctx, scratch, &a, &b, &out_pt[4], &out_pt[4 + depth], 1 << depth, secp256k1_bulletproof_abgh_callback, (void *) &abgh_data, commit);
    secp256k1_lr_generator_finalize(&abgh_data.lr_gen);

    /* Encode everything */
    secp256k1_scalar_get_b32(&proof[0], &t);
    secp256k1_scalar_get_b32(&proof[32], &taux);
    secp256k1_scalar_get_b32(&proof[64], &mu);
    secp256k1_scalar_get_b32(&proof[96], &a);
    secp256k1_scalar_get_b32(&proof[128], &b);
    secp256k1_bulletproof_serialize_points(&proof[160], out_pt, 4 + 2*depth);

    secp256k1_rfc6979_hmac_sha256_finalize(&rng);

    return 1;
}
#endif
