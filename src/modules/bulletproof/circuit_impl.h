/**********************************************************************
 * Copyright (c) 2018 Andrew Poelstra                                 *
 * Distributed under the MIT software license, see the accompanying   *
 * file COPYING or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#ifndef SECP256K1_MODULE_BULLETPROOF_CIRCUIT_IMPL
#define SECP256K1_MODULE_BULLETPROOF_CIRCUIT_IMPL

#include "modules/bulletproof/inner_product_impl.h"
#include "modules/bulletproof/util.h"
#include "group.h"

#include <stdlib.h>

typedef struct {
    secp256k1_rfc6979_hmac_sha256 rng;
    secp256k1_scalar y;
    secp256k1_scalar yinv;
    secp256k1_scalar z;
    secp256k1_scalar yn;
    secp256k1_scalar yinvn;
    secp256k1_scalar sl, sr;
    const secp256k1_scalar *al;
    const secp256k1_scalar *ar;
    const secp256k1_scalar *ao;
    const secp256k1_bulletproof_circuit *circ;
    size_t i;
} secp256k1_bulletproof_circuit_lr_generator;

static void secp256k1_bulletproof_circuit_lr_generator_update(secp256k1_bulletproof_circuit_lr_generator *generator);

static void secp256k1_bulletproof_circuit_lr_generator_init(secp256k1_bulletproof_circuit_lr_generator *generator, const secp256k1_rfc6979_hmac_sha256 *rng, const secp256k1_scalar *y, const secp256k1_scalar *z, const secp256k1_scalar *al, const secp256k1_scalar *ar, const secp256k1_scalar *ao, const secp256k1_bulletproof_circuit *circ) {
    memcpy(&generator->rng, rng, sizeof(*rng));
    generator->y = *y;
    generator->z = *z;
    secp256k1_scalar_set_int(&generator->yn, 1);
    secp256k1_scalar_set_int(&generator->yinvn, 1);
    secp256k1_scalar_inverse_var(&generator->yinv, y);
    generator->al = al;
    generator->ar = ar;
    generator->ao = ao;
    generator->circ = circ;
    generator->i = 0;
    secp256k1_bulletproof_genrand_pair(&generator->rng, &generator->sl, &generator->sr);
}

static void secp256k1_bulletproof_circuit_lr_generator_finalize(secp256k1_bulletproof_circuit_lr_generator *generator) {
    secp256k1_rfc6979_hmac_sha256_finalize(&generator->rng);
}

static void secp256k1_bulletproof_circuit_lr_generator_update(secp256k1_bulletproof_circuit_lr_generator *generator) {
    secp256k1_bulletproof_genrand_pair(&generator->rng, &generator->sl, &generator->sr);

    secp256k1_scalar_mul(&generator->yn, &generator->yn, &generator->y);
    secp256k1_scalar_mul(&generator->yinvn, &generator->yinvn, &generator->yinv);
    generator->i += 1;
}

static void secp256k1_bulletproof_circuit_lr_generate(const secp256k1_bulletproof_circuit_lr_generator *generator, secp256k1_scalar *lout, secp256k1_scalar *rout, const secp256k1_scalar *x) {
    secp256k1_scalar negone;
    secp256k1_scalar sl, sr;
    secp256k1_scalar x2, x3;
    secp256k1_scalar tmp;
    const size_t i = generator->i;

    secp256k1_scalar_set_int(&negone, 1);
    secp256k1_scalar_negate(&negone, &negone);

    secp256k1_scalar_sqr(&x2, x);
    secp256k1_scalar_mul(&x3, &x2, x);

    secp256k1_scalar_mul(&sl, &generator->sl, &x3);
    secp256k1_scalar_mul(&sr, &generator->sr, &x3);
    secp256k1_scalar_mul(&sr, &sr, &generator->yn);
    secp256k1_scalar_mul(lout, &generator->ao[i], x); /* l = a_O * x */
    secp256k1_scalar_add(lout, lout, &generator->al[i]); /* l = a_O * x + a_L */
    secp256k1_scalar_mul(rout, &generator->ar[i], x); /* r = a_R * X */
    secp256k1_scalar_add(rout, rout, &negone); /* r = a_R * X - 1 */
    secp256k1_scalar_mul(rout, rout, &generator->yn); /* r = y^n * a_R * x - y^n */

    secp256k1_scalar_mul(&tmp, &generator->circ->wr[i].cached_sum, &generator->yinvn);
    secp256k1_scalar_add(lout, lout, &tmp);
    /* ^  l = a_O * x + a_L + y^-n (z^Q . W_R)  */

    secp256k1_scalar_mul(&tmp, &generator->circ->wl[i].cached_sum, x);
    secp256k1_scalar_add(rout, rout, &tmp);

    secp256k1_scalar_add(rout, rout, &generator->circ->wo[i].cached_sum);
    /* ^  r = y^n * a_R * x - y^n + z^Q . (xW_L + W_O) */

    secp256k1_scalar_mul(lout, lout, x); /* l = a_O * x^2 + (a_L + y^-n (z^Q . W_R)) * x  */

    secp256k1_scalar_add(lout, lout, &sl); /* add s_L * x^3 */
    secp256k1_scalar_add(rout, rout, &sr); /* add s_R * x^3 */
}

typedef struct {
    secp256k1_scalar yinv;
    secp256k1_scalar x;
    secp256k1_scalar mul;
    secp256k1_scalar cache;
    secp256k1_bulletproof_circuit_lr_generator lr_gen;
    const secp256k1_ge *geng;
    const secp256k1_ge *genh;
} secp256k1_bulletproof_circuit_abgh_data;

static int secp256k1_bulletproof_circuit_abgh_callback(secp256k1_scalar *sc, secp256k1_ge *pt, size_t idx, void *data) {
    secp256k1_bulletproof_circuit_abgh_data *ctx = (secp256k1_bulletproof_circuit_abgh_data *) data;
    const int is_g = idx % 2 == 0;

    if (is_g) {
        secp256k1_bulletproof_circuit_lr_generate(&ctx->lr_gen, sc, &ctx->cache, &ctx->x);
        secp256k1_bulletproof_circuit_lr_generator_update(&ctx->lr_gen);
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

/* Proof format:
 *
 * Serialized scalars (32 bytes) t, tau_x, mu, a, b
 * Serialized points (bit array of parity followed by 32 bytes): A_I, A_O, S, T_1, T_3, T_4, T_5, T_6, [inner product proof points]
 */
static int secp256k1_bulletproof_relation66_prove_impl(const secp256k1_context *ctx, const secp256k1_ecmult_context *ecmult_ctx, secp256k1_scratch *scratch, unsigned char *proof, size_t *plen, const secp256k1_scalar *al, const secp256k1_scalar *ar, const secp256k1_scalar *ao, size_t na, const secp256k1_ge *commitp, const secp256k1_scalar *blinds, size_t nc, const secp256k1_ge *genp, secp256k1_bulletproof_circuit *circ, const secp256k1_ge *geng, const secp256k1_ge *genh, const unsigned char *nonce, const unsigned char *extra_commit, size_t extra_commit_len) {
    secp256k1_bulletproof_circuit_lr_generator lr_gen;
    secp256k1_bulletproof_circuit_abgh_data abgh_data;
    secp256k1_rfc6979_hmac_sha256 rng;
    secp256k1_sha256 sha256;
    const size_t depth = CTZ(circ->n_gates);
    unsigned char commit[32] = {0};
    secp256k1_scalar zero, one, onehalf, onethird, twothirds, fourthirds, eightthirds;
    secp256k1_scalar alpha, beta, rho, mu;
    secp256k1_scalar tau1, tau3, tau4, tau5, tau6, taux; /* tau2 missing on purpose */
    secp256k1_scalar t[7];  /* t[1..6] are coefficients; t[0] is the polynomial evaluated at x */
    secp256k1_scalar tauv;  /* <z, WV*gamma> term in eq (73) */
    secp256k1_scalar a, b;
    secp256k1_scalar x, xn, y, z;
    secp256k1_scalar tmp;
    secp256k1_gej aij, aoj, sj;
    secp256k1_ge tmpge;
    secp256k1_ge out_pt[8 + 2*SECP256K1_BULLETPROOF_MAX_DEPTH];
    int overflow;
    size_t i;

    if (na != circ->n_gates || nc != circ->n_commits) {
        return 0;
    }
    if (*plen < (13 + 2*depth) * 32 + (8 + 2*depth + 7) / 8) {
        return 0;
    } else {
        *plen = (13 + 2*depth) * 32 + (8 + 2*depth + 7) / 8;
    }

    /* Commit to all input data */
    secp256k1_bulletproof_update_commit_n(commit, commitp, nc);
    secp256k1_bulletproof_update_commit_n(commit, genp, 1);
    /* TODO commit to circuit */
    if (extra_commit != NULL) {
        secp256k1_sha256_initialize(&sha256);
        secp256k1_sha256_write(&sha256, commit, 32);
        secp256k1_sha256_write(&sha256, extra_commit, extra_commit_len);
        secp256k1_sha256_finalize(&sha256, commit);
    }

    /* Setup, generate randomness */
    secp256k1_scalar_set_int(&zero, 0);
    secp256k1_scalar_set_int(&one, 1);
    secp256k1_scalar_set_int(&tmp, 6);
    secp256k1_scalar_inverse_var(&tmp, &tmp);
    secp256k1_scalar_set_int(&onethird, 2);
    secp256k1_scalar_mul(&onethird, &onethird, &tmp);
    secp256k1_scalar_set_int(&onehalf, 3);
    secp256k1_scalar_mul(&onehalf, &onehalf, &tmp);
    secp256k1_scalar_add(&twothirds, &onethird, &onethird);
    secp256k1_scalar_add(&fourthirds, &twothirds, &twothirds);
    secp256k1_scalar_add(&eightthirds, &fourthirds, &fourthirds);

    secp256k1_rfc6979_hmac_sha256_initialize(&rng, nonce, 32); /* todo initialize from secret input */
    secp256k1_bulletproof_genrand_pair(&rng, &alpha, &beta);
    secp256k1_bulletproof_genrand_pair(&rng, &rho, &tau1);
    secp256k1_bulletproof_genrand_pair(&rng, &tau3, &tau4); /* t2 will be generated deterministically */
    secp256k1_bulletproof_genrand_pair(&rng, &tau5, &tau6);

    /* Compute A_I, A_O, S */
    lr_gen.rng = rng;

    secp256k1_bulletproof_vector_commit(&aij, al, geng, circ->n_gates, &alpha, NULL);
    secp256k1_ge_set_gej(&tmpge, &aij);
    secp256k1_bulletproof_vector_commit(&aij, ar, genh, circ->n_gates, NULL, NULL);
    secp256k1_gej_add_ge(&aij, &aij, &tmpge);

    secp256k1_bulletproof_vector_commit(&aoj, ao, geng, circ->n_gates, &beta, NULL);

    secp256k1_ecmult_const(&sj, &secp256k1_ge_const_g, &rho, 256);
    for (i = 0; i < circ->n_gates; i++) {
        secp256k1_scalar sl, sr;
        secp256k1_gej termj;
        secp256k1_ge term;

        secp256k1_bulletproof_genrand_pair(&lr_gen.rng, &sl, &sr);

        secp256k1_ecmult_const(&termj, &geng[i], &sl, 256);
        secp256k1_ge_set_gej(&term, &termj);
        secp256k1_gej_add_ge(&sj, &sj, &term);
        secp256k1_ecmult_const(&termj, &genh[i], &sr, 256);
        secp256k1_ge_set_gej(&term, &termj);
        secp256k1_gej_add_ge(&sj, &sj, &term);
    }
    secp256k1_rfc6979_hmac_sha256_finalize(&lr_gen.rng);

    /* get challenges y and z */
    secp256k1_ge_set_gej(&out_pt[0], &aij);
    secp256k1_ge_set_gej(&out_pt[1], &aoj);
    secp256k1_ge_set_gej(&out_pt[2], &sj);

    secp256k1_bulletproof_update_commit_n(commit, &out_pt[0], 3);
    secp256k1_scalar_set_b32(&y, commit, &overflow);
    if (overflow || secp256k1_scalar_is_zero(&y)) {
        return 0;
    }
    secp256k1_bulletproof_update_commit_n(commit, NULL, 0);
    secp256k1_scalar_set_b32(&z, commit, &overflow);
    if (overflow || secp256k1_scalar_is_zero(&z)) {
        return 0;
    }

    secp256k1_circuit_compress(circ, &z);

    /* Compute coefficients t[1..6] */

    /* Start by computing each entry of l and r, as
     *   l = l1 * X          + l2 * X^2 + l3 * X^3
     *   r = r0     + r1 * X            + r3 * X^3
     * and observe that
     *   t1 = <l1, r0>
     *   t2 = <l1, r1> + <l2, r0>
     *   t3 = <l2, r1> + <l3, r0>
     *   t4 = <l3, r1> + <l1, r3>
     *   t5 = <l2, r3>
     *   t6 = <l3, r3>
     * So we compute these terms and add them to t1,t3,etc as running sums.
     */

    for (i = 0; i < 6; i++) {
        secp256k1_scalar_clear(&t[i + 1]);
    }
    secp256k1_bulletproof_circuit_lr_generator_init(&lr_gen, &rng, &y, &z, al, ar, ao, circ);
    for (i = 0; i < circ->n_gates; i++) {
        secp256k1_scalar lone, rone;
        secp256k1_scalar lhalf, rhalf;
        secp256k1_scalar ltmp, rtmp;
        secp256k1_scalar l1, l3;     /* l coefficients -- l2 = a_O[i], l0 = 0 */
        secp256k1_scalar r0, r1, r3; /* r coefficients -- r2 = 0 */

        secp256k1_bulletproof_circuit_lr_generate(&lr_gen, &lone, &r0, &zero);
        secp256k1_bulletproof_circuit_lr_generate(&lr_gen, &lone, &rone, &one);
        secp256k1_bulletproof_circuit_lr_generate(&lr_gen, &lhalf, &rhalf, &onehalf);
        secp256k1_bulletproof_circuit_lr_generator_update(&lr_gen);

        secp256k1_scalar_add(&l1, &lone, &ao[i]); /* l1 = l(1) + l2 + l0 */
        secp256k1_scalar_add(&r1, &rone, &r0);
        secp256k1_scalar_mul(&l1, &l1, &onethird); /* l1 = 1/3 l(1) + 1/3 l2 + 1/3 l0 */
        secp256k1_scalar_mul(&r1, &r1, &onethird);
        secp256k1_scalar_add(&r1, &r1, &r0);
        secp256k1_scalar_add(&r1, &r1, &r0); /* l1 = 1/3 l(1) + 1/3 l2 + 7/3 l0 */
        secp256k1_scalar_negate(&l1, &l1); /* l1 = -1/3 l(1) - 1/3 l2 - 7/3 l0 */
        secp256k1_scalar_negate(&r1, &r1);

        secp256k1_scalar_mul(&ltmp, &lhalf, &eightthirds);
        secp256k1_scalar_mul(&rtmp, &rhalf, &eightthirds);
        secp256k1_scalar_add(&l1, &l1, &ltmp); /* l1 = -1/3 l(1) + 8/3 l(1/2) - 1/3 l2 - 7/3 l0 */
        secp256k1_scalar_add(&r1, &r1, &rtmp);

        secp256k1_scalar_mul(&l3, &ao[i], &twothirds); /* l3 = 2/3 l2 */
        secp256k1_scalar_add(&l3, &l3, &ltmp); /* l3 = 2/3 l2 + 8/3 l(1/2) */
        secp256k1_scalar_negate(&l3, &l3); /* l3 = -2/3 l2 - 8/3 l(1/2) */
        secp256k1_scalar_negate(&r3, &rtmp);

        secp256k1_scalar_mul(&ltmp, &lone, &fourthirds);
        secp256k1_scalar_add(&rtmp, &r0, &rone);
        secp256k1_scalar_mul(&rtmp, &rtmp, &fourthirds);
        secp256k1_scalar_add(&l3, &l3, &ltmp); /* l3 = -2/3 l2 - 8/3 l(1/2) + 4/3 l(1) + 4/3 l0 */
        secp256k1_scalar_add(&r3, &r3, &rtmp);

        /* Now that we have the individual coefficients, compute the dot product */
        secp256k1_scalar_mul(&ltmp, &l1, &r0);
        secp256k1_scalar_add(&t[1], &t[1], &ltmp);

        secp256k1_scalar_mul(&ltmp, &l1, &r1);
        secp256k1_scalar_add(&t[2], &t[2], &ltmp);
        secp256k1_scalar_mul(&ltmp, &ao[i], &r0);
        secp256k1_scalar_add(&t[2], &t[2], &ltmp);

        secp256k1_scalar_mul(&ltmp, &ao[i], &r1);
        secp256k1_scalar_add(&t[3], &t[3], &ltmp);
        secp256k1_scalar_mul(&ltmp, &l3, &r0);
        secp256k1_scalar_add(&t[3], &t[3], &ltmp);

        secp256k1_scalar_mul(&ltmp, &l3, &r1);
        secp256k1_scalar_add(&t[4], &t[4], &ltmp);
        secp256k1_scalar_mul(&ltmp, &l1, &r3);
        secp256k1_scalar_add(&t[4], &t[4], &ltmp);

        secp256k1_scalar_mul(&ltmp, &ao[i], &r3);
        secp256k1_scalar_add(&t[5], &t[5], &ltmp);

        secp256k1_scalar_mul(&ltmp, &l3, &r3);
        secp256k1_scalar_add(&t[6], &t[6], &ltmp);
    }
    secp256k1_bulletproof_circuit_lr_generator_finalize(&lr_gen);

    /* Compute T1, T3, T4, T5, T6 */
    secp256k1_bulletproof_vector_commit(&aij, &t[1], genp, 1, &tau1, NULL);
    secp256k1_ge_set_gej(&out_pt[3], &aij);

    secp256k1_bulletproof_vector_commit(&aij, &t[3], genp, 1, &tau3, NULL);
    secp256k1_ge_set_gej(&out_pt[4], &aij);

    secp256k1_bulletproof_vector_commit(&aij, &t[4], genp, 1, &tau4, NULL);
    secp256k1_ge_set_gej(&out_pt[5], &aij);

    secp256k1_bulletproof_vector_commit(&aij, &t[5], genp, 1, &tau5, NULL);
    secp256k1_ge_set_gej(&out_pt[6], &aij);

    secp256k1_bulletproof_vector_commit(&aij, &t[6], genp, 1, &tau6, NULL);
    secp256k1_ge_set_gej(&out_pt[7], &aij);

    /* Compute x, tau_x, mu and t */
    secp256k1_bulletproof_update_commit_n(commit, &out_pt[3], 5);
    secp256k1_scalar_set_b32(&x, commit, &overflow);
    if (overflow || secp256k1_scalar_is_zero(&x)) {
        return 0;
    }

    secp256k1_scalar_mul(&alpha, &alpha, &x);
    secp256k1_scalar_mul(&tau1, &tau1, &x);
    secp256k1_scalar_mul(&t[1], &t[1], &x);

    secp256k1_scalar_sqr(&xn, &x);
    secp256k1_scalar_mul(&beta, &beta, &xn);
    secp256k1_scalar_clear(&tauv);
    for (i = 0; i < circ->n_commits; i++) {
        secp256k1_scalar zwv;
        secp256k1_scalar_mul(&zwv, &circ->wv[i].cached_sum, &blinds[i]);
        secp256k1_scalar_add(&tauv, &tauv, &zwv);
    }
    secp256k1_scalar_mul(&tauv, &tauv, &xn);
    secp256k1_scalar_mul(&t[2], &t[2], &xn);

    secp256k1_scalar_mul(&xn, &xn, &x);
    secp256k1_scalar_mul(&rho, &rho, &xn);
    secp256k1_scalar_mul(&tau3, &tau3, &xn);
    secp256k1_scalar_mul(&t[3], &t[3], &xn);

    secp256k1_scalar_mul(&xn, &xn, &x);
    secp256k1_scalar_mul(&tau4, &tau4, &xn);
    secp256k1_scalar_mul(&t[4], &t[4], &xn);

    secp256k1_scalar_mul(&xn, &xn, &x);
    secp256k1_scalar_mul(&tau5, &tau5, &xn);
    secp256k1_scalar_mul(&t[5], &t[5], &xn);

    secp256k1_scalar_mul(&xn, &xn, &x);
    secp256k1_scalar_mul(&tau6, &tau6, &xn);
    secp256k1_scalar_mul(&t[6], &t[6], &xn);

    secp256k1_scalar_add(&taux, &tau1, &tauv);
    secp256k1_scalar_add(&taux, &taux, &tau3);
    secp256k1_scalar_add(&taux, &taux, &tau4);
    secp256k1_scalar_add(&taux, &taux, &tau5);
    secp256k1_scalar_add(&taux, &taux, &tau6);

    secp256k1_scalar_add(&t[0], &t[1], &t[2]);
    secp256k1_scalar_add(&t[0], &t[0], &t[3]);
    secp256k1_scalar_add(&t[0], &t[0], &t[4]);
    secp256k1_scalar_add(&t[0], &t[0], &t[5]);
    secp256k1_scalar_add(&t[0], &t[0], &t[6]);

#ifdef VERIFY
{
    secp256k1_scalar tcheck;
    secp256k1_scalar_clear(&tcheck);
    secp256k1_bulletproof_circuit_lr_generator_init(&lr_gen, &rng, &y, &z, al, ar, ao, circ);
    for (i = 0; i < circ->n_gates; i++) {
        secp256k1_scalar lx, rx;
        secp256k1_bulletproof_circuit_lr_generate(&lr_gen, &lx, &rx, &x);
        secp256k1_bulletproof_circuit_lr_generator_update(&lr_gen);
        secp256k1_scalar_mul(&lx, &lx, &rx);
        secp256k1_scalar_add(&tcheck, &tcheck, &lx);
    }
    secp256k1_bulletproof_circuit_lr_generator_finalize(&lr_gen);
    CHECK(secp256k1_scalar_eq(&t[0], &tcheck));
}
#endif

    secp256k1_scalar_add(&mu, &alpha, &beta);
    secp256k1_scalar_add(&mu, &mu, &rho);

    /* Negate taux and mu so verifier doesn't have to */
    secp256k1_scalar_negate(&mu, &mu);
    secp256k1_scalar_negate(&taux, &taux);

    /* Mix these scalars into the hash so the input to the inner product proof is fixed */
    secp256k1_sha256_initialize(&sha256);
    secp256k1_sha256_write(&sha256, commit, 32);
    secp256k1_scalar_get_b32(commit, &t[0]);
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
    secp256k1_bulletproof_circuit_lr_generator_init(&abgh_data.lr_gen, &rng, &y, &z, al, ar, ao, circ);
    if (secp256k1_bulletproof_inner_product_prove_impl(ctx, ecmult_ctx, scratch, &a, &b, &out_pt[8], &out_pt[8 + depth], 1ull << depth, secp256k1_bulletproof_circuit_abgh_callback, (void *) &abgh_data, commit) == 0) {
        return 0;
    }
    secp256k1_bulletproof_circuit_lr_generator_finalize(&abgh_data.lr_gen);

    /* Encode everything */
    secp256k1_scalar_get_b32(&proof[0], &t[0]);
    secp256k1_scalar_get_b32(&proof[32], &taux);
    secp256k1_scalar_get_b32(&proof[64], &mu);
    secp256k1_scalar_get_b32(&proof[96], &a);
    secp256k1_scalar_get_b32(&proof[128], &b);
    secp256k1_bulletproof_serialize_points(&proof[160], out_pt, 8 + 2*depth);

    secp256k1_rfc6979_hmac_sha256_finalize(&rng);

    return 1;
}

typedef struct  {
    secp256k1_scalar mul_by;
    secp256k1_scalar x;
    secp256k1_scalar yinv;
    secp256k1_scalar z;
    const secp256k1_bulletproof_circuit *circ;
    /* state tracking */
    size_t count;
    /* eq 83 */
    secp256k1_ge age[3];
    /* eq 82 */
    secp256k1_scalar randomizer82;
    secp256k1_ge tge[5];
    secp256k1_scalar t;
    const secp256k1_ge *genp;
    const secp256k1_ge *commits;
    size_t n_commits;
} secp256k1_bulletproof_circuit_vfy_ecmult_context;

static int secp256k1_bulletproof_circuit_vfy_callback(secp256k1_scalar *sc, secp256k1_ge *pt, secp256k1_scalar *randomizer, size_t idx, void *data) {
    secp256k1_bulletproof_circuit_vfy_ecmult_context *ctx = (secp256k1_bulletproof_circuit_vfy_ecmult_context *) data;

    if (idx < ctx->circ->n_gates) { /* Gi */
        /* mul_by tells the caller how to get H'; however note that the scalar we
         * return in `sc` will be multiplied by H, not H'. Update it before doing
         * this round. */
        secp256k1_scalar_mul(&ctx->mul_by, &ctx->mul_by, &ctx->yinv);

        secp256k1_scalar_mul(sc, &ctx->circ->wr[idx].cached_sum, &ctx->x);
        secp256k1_scalar_mul(sc, sc, &ctx->mul_by); /* mul_by is y^-n */
        secp256k1_scalar_mul(sc, sc, randomizer);
    } else if (idx < 2 * ctx->circ->n_gates) { /* Hi */
        secp256k1_scalar dot;
        idx -= ctx->circ->n_gates;
        /* TODO do this smarter */
        secp256k1_scalar_mul(&ctx->mul_by, &ctx->mul_by, &ctx->yinv);
if(idx == 0) {
secp256k1_scalar_set_int(&ctx->mul_by, 1);
}
        secp256k1_scalar_mul(sc, &ctx->circ->wl[idx].cached_sum, &ctx->x);
        secp256k1_scalar_add(sc, sc, &ctx->circ->wo[idx].cached_sum);
        secp256k1_scalar_mul(sc, sc, &ctx->mul_by);

        secp256k1_scalar_set_int(&dot, 1);
        secp256k1_scalar_negate(&dot, &dot);
        secp256k1_scalar_add(sc, sc, &dot);

        secp256k1_scalar_mul(sc, sc, randomizer);
    /* return a (scalar, point) pair to add to the multiexp */
    } else {
        switch(ctx->count) {
        /* g^(x^2(k + <z^Q, c>) - t) (82) */
        case 0: {
            size_t i;
            secp256k1_scalar yn;
            secp256k1_scalar tmp;
            secp256k1_scalar_set_int(&yn, 1);
            secp256k1_scalar_clear(sc);
            for (i = 0; i < ctx->circ->n_gates; i++) {
                secp256k1_scalar term;
                secp256k1_scalar_mul(&term, &ctx->circ->wr[i].cached_sum, &ctx->circ->wl[i].cached_sum);
                secp256k1_scalar_mul(&term, &term, &yn);
                secp256k1_scalar_add(sc, sc, &term);

                secp256k1_scalar_mul(&yn, &yn, &ctx->yinv);
            }
            secp256k1_scalar_add(sc, sc, &ctx->circ->cached_c_sum);

            secp256k1_scalar_sqr(&tmp, &ctx->x);
            secp256k1_scalar_mul(sc, sc, &tmp);
            secp256k1_scalar_negate(&tmp, &ctx->t);
            secp256k1_scalar_add(sc, sc, &tmp);

            secp256k1_scalar_mul(sc, sc, &ctx->randomizer82);
            *pt = *ctx->genp;
            break;
        }
        /* A_I^x (83) */
        case 1:
            *sc = ctx->x;
            *pt = ctx->age[0];
            break;
        /* A_O^(x^2) (83) */
        case 2:
            secp256k1_scalar_sqr(sc, &ctx->x);
            *pt = ctx->age[1];
            break;
        /* S^(x^3) (83) */
        case 3:
            secp256k1_scalar_sqr(sc, &ctx->x); /* TODO cache previous squaring */
            secp256k1_scalar_mul(sc, sc, &ctx->x);
            *pt = ctx->age[2];
            break;
        /* T_1^x (82) */
        case 4:
            secp256k1_scalar_mul(sc, &ctx->x, &ctx->randomizer82);
            *pt = ctx->tge[0];
            break;
        default:
            if (ctx->count < 9) {
                size_t i;
                secp256k1_scalar_mul(sc, &ctx->x, &ctx->randomizer82);
                for (i = 0; i < ctx->count - 3; i++) {
                    secp256k1_scalar_mul(sc, sc, &ctx->x);
                }
                *pt = ctx->tge[ctx->count - 4];
            } else {
                /* V^(x^2 . (z^Q . W_V)) (82) */
                VERIFY_CHECK(!"bulletproof: too many points added by circuit_verify_impl to inner_product_verify_impl");
            }
        }
        secp256k1_scalar_mul(sc, sc, randomizer);
        ctx->count++;
    }
    return 1;
}

static int secp256k1_bulletproof_relation66_verify_impl(const secp256k1_context *ctx, const secp256k1_ecmult_context *ecmult_ctx, secp256k1_scratch *scratch, const unsigned char **proof, size_t *plen, size_t n_proofs, const secp256k1_ge *commitp, size_t nc, const secp256k1_ge *genp, secp256k1_bulletproof_circuit **circ, const secp256k1_ge *geng, const secp256k1_ge *genh, const unsigned char *extra_commit, size_t extra_commit_len) {
    const size_t depth = CTZ(circ[0]->n_gates);
    secp256k1_bulletproof_circuit_vfy_ecmult_context ecmult_data[MAX_BATCH_QTY];
    secp256k1_bulletproof_innerproduct_context innp_ctx[MAX_BATCH_QTY];
    size_t i;

    /* sanity-check input */
    for (i = 0; i < n_proofs; i++) {
        if (plen[i] != (13 + 2*depth) * 32 + (8 + 2*depth + 7) / 8) {
            return 0;
        }
        if (depth > SECP256K1_BULLETPROOF_MAX_DEPTH || plen[i] > SECP256K1_BULLETPROOF_MAX_PROOF) {
            return 0;
        }
    }

    for (i = 0; i < n_proofs; i++) {
        secp256k1_sha256 sha256;
        unsigned char randomizer82[32] = {0};  /* randomizer for eq (82) so we can add it to eq (83) to save a separate multiexp */
        unsigned char commit[32] = {0};
        secp256k1_scalar taux, mu, a, b;
        secp256k1_scalar y;
        int overflow;

        /* Commit to all input data: pedersen commit, asset generator, extra_commit */
        secp256k1_bulletproof_update_commit_n(commit, commitp, nc);
        secp256k1_bulletproof_update_commit_n(commit, genp, 1);
        if (extra_commit != NULL) {
            secp256k1_sha256_initialize(&sha256);
            secp256k1_sha256_write(&sha256, commit, 32);
            secp256k1_sha256_write(&sha256, extra_commit, extra_commit_len);
            secp256k1_sha256_finalize(&sha256, commit);
        }

        /* Deserialize everything */
        secp256k1_bulletproof_deserialize_point(&ecmult_data[i].age[0], &proof[i][160], 0, 8 + 2*depth);
        secp256k1_bulletproof_deserialize_point(&ecmult_data[i].age[1], &proof[i][160], 1, 8 + 2*depth);
        secp256k1_bulletproof_deserialize_point(&ecmult_data[i].age[2], &proof[i][160], 2, 8 + 2*depth);
        secp256k1_bulletproof_deserialize_point(&ecmult_data[i].tge[0], &proof[i][160], 3, 8 + 2*depth);
        secp256k1_bulletproof_deserialize_point(&ecmult_data[i].tge[1], &proof[i][160], 4, 8 + 2*depth);
        secp256k1_bulletproof_deserialize_point(&ecmult_data[i].tge[2], &proof[i][160], 5, 8 + 2*depth);
        secp256k1_bulletproof_deserialize_point(&ecmult_data[i].tge[3], &proof[i][160], 6, 8 + 2*depth);
        secp256k1_bulletproof_deserialize_point(&ecmult_data[i].tge[4], &proof[i][160], 7, 8 + 2*depth);

        /* Compute y, z, x */
        secp256k1_bulletproof_update_commit_n(commit, ecmult_data[i].age, 3);
        secp256k1_scalar_set_b32(&y, commit, &overflow);
        if (overflow || secp256k1_scalar_is_zero(&y)) {
            return 0;
        }
        secp256k1_scalar_inverse_var(&ecmult_data[i].yinv, &y);  /* TODO batch this into another inverse */
        secp256k1_bulletproof_update_commit_n(commit, NULL, 0);
        secp256k1_scalar_set_b32(&ecmult_data[i].z, commit, &overflow);
        if (overflow || secp256k1_scalar_is_zero(&ecmult_data[i].z)) {
            return 0;
        }

        secp256k1_bulletproof_update_commit_n(commit, ecmult_data[i].tge, 5);
        secp256k1_scalar_set_b32(&ecmult_data[i].x, commit, &overflow);
        if (overflow || secp256k1_scalar_is_zero(&ecmult_data[i].x)) {
            return 0;
        }

        secp256k1_circuit_compress(circ[i], &ecmult_data[i].z);

        /* Extract scalars */
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

        /* Mix these scalars into the hash so the input to the inner product proof is fixed */
        secp256k1_sha256_initialize(&sha256);
        secp256k1_sha256_write(&sha256, commit, 32);
        secp256k1_sha256_write(&sha256, proof[i], 96);
        secp256k1_sha256_finalize(&sha256, commit);

        secp256k1_sha256_initialize(&sha256);
        secp256k1_sha256_write(&sha256, commit, 32);
        secp256k1_sha256_finalize(&sha256, randomizer82);
        secp256k1_scalar_set_b32(&ecmult_data[i].randomizer82, randomizer82, &overflow);
        if (overflow || secp256k1_scalar_is_zero(&ecmult_data[i].randomizer82)) {
            return 0;
        }

        /* compute exponent offsets */
        ecmult_data[i].circ = circ[i];
        ecmult_data[i].count = 0;
        ecmult_data[i].mul_by = y;

        ecmult_data[i].genp = genp;
        ecmult_data[i].commits = commitp;
        ecmult_data[i].n_commits = nc;

        secp256k1_scalar_mul(&taux, &taux, &ecmult_data[i].randomizer82);
        secp256k1_scalar_add(&mu, &mu, &taux);

        innp_ctx[i].serialized_points = &proof[i][160];
        innp_ctx[i].n_extra_ser_points = 8;
        innp_ctx[i].a = a;
        innp_ctx[i].b = b;
        innp_ctx[i].dot = ecmult_data[i].t;
        innp_ctx[i].p_offs = mu;
        memcpy(innp_ctx[i].commit, commit, 32);
        innp_ctx[i].rangeproof_cb = secp256k1_bulletproof_circuit_vfy_callback;
        innp_ctx[i].rangeproof_cb_data = (void *) &ecmult_data[i];
        innp_ctx[i].n_extra_rangeproof_points = 9 + nc;
    }
    return secp256k1_bulletproof_inner_product_verify_impl(ctx, ecmult_ctx, scratch, geng, genh, circ[0]->n_gates, innp_ctx, n_proofs);
}

#endif
