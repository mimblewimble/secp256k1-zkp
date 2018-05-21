/**********************************************************************
 * Copyright (c) 2017 Andrew Poelstra                                 *
 * Distributed under the MIT software license, see the accompanying   *
 * file COPYING or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#ifndef SECP256K1_MODULE_BULLETPROOF_TESTS
#define SECP256K1_MODULE_BULLETPROOF_TESTS

#include <string.h>

#include "group.h"
#include "scalar.h"
#include "testrand.h"
#include "util.h"

#include "modules/bulletproof/generators.h"

#include "include/secp256k1_bulletproof.h"

#define MAX_WIDTH (1ul << 20)
typedef struct {
    const secp256k1_scalar *a;
    const secp256k1_scalar *b;
    const secp256k1_ge *g;
    const secp256k1_ge *h;
    size_t n;
} test_bulletproof_ecmult_context;

static int test_bulletproof_ecmult_callback(secp256k1_scalar *sc, secp256k1_ge *pt, size_t idx, void *data) {
    test_bulletproof_ecmult_context *ecctx = (test_bulletproof_ecmult_context *) data;
    if (idx < ecctx->n) {
        *sc = ecctx->a[idx];
        *pt = ecctx->g[idx];
    } else {
        VERIFY_CHECK(idx < 2*ecctx->n);
        *sc = ecctx->b[idx - ecctx->n];
        *pt = ecctx->h[idx - ecctx->n];
    }
    return 1;
}

typedef struct {
    secp256k1_scalar offs;
    secp256k1_scalar ext_sc;
    secp256k1_scalar skew_sc;
    secp256k1_ge ext_pt;
    secp256k1_ge p;
    size_t n;
    int parity;
} test_bulletproof_offset_context;

static int test_bulletproof_offset_vfy_callback(secp256k1_scalar *sc, secp256k1_ge *pt, secp256k1_scalar *randomizer, size_t idx, void *data) {
    test_bulletproof_offset_context *ecctx = (test_bulletproof_offset_context *) data;
    secp256k1_scalar_set_int(&ecctx->offs, 1);
    if (idx < 2 * ecctx->n) {
        secp256k1_scalar idxsc;
        secp256k1_scalar_set_int(&idxsc, idx);
        secp256k1_scalar_mul(sc, &ecctx->skew_sc, &idxsc);
    } else {
        if (ecctx->parity) {
            *sc = ecctx->ext_sc;
            *pt = ecctx->ext_pt;
        } else {
            secp256k1_scalar_set_int(sc, 1);
            *pt = ecctx->p;
        }
    }
    secp256k1_scalar_mul(sc, sc, randomizer);
    ecctx->parity = !ecctx->parity;
    return 1;
}

typedef struct {
    const secp256k1_ge *geng;
    const secp256k1_ge *genh;
    const secp256k1_scalar *a_arr;
    const secp256k1_scalar *b_arr;
} secp256k1_bulletproof_ip_test_abgh_data;


static int secp256k1_bulletproof_ip_test_abgh_callback(secp256k1_scalar *sc, secp256k1_ge *pt, size_t idx, void *data) {
    secp256k1_bulletproof_ip_test_abgh_data *cbctx = (secp256k1_bulletproof_ip_test_abgh_data *) data;
    const int is_g = idx % 2 == 0;

    random_scalar_order(sc);
    if (is_g) {
        *sc = cbctx->a_arr[idx / 2];
        *pt = cbctx->geng[idx / 2];
    } else {
        *sc = cbctx->b_arr[idx / 2];
        *pt = cbctx->genh[idx / 2];
    }
    return 1;
}

void test_bulletproof_inner_product(size_t depth, const secp256k1_ge *geng, const secp256k1_ge *genh) {
    const secp256k1_scalar zero = SECP256K1_SCALAR_CONST(0,0,0,0,0,0,0,0);
    secp256k1_gej pj;
    secp256k1_gej tmpj, tmpj2;
    secp256k1_ge out_pt[128];
    secp256k1_scalar *a_arr = (secp256k1_scalar *)checked_malloc(&ctx->error_callback, MAX_WIDTH * sizeof(*a_arr));
    secp256k1_scalar *b_arr = (secp256k1_scalar *)checked_malloc(&ctx->error_callback, MAX_WIDTH * sizeof(*b_arr));
    unsigned char commit[32] = "hash of P, c, etc. all that jazz";
    unsigned char serialized_points[32 * 128];
    size_t j;
    test_bulletproof_offset_context offs_ctx;
    secp256k1_bulletproof_ip_test_abgh_data abgh_data;
    secp256k1_bulletproof_innerproduct_context innp_ctx;

    secp256k1_scratch *scratch = secp256k1_scratch_space_create(ctx, 1000000, 256 * MAX_WIDTH);

    CHECK(depth < SECP256K1_BULLETPROOF_MAX_DEPTH);
    CHECK(1u << depth <= MAX_WIDTH);

    for (j = 0; j < 1u << depth; j++) {
        random_scalar_order(&a_arr[j]);
        random_scalar_order(&b_arr[j]);
    }
    secp256k1_scalar_dot_product(&innp_ctx.dot, a_arr, b_arr, 1ul << depth);

    abgh_data.geng = geng;
    abgh_data.genh = genh;
    abgh_data.a_arr = a_arr;
    abgh_data.b_arr = b_arr;

    random_scalar_order(&innp_ctx.p_offs);
    random_group_element_test(&offs_ctx.ext_pt);
    random_scalar_order(&offs_ctx.ext_sc);
    secp256k1_scalar_clear(&offs_ctx.skew_sc);
    offs_ctx.n = 1 << depth;

    CHECK(secp256k1_bulletproof_inner_product_prove_impl(ctx, &ctx->ecmult_ctx, scratch, &innp_ctx.a, &innp_ctx.b, &out_pt[0], &out_pt[depth], 1ul << depth, secp256k1_bulletproof_ip_test_abgh_callback, (void *) &abgh_data, commit) == 1);

    innp_ctx.serialized_points = serialized_points;
    secp256k1_bulletproof_serialize_points(serialized_points, out_pt, 2 * depth);
    innp_ctx.n_extra_ser_points = 0;
    memcpy(innp_ctx.commit, commit, 32);
    innp_ctx.rangeproof_cb = NULL;
    innp_ctx.rangeproof_cb_data = NULL;
    innp_ctx.n_extra_rangeproof_points = 1;
    innp_ctx.rangeproof_cb = test_bulletproof_offset_vfy_callback;
    innp_ctx.rangeproof_cb_data = (void *) &offs_ctx;

    {
        test_bulletproof_ecmult_context ecmult_data;
        ecmult_data.n = 1 << depth;
        ecmult_data.a = a_arr;
        ecmult_data.b = b_arr;
        ecmult_data.g = geng;
        ecmult_data.h = genh;
        CHECK(secp256k1_ecmult_multi_var(&ctx->ecmult_ctx, scratch, &ctx->error_callback, &pj, &zero, test_bulletproof_ecmult_callback, (void*) &ecmult_data, 2 << depth));
    }

    /* skew P by a random amount and instruct the verifier to offset it */
    secp256k1_ecmult_gen(&ctx->ecmult_gen_ctx, &tmpj, &innp_ctx.p_offs);
    secp256k1_gej_add_var(&pj, &pj, &tmpj, NULL);
    secp256k1_ge_set_gej(&offs_ctx.p, &pj);

    /* wrong p_offs should fail */
    offs_ctx.parity = 0;
    CHECK(secp256k1_bulletproof_inner_product_verify_impl(ctx, &ctx->ecmult_ctx, scratch, geng, genh, 1 << depth, &innp_ctx, 1) == 0);

    secp256k1_scalar_negate(&innp_ctx.p_offs, &innp_ctx.p_offs);

    offs_ctx.parity = 0;
    CHECK(secp256k1_bulletproof_inner_product_verify_impl(ctx, &ctx->ecmult_ctx, scratch, geng, genh, 1 << depth, &innp_ctx, 1) == 1);
    /* check that verification did not trash anything */
    offs_ctx.parity = 0;
    CHECK(secp256k1_bulletproof_inner_product_verify_impl(ctx, &ctx->ecmult_ctx, scratch, geng, genh, 1 << depth, &innp_ctx, 1) == 1);
    /* check that adding a no-op rangeproof skew function doesn't break anything */
    offs_ctx.parity = 0;
    CHECK(secp256k1_bulletproof_inner_product_verify_impl(ctx, &ctx->ecmult_ctx, scratch, geng, genh, 1 << depth, &innp_ctx, 1) == 1);

    /* Offset P by some random point and then try to undo this in the verification */
    secp256k1_gej_set_ge(&tmpj2, &offs_ctx.ext_pt);
    secp256k1_ecmult(&ctx->ecmult_ctx, &tmpj, &tmpj2, &offs_ctx.ext_sc, &zero);
    secp256k1_gej_neg(&tmpj, &tmpj);
    secp256k1_gej_add_ge_var(&tmpj, &tmpj, &offs_ctx.p, NULL);
    secp256k1_ge_set_gej(&offs_ctx.p, &tmpj);
    offs_ctx.parity = 0;
    innp_ctx.n_extra_rangeproof_points = 2;
    CHECK(secp256k1_bulletproof_inner_product_verify_impl(ctx, &ctx->ecmult_ctx, scratch, geng, genh, 1 << depth, &innp_ctx, 1) == 1);

    /* Offset each basis by some random point and try to undo this in the verification */
    secp256k1_gej_set_infinity(&tmpj2);
    for (j = 0; j < 1u << depth; j++) {
        size_t k;
        /* Offset by k-times the kth G basis and (k+n)-times the kth H basis */
        for (k = 0; k < j; k++) {
            secp256k1_gej_add_ge_var(&tmpj2, &tmpj2, &geng[j], NULL);
            secp256k1_gej_add_ge_var(&tmpj2, &tmpj2, &genh[j], NULL);
        }
        for (k = 0; k < 1u << depth; k++) {
            secp256k1_gej_add_ge_var(&tmpj2, &tmpj2, &genh[j], NULL);
        }
    }
    random_scalar_order(&offs_ctx.skew_sc);
    secp256k1_ecmult(&ctx->ecmult_ctx, &tmpj, &tmpj2, &offs_ctx.skew_sc, &zero);
    secp256k1_gej_add_ge_var(&tmpj, &tmpj, &offs_ctx.p, NULL);
    secp256k1_ge_set_gej(&offs_ctx.p, &tmpj);
    secp256k1_scalar_negate(&offs_ctx.skew_sc, &offs_ctx.skew_sc);

    offs_ctx.parity = 0;
    CHECK(secp256k1_bulletproof_inner_product_verify_impl(ctx, &ctx->ecmult_ctx, scratch, geng, genh, 1 << depth, &innp_ctx, 1) == 1);

    /* Try to validate the same proof twice */
{
    test_bulletproof_offset_context offs_ctxs[2];
    secp256k1_bulletproof_innerproduct_context innp_ctxs[2];
    offs_ctx.parity = 1;  /* set parity to 1 so the common point will be returned first, as required by the multi-proof verifier */
    memcpy(&innp_ctxs[0], &innp_ctx, sizeof(innp_ctx));
    memcpy(&innp_ctxs[1], &innp_ctx, sizeof(innp_ctx));
    memcpy(&offs_ctxs[0], &offs_ctx, sizeof(offs_ctx));
    memcpy(&offs_ctxs[1], &offs_ctx, sizeof(offs_ctx));
    innp_ctxs[0].rangeproof_cb_data = &offs_ctxs[0];
    innp_ctxs[1].rangeproof_cb_data = &offs_ctxs[1];
    CHECK(secp256k1_bulletproof_inner_product_verify_impl(ctx, &ctx->ecmult_ctx, scratch, geng, genh, 1 << depth, innp_ctxs, 2) == 1);
}

    free(a_arr);
    free(b_arr);
    secp256k1_scratch_destroy(scratch);
}

void test_bulletproof_rangeproof(size_t nbits, size_t expected_size, const secp256k1_ge *geng, const secp256k1_ge *genh) {
    secp256k1_scalar blind;
    unsigned char proof[1024];
    const unsigned char *proof_ptr[2];
    size_t plen[2] = { sizeof(proof), sizeof(proof) };
    uint64_t v = 123456;
    secp256k1_gej commitj;
    secp256k1_ge commitp;
    const secp256k1_ge *commitp_ptr[2];
    secp256k1_ge genp[2];
    secp256k1_scalar vs;
    secp256k1_gej altgen;
    unsigned char nonce[32] = "my kingdom for some randomness!!";

    secp256k1_scratch *scratch = secp256k1_scratch_space_create(ctx, 1000000, 10000000);

    if (v >> nbits > 0) {
        v = 0;
    }

    proof_ptr[0] = proof_ptr[1] = proof;

    secp256k1_gej_set_ge(&altgen, &secp256k1_ge_const_g2);
    random_scalar_order(&blind);
    secp256k1_scalar_set_u64(&vs, v);
    secp256k1_ecmult(&ctx->ecmult_ctx, &commitj, &altgen, &vs, &blind);
    secp256k1_ge_set_gej(&commitp, &commitj);
    commitp_ptr[0] = commitp_ptr[1] = &commitp;

    genp[0] = genp[1] = secp256k1_ge_const_g2;

    CHECK(secp256k1_bulletproof_rangeproof_prove_impl(ctx, &ctx->ecmult_gen_ctx, &ctx->ecmult_ctx, scratch, proof, &plen[0], nbits, &v, &blind, &commitp, 1, &secp256k1_ge_const_g2, geng, genh, nonce, NULL, 0, NULL) == 1);
    CHECK(plen[0] == expected_size);
    plen[1] = plen[0];
    CHECK(secp256k1_bulletproof_rangeproof_verify_impl(ctx, &ctx->ecmult_ctx, scratch, proof_ptr, plen, 2, nbits, NULL, NULL, NULL, NULL, commitp_ptr, 1, genp, geng, genh, NULL, 0) == 1);

    secp256k1_scratch_destroy(scratch);
}

void test_bulletproof_rangeproof_message(size_t nbits, size_t expected_size, secp256k1_ge *geng) {
    secp256k1_scalar blind;
    unsigned char proof[1024];
    const unsigned char *blind_ptr[1];
    size_t plen[1] = { sizeof(proof) };
    uint64_t v = 123456;
    secp256k1_ge commitp;
    secp256k1_gej commitj;
    secp256k1_generator tmp_gen;
    unsigned char message[64] = "inserting a 64 byte message into a rangeproof to recover later\n";
    unsigned char recovered_message[64];
    unsigned char blind_bytes[32];
    secp256k1_scalar vs;
    secp256k1_gej altgen;
    secp256k1_pedersen_commitment p_commit;
 
    secp256k1_scratch *scratch = secp256k1_scratch_space_create(ctx, 1000000, 10000000);

    if (v >> nbits > 0) {
        v = 0;
    }

    secp256k1_gej_set_ge(&altgen, &secp256k1_ge_const_g2);
    random_scalar_order(&blind);
    secp256k1_scalar_set_u64(&vs, v);
    secp256k1_ecmult(&ctx->ecmult_ctx, &commitj, &altgen, &vs, &blind);
    secp256k1_ge_set_gej(&commitp, &commitj);
    secp256k1_scalar_get_b32(blind_bytes, &blind);
    secp256k1_pedersen_commitment_save(&p_commit, &commitp);

    blind_ptr[0] = blind_bytes;

    secp256k1_generator_save(&tmp_gen, geng);
    CHECK(secp256k1_bulletproof_rangeproof_prove(ctx, scratch, proof, &plen[0], &v, blind_ptr, 1, &tmp_gen, nbits, blind_bytes, NULL, 0, message) == 1);
    secp256k1_scratch_destroy(scratch);
    CHECK(plen[0] == expected_size);
    secp256k1_pedersen_commitment_load(&commitp, &p_commit);
    CHECK(secp256k1_bulletproof_rangeproof_unwind_message(ctx, proof, plen[0], &p_commit, nbits, &tmp_gen, NULL, 0, blind_bytes, blind_bytes, recovered_message) == 1);
    CHECK(strcmp((const char*) message, (const char *) recovered_message) == 0);

}

void test_bulletproof_rangeproof_aggregate(size_t nbits, size_t n_commits, size_t expected_size, const secp256k1_ge *geng, const secp256k1_ge *genh) {
    unsigned char proof[1024];
    const unsigned char *proof_ptr = proof;
    size_t plen = sizeof(proof);
    secp256k1_scalar *blind = (secp256k1_scalar *)checked_malloc(&ctx->error_callback, n_commits * sizeof(*blind));
    uint64_t *v = (uint64_t *)checked_malloc(&ctx->error_callback, n_commits * sizeof(*v));
    secp256k1_ge *commitp = (secp256k1_ge *)checked_malloc(&ctx->error_callback, n_commits * sizeof(*commitp));
    const secp256k1_ge *constptr = commitp;
    secp256k1_ge genp;
    unsigned char commit[32] = {0};
    unsigned char nonce[32] = "mary, mary quite contrary how do";
    size_t i;

    secp256k1_scratch *scratch = secp256k1_scratch_space_create(ctx, 1000000, 10000000);

    genp = secp256k1_ge_const_g2;
    for (i = 0; i < n_commits; i++) {
        secp256k1_scalar vs;
        secp256k1_gej commitj;
        secp256k1_gej genpj;

        v[i] = 223 * i; /* dice-roll random # */
        if (v[i] >> nbits > 0) {
            v[i] = 0;
        }
        secp256k1_scalar_set_u64(&vs, v[i]);
        random_scalar_order(&blind[i]);
        secp256k1_gej_set_ge(&genpj, &genp);
        secp256k1_ecmult(&ctx->ecmult_ctx, &commitj, &genpj, &vs, &blind[i]);
        secp256k1_ge_set_gej(&commitp[i], &commitj);

        secp256k1_bulletproof_update_commit(commit, &commitp[i], &genp);
    }

    CHECK(secp256k1_bulletproof_rangeproof_prove_impl(ctx, &ctx->ecmult_gen_ctx, &ctx->ecmult_ctx, scratch, proof, &plen, nbits, v, blind, commitp, n_commits, &genp, geng, genh, nonce, NULL, 0, NULL) == 1);
    CHECK(plen == expected_size);
    CHECK(secp256k1_bulletproof_rangeproof_verify_impl(ctx, &ctx->ecmult_ctx, scratch, &proof_ptr, &plen, 1, nbits, NULL, NULL, NULL, NULL, &constptr, n_commits, &genp, geng, genh, NULL, 0) == 1);

    secp256k1_scratch_destroy(scratch);
    free(commitp);
    free(v);
    free(blind);
}

void test_bulletproof_circuit(const secp256k1_ge *geng, const secp256k1_ge *genh) {
    unsigned char proof[2000];
    const unsigned char nonce[32] = "ive got a bit won't tell u which";
    const unsigned char *proof_ptr = proof;
    size_t plen = sizeof(proof);
    secp256k1_scalar one;
    secp256k1_scalar al[2];
    secp256k1_scalar ar[2];
    secp256k1_scalar ao[2];
    secp256k1_scratch *scratch = secp256k1_scratch_space_create(ctx, 1000000, 100000000);
#include "circuits/jubjub-12.circuit"
#include "circuits/jubjub-12.assn"

    const char inv_17_19_circ[] = "2,0,4; L0 = 17; 2*L1 - L0 = 21; O0 = 1; O1 = 1;";
    secp256k1_bulletproof_circuit *simple = secp256k1_parse_circuit(ctx, inv_17_19_circ);
    secp256k1_bulletproof_circuit *incl = secp256k1_parse_circuit(ctx, incl_desc);

secp256k1_scalar challenge;
secp256k1_scalar answer;

    CHECK (simple != NULL);
secp256k1_scalar_set_int(&challenge, 17);
secp256k1_scalar_inverse(&answer, &challenge);

    secp256k1_scalar_set_int(&one, 1);

    /* Try to prove with input 0, 1, 0 */
    al[0] = al[1] = challenge;
    ar[0] = ar[1] = answer;
    ao[0] = ao[1] = one;

secp256k1_scalar_set_int(&challenge, 19);
secp256k1_scalar_inverse(&answer, &challenge);
    al[1] = challenge;
    ar[1] = answer;

    CHECK(secp256k1_bulletproof_relation66_prove_impl(
        ctx,
        &ctx->ecmult_ctx,
        scratch,
        proof, &plen,
        al, ar, ao, 2,
        NULL, NULL, 0,
        &secp256k1_ge_const_g2,
        simple,
        geng, genh,
        nonce,
        NULL, 0
    ));

    CHECK(secp256k1_bulletproof_relation66_verify_impl(
        ctx,
        &ctx->ecmult_ctx,
        scratch,
        &proof_ptr, &plen, 1,
        NULL, 0,
        &secp256k1_ge_const_g2,
        &simple,
        geng, genh,
        NULL, 0
    ));

    plen = 2000;
    CHECK(secp256k1_bulletproof_relation66_prove_impl(
        ctx,
        &ctx->ecmult_ctx,
        scratch,
        proof, &plen,
        incl_al, incl_ar, incl_ao, incl->n_gates,
        NULL, NULL, 0,
        &secp256k1_ge_const_g2,
        incl,
        geng, genh,
        nonce,
        NULL, 0
    ));

    CHECK(secp256k1_bulletproof_relation66_verify_impl(
        ctx,
        &ctx->ecmult_ctx,
        scratch,
        &proof_ptr, &plen, 1,
        NULL, 0,
        &secp256k1_ge_const_g2,
        &incl,
        geng, genh,
        NULL, 0
    ));

    secp256k1_circuit_destroy(ctx, simple);
    secp256k1_circuit_destroy(ctx, incl);
    secp256k1_scratch_destroy(scratch);
}

void run_bulletproof_tests(void) {
    size_t i;

    /* Make a ton of generators */
    size_t n_gens = 16384;
    secp256k1_ge *geng = (secp256k1_ge *)checked_malloc(&ctx->error_callback, sizeof(secp256k1_ge) * n_gens);
    secp256k1_ge *genh = (secp256k1_ge *)checked_malloc(&ctx->error_callback, sizeof(secp256k1_ge) * n_gens);
    for (i = 0; i < n_gens; i++) {
       secp256k1_generator tmpgen;
       unsigned char commit[32] = { 0 };
       commit[0] = i;
       commit[1] = i >> 8;
       commit[2] = i >> 16;
       commit[3] = i >> 24;

       commit[31] = 'G';
       CHECK(secp256k1_generator_generate(ctx, &tmpgen, commit));
       secp256k1_generator_load(&geng[i], &tmpgen);
       commit[31] = 'H';
       CHECK(secp256k1_generator_generate(ctx, &tmpgen, commit));
       secp256k1_generator_load(&genh[i], &tmpgen);
    }

    test_bulletproof_rangeproof(64, 674, geng, genh);
    test_bulletproof_rangeproof_message(64, 674, geng);
#if 0
    /* sanity checks */
    test_bulletproof_inner_product(0, geng, genh);
    test_bulletproof_inner_product(1, geng, genh);
    test_bulletproof_inner_product(2, geng, genh);
    for (i = 0; i < (size_t) count; i++) {
        test_bulletproof_inner_product(5, geng, genh);
        test_bulletproof_inner_product(6, geng, genh);
    }
    test_bulletproof_inner_product(10, geng, genh);

    test_bulletproof_rangeproof(1, 289, geng, genh);
    test_bulletproof_rangeproof(2, 353, geng, genh);
    test_bulletproof_rangeproof(16, 546, geng, genh);
    test_bulletproof_rangeproof(32, 610, geng, genh);
    test_bulletproof_rangeproof(64, 674, geng, genh);

    test_bulletproof_rangeproof_aggregate(64, 1, 674, geng, genh);
    test_bulletproof_rangeproof_aggregate(8, 2, 546, geng, genh);
    test_bulletproof_rangeproof_aggregate(8, 4, 610, geng, genh);

#endif
    test_bulletproof_circuit(geng, genh);

    free(geng);
    free(genh);
}
#undef MAX_WIDTH

#endif
