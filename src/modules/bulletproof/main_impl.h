/**********************************************************************
 * Copyright (c) 2017 Andrew Poelstra                                 *
 * Distributed under the MIT software license, see the accompanying   *
 * file COPYING or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#ifndef SECP256K1_MODULE_BULLETPROOF_MAIN_IMPL
#define SECP256K1_MODULE_BULLETPROOF_MAIN_IMPL

#include "group.h"
#include "scalar.h"

#include "modules/rangeproof/main_impl.h"
#include "modules/rangeproof/pedersen_impl.h"

#include "modules/bulletproof/generators.h"
#include "modules/bulletproof/parser_impl.h"
#include "modules/bulletproof/inner_product_impl.h"
#include "modules/bulletproof/circuit_impl.h"
#include "modules/bulletproof/rangeproof_impl.h"
#include "modules/bulletproof/util.h"

#define MAX_WIDTH (1ul << 20)

int secp256k1_bulletproof_rangeproof_verify(const secp256k1_context* ctx, secp256k1_scratch_space *scratch, const unsigned char *proof, size_t plen,
 const secp256k1_pedersen_commitment* commit, size_t n_commits, size_t nbits, const secp256k1_generator* gen, const unsigned char *extra_commit, size_t extra_commit_len) {
    size_t i;
    secp256k1_ge genp;
    secp256k1_ge commitp[100];
    const secp256k1_ge *commitp_ptr = commitp;

    VERIFY_CHECK(ctx != NULL);
    ARG_CHECK(scratch != NULL);
    ARG_CHECK(commit != NULL);
    ARG_CHECK(gen != NULL);
    ARG_CHECK(n_commits <= sizeof(commitp) / sizeof(commitp[0]));
    ARG_CHECK(extra_commit != NULL || extra_commit_len == 0);
    ARG_CHECK(secp256k1_ecmult_context_is_built(&ctx->ecmult_ctx));

    secp256k1_generator_load(&genp, gen);
    for (i = 0; i < n_commits; i++) {
        secp256k1_pedersen_commitment_load(&commitp[i], &commit[i]);
    }

    return secp256k1_bulletproof_rangeproof_verify_impl(ctx, &ctx->ecmult_ctx, scratch, &proof, &plen, 1, nbits, NULL, NULL, NULL, NULL, &commitp_ptr, n_commits, &genp, &secp256k1_ge_const_gi[0], &secp256k1_ge_const_gi[64], extra_commit, extra_commit_len);
}

SECP256K1_API int secp256k1_bulletproof_rangeproof_unwind_message(
    const secp256k1_context* ctx,
    const unsigned char *proof, 
    size_t plen,
    const secp256k1_pedersen_commitment* commit,
    size_t nbits, 
    const secp256k1_generator* gen,
    const unsigned char *extra_commit, 
    size_t extra_commit_len, 
    const unsigned char* nonce,
    unsigned char* message) {

    secp256k1_ge genp;
    secp256k1_ge commitp[1];
    const secp256k1_ge *commitp_ptr = commitp;
    /* Recover tau_1 and mu that should have been XORed with message, using nonce+commit to seed the rng generator */
    unsigned char mu[32];
    unsigned char taux[32];
    unsigned char x[32];
    unsigned char z[32];
    unsigned char commit_bytes[32] = {0};
    int returnval;
    secp256k1_sha256 sha256;
    secp256k1_scratch_space *scratch = secp256k1_scratch_space_create(ctx, 1000000, 256 * MAX_WIDTH);

    VERIFY_CHECK(ctx != NULL);
    ARG_CHECK(commit != NULL);
    ARG_CHECK(nonce != NULL);
    ARG_CHECK(gen != NULL);
    ARG_CHECK(extra_commit != NULL || extra_commit_len == 0);
    ARG_CHECK(secp256k1_ecmult_context_is_built(&ctx->ecmult_ctx));

    secp256k1_generator_load(&genp, gen);
    secp256k1_pedersen_commitment_load(&commitp[0], &commit[0]);
    secp256k1_bulletproof_update_commit(commit_bytes, &commitp[0], &genp);

    returnval = secp256k1_bulletproof_rangeproof_verify_impl(ctx, &ctx->ecmult_ctx, scratch, &proof, &plen, 1, nbits, mu, taux, x, z, &commitp_ptr, 1, &genp, &secp256k1_ge_const_gi[0], &secp256k1_ge_const_gi[64], extra_commit, extra_commit_len);
    if (!returnval) return 0;

    secp256k1_scratch_space_destroy(scratch);

    if (extra_commit_len > 0) {
        secp256k1_sha256_initialize(&sha256);
        secp256k1_sha256_write(&sha256, commit_bytes, 32);
        secp256k1_sha256_write(&sha256, extra_commit, extra_commit_len);
        secp256k1_sha256_finalize(&sha256, commit_bytes);
    }
    /* now it's verified with tau1, mu, x and z recovered, which we need to get the XORed message back */
    secp256k1_bulletproof_unwind_message(commit_bytes, nonce, mu, taux, x, z, message);

    return returnval;
}

SECP256K1_API int secp256k1_bulletproof_rangeproof_verify_single_w_scratch(
    const secp256k1_context* ctx,
    const unsigned char* proof,
    size_t plen,
    const secp256k1_pedersen_commitment* commit,
    size_t nbits,
    const secp256k1_generator* gen,
    const unsigned char* extra_commit,
    size_t extra_commit_len
) {
    int returnval;
    secp256k1_scratch_space *scratch = secp256k1_scratch_space_create(ctx, 1000000, 256 * MAX_WIDTH);
    returnval=secp256k1_bulletproof_rangeproof_verify(ctx, scratch, proof, plen, commit, 1, nbits, gen, extra_commit, extra_commit_len);
    secp256k1_scratch_space_destroy(scratch);
    return returnval;
}


int secp256k1_bulletproof_rangeproof_verify_multi(const secp256k1_context* ctx, secp256k1_scratch_space *scratch, const unsigned char *proof, size_t plen, size_t n_proofs, const secp256k1_pedersen_commitment* commit, size_t n_commits, size_t nbits, const secp256k1_generator* gen, const unsigned char *extra_commit, size_t extra_commit_len) {
    const unsigned char *proof_ptr[MAX_BATCH_QTY];
    size_t plens[MAX_BATCH_QTY];
    secp256k1_ge genp;
    secp256k1_ge commitp[100];
    const secp256k1_ge *commitp_ptr[MAX_BATCH_QTY];
    size_t i;

    VERIFY_CHECK(ctx != NULL);
    ARG_CHECK(scratch != NULL);
    ARG_CHECK(commit != NULL);
    ARG_CHECK(gen != NULL);
    ARG_CHECK(extra_commit != NULL || extra_commit_len == 0);
    ARG_CHECK(secp256k1_ecmult_context_is_built(&ctx->ecmult_ctx));

    secp256k1_generator_load(&genp, gen);
for (i = 0; i < n_commits; i++) {
    secp256k1_pedersen_commitment_load(&commitp[i], &commit[i]);
}
for (i = 0; i < 100; i++) {
    proof_ptr[i] = proof;
    commitp_ptr[i] = commitp;
    plens[i] = plen;
}

    return secp256k1_bulletproof_rangeproof_verify_impl(ctx, &ctx->ecmult_ctx, scratch, proof_ptr, plens, n_proofs, nbits, NULL, NULL, NULL, NULL, commitp_ptr, n_commits, &genp, &secp256k1_ge_const_gi[0], &secp256k1_ge_const_gi[64], extra_commit, extra_commit_len);
}

int secp256k1_bulletproof_rangeproof_prove(const secp256k1_context* ctx, secp256k1_scratch_space *scratch, unsigned char *proof, size_t *plen, uint64_t *value, const unsigned char **blind, size_t n_commits,
 const secp256k1_generator* gen, size_t nbits, const unsigned char *nonce, const unsigned char *extra_commit, size_t extra_commit_len, const unsigned char* message) {
    secp256k1_ge commitp[100]; /* TODO choose a sane limit */
    secp256k1_scalar blinds[100];
    secp256k1_ge genp;
    size_t i;

    VERIFY_CHECK(ctx != NULL);
    ARG_CHECK(scratch != NULL);
    ARG_CHECK(proof != NULL);
    ARG_CHECK(plen != NULL);
    ARG_CHECK(blind != NULL);
    ARG_CHECK(gen != NULL);
    ARG_CHECK(nonce != NULL);
    ARG_CHECK(n_commits <= sizeof(commitp) / sizeof(commitp[0]));
    ARG_CHECK(nbits <= 64);
    if (nbits < 64) {
        for (i = 0; i < n_commits; i++) {
            ARG_CHECK(value[i] < (1ull << nbits));
            ARG_CHECK(blind[i] != NULL);
        }
    }
    ARG_CHECK(extra_commit != NULL || extra_commit_len == 0);
    ARG_CHECK(secp256k1_ecmult_context_is_built(&ctx->ecmult_ctx));
    ARG_CHECK(secp256k1_ecmult_gen_context_is_built(&ctx->ecmult_gen_ctx));

    secp256k1_generator_load(&genp, gen);

    for (i = 0; i < n_commits; i++) {
        int overflow;
        secp256k1_gej commitj;
        secp256k1_scalar_set_b32(&blinds[i], blind[i], &overflow);
        if (overflow || secp256k1_scalar_is_zero(&blinds[i])) {
            return 0;
        }
        secp256k1_pedersen_ecmult(&ctx->ecmult_gen_ctx, &commitj, &blinds[i], value[i], &genp);
        secp256k1_ge_set_gej(&commitp[i], &commitj);
    }

    return secp256k1_bulletproof_rangeproof_prove_impl(ctx, &ctx->ecmult_gen_ctx, &ctx->ecmult_ctx, scratch,
        proof, plen, nbits, value, blinds, commitp, n_commits, &genp, &secp256k1_ge_const_gi[0], &secp256k1_ge_const_gi[64], nonce, extra_commit, extra_commit_len, message);
}

SECP256K1_API int secp256k1_bulletproof_rangeproof_prove_single_w_scratch(
    const secp256k1_context* ctx,
    unsigned char* proof,
    size_t* plen,
    uint64_t value,
    const unsigned char* blind,
    const secp256k1_generator* gen,
    size_t nbits,
    const unsigned char* nonce,
    const unsigned char* extra_commit,
    size_t extra_commit_len,
    const unsigned char* message){

    int returnval;
    secp256k1_scratch_space *scratch = secp256k1_scratch_space_create(ctx, 1000000, 256 * MAX_WIDTH);
    const unsigned char* blinds[1];
    blinds[0] = blind;
    returnval=secp256k1_bulletproof_rangeproof_prove(ctx, scratch, proof, plen, &value, blinds, 1, gen, nbits, nonce, extra_commit, extra_commit_len, message);
    secp256k1_scratch_space_destroy(scratch);
    return returnval;

}


secp256k1_bulletproof_circuit *secp256k1_circuit_parse(const secp256k1_context *ctx, const char *description) {
    VERIFY_CHECK(ctx != NULL);
    ARG_CHECK(description != NULL);
    return secp256k1_parse_circuit(ctx, description);
}

void secp256k1_circuit_destroy(const secp256k1_context *ctx, secp256k1_bulletproof_circuit *circ) {
    VERIFY_CHECK(ctx != NULL);
    secp256k1_circuit_destroy_impl(circ);
}

int secp256k1_bulletproof_circuit_prove(const secp256k1_context* ctx, secp256k1_scratch_space *scratch, unsigned char *proof, size_t *plen, secp256k1_bulletproof_circuit *circ, unsigned char *nonce) {
#include "circuits/jubjub-3072.assn"

    VERIFY_CHECK(ctx != NULL);
    ARG_CHECK(scratch != NULL);
    ARG_CHECK(proof != NULL);
    ARG_CHECK(plen != NULL);
    ARG_CHECK(circ != NULL);
    ARG_CHECK(nonce != NULL);

    return secp256k1_bulletproof_relation66_prove_impl(
        ctx,
        &ctx->ecmult_ctx,
        scratch,
        proof, plen,
        incl_al, incl_ar, incl_ao, circ->n_gates,
        NULL, NULL, 0,
        &secp256k1_ge_const_g2,
        circ,
        &secp256k1_ge_const_gi[0], &secp256k1_ge_const_gi[32768],
        nonce,
        NULL, 0
    );
}

int secp256k1_bulletproof_circuit_verify(const secp256k1_context* ctx, secp256k1_scratch_space *scratch, const unsigned char *proof, size_t plen, secp256k1_bulletproof_circuit *circ) {
    VERIFY_CHECK(ctx != NULL);
    ARG_CHECK(scratch != NULL);
    ARG_CHECK(proof != NULL);
    ARG_CHECK(circ != NULL);

    return secp256k1_bulletproof_relation66_verify_impl(
        ctx,
        &ctx->ecmult_ctx,
        scratch,
        &proof, &plen, 1,
        NULL, 0,
        &secp256k1_ge_const_g2,
        &circ,
        &secp256k1_ge_const_gi[0], &secp256k1_ge_const_gi[32768],
        NULL, 0
    );
}

int secp256k1_bulletproof_circuit_verify_multi(const secp256k1_context* ctx, secp256k1_scratch_space *scratch, const unsigned char *proof, size_t plen, size_t n_proofs, secp256k1_bulletproof_circuit **circ) {
    const unsigned char *proof_ptr[MAX_BATCH_QTY];
    size_t plens[MAX_BATCH_QTY];
    size_t i;

    VERIFY_CHECK(ctx != NULL);
    ARG_CHECK(scratch != NULL);
    ARG_CHECK(proof != NULL);
    ARG_CHECK(circ != NULL);

    for (i = 0; i < n_proofs; i++) {
        proof_ptr[i] = proof;
        plens[i] = plen;
    }

    return secp256k1_bulletproof_relation66_verify_impl(
        ctx,
        &ctx->ecmult_ctx,
        scratch,
        proof_ptr, plens, n_proofs,
        NULL, 0,
        &secp256k1_ge_const_g2,
        circ,
        &secp256k1_ge_const_gi[0], &secp256k1_ge_const_gi[32768],
        NULL, 0
    );
}

#endif
