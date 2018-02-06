#ifndef _SECP256K1_BULLETPROOF_
# define _SECP256K1_BULLETPROOF_

# include "secp256k1.h"
# include "secp256k1_generator.h"
# include "secp256k1_rangeproof.h"

# ifdef __cplusplus
extern "C" {
# endif

#define MAX_BATCH_QTY	100

/* Maximum depth of 31 lets us validate an aggregate of 2^25 64-bit proofs */
#define SECP256K1_BULLETPROOF_MAX_DEPTH 60

/* Size of a hypothetical 31-depth rangeproof, in bytes */
#define SECP256K1_BULLETPROOF_MAX_PROOF (160 + 66*32 + 7)

typedef struct secp256k1_bulletproof_circuit secp256k1_bulletproof_circuit;

SECP256K1_API int secp256k1_bulletproof_rangeproof_verify(
    const secp256k1_context* ctx,
    secp256k1_scratch_space* scratch,
    const unsigned char* proof,
    size_t plen,
    const secp256k1_pedersen_commitment* commit,
    size_t n_commits,
    size_t nbits,
    const secp256k1_generator* gen,
    const unsigned char* extra_commit,
    size_t extra_commit_len
) SECP256K1_ARG_NONNULL(1) SECP256K1_ARG_NONNULL(2) SECP256K1_ARG_NONNULL(3) SECP256K1_ARG_NONNULL(5) SECP256K1_ARG_NONNULL(8);

SECP256K1_API int secp256k1_bulletproof_rangeproof_verify_single_w_scratch(
    const secp256k1_context* ctx,
    const unsigned char* proof,
    size_t plen,
    const secp256k1_pedersen_commitment* commit,
    size_t nbits,
    const secp256k1_generator* gen,
    const unsigned char* extra_commit,
    size_t extra_commit_len
) SECP256K1_ARG_NONNULL(1) SECP256K1_ARG_NONNULL(2) SECP256K1_ARG_NONNULL(4) SECP256K1_ARG_NONNULL(6);


SECP256K1_API int secp256k1_bulletproof_rangeproof_verify_multi(
    const secp256k1_context* ctx,
    secp256k1_scratch_space* scratch,
    const unsigned char* proof,
    size_t plen,
    size_t n_proofs,
    const secp256k1_pedersen_commitment* commit,
    size_t n_commits,
    size_t nbits,
    const secp256k1_generator* gen,
    const unsigned char* extra_commit,
    size_t extra_commit_len
) SECP256K1_ARG_NONNULL(1) SECP256K1_ARG_NONNULL(2) SECP256K1_ARG_NONNULL(3) SECP256K1_ARG_NONNULL(6) SECP256K1_ARG_NONNULL(9);


SECP256K1_API int secp256k1_bulletproof_rangeproof_prove(
    const secp256k1_context* ctx,
    secp256k1_scratch_space* scratch,
    unsigned char* proof,
    size_t* plen,
    uint64_t *value,
    const unsigned char** blind,
    size_t n_commits,
    const secp256k1_generator* gen,
    size_t nbits,
    const unsigned char* nonce,
    const unsigned char* extra_commit,
    size_t extra_commit_len
) SECP256K1_ARG_NONNULL(1) SECP256K1_ARG_NONNULL(2) SECP256K1_ARG_NONNULL(3) SECP256K1_ARG_NONNULL(4) SECP256K1_ARG_NONNULL(5) SECP256K1_ARG_NONNULL(6) SECP256K1_ARG_NONNULL(8) SECP256K1_ARG_NONNULL(10);

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
    size_t extra_commit_len
) SECP256K1_ARG_NONNULL(1) SECP256K1_ARG_NONNULL(2) SECP256K1_ARG_NONNULL(3) SECP256K1_ARG_NONNULL(5) SECP256K1_ARG_NONNULL(6) SECP256K1_ARG_NONNULL(8);

/* circuit stuff */
SECP256K1_API secp256k1_bulletproof_circuit *secp256k1_circuit_parse(
    const secp256k1_context *ctx,
    const char *description
) SECP256K1_ARG_NONNULL(1) SECP256K1_ARG_NONNULL(2);

SECP256K1_API void secp256k1_circuit_destroy(
    const secp256k1_context *ctx,
    secp256k1_bulletproof_circuit *circ
) SECP256K1_ARG_NONNULL(1);

SECP256K1_API int secp256k1_bulletproof_circuit_prove(
    const secp256k1_context* ctx,
    secp256k1_scratch_space *scratch,
    unsigned char *proof,
    size_t *plen,
    secp256k1_bulletproof_circuit *circ,
    unsigned char *nonce
) SECP256K1_ARG_NONNULL(1) SECP256K1_ARG_NONNULL(2) SECP256K1_ARG_NONNULL(3) SECP256K1_ARG_NONNULL(4) SECP256K1_ARG_NONNULL(5) SECP256K1_ARG_NONNULL(6);

SECP256K1_API int secp256k1_bulletproof_circuit_verify(
    const secp256k1_context* ctx,
    secp256k1_scratch_space *scratch,
    const unsigned char *proof,
    size_t plen,
    secp256k1_bulletproof_circuit *circ
) SECP256K1_ARG_NONNULL(1) SECP256K1_ARG_NONNULL(2) SECP256K1_ARG_NONNULL(3) SECP256K1_ARG_NONNULL(5);

SECP256K1_API int secp256k1_bulletproof_circuit_verify_multi(
    const secp256k1_context* ctx,
    secp256k1_scratch_space *scratch,
    const unsigned char *proof,
    size_t plen,
    size_t n_proofs,
    secp256k1_bulletproof_circuit **circ
) SECP256K1_ARG_NONNULL(1) SECP256K1_ARG_NONNULL(2) SECP256K1_ARG_NONNULL(3) SECP256K1_ARG_NONNULL(6);


# ifdef __cplusplus
}
# endif

#endif
