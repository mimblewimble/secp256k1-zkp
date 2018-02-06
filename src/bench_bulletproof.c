/**********************************************************************
 * Copyright (c) 2017 Andrew Poelstra                                 *
 * Distributed under the MIT software license, see the accompanying   *
 * file COPYING or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#include <stdint.h>

#include "include/secp256k1_generator.h"
#include "include/secp256k1_bulletproof.h"
#include "include/secp256k1_rangeproof.h"
#include "util.h"
#include "bench.h"

#define CIRC "jubjub-3072"
#define FILENAME "src/modules/bulletproof/circuits/jubjub-3072.circuit"
#define BENCHPROVE "bulletproof_" CIRC "_prove,"
#define BENCHVFY "bulletproof_" CIRC "_verify,"
#define BENCHVFY2 "bulletproof_" CIRC "_verify+1,"

typedef struct {
    secp256k1_context *ctx;
    secp256k1_scratch_space *scratch;
    secp256k1_generator altgen;
    secp256k1_pedersen_commitment *commit;
    secp256k1_bulletproof_circuit *circ;
    secp256k1_bulletproof_circuit *circ2;
    size_t *value;
    unsigned char nonce[32];
    const unsigned char **blind;
    size_t nbits;
    size_t n_commits;
    size_t n_proofs;
    unsigned char proof[2000];
    size_t plen;
} bench_bulletproof_t;

static void bench_bulletproof_setup(void* arg) {
    bench_bulletproof_t *data = (bench_bulletproof_t*)arg;
    size_t i;

    const unsigned char nonce[32] = "my kingdom for some randomness!!";
    const unsigned char genbd[32] = "yet more blinding, for the asset";
    unsigned char blind[32] = "and my kingdom too for a blinder";

    memcpy(data->nonce, nonce, 32);
    data->commit = (secp256k1_pedersen_commitment *)malloc(data->n_commits * sizeof(*data->commit));
    data->blind = (const unsigned char **)malloc(data->n_commits * sizeof(*data->commit));
    data->value = (size_t *)malloc(data->n_commits * sizeof(*data->commit));

    CHECK(secp256k1_generator_generate(data->ctx, &data->altgen, genbd));
    for (i = 0; i < data->n_commits; i++) {
        data->blind[i] = malloc(32);
        blind[0] = i;
        blind[1] = i >> 8;
        memcpy((unsigned char*) data->blind[i], blind, 32);
        data->value[i] = i * 17;
        CHECK(secp256k1_pedersen_commit(data->ctx, &data->commit[i], data->blind[i], data->value[i], &data->altgen));
    }

    data->plen = sizeof(data->proof);
    CHECK(secp256k1_bulletproof_rangeproof_prove(data->ctx, data->scratch, data->proof, &data->plen, data->value, data->blind, data->n_commits, &data->altgen, data->nbits, data->nonce, NULL, 0) == 1);
    CHECK(secp256k1_bulletproof_rangeproof_verify(data->ctx, data->scratch, data->proof, data->plen, data->commit, data->n_commits, data->nbits, &data->altgen, NULL, 0) == 1);
}

static void bench_bulletproof_teardown(void* arg) {
    bench_bulletproof_t *data = (bench_bulletproof_t*)arg;
    size_t i;

    if (data->blind != NULL) {
        for (i = 0; i < data->n_commits; i++) {
            free((unsigned char*) data->blind[i]);
        }
    }
    free(data->blind);
    free(data->value);
    free(data->commit);
}

static void bench_bulletproof_circuit_setup(void* arg) {
    bench_bulletproof_t *data = (bench_bulletproof_t*)arg;
    secp256k1_bulletproof_circuit *circ_ptr[2];
    data->blind = NULL;
    data->value = NULL;
    data->commit = NULL;
    data->plen = sizeof(data->proof);
    circ_ptr[0] = data->circ;
    circ_ptr[1] = data->circ2;
    CHECK(secp256k1_bulletproof_circuit_prove(data->ctx, data->scratch, data->proof, &data->plen, data->circ, data->nonce) == 1);
    CHECK(secp256k1_bulletproof_circuit_verify(data->ctx, data->scratch, data->proof, data->plen, data->circ) == 1);
    CHECK(secp256k1_bulletproof_circuit_verify_multi(data->ctx, data->scratch, data->proof, data->plen, 2, circ_ptr) == 1);
}

static void bench_bulletproof_jubjub_prove(void* arg) {
    bench_bulletproof_t *data = (bench_bulletproof_t*)arg;
    size_t i;
    for (i = 0; i < 5; i++) {
        CHECK(secp256k1_bulletproof_circuit_prove(data->ctx, data->scratch, data->proof, &data->plen, data->circ, data->nonce) == 1);
    }
}

static void bench_bulletproof_jubjub_verify(void* arg) {
    bench_bulletproof_t *data = (bench_bulletproof_t*)arg;
    size_t i;
    for (i = 0; i < 10; i++) {
        CHECK(secp256k1_bulletproof_circuit_verify(data->ctx, data->scratch, data->proof, data->plen, data->circ) == 1);
    }
}

static void bench_bulletproof_jubjub_verify_plus1(void* arg) {
    bench_bulletproof_t *data = (bench_bulletproof_t*)arg;
    secp256k1_bulletproof_circuit *circ_ptr[2];
    size_t i;
    circ_ptr[0] = data->circ;
    circ_ptr[1] = data->circ2;
    for (i = 0; i < 10; i++) {
        CHECK(secp256k1_bulletproof_circuit_verify_multi(data->ctx, data->scratch, data->proof, data->plen, 2, circ_ptr) == 1);
    }
}


static void bench_bulletproof_prove(void* arg) {
    bench_bulletproof_t *data = (bench_bulletproof_t*)arg;
    size_t i;
    for (i = 0; i < 25; i++) {
        CHECK(secp256k1_bulletproof_rangeproof_prove(data->ctx, data->scratch, data->proof, &data->plen, data->value, data->blind, data->n_commits, &data->altgen, data->nbits, data->nonce, NULL, 0) == 1);
    }
}

static void bench_bulletproof_verify(void* arg) {
    int i;
    bench_bulletproof_t *data = (bench_bulletproof_t*)arg;

    for (i = 0; i < 100; i++) {
        CHECK(secp256k1_bulletproof_rangeproof_verify_multi(data->ctx, data->scratch, data->proof, data->plen, data->n_proofs, data->commit, data->n_commits, data->nbits, &data->altgen, NULL, 0) == 1);
    }
}

static void run_test(bench_bulletproof_t *data, size_t nbits, size_t n_commits) {
    char str[64];

    data->nbits = nbits;
    data->n_commits = n_commits;

    sprintf(str, "bulletproof_prove, %i, %i, 0, ", (int)nbits, (int) n_commits);
    run_benchmark(str, bench_bulletproof_prove, bench_bulletproof_setup, bench_bulletproof_teardown, (void *)data, 10, 25);

    data->n_proofs = 1;
    sprintf(str, "bulletproof_verify, %i, %i, 0, ", (int)nbits, (int) n_commits);
    run_benchmark(str, bench_bulletproof_verify, bench_bulletproof_setup, bench_bulletproof_teardown, (void *)data, 10, 100);

    data->n_proofs = 2;
    sprintf(str, "bulletproof_verify, %i, %i, 1, ", (int)nbits, (int) n_commits);
    run_benchmark(str, bench_bulletproof_verify, bench_bulletproof_setup, bench_bulletproof_teardown, (void *)data, 10, 100);

    data->n_proofs = 5;
    sprintf(str, "bulletproof_verify, %i, %i, 4, ", (int)nbits, (int) n_commits);
    run_benchmark(str, bench_bulletproof_verify, bench_bulletproof_setup, bench_bulletproof_teardown, (void *)data, 10, 100);

    data->n_proofs = 10;
    sprintf(str, "bulletproof_verify, %i, %i, 9, ", (int)nbits, (int) n_commits);
    run_benchmark(str, bench_bulletproof_verify, bench_bulletproof_setup, bench_bulletproof_teardown, (void *)data, 10, 100);

    data->n_proofs = 20;
    sprintf(str, "bulletproof_verify, %i, %i, 19, ", (int)nbits, (int) n_commits);
    run_benchmark(str, bench_bulletproof_verify, bench_bulletproof_setup, bench_bulletproof_teardown, (void *)data, 10, 100);

    data->n_proofs = 30;
    sprintf(str, "bulletproof_verify, %i, %i, 29, ", (int)nbits, (int) n_commits);
    run_benchmark(str, bench_bulletproof_verify, bench_bulletproof_setup, bench_bulletproof_teardown, (void *)data, 10, 100);

    data->n_proofs = 40;
    sprintf(str, "bulletproof_verify, %i, %i, 39, ", (int)nbits, (int) n_commits);
    run_benchmark(str, bench_bulletproof_verify, bench_bulletproof_setup, bench_bulletproof_teardown, (void *)data, 10, 100);

    data->n_proofs = 50;
    sprintf(str, "bulletproof_verify, %i, %i, 49, ", (int)nbits, (int) n_commits);
    run_benchmark(str, bench_bulletproof_verify, bench_bulletproof_setup, bench_bulletproof_teardown, (void *)data, 10, 100);

    data->n_proofs = 75;
    sprintf(str, "bulletproof_verify, %i, %i, 74, ", (int)nbits, (int) n_commits);
    run_benchmark(str, bench_bulletproof_verify, bench_bulletproof_setup, bench_bulletproof_teardown, (void *)data, 10, 100);

    data->n_proofs = 100;
    sprintf(str, "bulletproof_verify, %i, %i, 99, ", (int)nbits, (int) n_commits);
    run_benchmark(str, bench_bulletproof_verify, bench_bulletproof_setup, bench_bulletproof_teardown, (void *)data, 10, 100);
}

int main(void) {
    bench_bulletproof_t data;
#include FILENAME

    data.ctx = secp256k1_context_create(SECP256K1_CONTEXT_SIGN | SECP256K1_CONTEXT_VERIFY);
    data.scratch = secp256k1_scratch_space_create(data.ctx, 10000000, 10000000);  /* 10M should be waay overkill */
    data.circ = secp256k1_circuit_parse(data.ctx, incl_desc);
    data.circ2 = secp256k1_circuit_parse(data.ctx, incl_desc);

#if 0
    run_benchmark(BENCHPROVE, bench_bulletproof_jubjub_prove, bench_bulletproof_circuit_setup, bench_bulletproof_teardown, (void *)&data, 5, 5);
    run_benchmark(BENCHVFY, bench_bulletproof_jubjub_verify, bench_bulletproof_circuit_setup, bench_bulletproof_teardown, (void *)&data, 100, 5);
    run_benchmark(BENCHVFY2, bench_bulletproof_jubjub_verify_plus1, bench_bulletproof_circuit_setup, bench_bulletproof_teardown, (void *)&data, 10, 5);
#endif

    run_test(&data, 8, 1);
    run_test(&data, 16, 1);
    run_test(&data, 32, 1);

    run_test(&data, 64, 1);
    run_test(&data, 64, 2);
    run_test(&data, 64, 4);
    run_test(&data, 64, 8);
    run_test(&data, 64, 16);
    run_test(&data, 64, 32);
    run_test(&data, 64, 64);
    run_test(&data, 64, 128);
    run_test(&data, 64, 256);
    run_test(&data, 64, 512);
    run_test(&data, 64, 1024);
    run_test(&data, 64, 2048);
    run_test(&data, 64, 4096);
    run_test(&data, 64, 8192);
    run_test(&data, 64, 16384);
    run_test(&data, 64, 32768);

    secp256k1_circuit_destroy(data.ctx, data.circ);
    secp256k1_circuit_destroy(data.ctx, data.circ2);
    secp256k1_scratch_space_destroy(data.scratch);
    secp256k1_context_destroy(data.ctx);
    return 0;
}
