/**********************************************************************
 * Copyright (c) 2017 Andrew Poelstra                                 *
 * Distributed under the MIT software license, see the accompanying   *
 * file COPYING or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#ifndef _SECP256K1_MODULE_AGGSIG_TESTS_
#define _SECP256K1_MODULE_AGGSIG_TESTS_

#include "secp256k1_aggsig.h"

void test_aggsig_api(void) {
    /* Setup contexts that just count errors */
    secp256k1_context *none = secp256k1_context_create(SECP256K1_CONTEXT_NONE);
    secp256k1_context *sign = secp256k1_context_create(SECP256K1_CONTEXT_SIGN);
    secp256k1_context *vrfy = secp256k1_context_create(SECP256K1_CONTEXT_VERIFY);
    secp256k1_context *both = secp256k1_context_create(SECP256K1_CONTEXT_SIGN | SECP256K1_CONTEXT_VERIFY);
    secp256k1_scratch_space *scratch = secp256k1_scratch_space_create(none, 1024*4096);
    secp256k1_scalar tmp_s;
    unsigned char seckeys[5][32];
    secp256k1_pubkey pubkeys[5];
    secp256k1_aggsig_partial_signature partials[5];
    unsigned char sub_result[64];
    unsigned char sub_result_alt[64];
    int sub_return_val;
    int sub_check;
    secp256k1_aggsig_context *aggctx;
    unsigned char seed[32] = { 1, 2, 3, 4, 0 };
    unsigned char sig[64];
    unsigned char sig2[64];
    unsigned char* sigs[2];
    unsigned char combined_sig[64];
    unsigned char sec_nonces[5][32];
    secp256k1_pubkey pub_nonces[5];
    unsigned char orig_sig;
    unsigned char orig_msg;
    unsigned char *msg = seed;  /* shh ;) */
    const secp256k1_pubkey* pubkey_combiner[2];
    secp256k1_pubkey combiner_sum;
    secp256k1_pubkey combiner_sum_2;
    int32_t ecount = 0;

    size_t i;
    size_t j;

    const unsigned char no_qr_const_1[] = {
        0xCA,0x69,0xE3,0x16,0xC1,0xF5,0x4E,0xFE,
        0x6C,0x23,0xFF,0x8F,0x1D,0x29,0xB9,0xB3,
        0xE2,0x22,0xA6,0x79,0x9D,0x3E,0xB1,0x1C,
        0xF5,0xC,0x8A,0x5,0x26,0x47,0xE8,0xDF
    };

    const unsigned char no_qr_const_2[] = {
        0xA1,0xB9,0xFF,0xD9,0x7E,0xC9,0x11,0xD,
        0x53,0xA2,0xCC,0x88,0x63,0x9A,0x31,0x9C,
        0x84,0x36,0xE5,0x0,0xB4,0x1C,0x18,0x2A,
        0xC4,0x9F,0x2,0x68,0x64,0x49,0x3F,0xC7,
    };

    secp256k1_context_set_error_callback(none, counting_illegal_callback_fn, &ecount);
    secp256k1_context_set_error_callback(sign, counting_illegal_callback_fn, &ecount);
    secp256k1_context_set_error_callback(vrfy, counting_illegal_callback_fn, &ecount);
    secp256k1_context_set_error_callback(both, counting_illegal_callback_fn, &ecount);
    secp256k1_context_set_illegal_callback(none, counting_illegal_callback_fn, &ecount);
    secp256k1_context_set_illegal_callback(sign, counting_illegal_callback_fn, &ecount);
    secp256k1_context_set_illegal_callback(vrfy, counting_illegal_callback_fn, &ecount);
    secp256k1_context_set_illegal_callback(both, counting_illegal_callback_fn, &ecount);

    for (i = 0; i < 5; i++) {
        random_scalar_order_test(&tmp_s);
        secp256k1_scalar_get_b32(seckeys[i], &tmp_s);
        CHECK(secp256k1_ec_pubkey_create(ctx, &pubkeys[i], seckeys[i]) == 1);
    }

    aggctx = secp256k1_aggsig_context_create(none, pubkeys, 5, seed);
    CHECK(ecount == 0);
    CHECK(aggctx != NULL);
    secp256k1_aggsig_context_destroy(NULL);  /* harmless */
    secp256k1_aggsig_context_destroy(aggctx);

    aggctx = secp256k1_aggsig_context_create(none, pubkeys, 0, seed);
    CHECK(ecount == 0);
    CHECK(aggctx != NULL);
    secp256k1_aggsig_context_destroy(aggctx);

    aggctx = secp256k1_aggsig_context_create(none, pubkeys, 0, NULL);
    CHECK(ecount == 1);
    CHECK(aggctx == NULL);
    aggctx = secp256k1_aggsig_context_create(none, NULL, 0, seed);
    CHECK(ecount == 2);
    CHECK(aggctx == NULL);

    aggctx = secp256k1_aggsig_context_create(none, pubkeys, 5, seed);

    CHECK(!secp256k1_aggsig_generate_nonce(none, aggctx, 0));
    CHECK(ecount == 3);
    CHECK(secp256k1_aggsig_generate_nonce(sign, aggctx, 0));
    CHECK(ecount == 3);
    CHECK(!secp256k1_aggsig_generate_nonce(vrfy, aggctx, 0));
    CHECK(ecount == 4);
    CHECK(!secp256k1_aggsig_generate_nonce(both, aggctx, 0));  /* double-generate, not API error */
    CHECK(ecount == 4);
    CHECK(secp256k1_aggsig_generate_nonce(both, aggctx, 1));
    CHECK(ecount == 4);
    CHECK(!secp256k1_aggsig_generate_nonce(both, NULL, 2));
    CHECK(ecount == 5);
    CHECK(!secp256k1_aggsig_generate_nonce(both, aggctx, 5));   /* out of range, API error */
    CHECK(ecount == 6);

    CHECK(!secp256k1_aggsig_partial_sign(both, aggctx, &partials[0], msg, seckeys[0], 0));  /* not all nonces generated, not API error */
    CHECK(secp256k1_aggsig_generate_nonce(both, aggctx, 2));
    CHECK(secp256k1_aggsig_generate_nonce(both, aggctx, 3));
    CHECK(secp256k1_aggsig_generate_nonce(both, aggctx, 4));
    CHECK(secp256k1_aggsig_partial_sign(both, aggctx, &partials[0], msg, seckeys[0], 0));
    CHECK(!secp256k1_aggsig_partial_sign(both, aggctx, &partials[0], msg, seckeys[0], 0));  /* double sign, not API error */
    CHECK(ecount == 6);

    CHECK(!secp256k1_aggsig_partial_sign(none, aggctx, &partials[1], msg, seckeys[1], 1));
    CHECK(ecount == 7);
    CHECK(!secp256k1_aggsig_partial_sign(vrfy, aggctx, &partials[1], msg, seckeys[1], 1));
    CHECK(ecount == 8);
    CHECK(secp256k1_aggsig_partial_sign(sign, aggctx, &partials[1], msg, seckeys[1], 1));
    CHECK(ecount == 8);
    CHECK(!secp256k1_aggsig_partial_sign(sign, aggctx, NULL, msg, seckeys[2], 2));
    CHECK(ecount == 9);
    CHECK(!secp256k1_aggsig_partial_sign(sign, aggctx, &partials[2], NULL, seckeys[2], 2));
    CHECK(ecount == 10);
    CHECK(!secp256k1_aggsig_partial_sign(sign, aggctx, &partials[2], msg, NULL, 2));
    CHECK(ecount == 11);
    CHECK(!secp256k1_aggsig_partial_sign(sign, aggctx, &partials[2], msg, seckeys[2], 5));  /* out of range, API error */
    CHECK(ecount == 12);
    CHECK(secp256k1_aggsig_partial_sign(both, aggctx, &partials[2], msg, seckeys[2], 2));
    CHECK(secp256k1_aggsig_partial_sign(both, aggctx, &partials[3], msg, seckeys[3], 3));
    CHECK(secp256k1_aggsig_partial_sign(both, aggctx, &partials[4], msg, seckeys[4], 4));
    CHECK(ecount == 12);

    CHECK(secp256k1_aggsig_combine_signatures(none, aggctx, sig, partials, 5));
    CHECK(!secp256k1_aggsig_combine_signatures(none, aggctx, sig, partials, 4)); /* wrong sig count, not API error (should it be?)  */
    CHECK(!secp256k1_aggsig_combine_signatures(none, aggctx, sig, partials, 0));
    CHECK(ecount == 12);
    CHECK(!secp256k1_aggsig_combine_signatures(none, NULL, sig, partials, 5));
    CHECK(ecount == 13);
    CHECK(!secp256k1_aggsig_combine_signatures(none, aggctx, NULL, partials, 5));
    CHECK(ecount == 14);
    CHECK(!secp256k1_aggsig_combine_signatures(none, aggctx, sig, NULL, 5));
    CHECK(ecount == 15);

    memset(sig, 0, sizeof(sig));
    CHECK(!secp256k1_aggsig_verify(vrfy, scratch, sig, msg, pubkeys, 5));
    CHECK(secp256k1_aggsig_combine_signatures(none, aggctx, sig, partials, 5));
    CHECK(!secp256k1_aggsig_verify(vrfy, scratch, sig, msg, pubkeys, 4));
    CHECK(!secp256k1_aggsig_verify(vrfy, scratch, sig, msg, pubkeys, 0));
    CHECK(secp256k1_aggsig_verify(vrfy, scratch, sig, msg, pubkeys, 5));
    CHECK(ecount == 15);

    CHECK(!secp256k1_aggsig_verify(none, scratch, sig, msg, pubkeys, 5));
    CHECK(ecount == 16);

    CHECK(!secp256k1_aggsig_verify(vrfy, NULL, sig, msg, pubkeys, 5));
    CHECK(ecount == 17);
    CHECK(!secp256k1_aggsig_verify(vrfy, scratch, NULL, msg, pubkeys, 5));
    CHECK(ecount == 18);
    CHECK(!secp256k1_aggsig_verify(vrfy, scratch, sig, NULL, pubkeys, 5));
    CHECK(ecount == 19);
    CHECK(!secp256k1_aggsig_verify(vrfy, scratch, sig, msg, NULL, 5));
    CHECK(ecount == 20);

    /* Test single api */
    memset(sig, 0, sizeof(sig));
    CHECK(secp256k1_aggsig_sign_single(sign, sig, msg, seckeys[0], NULL, NULL, NULL, NULL, NULL, seed));
    CHECK(ecount == 20);
    CHECK(secp256k1_aggsig_verify_single(vrfy, sig, msg, NULL, &pubkeys[0], NULL, NULL, 0));
    CHECK(!secp256k1_aggsig_verify_single(vrfy, sig, msg, NULL, &pubkeys[1], NULL, NULL, 0));
    CHECK(!secp256k1_aggsig_verify_single(vrfy, sig, msg, NULL, &pubkeys[0], &pubkeys[1], NULL, 0));
    orig_sig=sig[0];
    sig[0]=99;
    CHECK(!secp256k1_aggsig_verify_single(vrfy, sig, msg, NULL, &pubkeys[0], NULL, NULL, 0));
    sig[0]=orig_sig;
    CHECK(secp256k1_aggsig_verify_single(vrfy, sig, msg, NULL, &pubkeys[0], NULL, NULL, 0));
    orig_msg=msg[0];
    msg[0]=99;
    CHECK(!secp256k1_aggsig_verify_single(vrfy, sig, msg, NULL, &pubkeys[0], NULL, NULL, 0));
    msg[0]=orig_msg;

    /* Test single api with pubkey in e */
    memset(sig, 0, sizeof(sig));
    CHECK(secp256k1_aggsig_sign_single(sign, sig, msg, seckeys[0], NULL, NULL, NULL, NULL, &pubkeys[2], seed));
    CHECK(secp256k1_aggsig_verify_single(vrfy, sig, msg, NULL, &pubkeys[0], &pubkeys[2], NULL, 0));
    CHECK(!secp256k1_aggsig_verify_single(vrfy, sig, msg, NULL, &pubkeys[0], &pubkeys[3], NULL, 0));
    CHECK(!secp256k1_aggsig_verify_single(vrfy, sig, msg, NULL, &pubkeys[0], NULL, NULL, 0));
    CHECK(!secp256k1_aggsig_verify_single(vrfy, sig, msg, NULL, &pubkeys[1], &pubkeys[2], NULL, 0));
    CHECK(!secp256k1_aggsig_verify_single(vrfy, sig, msg, NULL, &pubkeys[1], NULL, NULL, 0));
    orig_sig=sig[0];
    sig[0]=99;
    CHECK(!secp256k1_aggsig_verify_single(vrfy, sig, msg, NULL, &pubkeys[0], &pubkeys[2], NULL, 0));
    sig[0]=orig_sig;
    CHECK(secp256k1_aggsig_verify_single(vrfy, sig, msg, NULL, &pubkeys[0], &pubkeys[2], NULL, 0));
    msg[0]=99;
    CHECK(!secp256k1_aggsig_verify_single(vrfy, sig, msg, NULL, &pubkeys[0], &pubkeys[2], NULL, 0));

    /* Overriding sec nonce */
    memset(sig, 0, sizeof(sig));
    CHECK(secp256k1_aggsig_sign_single(sign, sig, msg, seckeys[0], seckeys[1], NULL, NULL, NULL, NULL, seed));
    CHECK(secp256k1_aggsig_verify_single(vrfy, sig, msg, NULL, &pubkeys[0], NULL, NULL, 0));

    /* Overriding sec nonce and pub nonce encoded in e */
    memset(sig, 0, sizeof(sig));
    CHECK(secp256k1_aggsig_sign_single(sign, sig, msg, seckeys[0], seckeys[1], NULL, &pubkeys[3], NULL, NULL, seed));
    CHECK(secp256k1_aggsig_verify_single(vrfy, sig, msg, &pubkeys[3], &pubkeys[0], NULL, NULL, 0));

    /* Add extra key to the signature */
    memset(sig, 0, sizeof(sig));
    CHECK(secp256k1_aggsig_sign_single(sign, sig, msg, seckeys[0], seckeys[1], seckeys[2], &pubkeys[3], &pubkeys[3], &pubkeys[4], seed));
    /* Check that it doesn't verify without the extra key */
    CHECK(!secp256k1_aggsig_verify_single(vrfy, sig, msg, &pubkeys[3], &pubkeys[0], &pubkeys[4], NULL, 1));
    /* And that it does with */
    CHECK(secp256k1_aggsig_verify_single(vrfy, sig, msg, &pubkeys[3], &pubkeys[0], &pubkeys[4], &pubkeys[2], 1));

    /* Testing aggsig exchange algorithm for Grin */
    /* ****************************************** */

    for (i=0;i<20;i++){
        memset(sig, 0, sizeof(sig));
        memset(sig, 0, sizeof(sig2));
        memset(sig, 0, sizeof(combined_sig));

        /* Create a couple of nonces */
        /* Randomise seed to make it more interesting */
        random_scalar_order_test(&tmp_s);
        secp256k1_scalar_get_b32(seed, &tmp_s);
        CHECK(secp256k1_aggsig_export_secnonce_single(sign, sec_nonces[0], seed));
        random_scalar_order_test(&tmp_s);
        secp256k1_scalar_get_b32(seed, &tmp_s);
        CHECK(secp256k1_aggsig_export_secnonce_single(sign, sec_nonces[1], seed));

        for (j = 0; j < 2; j++) {
            CHECK(secp256k1_ec_pubkey_create(ctx, &pub_nonces[j], sec_nonces[j]) == 1);
        }

        /* Randomize keys */
        for (j = 0; j < 2; j++) {
            random_scalar_order_test(&tmp_s);
            secp256k1_scalar_get_b32(seckeys[j], &tmp_s);
            CHECK(secp256k1_ec_pubkey_create(ctx, &pubkeys[j], seckeys[j]) == 1);
        }

        /* Combine pubnonces */
        pubkey_combiner[0]=&pub_nonces[0];
        pubkey_combiner[1]=&pub_nonces[1];
        CHECK(secp256k1_ec_pubkey_combine(ctx, &combiner_sum, pubkey_combiner, 2) == 1);

        /* Combine pubkeys */
        pubkey_combiner[0]=&pubkeys[0];
        pubkey_combiner[1]=&pubkeys[1];
        CHECK(secp256k1_ec_pubkey_combine(ctx, &combiner_sum_2, pubkey_combiner, 2) == 1);

        /* Create 2 partial signatures (Sender, Receiver)*/
        CHECK(secp256k1_aggsig_sign_single(sign, sig, msg, seckeys[0], sec_nonces[0], NULL, &combiner_sum, &combiner_sum, &combiner_sum_2, seed));

        /* Receiver verifies sender's Sig and signs */
        CHECK(secp256k1_aggsig_verify_single(vrfy, sig, msg, &combiner_sum, &pubkeys[0], &combiner_sum_2, NULL, 1));
        CHECK(secp256k1_aggsig_sign_single(sign, sig2, msg, seckeys[1], sec_nonces[1], NULL, &combiner_sum, &combiner_sum, &combiner_sum_2, seed));
        /* sender verifies receiver's Sig then creates final combined sig */
        CHECK(secp256k1_aggsig_verify_single(vrfy, sig2, msg, &combiner_sum, &pubkeys[1], &combiner_sum_2, NULL, 1));

        sigs[0] = sig;
        sigs[1] = sig2;
        /* Add 2 sigs and nonces */
        CHECK(secp256k1_aggsig_add_signatures_single(sign, combined_sig, (const unsigned char **) sigs, 2, &combiner_sum));

        /* Ensure added sigs verify properly (with and without providing nonce_sum */
        CHECK(secp256k1_aggsig_verify_single(vrfy, combined_sig, msg, &combiner_sum, &combiner_sum_2, &combiner_sum_2, NULL, 0));
        CHECK(secp256k1_aggsig_verify_single(vrfy, combined_sig, msg, NULL, &combiner_sum_2, &combiner_sum_2, NULL, 0));

        /* And anything else doesn't */
        CHECK(!secp256k1_aggsig_verify_single(vrfy, combined_sig, msg, &combiner_sum, &combiner_sum_2, NULL, NULL, 0));
        CHECK(!secp256k1_aggsig_verify_single(vrfy, combined_sig, msg, NULL, &pub_nonces[1], NULL, NULL, 0));
        CHECK(!secp256k1_aggsig_verify_single(vrfy, combined_sig, msg, NULL, &pub_nonces[1], &combiner_sum_2, NULL, 0));
        CHECK(!secp256k1_aggsig_verify_single(vrfy, combined_sig, msg, &pub_nonces[0], &combiner_sum_2, NULL, NULL, 0));
        CHECK(!secp256k1_aggsig_verify_single(vrfy, combined_sig, msg, &pub_nonces[0], &combiner_sum_2, &combiner_sum_2, NULL, 0));
        msg[0]=1;
        msg[1]=2;
        msg[2]=3;
        CHECK(!secp256k1_aggsig_verify_single(vrfy, combined_sig, msg, NULL, &combiner_sum_2, NULL, NULL, 0));
        CHECK(!secp256k1_aggsig_verify_single(vrfy, combined_sig, msg, NULL, &combiner_sum_2, &combiner_sum_2, NULL, 0));

        /* Test subtracting partial sig from the completed sig (as used in early payment proofs)*/
        sub_return_val = secp256k1_aggsig_subtract_partial_signature(none, sub_result, sub_result_alt, combined_sig, sigs[0]);
        
        /* Check non-0 return */
        CHECK(sub_return_val);

        /* And this will never be -1 */
        CHECK(sub_return_val != -1);

        /* If there's only one possible result, it should be the only one that needs to be tested */
        if (sub_return_val == 1) {
            for (j = 0; j < 64; j++){
                CHECK(sub_result[j] == sigs[1][j]);
            }
        }

        /* if there are two possible results, both potentially need to be checked */
        if (sub_return_val == 2) {
            sub_check = 1;
            for (j = 0; j < 64; j++){
                if (sub_result[j] != sigs[1][j]) {
                    sub_check = 0;
                }
            }
            if (!sub_check) {
                sub_check = 1;
                for (j = 0; j < 64; j++){
                    if (sub_result_alt[j] != sigs[1][j]) {
                        sub_check = 0;
                    }
                }
            }
            CHECK(sub_check);
        }
    }
    /* check subtraction function returns -1 for values known to result in no possiblity being a quadratic residue */
    memcpy(combined_sig, no_qr_const_1, 32);
    memcpy(sigs[0], no_qr_const_2, 32);
    sub_return_val = secp256k1_aggsig_subtract_partial_signature(none, sub_result, sub_result_alt, combined_sig, sigs[0]);
    CHECK(sub_return_val == -1);

    /*** End aggsig for Grin exchange test ***/

    /* cleanup */
    secp256k1_aggsig_context_destroy(aggctx);
    secp256k1_scratch_space_destroy(scratch);
    secp256k1_context_destroy(none);
    secp256k1_context_destroy(sign);
    secp256k1_context_destroy(vrfy);
    secp256k1_context_destroy(both);
}

#define N_KEYS 200
void test_aggsig_onesigner(void) {
    secp256k1_pubkey pubkeys[N_KEYS];
    unsigned char seckeys[N_KEYS][32];
    secp256k1_aggsig_partial_signature partials[N_KEYS];
    const size_t n_pubkeys = sizeof(pubkeys) / sizeof(pubkeys[0]);
    secp256k1_scalar tmp_s;
    size_t i;
    size_t n_signers[] = { 1, 2, N_KEYS / 5, N_KEYS - 1, N_KEYS };
    const size_t n_n_signers = sizeof(n_signers) / sizeof(n_signers[0]);
    secp256k1_scratch_space *scratch = secp256k1_scratch_space_create(ctx, 1024*4096);

    unsigned char msg[32];

    random_scalar_order_test(&tmp_s);
    secp256k1_scalar_get_b32(msg, &tmp_s);

    for (i = 0; i < n_pubkeys; i++) {
        random_scalar_order_test(&tmp_s);
        secp256k1_scalar_get_b32(seckeys[i], &tmp_s);
        CHECK(secp256k1_ec_pubkey_create(ctx, &pubkeys[i], seckeys[i]) == 1);
    }

    for (i = 0; i < n_n_signers; i++) {
        size_t j;
        unsigned char seed[32];
        unsigned char sig[64];
        secp256k1_aggsig_context *aggctx;

        random_scalar_order_test(&tmp_s);
        secp256k1_scalar_get_b32(seed, &tmp_s);
        aggctx = secp256k1_aggsig_context_create(ctx, pubkeys, n_signers[i], seed);

        /* all nonces must be generated before signing */
        for (j = 0; j < n_signers[i]; j++) {
            CHECK(secp256k1_aggsig_generate_nonce(ctx, aggctx, j));
        }
        for (j = 0; j < n_signers[i]; j++) {
            CHECK(secp256k1_aggsig_partial_sign(ctx, aggctx, &partials[j], msg, seckeys[j], j));
        }
        CHECK(secp256k1_aggsig_combine_signatures(ctx, aggctx, sig, partials, n_signers[i]));
        CHECK(secp256k1_aggsig_verify(ctx, scratch, sig, msg, pubkeys, n_signers[i]));
        /* Make sure verification with 0 pubkeys fails without Bad Things happenings */
        CHECK(!secp256k1_aggsig_verify(ctx, scratch, sig, msg, pubkeys, 0));

        secp256k1_aggsig_context_destroy(aggctx);
    }

    secp256k1_scratch_space_destroy(scratch);
}
#undef N_KEYS

void run_aggsig_tests(void) {
    test_aggsig_api();
    test_aggsig_onesigner();
}

#endif
