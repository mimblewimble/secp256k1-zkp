/**********************************************************************
 * Copyright (c) 2017 Andrew Poelstra                                 *
 * Distributed under the MIT software license, see the accompanying   *
 * file COPYING or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#ifndef SECP256K1_MODULE_BULLETPROOF_UTIL
#define SECP256K1_MODULE_BULLETPROOF_UTIL

SECP256K1_INLINE static size_t secp256k1_ceil_lg(size_t n) {
    VERIFY_CHECK(n > 0);
    switch (n) {
    case 1: return 0;
    case 2: return 1;
    case 3: return 2;
    case 4: return 2;
    case 5: return 3;
    case 6: return 3;
    case 7: return 3;
    case 8: return 3;
    default: {
        size_t i = 0;
        n--;
        while (n > 0) {
            n /= 2;
            i++;
        }
        return i;
    }
    }
}

SECP256K1_INLINE static void secp256k1_bulletproof_genrand_pair(secp256k1_rfc6979_hmac_sha256 *rng, secp256k1_scalar *out1, secp256k1_scalar *out2) {
    unsigned char tmp[32];
    int overflow;

    secp256k1_rfc6979_hmac_sha256_generate(rng, tmp, 32);
    do {
        secp256k1_rfc6979_hmac_sha256_generate(rng, tmp, 32);
        secp256k1_scalar_set_b32(out1, tmp, &overflow);
    } while (overflow || secp256k1_scalar_is_zero(out1));

    secp256k1_rfc6979_hmac_sha256_generate(rng, tmp, 32);
    do {
        secp256k1_rfc6979_hmac_sha256_generate(rng, tmp, 32);
        secp256k1_scalar_set_b32(out2, tmp, &overflow);
    } while (overflow || secp256k1_scalar_is_zero(out2));
}

SECP256K1_INLINE static void secp256k1_bulletproof_serialize_points(unsigned char *out, secp256k1_ge *pt, size_t n) {
    const size_t bitveclen = (n + 7) / 8;
    size_t i;

    memset(out, 0, bitveclen);
    for (i = 0; i < n; i++) {
        secp256k1_fe pointx;
        pointx = pt[i].x;
        secp256k1_fe_normalize(&pointx);
        secp256k1_fe_get_b32(&out[bitveclen + i*32], &pointx);
        if (!secp256k1_fe_is_quad_var(&pt[i].y)) {
            out[i/8] |= (1ull << (i % 8));
        }
    }
}

SECP256K1_INLINE static void secp256k1_bulletproof_deserialize_point(secp256k1_ge *pt, const unsigned char *data, size_t i, size_t n) {
    const size_t bitveclen = (n + 7) / 8;
    const size_t offset = bitveclen + i*32;
    secp256k1_fe fe;

    secp256k1_fe_set_b32(&fe, &data[offset]);
    secp256k1_ge_set_xquad(pt, &fe);
    if (data[i / 8] & (1 << (i % 8))) {
        secp256k1_ge_neg(pt, pt);
    }
}

static void secp256k1_bulletproof_update_commit(unsigned char *commit, const secp256k1_ge *lpt, const secp256k1_ge *rpt) {
    secp256k1_fe pointx;
    secp256k1_sha256 sha256;
    unsigned char lrparity;
    lrparity = (!secp256k1_fe_is_quad_var(&lpt->y) << 1) + !secp256k1_fe_is_quad_var(&rpt->y);
    secp256k1_sha256_initialize(&sha256);
    secp256k1_sha256_write(&sha256, commit, 32);
    secp256k1_sha256_write(&sha256, &lrparity, 1);
    pointx = lpt->x;
    secp256k1_fe_normalize(&pointx);
    secp256k1_fe_get_b32(commit, &pointx);
    secp256k1_sha256_write(&sha256, commit, 32);
    pointx = rpt->x;
    secp256k1_fe_normalize(&pointx);
    secp256k1_fe_get_b32(commit, &pointx);
    secp256k1_sha256_write(&sha256, commit, 32);
    secp256k1_sha256_finalize(&sha256, commit);
}

static void secp256k1_bulletproof_update_commit_n(unsigned char *commit, const secp256k1_ge *pt, size_t n) {
    secp256k1_sha256 sha256;
    unsigned char lrparity = 0;
    size_t i;

    VERIFY_CHECK(n < 8);

    for (i = 0; i < n; i++) {
        lrparity |= secp256k1_fe_is_quad_var(&pt[i].y) << i;
    }

    secp256k1_sha256_initialize(&sha256);
    secp256k1_sha256_write(&sha256, commit, 32);
    secp256k1_sha256_write(&sha256, &lrparity, 1);
    for (i = 0; i < n; i++) {
        secp256k1_fe pointx;
        pointx = pt[i].x;
        secp256k1_fe_normalize(&pointx);
        secp256k1_fe_get_b32(commit, &pointx);
        secp256k1_sha256_write(&sha256, commit, 32);
    }
    secp256k1_sha256_finalize(&sha256, commit);
}

/* Convenience function to compute blind*G + sum_i (s[i] * gen[i])
 * If G is passed as NULL, use the standard generator. While in the
 * standard-generator case we could use ecmult_gen rather than
 * ecmult_const, we don't. This function is only used during proof
 * generation so performance is not critical.
 *
 * If `blind` is NULL it is treated as zero.
 *
 * This function is not constant-time with respect to the NULLness
 * of its inputs. NULLness should never be correlated with secret data.
 */
static void secp256k1_bulletproof_vector_commit(secp256k1_gej *r, const secp256k1_scalar *s, const secp256k1_ge *gen, size_t n, const secp256k1_scalar *blind, const secp256k1_ge *g) {
    secp256k1_scalar zero;
    secp256k1_ge rge;

    if (g == NULL) {
        g = &secp256k1_ge_const_g;
    }
    if (blind == NULL) {
        secp256k1_scalar_clear(&zero);
        blind = &zero;
    }
 
    /* Multiply by blinding factor */
    secp256k1_ecmult_const(r, g, blind, 256);

    /* Do the non-blind sum, going through contortions to avoid adding infinities */
    while (n--) {
        int inf;
        secp256k1_ge tmpge;
        secp256k1_ge negg;
        secp256k1_gej tmpj;

        /* Add G, undoing it if this causes rge == infinity */
        secp256k1_ge_set_gej(&tmpge, r);
        secp256k1_gej_add_ge(r, r, g);
        secp256k1_ge_set_gej(&rge, r);

        inf = secp256k1_ge_is_infinity(&rge);
        secp256k1_fe_cmov(&rge.x, &tmpge.x, inf);
        secp256k1_fe_cmov(&rge.y, &tmpge.y, inf);
        rge.infinity = 0;

        /* Add the next term to our now-guaranteed-noninfinite R */
        secp256k1_ecmult_const(&tmpj, &gen[n], &s[n], 256);
        secp256k1_gej_add_ge(r, &tmpj, &rge); /* here tmpj may be infinite but tmpge won't be */

        /* Subtract G, undoing it if we undid the addition above */
        secp256k1_ge_neg(&negg, g);
        secp256k1_ge_set_gej(&tmpge, r);
        secp256k1_gej_add_ge(r, r, &negg);
        secp256k1_ge_set_gej(&rge, r);

        secp256k1_fe_cmov(&rge.x, &tmpge.x, inf);
        secp256k1_fe_cmov(&rge.y, &tmpge.y, inf);
        rge.infinity = rge.infinity * (1 - inf) + tmpge.infinity * inf;
    }
}

#endif
