/**********************************************************************
 * Copyright (c) 2018 Andrew Poelstra                                 *
 * Distributed under the MIT software license, see the accompanying   *
 * file COPYING or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#ifndef SECP256K1_MODULE_BULLETPROOF_PARSER_IMPL
#define SECP256K1_MODULE_BULLETPROOF_PARSER_IMPL

#include <ctype.h>

typedef struct {
    size_t idx;
    int special;
    secp256k1_scalar scal;
} secp256k1_bulletproof_wmatrix_entry;

typedef struct {
    size_t size;
    secp256k1_scalar cached_sum;
    secp256k1_bulletproof_wmatrix_entry *entry;
} secp256k1_bulletproof_wmatrix_row;

struct secp256k1_bulletproof_circuit {
    size_t n_gates; /* n */
    size_t n_commits; /* m */
    size_t n_constraints; /* Q */
    secp256k1_bulletproof_wmatrix_row *wl;
    secp256k1_bulletproof_wmatrix_row *wr;
    secp256k1_bulletproof_wmatrix_row *wo;
    secp256k1_bulletproof_wmatrix_row *wv;
    secp256k1_scalar cached_c_sum;
    secp256k1_scalar *c;

    secp256k1_scalar *scratch;
};

void secp256k1_parse_scalar(secp256k1_scalar *r, const char *c, const char **end) {
    int neg = 0;
    int null = 1;
    while (isspace(*c)) {
        c++;
    }
    if (*c == '-') {
        neg = 1;
    }
    if (*c == '-' || *c == '+') {
        c++;
    }
    while (isspace(*c)) {
        c++;
    }
    secp256k1_scalar_clear(r);
    while (isdigit(*c)) {
        secp256k1_scalar digit;
        secp256k1_scalar_set_int(&digit, 10);
        secp256k1_scalar_mul(r, r, &digit);
        secp256k1_scalar_set_int(&digit, *c - '0');
        secp256k1_scalar_add(r, r, &digit);
        null = 0;
        c++;
    }
    /* interpret empty string as 1 */
    if (null == 1) {
        secp256k1_scalar_set_int(r, 1);
    }
    while (*c && *c != ',' && *c != ';' && *c != '=' && *c != 'L' && *c != 'R' && *c != 'O') {
        c++;
    }
    if (neg) {
        secp256k1_scalar_negate(r, r);
    }
    if (end != NULL) {
        *end = c;
    }
}

static secp256k1_bulletproof_circuit *secp256k1_parse_circuit(const secp256k1_context *ctx, const char *c) {
    size_t i;
    int chars_read;
    int n_gates;
    int n_commits;
    int n_constraints;
    secp256k1_bulletproof_circuit *ret = (secp256k1_bulletproof_circuit*)checked_malloc(&ctx->error_callback, sizeof(*ret));

    if (sscanf(c, "%d,%d,%d; %n", &n_gates, &n_commits, &n_constraints, &chars_read) != 3) {
        free (ret);
        return NULL;
    }
    c += chars_read;

    ret->n_gates = n_gates;
    ret->n_commits = n_commits;
    ret->n_constraints = n_constraints;
    ret->wl = (secp256k1_bulletproof_wmatrix_row *)checked_malloc(&ctx->error_callback, ret->n_gates * sizeof(*ret->wl));
    ret->wr = (secp256k1_bulletproof_wmatrix_row *)checked_malloc(&ctx->error_callback, ret->n_gates * sizeof(*ret->wr));
    ret->wo = (secp256k1_bulletproof_wmatrix_row *)checked_malloc(&ctx->error_callback, ret->n_gates * sizeof(*ret->wo));
    ret->wv = (secp256k1_bulletproof_wmatrix_row *)checked_malloc(&ctx->error_callback, ret->n_commits * sizeof(*ret->wv));
    ret->c = (secp256k1_scalar *)checked_malloc(&ctx->error_callback, ret->n_constraints * sizeof(*ret->wl));

    ret->scratch = (secp256k1_scalar *)checked_malloc(&ctx->error_callback, ret->n_constraints * sizeof(*ret->scratch));

    memset(ret->wl, 0, ret->n_gates * sizeof(*ret->wl));
    memset(ret->wr, 0, ret->n_gates * sizeof(*ret->wr));
    memset(ret->wo, 0, ret->n_gates * sizeof(*ret->wo));
    memset(ret->wv, 0, ret->n_commits * sizeof(*ret->wv));
    memset(ret->c, 0, ret->n_constraints * sizeof(*ret->c));

    for (i = 0; i < ret->n_constraints; i++) {
        int index;
        size_t j;

        j = 0;
        while (*c && *c != '=') {
            secp256k1_bulletproof_wmatrix_row *w;
            secp256k1_bulletproof_wmatrix_row *row;
            secp256k1_bulletproof_wmatrix_entry *entry;
            secp256k1_scalar mul;

            secp256k1_parse_scalar(&mul, c, &c);
            switch (*c) {
            case 'L':
                w = ret->wl;
                break;
            case 'R':
                w = ret->wr;
                break;
            case 'O':
                w = ret->wo;
                break;
            case 'V':
                w = ret->wv;
                break;
            default:
                secp256k1_circuit_destroy(ctx, ret);
                return NULL;
            }
            c++;
            if (sscanf(c, "%d %n", &index, &chars_read) != 1) {
                secp256k1_circuit_destroy(ctx, ret);
        free (ret);
                return NULL;
            }
            if ((w != ret->wv && index >= n_gates) || (w == ret->wv && index >= n_commits)) {
                secp256k1_circuit_destroy(ctx, ret);
                return NULL;
            }
            row = &w[index];

            row->size++;
            row->entry = checked_realloc(&ctx->error_callback, row->entry, row->size * sizeof(*row->entry));
            entry = &row->entry[row->size - 1];
            entry->idx = i;
            entry->special = 0;
            if (secp256k1_scalar_is_one(&mul)) {
                entry->special = 1;
            }
            secp256k1_scalar_negate(&mul, &mul);
            if (secp256k1_scalar_is_one(&mul)) {
                entry->special = -1;
            }
            secp256k1_scalar_negate(&mul, &mul);

            c += chars_read;
            entry->idx = i;
            entry->scal = mul;
            j++;
        }
        if (*c == '=') {
            c++;
            secp256k1_parse_scalar(&ret->c[i], c, &c);
            if (*c != ';') {
                secp256k1_circuit_destroy(ctx, ret);
                return NULL;
            }
            c++;
        } else {
            secp256k1_circuit_destroy(ctx, ret);
            return NULL;
        }
    }

    return ret;
}

static void secp256k1_wmatrix_row_compress(secp256k1_bulletproof_wmatrix_row *row, const secp256k1_scalar *zn) {
    size_t j;
    secp256k1_scalar_clear(&row->cached_sum);
    for (j = 0; j < row->size; j++) {
        secp256k1_scalar term;
        switch (row->entry[j].special) {
        case -1:
            secp256k1_scalar_negate(&term, &zn[row->entry[j].idx]);
            break;
        case 1:
            term = zn[row->entry[j].idx];
            break;
        default:
            secp256k1_scalar_mul(&term, &row->entry[j].scal, &zn[row->entry[j].idx]);
            break;
        }
        secp256k1_scalar_add(&row->cached_sum, &row->cached_sum, &term);
    }
}

static void secp256k1_circuit_compress(secp256k1_bulletproof_circuit *circ, const secp256k1_scalar *z) {
    size_t i;

    circ->scratch[0] = *z;
    for (i = 1; i < circ->n_constraints; i++) {
        secp256k1_scalar_mul(&circ->scratch[i], &circ->scratch[i-1], z);
    }

    for (i = 0; i < circ->n_gates; i++) {
        secp256k1_wmatrix_row_compress(&circ->wl[i], circ->scratch);
        secp256k1_wmatrix_row_compress(&circ->wr[i], circ->scratch);
        secp256k1_wmatrix_row_compress(&circ->wo[i], circ->scratch);
    }

    for (i = 0; i < circ->n_commits; i++) {
        secp256k1_wmatrix_row_compress(&circ->wv[i], circ->scratch);
    }

    secp256k1_scalar_clear(&circ->cached_c_sum);
    for (i = 0; i < circ->n_constraints; i++) {
        secp256k1_scalar term;
        secp256k1_scalar_mul(&term, &circ->c[i], &circ->scratch[i]);
        secp256k1_scalar_add(&circ->cached_c_sum, &circ->cached_c_sum, &term);
    }
}

static void secp256k1_circuit_destroy_impl(secp256k1_bulletproof_circuit *circ) {
    if (circ != NULL) {
        free(circ->wl);
        free(circ->wr);
        free(circ->wo);
        free(circ->wv);
        free(circ->c);
        free(circ);
    }
}

#endif
