/***********************************************************************
 * Copyright (c) 2017 Gregory Maxwell                                  *
 * Distributed under the MIT software license, see the accompanying    *
 * file COPYING or http://www.opensource.org/licenses/mit-license.php. *
 ***********************************************************************/

#ifndef _SECP256K1_SWITCH_IMPL_H_
#define _SECP256K1_SWITCH_IMPL_H_

#include "switch.h"

/* sec * G3 */
SECP256K1_INLINE static void secp256k1_switch_ecmult(secp256k1_gej *rj,
                                                     const secp256k1_scalar *sec,
                                                     const secp256k1_ge* genp) {
    secp256k1_ecmult_const(rj, genp, sec, 64);
}

#endif
