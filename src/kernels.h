#ifndef KERNELS_H
#define KERNELS_H

#include <adf.h>

/* Tailles des buffers, en échantillons int32 (PLIO 32 bits : 1 échantillon par ligne).
 *
 * Entrée  : 507 mots = 16200 bits empaquetés MSB first (le dernier mot n'en
 *           contient que 8, le reste est à 0).
 * Sortie  : 507 mots corrigés, mot 507 = nombre d'erreurs corrigées,
 *           mot 508 = flag Fail.                                           */
#define BCH_IN_SAMPLES   507
#define BCH_OUT_SAMPLES  509

void bch_decode(adf::input_buffer<int32> & __restrict in,
                adf::output_buffer<int32> & __restrict out);

#endif