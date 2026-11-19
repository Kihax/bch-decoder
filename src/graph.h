#ifndef GRAPH_H
#define GRAPH_H

#include <adf.h>
#include "kernels.h"

using namespace adf;

class bch_graph : public graph {
private:
    kernel k_bch;

public:
    input_plio  plio_in;
    output_plio plio_out;

    bch_graph() {
        /* kernel */
        k_bch = kernel::create(bch_decode);
        source(k_bch) = "vector/kernels.cc";
        runtime<ratio>(k_bch) = 0.9;   /* kernel lourd : une tuile pour lui seul */

        /* E/S PLIO 32 bits : 2 échantillons int16 par ligne des fichiers */
        plio_in  = input_plio::create("plio_in",   plio_32_bits, "data/input.txt");
        plio_out = output_plio::create("plio_out", plio_32_bits, "data/output.txt");

        /* connexions : une itération du graphe = une trame */
        connect(plio_in.out[0], k_bch.in[0]);
        dimensions(k_bch.in[0]) = { BCH_IN_SAMPLES };

        connect(k_bch.out[0], plio_out.in[0]);
        dimensions(k_bch.out[0]) = { BCH_OUT_SAMPLES };
    }
};

#endif