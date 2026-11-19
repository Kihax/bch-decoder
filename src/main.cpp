#include "graph.h"

/* Doit correspondre au nombre de trames de data/input.txt (option -f du générateur) */
#ifndef NB_FRAMES
#define NB_FRAMES 4
#endif

bch_graph bch_top;

#if defined(__AIESIM__) || defined(__X86SIM__)
int main(int argc, char **argv)
{
    bch_top.init();
    bch_top.run(NB_FRAMES);
    bch_top.end();
    return 0;
}
#endif
