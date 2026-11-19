#include <adf.h>
#include "kernels.h"

/* ---------------- Paramètres DVB-S2 normal frame, t = 12 ---------------- */
#define BCH_M       14
#define BCH_T       12
#define BCH_N       16200                 /* code raccourci (n < 2^m - 1)    */
#define GF_ORDER    ((1 << BCH_M) - 1)    /* 16383                           */
#define GF_POLY     0x402B                /* x^14 + x^5 + x^3 + x + 1        */
#define NB_WORDS    ((BCH_N + 31) / 32)   /* 507 mots de 32 bits             */
#define LAST_BITS   (BCH_N - 32 * (NB_WORDS - 1))   /* 8 bits utiles         */
#define FAIL_NOCOR  0

/* Tables de multiplication par constante : les éléments de GF(2^14) tiennent
 * sur 14 bits, on les garde en int16 (6 Ko au lieu de 12 Ko). */
static int16 tabL[BCH_T][128];
static int16 tabH[BCH_T][128];

static int S[2 * BCH_T];
static int sigma[BCH_T + 1];
static int Bp[BCH_T + 1], B0[BCH_T + 1], Bs[BCH_T + 1];
static int reg[BCH_T + 1];
static int racines[BCH_T];

/* ---------------------------- Arithmétique GF ---------------------------- */
static inline int gf_mulx(int a)
{
    return (a << 1) ^ (GF_POLY & -((a >> (BCH_M - 1)) & 1));
}

static int gf_mul(int a, int b)
{
    int r = 0;
    for (int k = 0; k < BCH_M; k++) {
        r ^= a & -(b & 1);
        b >>= 1;
        a = gf_mulx(a);
    }
    return r;
}

static int gf_pow(int a, int e)
{
    int r = 1;
    while (e) {
        if (e & 1) r = gf_mul(r, a);
        a = gf_mul(a, a);
        e >>= 1;
    }
    return r;
}

static void build_table(int idx, int beta)
{
    int basis[BCH_M];
    basis[0] = beta;
    for (int b = 1; b < BCH_M; b++) basis[b] = gf_mulx(basis[b - 1]);

    tabL[idx][0] = 0;
    tabH[idx][0] = 0;
    for (int b = 0; b < 7; b++) {
        for (int x = 0; x < (1 << b); x++) {
            tabL[idx][(1 << b) | x] = tabL[idx][x] ^ basis[b];
            tabH[idx][(1 << b) | x] = tabH[idx][x] ^ basis[7 + b];
        }
    }
}

static inline int mul_const(int idx, int v)
{
    return tabL[idx][v & 127] ^ tabH[idx][v >> 7];
}

static int deg_poly(const int *p, int len)
{
    int c = len - 1;
    while (p[c] == 0 && c != 0) c--;
    return c;
}

/* ------------------------------- Syndromes ------------------------------- */
/* Retourne 1 si tous les syndromes sont nuls. */
static int compute_syndromes(const int32 *__restrict y)
{
    int nonzero = 0;

    printf("bch decode : compute syndromes \n");
    for (int s = 0; s < BCH_T; s++)
        build_table(s, gf_pow(2, 2 * s + 1));
    printf("bch decode : table built \n");

    printf("bch decode : compute syndromes \n");
    for (int s = 0; s < BCH_T; s++) {
        int acc = 0;
        for (int w = 0; w < NB_WORDS; w++) {
            unsigned int word = (unsigned int)y[w];
            int nb = (w == NB_WORDS - 1) ? LAST_BITS : 32;
            for (int b = 0; b < nb; b++) {
                int bit = (int)(word >> 31);        /* bit de poids fort d'abord */
                word <<= 1;
                acc = mul_const(s, acc) ^ bit;      /* Horner */
            }
        }
        int beta = gf_pow(2, 2 * s + 1);
        S[2 * s] = gf_mul(acc, gf_pow(beta, GF_ORDER - BCH_N));
    }

    /* syndromes pairs : S_{2j} = S_j^2 */
    for (int k = 1; k < 2 * BCH_T; k += 2)
        S[k] = gf_mul(S[k / 2], S[k / 2]);

    for (int k = 0; k < 2 * BCH_T; k++)
        if (S[k] != 0) nonzero = 1;

    return !nonzero;
}

/* ------------------- Berlekamp-Massey sans inversion --------------------- */
static void key_equation_bm(void)
{
    int d = 0, gamma = 1, delta, r, degB;

    for (int j = 0; j <= BCH_T; j++) { Bp[j] = 0; B0[j] = 0; Bs[j] = 0; }
    Bp[0] = 1;
    B0[1] = 1;

    for (int i = 0; i < BCH_T - 1; i++) {
        r = 2 * i;
        degB = deg_poly(Bp, BCH_T + 1);

        delta = 0;
        for (int j = 0; j <= degB; j++)
            if (r >= j) delta ^= gf_mul(Bp[j], S[r - j]);

        for (int j = 0; j <= BCH_T; j++) {
            Bs[j] = Bp[j];
            Bp[j] = gf_mul(gamma, Bp[j]) ^ gf_mul(delta, B0[j]);
        }

        for (int j = BCH_T; j >= 2; j--) B0[j] = B0[j - 2];
        B0[0] = 0;
        B0[1] = 0;

        if ((delta != 0) && (2 * d <= r)) {
            for (int j = BCH_T; j >= 2; j--) B0[j] = Bs[j - 2];
            d = r + 1 - d;
            gamma = delta;
        }
    }

    r = 2 * BCH_T - 2;
    degB = deg_poly(Bp, BCH_T + 1);
    delta = 0;
    for (int j = 0; j <= degB; j++)
        if (r >= j) delta ^= gf_mul(Bp[j], S[r - j]);

    for (int j = 0; j <= BCH_T; j++)
        sigma[j] = gf_mul(delta, B0[j]) ^ gf_mul(gamma, Bp[j]);
}

/* ------------------------------ Chien search ----------------------------- */
static int chien_search(int *deg_out)
{
    int deg = deg_poly(sigma, BCH_T + 1);
    int found = 0;

    for (int j = 1; j <= deg; j++) {
        build_table(j - 1, gf_pow(2, GF_ORDER - j));
        reg[j] = sigma[deg - j];
    }

    for (int i = 0; i < BCH_N; i++) {
        int val = 0;
        for (int j = 1; j <= deg; j++) {
            reg[j] = mul_const(j - 1, reg[j]);
            val ^= reg[j];
        }
        if (val == sigma[deg]) {
            if (found < BCH_T) racines[found] = i;
            found++;
        }
    }
    *deg_out = deg;
    return found;
}

/* --------------------------------- Kernel -------------------------------- */
/* Bit Y[i] du mot reçu = bit (31 - i%32) du mot i/32 (MSB en premier). */
void bch_decode(adf::input_buffer<int32> & __restrict in,
                adf::output_buffer<int32> & __restrict out)
{
    printf("bch decode : start \n");
    const int32 *__restrict in_t  = in.data();
    int32       *__restrict out_t = out.data();
    printf("bch decode : load");

    int nb_roots = 0, fail = 0;

    for (unsigned i = 0; i < NB_WORDS; i++)
        out_t[i] = in_t[i];
     printf("bch decode : copy done");

    if (!compute_syndromes(in_t)) {
        printf("bch decode : compute syndromes done");
        int deg;
        key_equation_bm();
        printf("bch decode : key equation done");
        int found = chien_search(&deg);
        printf("bch decode : chien search done");

        if (found == 0) {
            nb_roots = 0;
            fail = 1;
        } else if (found != deg) {
            nb_roots = (FAIL_NOCOR == 0) ? found : 0;
            fail = 1;
        } else {
            nb_roots = found;
            fail = 0;
        }

        int ncor = (nb_roots > BCH_T) ? BCH_T : nb_roots;
        for (int i = 0; i < ncor; i++) {
            int pos = racines[i];
            out_t[pos >> 5] ^= (int32)(0x80000000u >> (pos & 31));
        }
    }

    out_t[NB_WORDS]     = nb_roots;
    out_t[NB_WORDS + 1] = fail;
}