#include <adf.h>
#include <aie_api/aie.hpp>
#include "../kernels.h"

/* ---------------- Paramètres DVB-S2 normal frame, t = 12 ---------------- */
#define BCH_M       14
#define BCH_T       12
#define BCH_N       16200
#define GF_ORDER    ((1 << BCH_M) - 1)    /* 16383                           */
#define GF_POLY     0x402B                /* x^14 + x^5 + x^3 + x + 1        */
#define NB_WORDS    ((BCH_N + 31) / 32)   /* 507 mots de 32 bits             */
#define LAST_BITS   (BCH_N - 32 * (NB_WORDS - 1))      /* 8 bits utiles      */
#define PAD_BITS    (32 * NB_WORDS - BCH_N)            /* 24 bits de bourrage*/
#define CORR_EXP    (GF_ORDER - BCH_N - PAD_BITS)      /* 159                */
#define FAIL_NOCOR  0

#define SL 16     /* lanes int16 pour les syndromes (12 utilisées)           */
#define CL 32     /* lanes int16 pour Chien (puissance de 2 ; 16 si mémoire  */
                  /* trop juste : la table de bases est 12*14*CL*2 octets)   */

template <unsigned N> using vec = aie::vector<int16, N>;

/* ====================== Arithmétique GF(2^14) scalaire ===================== */
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

static int deg_poly(const int *p, int len)
{
    int c = len - 1;
    while (p[c] == 0 && c != 0) c--;
    return c;
}

/* ====================== Arithmétique GF(2^14) vectorielle ================== */
/* a * alpha, élément par élément */
template <unsigned N>
static inline vec<N> vmulx(vec<N> a)
{
    vec<N> a2 = aie::add(a, a);                                   /* a << 1 */
    aie::mask<N> ov = aie::neq(
        aie::bit_and(a, aie::broadcast<int16, N>(1 << (BCH_M - 1))),
        aie::zeros<int16, N>());                                  /* bit 13 */
    return aie::select(a2, aie::bit_xor(a2, aie::broadcast<int16, N>(GF_POLY)), ov);
}

/* bas[k] = beta * alpha^k (par lane). Rend a * beta :
 * a est linéaire sur GF(2), donc a*beta = XOR des bas[k] pour les bits k de a */
template <unsigned N>
static inline vec<N> vmul_basis(vec<N> a, const int16 (*bas)[N])
{
    vec<N> r = aie::zeros<int16, N>();
    vec<N> z = aie::zeros<int16, N>();
    for (int k = 0; k < BCH_M; k++) {
        aie::mask<N> m = aie::neq(aie::bit_and(a, aie::broadcast<int16, N>(1 << k)), z);
        r = aie::select(r, aie::bit_xor(r, aie::load_v<N>(bas[k])), m);
    }
    return r;
}

template <unsigned N>
static void build_basis(vec<N> beta, int16 (*bas)[N])
{
    for (int k = 0; k < BCH_M; k++) {
        aie::store_v(bas[k], beta);
        beta = vmulx<N>(beta);
    }
}

/* ============================ Zones de travail ============================= */
alignas(aie::vector_decl_align) static int16 basS[BCH_M][SL];
alignas(aie::vector_decl_align) static int16 pool[BCH_T * BCH_M * CL];  /* c[i] (syndromes) puis bases (Chien) */
alignas(aie::vector_decl_align) static int16 a_beta[SL];
alignas(aie::vector_decl_align) static int16 a_b32[SL];
alignas(aie::vector_decl_align) static int16 a_corr[SL];
alignas(aie::vector_decl_align) static int16 a_raw[SL];
alignas(aie::vector_decl_align) static int16 regm[BCH_T][CL];
alignas(aie::vector_decl_align) static int16 tmpv[CL];

static int S[2 * BCH_T];
static int sigma[BCH_T + 1];
static int Bp[BCH_T + 1], B0[BCH_T + 1], Bs[BCH_T + 1];
static int racines[BCH_T];

/* ------------------------------- Syndromes ------------------------------- */
/* Vectorisé sur les 12 syndromes impairs S_1, S_3, ..., S_23 (une lane chacun).
 * Horner par mots de 32 bits :
 *     acc <- acc * beta^32  XOR  V(mot),   V(mot) = sum_i bit_i * beta^(31-i)
 * Le dernier mot est complété par 24 zéros : acc = R(beta) * beta^24, d'où la
 * correction finale par beta^(2^m - 1 - n - 24).
 * Retourne 1 si tous les syndromes sont nuls.                                */
static int compute_syndromes(const int32 *__restrict y)
{
    printf("bch decode : compute syndromes \n");
    for (int l = 0; l < SL; l++) {
        int b = (l < BCH_T) ? gf_pow(2, 2 * l + 1) : 1;
        a_beta[l] = (int16)b;
        a_b32[l]  = (int16)gf_pow(b, 32);
        a_corr[l] = (int16)gf_pow(b, CORR_EXP);
    }
    printf("bch decode : basis built \n");

    /* c[i] = beta^(31-i) */
    int16 (*cm)[SL] = (int16 (*)[SL])pool;
    build_basis<SL>(aie::load_v<SL>(a_beta), basS);
    vec<SL> c = aie::broadcast<int16, SL>(1);
    aie::store_v(cm[31], c);
    for (int i = 30; i >= 0; i--) {
        c = vmul_basis<SL>(c, basS);
        aie::store_v(cm[i], c);
    }
    printf("bch decode : syndromes computed \n");

    /* Horner */
    build_basis<SL>(aie::load_v<SL>(a_b32), basS);
    vec<SL> acc = aie::zeros<int16, SL>();
    for (int w = 0; w < NB_WORDS; w++) {
        unsigned int word = (unsigned int)y[w];
        if (w == NB_WORDS - 1)
            word &= (0xFFFFFFFFu << (32 - LAST_BITS));      /* bits hors trame = 0 */

        vec<SL> V = aie::zeros<int16, SL>();
        for (int i = 0; i < 32; i++) {
            int16 mk = (int16)(-(int)(word >> 31));         /* 0 ou -1 */
            word <<= 1;
            V = aie::bit_xor(V, aie::bit_and(aie::load_v<SL>(cm[i]),
                                             aie::broadcast<int16, SL>(mk)));
        }
        acc = aie::bit_xor(vmul_basis<SL>(acc, basS), V);
    }
    printf("bch decode : Horner done \n");

    /* correction du raccourcissement */
    build_basis<SL>(aie::load_v<SL>(a_corr), basS);
    aie::store_v(a_raw, vmul_basis<SL>(acc, basS));

    for (int s = 0; s < BCH_T; s++) S[2 * s] = a_raw[s];
    for (int k = 1; k < 2 * BCH_T; k += 2)                  /* S_2j = S_j^2 */
        S[k] = gf_mul(S[k / 2], S[k / 2]);

    int nonzero = 0;
    for (int k = 0; k < 2 * BCH_T; k++)
        if (S[k] != 0) nonzero = 1;
    return !nonzero;
}

/* ------------------- Berlekamp-Massey sans inversion --------------------- */
/* Reste scalaire : ~440 multiplications GF, négligeable devant Chien. */
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
/* CL positions évaluées en parallèle (une lane par position i = q*CL + l).
 * Terme j : reg_j = sigma[deg-j] * alpha^(-j(i+1)) ; entre deux blocs, chaque
 * terme est multiplié par alpha^(-j*CL) (constante -> table de bases).
 * Racine en i  <=>  sigma[deg] XOR (somme des termes) == 0.                  */
static int chien_search(int *deg_out)
{
    int deg = deg_poly(sigma, BCH_T + 1);
    int found = 0;
    int16 (*bc)[BCH_M][CL] = (int16 (*)[BCH_M][CL])pool;

    int g1 = gf_pow(2, GF_ORDER - 1);          /* alpha^-1 */
    int g  = 1;
    for (int j = 1; j <= deg; j++) {
        g = gf_mul(g, g1);                     /* alpha^-j */

        int t = sigma[deg - j];
        for (int l = 0; l < CL; l++) {         /* valeur initiale, lane l */
            t = gf_mul(t, g);                  /* sigma * alpha^(-j(l+1)) */
            regm[j - 1][l] = (int16)t;
        }

        int gam = g;                           /* alpha^(-j*CL) */
        for (int s = 1; s < CL; s <<= 1) gam = gf_mul(gam, gam);
        build_basis<CL>(aie::broadcast<int16, CL>((int16)gam), bc[j - 1]);
    }

    const int nblk = (BCH_N + CL - 1) / CL;
    const vec<CL> sdeg = aie::broadcast<int16, CL>((int16)sigma[deg]);

    for (int q = 0; q < nblk; q++) {
        vec<CL> acc = sdeg;
        for (int j = 0; j < deg; j++) {
            vec<CL> r = aie::load_v<CL>(regm[j]);
            acc = aie::bit_xor(acc, r);
            aie::store_v(regm[j], vmul_basis<CL>(r, bc[j]));
        }

        if (aie::reduce_min(acc) == 0) {       /* au moins une racine dans le bloc */
            aie::store_v(tmpv, acc);
            for (int l = 0; l < CL; l++) {
                int i = q * CL + l;
                if (i < BCH_N && tmpv[l] == 0) {
                    if (found < BCH_T) racines[found] = i;
                    found++;
                }
            }
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
    printf("bch decode : copy");

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
        printf("bch decode : correction done");
    }

    out_t[NB_WORDS]     = nb_roots;
    out_t[NB_WORDS + 1] = fail;
}
