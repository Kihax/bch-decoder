#!/usr/bin/env python3
"""
Génère des trames BCH encodées (DVB-S2 normal frame, GF(2^14), t=12, n=16200,
k=16032) au format attendu par le kernel AIE `bch_decode` :

  - un mot reçu = 16200 bits, empaquetés MSB en premier dans des mots int32
    (bit Y[i] -> bit (31 - i%32) du mot i//32), soit 507 mots par trame ;
  - la trame est complétée avec des zéros jusqu'à un multiple de
    --samples-per-line, puis écrite en décimal signé (int32).

Fichiers produits :
  input.txt   trames (éventuellement bruitées) à fournir à l'AIE
  golden.txt  trames sans erreur (référence), même format
  info.txt    positions des erreurs injectées, trame par trame

Encodeur identique à BCHencode() de bch.c (systématique, mot[0] = x^(n-1)).
"""
import argparse
import random

M = 14
T = 12
N = 16200
R = M * T                       # 168 bits de redondance
K = N - R                       # 16032
NB_WORDS = (N + 31) // 32       # 507
GF_POLY = 0x402B                # x^14 + x^5 + x^3 + x + 1

# Polynômes (coefficient de x^i en position i), extraits de bch.c
POLY14 = [1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1]
POLY14_A = [
    [1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1],   # a3
    [1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1],   # a5
    [1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1],   # a7
    [1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1],   # a9
    [1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1],   # a11
    [1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1],   # a13
    [1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1],   # a15
    [1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1],   # a17
    [1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1],   # a19
    [1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1],   # a21
    [1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1],   # a23
]


def to_int(coeffs):
    return sum(c << i for i, c in enumerate(coeffs))


def clmul(a, b):
    """Produit de polynômes binaires (carry-less)."""
    r = 0
    while b:
        if b & 1:
            r ^= a
        a <<= 1
        b >>= 1
    return r


def generator_poly():
    g = to_int(POLY14)
    for p in POLY14_A:
        g = clmul(g, to_int(p))
    assert g.bit_length() - 1 == R, "degré du polynôme générateur inattendu"
    return g


# ----------------------------- arithmétique GF ------------------------------
def gf_mul(a, b):
    r = 0
    for _ in range(M):
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a >> M:
            a ^= GF_POLY
    return r


def selftest_generator(g):
    """g(x) doit s'annuler en alpha^1 ... alpha^(2t)."""
    alpha_j = 1
    for j in range(1, 2 * T + 1):
        alpha_j = gf_mul(alpha_j, 2)
        acc = 0
        for i in range(g.bit_length() - 1, -1, -1):      # Horner
            acc = gf_mul(acc, alpha_j) ^ ((g >> i) & 1)
        if acc != 0:
            raise SystemExit(f"selftest: g(alpha^{j}) != 0")
    print("selftest OK : g(x) s'annule en alpha^1..alpha^%d" % (2 * T))


# ------------------------------- encodage -----------------------------------
def bch_encode(msg_bits, g):
    """msg_bits : entier de K bits (bit K-1 = premier bit émis).
    Retourne le mot de code (entier de N bits, bit N-1 = premier bit émis)."""
    mask = (1 << R) - 1
    gl = g & mask
    v = 0
    for i in range(K - 1, -1, -1):
        fb = ((v >> (R - 1)) & 1) ^ ((msg_bits >> i) & 1)
        v = (v << 1) & mask
        if fb:
            v ^= gl
    return (msg_bits << R) | v


def check_codeword(cw, g):
    """Reste de cw(x) mod g(x) == 0 ?"""
    dg = R
    for i in range(N - 1, dg - 1, -1):
        if (cw >> i) & 1:
            cw ^= g << (i - dg)
    return cw == 0


# ------------------------------ empaquetage ---------------------------------
def pack_words(cw):
    """N bits (MSB = premier émis) -> NB_WORDS mots de 32 bits, MSB first."""
    padded = cw << (NB_WORDS * 32 - N)
    return [(padded >> (32 * (NB_WORDS - 1 - w))) & 0xFFFFFFFF for w in range(NB_WORDS)]


def to_int32(w):
    return w - (1 << 32) if w >= (1 << 31) else w


def write_frames(path, frames, spl):
    with open(path, "w") as f:
        for words in frames:
            words = list(words)
            while len(words) % spl:
                words.append(0)
            for i in range(0, len(words), spl):
                f.write(" ".join(str(to_int32(w)) for w in words[i:i + spl]) + "\n")


# ---------------------------------- main -------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-f", "--frames", type=int, default=4, help="nombre de trames")
    ap.add_argument("-e", "--errors", type=int, default=0,
                    help="erreurs binaires injectées par trame (t=12 corrigeables)")
    ap.add_argument("-s", "--seed", type=int, default=1)
    ap.add_argument("--samples-per-line", type=int, default=1,
                    help="échantillons int32 par ligne : 1 (PLIO 32 b), 2 (64 b), 4 (128 b)")
    ap.add_argument("--check", action="store_true",
                    help="vérifie chaque mot de code (reste mod g == 0)")
    ap.add_argument("-o", "--output", default="input.txt")
    ap.add_argument("--golden", default="golden.txt")
    ap.add_argument("--info", default="info.txt")
    a = ap.parse_args()

    rng = random.Random(a.seed)
    g = generator_poly()
    selftest_generator(g)

    noisy, clean, info = [], [], []
    for fr in range(a.frames):
        msg = rng.getrandbits(K)
        cw = bch_encode(msg, g)
        if a.check and not check_codeword(cw, g):
            raise SystemExit(f"trame {fr} : mot de code invalide")

        pos = sorted(rng.sample(range(N), a.errors))
        rx = cw
        for p in pos:                       # position p = indice Y[p], bit N-1-p
            rx ^= 1 << (N - 1 - p)

        clean.append(pack_words(cw))
        noisy.append(pack_words(rx))
        info.append(pos)
        print(f"trame {fr} : {a.errors} erreur(s) injectée(s)")

    write_frames(a.output, noisy, a.samples_per_line)
    write_frames(a.golden, clean, a.samples_per_line)
    with open(a.info, "w") as f:
        for fr, pos in enumerate(info):
            f.write(f"frame {fr}: {' '.join(map(str, pos))}\n")

    print(f"{a.frames} trame(s) x {NB_WORDS} mots -> {a.output} "
          f"({a.samples_per_line} échantillon(s)/ligne)")


if __name__ == "__main__":
    main()