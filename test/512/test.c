/* ChainHash-512 tests.
 * 1. test/512/vectors.txt (independent Python oracle, key[i] = (73i+11) & 255): the reference and every
 *    available backend, one-shot and streaming; the same key given as words.
 * 2. Every backend and streaming (random splits) against chainhash512_reference on random inputs: every length
 *    up to 4200, the last-region shapes q = 1..8 in regions 1-3, random lengths up to MAXLEN, random
 *    misalignment, random keys.
 * usage: test [NRAND 16000] [MAXLEN 300000] [vectors test/512/vectors.txt] */
#include <stdio.h>
#include <stdlib.h>
#include "chainhash512.h"
static uint64_t rs = 0x5EED;
static uint64_t rnd(void) { uint64_t z = (rs += 0x9e3779b97f4a7c15ULL); z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL; z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; return z ^ (z >> 31); }
static chainhash512_key K[4], KW;
static chainhash512_stream S;
static const char *const NAME[4] = {"portable", "pclmul-sse4.1", "avx512-vpclmulqdq", "neon-pmull"};
int main(int argc, char **argv) {
    long nrand = argc > 1 ? atol(argv[1]) : 16000, it; size_t maxlen = argc > 2 ? (size_t)atol(argv[2]) : 300000, i, cap;
    const char *vec = argc > 3 ? argv[3] : "test/512/vectors.txt";
    uint8_t raw[CHAINHASH512_KEY_BYTES], *buf, *msg; uint64_t w[CHAINHASH512_KEY_WORDS];
    int det = chainhash512_backend(), bk[4], nb = 0, b, nv = 0, bad = 0, nL = 0, r, qq; long tot = 0, rbad = 0, sbad = 0, stot = 0;
    FILE *f; char line[512]; size_t Ls[64];
    bk[nb++] = CH512_PORTABLE;
    if (det == CH512_AVX512) bk[nb++] = CH512_PCLMUL;
    if (det != CH512_PORTABLE) bk[nb++] = det;
    printf("ChainHash-512; detected backend: %s\n", NAME[det]);
    /* 1. vectors */
    for (i = 0; i < CHAINHASH512_KEY_BYTES; i++) raw[i] = (uint8_t)((i * 73 + 11) & 255);
    for (b = 0; b < nb; b++) { chainhash512_key_from_bytes_with_backend(&K[b], raw, bk[b]); if (chainhash512_key_backend(&K[b]) != bk[b]) { printf("backend %d not set\n", bk[b]); bad++; } }
    for (i = 0; i < CHAINHASH512_KEY_WORDS; i++) { int j; w[i] = 0; for (j = 7; j >= 0; j--) w[i] = (w[i] << 8) | raw[8 * i + j]; }
    chainhash512_key_from_words(&KW, w);
    cap = 400000; msg = (uint8_t *)malloc(cap); for (i = 0; i < cap; i++) msg[i] = (uint8_t)((i * 137 + 29) & 255);
    if (!(f = fopen(vec, "r"))) { printf("missing %s\n", vec); return 2; }
    while (fgets(line, sizeof line, f)) {
        size_t n, o, step; char hex[200]; uint8_t want[64], got[64]; int j;
        if (line[0] == '#' || sscanf(line, "%zu %199s", &n, hex) != 2) continue;
        for (j = 0; j < 64; j++) { unsigned x; sscanf(hex + 2 * j, "%2x", &x); want[j] = (uint8_t)x; }
        chainhash512_reference(&K[0], msg, n, got); bad += memcmp(got, want, 64) != 0; nv++;
        chainhash512(&KW, msg, n, got); bad += memcmp(got, want, 64) != 0;
        for (b = 0; b < nb; b++) {
            chainhash512(&K[b], msg, n, got); bad += memcmp(got, want, 64) != 0;
            chainhash512_init(&S, &K[b]); o = 0; step = 777;
            while (o < n) { size_t m = n - o < step ? n - o : step; chainhash512_update(&S, msg + o, m); o += m; step = step * 3 % 20011 + 1; }
            chainhash512_final(&S, got); bad += memcmp(got, want, 64) != 0;
        }
    }
    fclose(f);
    printf("vectors: %d lengths x (reference + key_from_words + %d backends x {one-shot, stream}): %d mismatches -> %s\n", nv, nb, bad, bad || !nv ? "FAIL" : "PASS");
    if (!nv) bad++;
    /* 2. random, against the reference */
    for (i = 0; i < CHAINHASH512_KEY_BYTES; i++) raw[i] = (uint8_t)rnd();
    for (b = 0; b < nb; b++) chainhash512_key_from_bytes_with_backend(&K[b], raw, bk[b]);
    cap = maxlen > 3 * PH1_R * PH1_BLOCK ? maxlen : 3 * PH1_R * PH1_BLOCK; buf = (uint8_t *)malloc(cap + 64); for (i = 0; i < cap + 64; i++) buf[i] = (uint8_t)rnd();
    for (r = 0; r <= 2; r++) for (qq = 1; qq <= 8; qq++) { size_t base = (size_t)r * 8 * PH1_BLOCK + (size_t)(qq - 1) * PH1_BLOCK; Ls[nL++] = base + 1; Ls[nL++] = base + PH1_BLOCK; }
    for (it = -(long)nL - 4201; it < nrand; it++) {
        size_t n; const uint8_t *p; uint8_t ref[64], got[64];
        if (it < -(long)nL) n = (size_t)(it + nL + 4201); else if (it < 0) n = Ls[it + nL];
        else { n = (size_t)(rnd() % (maxlen + 1)); if (it % 3 == 0) n = (rnd() % 24) * PH1_BLOCK + (rnd() % 5) * PH1_CHUNK + (rnd() % 300) - 150; if (n > maxlen) n = maxlen; }
        if (it >= 0 && n > 60000 && (it % 8)) n %= 60000;                       /* keep the bit-serial reference affordable */
        p = buf + (rnd() % 33); chainhash512_reference(&K[0], p, n, ref);
        for (b = 1; b < nb; b++) { chainhash512(&K[b], p, n, got); tot++; if (memcmp(ref, got, 64)) { if (rbad < 5) printf("MISMATCH %s n=%zu\n", NAME[bk[b]], n); rbad++; } }
        if ((it & 7) == 0) for (b = 0; b < nb; b++) { size_t o = 0; chainhash512_init(&S, &K[b]);
            while (o < n) { size_t m = (rnd() & 3) == 0 ? (size_t)(rnd() % (3 * PH1_BLOCK)) : (size_t)(rnd() % 700); if (m > n - o) m = n - o; chainhash512_update(&S, p + o, m); o += m; }
            chainhash512_final(&S, got); stot++; if (memcmp(ref, got, 64)) { if (sbad < 5) printf("STREAM MISMATCH %s n=%zu\n", NAME[bk[b]], n); sbad++; } }
        if (it >= 0 && it % 1000 == 0) { for (i = 0; i < CHAINHASH512_KEY_BYTES; i++) raw[i] = (uint8_t)rnd(); for (b = 0; b < nb; b++) chainhash512_key_from_bytes_with_backend(&K[b], raw, bk[b]); }
    }
    chainhash512_key_from_seed(&K[0], 42);
    printf("random: backends vs reference %ld hashes, %ld mismatches; streaming vs reference %ld, %ld mismatches -> %s\n", tot, rbad, stot, sbad, (rbad || sbad) ? "FAIL" : "PASS");
    printf("ChainHash-512 ALL %s\n", (bad || rbad || sbad) ? "FAIL" : "PASS");
    free(buf); free(msg);
    return (bad || rbad || sbad) != 0;
}
