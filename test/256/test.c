/* ChainHash-256 tests: every backend available on this machine, one-shot and streaming (8 split modes),
 * against the definition chainhash256_reference (bit-serial; with -DCH256_HW_CLMUL its 64x64 carry-less
 * product uses the hardware instruction, same arithmetic, for speed).
 *   fast backends: every length 0..LALL (misaligned by n%13), NRAND random inputs (fresh key every 400,
 *   0..63-byte misalignment, lengths up to 4 regions), NSTREAM streaming runs, large inputs;
 *   portable: 0..300 and 100 streaming runs.
 *   vectors: every line "64 key n digest" of the vectors file (message m[j] = j*137+29): key = a decimal seed
 *   (chainhash256_key_from_seed) or R = the raw 288-byte key b[i] = (73 i + 11) mod 256, given both as bytes
 *   (chainhash256_key_from_bytes) and as words (chainhash256_key_from_words).
 * usage: test [NRAND 20000] [NSTREAM 4000] [LALL 41024] [vectors test/256/vectors.txt] */
#include <stdio.h>
#include <stdlib.h>
#include "chainhash256.h"
static uint64_t rs = 7;
static uint64_t splitmix(uint64_t *st) {
    uint64_t z = (*st += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL; return z ^ (z >> 31);
}
static uint64_t rnd(void) { return splitmix(&rs); }
static chainhash256_key K, F;
static chainhash256_stream S;
static void stream_hash(const chainhash256_key *k, const uint8_t *m, size_t n, uint8_t o[32], int mode) {
    size_t i = 0; chainhash256_init(&S, k);
    while (i < n) { size_t c;
        switch (mode) { case 0: c = 1 + rnd() % 7; break; case 1: c = PH_REGION; break; case 2: c = PH_REGION - 1; break;
                        case 3: c = PH_REGION + 1; break; case 4: c = PH_BLOCK; break; default: c = 1 + rnd() % (rnd() % 2 ? 64 : 3 * PH_REGION); }
        if (c > n - i) c = n - i;
        chainhash256_update(&S, m + i, c); i += c;
        if (mode >= 4 && rnd() % 7 == 0) chainhash256_update(&S, m + i, 0); }
    chainhash256_final(&S, o);
}
static const char *const NAME[4] = {"portable", "pclmul-sse4.1", "avx512-vpclmulqdq", "neon-pmull"};
int main(int argc, char **argv) {
    int NRAND = argc > 1 ? atoi(argv[1]) : 20000, NSTREAM = argc > 2 ? atoi(argv[2]) : 4000;
    size_t LMAX = argc > 3 ? (size_t)atol(argv[3]) : (size_t)(PH_REGION + 2 * PH_BLOCK + 64);
    const char *vec = argc > 4 ? argv[4] : "test/256/vectors.txt";
    size_t maxn = 4 * PH_REGION + PH_BLOCK + 200, i; if (maxn < 262145 + 64) maxn = 262145 + 64;
    uint8_t *buf = (uint8_t *)malloc(maxn + 64), *mm = (uint8_t *)malloc(maxn);
    int best = chainhash256_backend(), list[3], nb = 0, fails = 0, bi;
    for (i = 0; i < maxn + 64; i++) buf[i] = (uint8_t)rnd();
    for (i = 0; i < maxn; i++) mm[i] = (uint8_t)(i * 137 + 29);
    list[nb++] = CH256_PORTABLE;
    if (best == CH256_AVX512) { list[nb++] = CH256_PCLMUL; list[nb++] = CH256_AVX512; } else if (best != CH256_PORTABLE) list[nb++] = best;
    for (bi = 0; bi < nb; bi++) {
        int be = list[bi], fast = be != CH256_PORTABLE, f = 0, c = 0, sf = 0, sc = 0, t; uint8_t r[32], o[32];
        size_t LALL = fast ? LMAX : 300, n;
        if (!chainhash256_has_backend(be)) { printf("  has_backend(%d) = 0\n", be); fails++; continue; }
        chainhash256_key_from_seed_with_backend(&K, 777, be);
        if (chainhash256_key_backend(&K) != be) { printf("  key backend %d != %d\n", chainhash256_key_backend(&K), be); fails++; }
        for (n = 0; n <= LALL; n++) { chainhash256_reference(&K, buf + (n % 13), n, r); chainhash256(&K, buf + (n % 13), n, o); c++;
            if (memcmp(r, o, 32)) { f++; if (f < 5) printf("  one-shot mismatch %s n=%zu\n", NAME[be], n); } }
        for (t = 0; t < (fast ? NRAND : 0); t++) { size_t off;
            if (t % 400 == 0) chainhash256_key_from_seed_with_backend(&K, rnd(), be);
            n = (t % 4 == 0) ? (size_t)(rnd() % (4 * PH_REGION)) : (size_t)(rnd() % (2 * PH_REGION)); off = rnd() % 64;
            chainhash256_reference(&K, buf + off, n, r); chainhash256(&K, buf + off, n, o); c++;
            if (memcmp(r, o, 32)) { f++; if (f < 5) printf("  one-shot mismatch %s n=%zu off=%zu\n", NAME[be], n, off); } }
        chainhash256_key_from_seed_with_backend(&K, 4242, be);
        for (t = 0; t < (fast ? NSTREAM : 100); t++) { size_t off; int mode = t % 8;
            n = (t % 5 == 0) ? (size_t)(rnd() % (4 * PH_REGION)) : (size_t)(rnd() % (2 * PH_REGION + PH_BLOCK)); off = rnd() % 64;
            if (!fast) n %= 2 * PH_REGION + 100;
            chainhash256_reference(&K, buf + off, n, r); stream_hash(&K, buf + off, n, o, mode); sc++;
            if (memcmp(r, o, 32)) { sf++; if (sf < 5) printf("  stream mismatch %s n=%zu mode=%d\n", NAME[be], n, mode); } }
        { size_t big[] = {PH_REGION, PH_REGION + 1, 2 * PH_REGION, 3 * PH_REGION - 1, 4 * PH_REGION, 4 * PH_REGION + PH_BLOCK + 100}; int j;
          for (j = 0; j < (fast ? 6 : 2); j++) { chainhash256_reference(&K, buf + 1, big[j], r); chainhash256(&K, buf + 1, big[j], o); c++; if (memcmp(r, o, 32)) f++;
              stream_hash(&K, buf + 1, big[j], o, 7); sc++; if (memcmp(r, o, 32)) sf++; } }
        printf("backend %-18s one-shot == reference: %s (%d)   streaming == reference: %s (%d)\n", NAME[be], f ? "FAIL" : "PASS", c, sf ? "FAIL" : "PASS", sc);
        fails += f + sf;
    }
    { FILE *fp = fopen(vec, "r"); int nv = 0, bad = 0, m; char key[24], hex[80], h2[80]; size_t n;
      static uint8_t kb[CHAINHASH256_KEY_BYTES]; static uint64_t kw[CHAINHASH256_KEY_WORDS]; static chainhash256_key B, W;
      int e;
      for (e = 0; e < CHAINHASH256_KEY_BYTES; e++) kb[e] = (uint8_t)(73 * e + 11);
      for (e = 0; e < CHAINHASH256_KEY_WORDS; e++) { int j; kw[e] = 0; for (j = 0; j < 8; j++) kw[e] |= (uint64_t)kb[8 * e + j] << (8 * j); }
      chainhash256_key_from_bytes(&B, kb); chainhash256_key_from_words(&W, kw);
      if (!fp) { printf("  cannot open %s\n", vec); fails++; }
      else { while (fscanf(fp, "%d %20s %zu %70s", &m, key, &n, hex) == 4) { uint8_t o[32]; int j;
                 if (m != PH_M) continue;
                 if (strcmp(key, "R")) { chainhash256_key_from_seed(&F, strtoull(key, NULL, 10)); chainhash256(&F, mm, n, o); nv++;
                     for (j = 0; j < 32; j++) snprintf(h2 + 2 * j, 3, "%02x", o[j]);
                     if (strcmp(h2, hex)) { bad++; if (bad < 4) printf("  vector mismatch seed=%s n=%zu\n", key, n); } }
                 else { nv++;
                     chainhash256(&B, mm, n, o); for (j = 0; j < 32; j++) snprintf(h2 + 2 * j, 3, "%02x", o[j]);
                     if (strcmp(h2, hex)) { bad++; printf("  vector mismatch (raw key from bytes) n=%zu\n", n); }
                     chainhash256(&W, mm, n, o); for (j = 0; j < 32; j++) snprintf(h2 + 2 * j, 3, "%02x", o[j]);
                     if (strcmp(h2, hex)) { bad++; printf("  vector mismatch (raw key from words) n=%zu\n", n); } } }
             fclose(fp);
             printf("vectors (%s backend, %s): %d lines, %d mismatches (seed keys; raw key R from bytes and from words)\n", NAME[chainhash256_backend()], vec, nv, bad);
             if (bad || !nv) fails++; } }
    printf("ChainHash-256 ALL %s\n", fails ? "FAIL" : "PASS");
    free(buf); free(mm);
    return fails != 0;
}
