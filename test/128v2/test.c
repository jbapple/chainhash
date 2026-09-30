/* ChainHash-128 v2 tests.
 * 1. test/128v2/vectors.txt (independent Python oracle, key from bytes) on the reference
 *    chainhash128v2_reference and on every available backend, one-shot and streaming.
 * 2. Random inputs against the reference: one-shot and random-split streaming on every backend,
 *    random misalignment, fresh keys every 50 inputs, including all-0x00 and all-0xff keys.
 * 3. chainhash128v2_key_from_words gives the key of the same bytes.
 * usage: test [random inputs (20000)] [vectors file (test/128v2/vectors.txt)] */
#include <stdio.h>
#include <stdlib.h>
#include "chainhash128v2.h"
static uint64_t rs = 0x9e3779b97f4a7c15ULL;
static uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static chainhash128v2_key K, K2;
static chainhash128v2_stream S;
int main(int argc, char **argv) {
    int N = argc > 1 ? atoi(argv[1]) : 20000, backs[5], nb = 0, bi, bad = 0, nv = 0, t;
    size_t cap = 3 * CH128P_REGION + 1024 > 40064 ? 3 * CH128P_REGION + 1024 : 40064, i;
    uint8_t *buf = (uint8_t *)malloc(cap + 64), key[CHAINHASH128V2_KEY_BYTES], *m; FILE *f; char line[256];
    const char *fn = argc > 2 ? argv[2] : "test/128v2/vectors.txt";
    ch128_word w[CHAINHASH128V2_KEY_WORDS];
    for (bi = 0; bi < 5; bi++) if (chainhash128v2_has_backend(bi)) backs[nb++] = bi;
    printf("ChainHash-128 v2; dispatched backend %d; available:", chainhash128v2_backend());
    for (bi = 0; bi < nb; bi++) printf(" %d", backs[bi]);
    printf("\n");
    for (i = 0; i < CHAINHASH128V2_KEY_BYTES; i++) key[i] = (uint8_t)(i * 73 + 11);
    chainhash128v2_key_from_bytes(&K, key);
    for (i = 0; i < CHAINHASH128V2_KEY_WORDS; i++) w[i] = ch128_load(key + 16 * i);
    chainhash128v2_key_from_words(&K2, w);
    for (i = 0; i < cap; i++) buf[i] = (uint8_t)(i * 137 + 29);
    if (!(f = fopen(fn, "r"))) { printf("missing %s\n", fn); return 1; }
    while (fgets(line, sizeof line, f)) {
        size_t L; char hex[40]; uint8_t want[16], got[16]; int j;
        if (line[0] == '#' || sscanf(line, "%zu %32s", &L, hex) != 2) continue;
        for (j = 0; j < 16; j++) { unsigned x; sscanf(hex + 2 * j, "%2x", &x); want[j] = (uint8_t)x; }
        nv++;
        chainhash128_store(got, chainhash128v2_reference(&K, buf, L));
        if (memcmp(got, want, 16)) { bad++; printf("VECTOR MISMATCH reference len=%zu\n", L); }
        chainhash128_store(got, chainhash128v2(&K2, buf, L));
        if (memcmp(got, want, 16)) { bad++; printf("VECTOR MISMATCH key_from_words len=%zu\n", L); }
        for (bi = 0; bi < nb; bi++) {
            chainhash128_store(got, chainhash128v2_with_backend(&K, buf, L, backs[bi]));
            if (memcmp(got, want, 16)) { bad++; printf("VECTOR MISMATCH backend %d len=%zu\n", backs[bi], L); }
            chainhash128v2_init(&S, &K, backs[bi]); chainhash128v2_update(&S, buf, L / 3);
            chainhash128v2_update(&S, buf + L / 3, L - L / 3); chainhash128_store(got, chainhash128v2_final(&S));
            if (memcmp(got, want, 16)) { bad++; printf("VECTOR MISMATCH stream backend %d len=%zu\n", backs[bi], L); }
        }
    }
    fclose(f);
    printf("vectors (%s): %d lengths x (reference + key_from_words + %d backends x {one-shot, stream}): %d mismatches\n", fn, nv, nb, bad);
    if (!nv) bad++;
    for (t = 0; t < N; t++) {
        size_t L = (rnd() & 1) ? (size_t)(rnd() % 1200) : (size_t)(rnd() % cap), off = (size_t)(rnd() % 16), j; ch128_word r;
        if (t % 50 == 0) {
            for (j = 0; j < CHAINHASH128V2_KEY_BYTES; j++) key[j] = (uint8_t)rnd();
            if (t % 500 == 0) memset(key, (t / 500) & 1 ? 0xff : 0, sizeof key);
            chainhash128v2_key_from_bytes(&K, key);
        }
        m = buf + off; for (j = 0; j < L; j++) m[j] = (uint8_t)rnd();
        r = chainhash128v2_reference(&K, m, L);
        for (bi = 0; bi < nb; bi++) { size_t pos = 0;
            if (!ch128_equal(chainhash128v2_with_backend(&K, m, L, backs[bi]), r)) { if (bad++ < 10) printf("MISMATCH backend %d len=%zu\n", backs[bi], L); }
            chainhash128v2_init(&S, &K, backs[bi]);
            while (pos < L) { size_t c = (size_t)(rnd() % (L - pos + 1)); if (rnd() & 1) c = c % 64; chainhash128v2_update(&S, m + pos, c); pos += c; }
            if (!ch128_equal(chainhash128v2_final(&S), r)) { if (bad++ < 10) printf("MISMATCH stream backend %d len=%zu\n", backs[bi], L); }
        }
    }
    chainhash128v2_key_from_seed(&K, 42);
    if (!ch128_equal(chainhash128v2(&K, buf, 5000), chainhash128v2_reference(&K, buf, 5000))) bad++;
    printf("random: %d inputs (incl. all-0x00/0xff keys) x %d backends x {one-shot, stream}\n%s (%d mismatches)\n", N, nb, bad ? "FAIL" : "PASS", bad);
    free(buf);
    return bad != 0;
}
