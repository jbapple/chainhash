/* ph_ref.h -- the ChainHash-256 specification as C (bit-serial reference), at any PH_M (16/32/64 pairs per
 * block). Extracted verbatim from include/chainhash256.h, where PH_M is fixed at 64; used by the certificates
 * (cert_ph.c, cert_2l.c, cert_pk.c), which check it at several block sizes. */
#ifndef PH_V2
#define PH_V2 1
#endif
/* The specification as C (bit-serial reference), docs/SPEC-256.md.
 *   field   L = GF(2)[x]/(f), f = x^256 + x^10 + x^5 + x^2 + 1 (irreducible: test/256/ph_field.py);
 *           element = 4 little-endian 64-bit limbs (limb 0 = coefficients x^0..x^63)
 *   block   PH_M pairs (x_i, y_i) of 32-byte elements; c = sum_i (x_i + k_i)(y_i + l_i) with the power-key masks
 *           k_i = s^(2i-1), l_i = s^(2i) (ph2_derive from the 288-byte raw key ph2_raw)
 *   layout  regions of 8 interleaved blocks and zero-padded tail blocks, as in docs/SPEC-256.md section 3
 *   outer   ph2l_hash_ref below (the two-level stage, normative); ph_hash_ref is the earlier outer stage (a Horner
 *           chain in Y over the block values), kept for comparisons and not used by the hash
 * The multiply here is schoolbook over limbs with ch_clmul64 (bit-serial unless CH256_HW_CLMUL) and a bit-level
 * reduction; it never uses the Karatsuba evaluation the kernels use. */
#ifndef PH_REF_H
#define PH_REF_H
#include <stdint.h>
#include <string.h>
/* 64x64 -> 128 carry-less product: bit-serial unless CH256_HW_CLMUL (then PMULL / PCLMULQDQ); same result */
#if defined(CH256_HW_CLMUL) && defined(__aarch64__) && (defined(__ARM_FEATURE_AES) || defined(__ARM_FEATURE_CRYPTO))
#include <arm_neon.h>
static inline void ch_clmul64(uint64_t a, uint64_t b, uint64_t *lo, uint64_t *hi){
    uint64x2_t v = vreinterpretq_u64_p128(vmull_p64((poly64_t)a,(poly64_t)b));
    *lo=vgetq_lane_u64(v,0); *hi=vgetq_lane_u64(v,1);
}
#elif defined(CH256_HW_CLMUL) && defined(__PCLMUL__)
#include <immintrin.h>
static inline void ch_clmul64(uint64_t a, uint64_t b, uint64_t *lo, uint64_t *hi){
    __m128i r=_mm_clmulepi64_si128(_mm_set_epi64x(0,a),_mm_set_epi64x(0,b),0x00);
    *lo=(uint64_t)_mm_cvtsi128_si64(r); *hi=(uint64_t)_mm_extract_epi64(r,1);
}
#else
static inline void ch_clmul64(uint64_t a, uint64_t b, uint64_t *lo, uint64_t *hi){
    uint64_t l=0,h=0;
    for(int i=0;i<64;i++){
        if((b>>i)&1u){
            if(i==0){ l^=a; }
            else { l^=a<<i; h^=a>>(64-i); }
        }
    }
    *lo=l; *hi=h;
}
#endif
static inline uint64_t ch_splitmix(uint64_t *st){
    uint64_t z=(*st += 0x9E3779B97F4A7C15ULL);
    z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL;
    z=(z^(z>>27))*0x94D049BB133111EBULL;
    return z^(z>>31);
}
static inline uint64_t ch_ld64(const uint8_t *p){ uint64_t v; memcpy(&v,p,8); return v; }
#include <stdlib.h>
#ifndef PH_M
#define PH_M 16
#endif
#if PH_M % 8
#error "PH_M must be a multiple of 8"
#endif
#define PH_BLOCK (64*PH_M)
#define PH_REGION (8*PH_BLOCK)
typedef struct { uint64_t w[4]; } ph_el;
static inline ph_el ph_xor(ph_el a,ph_el b){ ph_el r; for(int i=0;i<4;i++) r.w[i]=a.w[i]^b.w[i]; return r; }
/* reduce a 512-bit product r[0..7] mod f, bit by bit from the top */
static inline ph_el ph_reduce512(const uint64_t r0[8]){
    uint64_t r[8]; memcpy(r,r0,sizeof r);
    for(int i=511;i>=256;i--) if((r[i>>6]>>(i&63))&1){ r[i>>6]^=1ULL<<(i&63);
        int t[4]={i-256+10,i-256+5,i-256+2,i-256}; for(int j=0;j<4;j++) r[t[j]>>6]^=1ULL<<(t[j]&63); }
    ph_el e; memcpy(e.w,r,32); return e; }
static inline ph_el ph_mul(ph_el a,ph_el b){
    uint64_t r[8]={0};
    for(int i=0;i<4;i++) for(int j=0;j<4;j++){ uint64_t lo,hi; ch_clmul64(a.w[i],b.w[j],&lo,&hi); r[i+j]^=lo; r[i+j+1]^=hi; }
    return ph_reduce512(r); }
static inline ph_el ph_addint(ph_el v,ph_el t){ ph_el r; unsigned __int128 c=0;
    for(int i=0;i<4;i++){ unsigned __int128 s=(unsigned __int128)v.w[i]+t.w[i]+c; r.w[i]=(uint64_t)s; c=s>>64; } return r; }
static inline ph_el ph_pow(ph_el b,uint64_t e){ ph_el r={{1,0,0,0}}; while(e){ if(e&1) r=ph_mul(r,b); b=ph_mul(b,b); e>>=1; } return r; }
typedef struct { ph_el k[PH_M], l[PH_M]; ph_el Y, c[5], tau; ph_el z; } ph_key;   /* z: two-level (2L) region-chain key */
static inline void ph_fill(ph_el*e,uint64_t*st){ for(int i=0;i<4;i++) e->w[i]=ch_splitmix(st); }
/* ---- PH-256 v2: POWER KEY.  Raw key = s | y | z | t | c0..c4 (9 elements, 288 bytes, each 32 bytes little-endian,
 * limb 0 first); the mask table is derived: k_j = s^(2j+1), l_j = s^(2j+2) for pair j = 0..PH_M-1 (pairs 1..m of the
 * spec: x_i masked by s^(2i-1), y_i by s^(2i)).  Y = y, tau = t.  All 9 elements independent and uniform. ---- */
typedef struct { ph_el s, y, z, t, c[5]; } ph2_raw;
static inline void ph2_derive_with(const ph2_raw*R,ph_key*K,ph_el (*mul)(ph_el,ph_el)){
    ph_el p=R->s;                                                          /* s^1 */
    for(int j=0;j<PH_M;j++){ K->k[j]=p; p=mul(p,R->s); K->l[j]=p; p=mul(p,R->s); }
    K->Y=R->y; K->z=R->z; K->tau=R->t; for(int i=0;i<5;i++) K->c[i]=R->c[i]; }
static inline ph_el ph_mul(ph_el a,ph_el b);
static inline void ph2_derive(const ph2_raw*R,ph_key*K){ ph2_derive_with(R,K,ph_mul); }
static inline void ph2_raw_from_bytes(const uint8_t b[288],ph2_raw*R){ ph_el*e=&R->s;
    for(int i=0;i<9;i++) for(int l=0;l<4;l++) e[i].w[l]=ch_ld64(b+32*i+8*l); }
static inline void ph2_raw_from_seed(ph2_raw*R,uint64_t seed){          /* tests/vectors only */
    uint64_t st=seed^0x5048323536763200ULL;                                /* "PH256v2" */
    ph_el*e=&R->s; for(int i=0;i<9;i++) ph_fill(&e[i],&st);
    if(!(R->y.w[0]|R->y.w[1]|R->y.w[2]|R->y.w[3])) R->y.w[0]=1; }
static inline void ph_key_from_seed(ph_key*K,uint64_t seed){        /* tests/benchmarks only */
#ifdef PH_V2
    ph2_raw R; ph2_raw_from_seed(&R,seed); ph2_derive(&R,K);
#else
        /* tests/benchmarks only */
    uint64_t st=seed^0x5048323536000000ULL^(uint64_t)PH_M;            /* "PH256" | m */
    for(int i=0;i<PH_M;i++){ ph_fill(&K->k[i],&st); ph_fill(&K->l[i],&st); }
    ph_fill(&K->Y,&st); for(int j=0;j<5;j++) ph_fill(&K->c[j],&st); ph_fill(&K->tau,&st);
    if(!(K->Y.w[0]|K->Y.w[1]|K->Y.w[2]|K->Y.w[3])) K->Y.w[0]=1;
    ph_fill(&K->z,&st);
#endif
    }

static inline ph_el ph_block(const ph_el*x,const ph_el*y,const ph_key*K){
    ph_el c={{0,0,0,0}}; for(int i=0;i<PH_M;i++) c=ph_xor(c,ph_mul(ph_xor(x[i],K->k[i]),ph_xor(y[i],K->l[i]))); return c; }
static inline ph_el ph_ldr(const uint8_t*p,size_t stride){ ph_el e; for(int l=0;l<4;l++) e.w[l]=ch_ld64(p+stride*l); return e; }
/* block b of region R / tail block blk (zero padded): the words of pair i */
static inline void ph_region_block(const uint8_t*R,int b,ph_el*x,ph_el*y){
    for(int i=0;i<PH_M;i++){ x[i]=ph_ldr(R+512*i+8*b,64); y[i]=ph_ldr(R+512*i+256+8*b,64); } }
static inline void ph_tail_block(const uint8_t*blk,ph_el*x,ph_el*y){
    for(int i=0;i<PH_M;i++){ x[i]=ph_ldr(blk+512*(i/8)+8*(i%8),64); y[i]=ph_ldr(blk+512*(i/8)+256+8*(i%8),64); } }
static inline size_t ph_nblocks(size_t n){ size_t rem=n%PH_REGION; return (n/PH_REGION)*8+(rem+PH_BLOCK-1)/PH_BLOCK; }
static inline ph_el ph_finish(ph_el v,const ph_key*K){
    ph_el X=ph_addint(v,K->tau), G=ph_mul(X,X);
    ph_el t=ph_mul(ph_xor(G,K->c[0]),ph_xor(ph_xor(X,G),K->c[1]));
    return ph_xor(ph_mul(ph_xor(X,K->c[2]),ph_xor(t,K->c[3])),K->c[4]); }
static inline void ph_store(ph_el d,uint8_t out[32]){ for(int i=0;i<4;i++) for(int b=0;b<8;b++) out[8*i+b]=(uint8_t)(d.w[i]>>(8*b)); }
static inline void ph_hash_ref(const ph_key*K,const uint8_t*m,size_t n,uint8_t out[32]){
    ph_el v={{(uint64_t)n,0,0,0}}, x[PH_M], y[PH_M]; size_t nfull=n/PH_REGION, rem=n-nfull*PH_REGION;
    for(size_t r=0;r<nfull;r++) for(int b=0;b<8;b++){ ph_region_block(m+r*PH_REGION,b,x,y); v=ph_xor(ph_mul(v,K->Y),ph_block(x,y,K)); }
    const uint8_t*T=m+nfull*PH_REGION;
    for(size_t b=0;b*PH_BLOCK<rem;b++){ uint8_t blk[PH_BLOCK]; size_t av=rem-b*PH_BLOCK; if(av>PH_BLOCK) av=PH_BLOCK;
        memcpy(blk,T+b*PH_BLOCK,av); memset(blk+av,0,PH_BLOCK-av); ph_tail_block(blk,x,y); v=ph_xor(ph_mul(v,K->Y),ph_block(x,y,K)); }
    ph_store(ph_finish(v,K),out); }
/* ======================= PH-256 v1: two-level (2L) outer stage — the NORMATIVE definition =======================
 * Block values b_1..b_m (level 1 as above) in message order: blocks 0..7 of each full region, then the tail blocks;
 * m = ph_nblocks(n), except n = 0: m = 1 and b_1 = 0.  Regions of R = 8 values; m' = ceil(m/8).
 * Region value (q values a_1..a_q, h = ceil(q/2), f = floor(q/2); y = K->Y):
 *     c = sum_{i=1}^{f} (a_i + y^(2i-1)) (a_{h+i} + y^(2i))  +  [q odd] a_h
 * Outer: V = n z^{m'} + sum_rho c_rho z^{m'-rho} (Horner in the independent key z, length leading);
 *     then X = V +_Z tau and the unchanged quintic (ph_finish).  For n <= one block: V = n z + b_1. */
static inline ph_el ph2l_region(const ph_el*a,int q,const ph_el*Y){
    int h=(q+1)/2, f=q/2; ph_el c={{0,0,0,0}}, yp[9]; yp[0]=(ph_el){{1,0,0,0}}; for(int e=1;e<=8;e++) yp[e]=ph_mul(yp[e-1],*Y);
    for(int i=1;i<=f;i++) c=ph_xor(c,ph_mul(ph_xor(a[i-1],yp[2*i-1]),ph_xor(a[h+i-1],yp[2*i])));
    if(q&1) c=ph_xor(c,a[h-1]);
    return c; }
/* all block values of the message (b[0..m-1]); returns m */
static inline size_t ph2l_blocks(const ph_key*K,const uint8_t*m,size_t n,ph_el*b){
    if(!n){ b[0]=(ph_el){{0,0,0,0}}; return 1; }
    ph_el x[PH_M], y[PH_M]; size_t nfull=n/PH_REGION, rem=n-nfull*PH_REGION, t=0;
    for(size_t r=0;r<nfull;r++) for(int bb=0;bb<8;bb++){ ph_region_block(m+r*PH_REGION,bb,x,y); b[t++]=ph_block(x,y,K); }
    const uint8_t*T=m+nfull*PH_REGION;
    for(size_t bb=0;bb*PH_BLOCK<rem;bb++){ uint8_t blk[PH_BLOCK]; size_t av=rem-bb*PH_BLOCK; if(av>PH_BLOCK) av=PH_BLOCK;
        memcpy(blk,T+bb*PH_BLOCK,av); memset(blk+av,0,PH_BLOCK-av); ph_tail_block(blk,x,y); b[t++]=ph_block(x,y,K); }
    return t; }
/* V from the block values (the outer stage by itself: used by the 2L certificate) */
static inline ph_el ph2l_outer(const ph_el*b,size_t m,uint64_t n,const ph_el*Y,const ph_el*z){
    size_t mp=(m+7)/8; ph_el V={{n,0,0,0}};
    for(size_t r=0;r<mp;r++){ int q=(int)((m-8*r)<8? m-8*r : 8); V=ph_xor(ph_mul(V,*z),ph2l_region(b+8*r,q,Y)); }
    return V; }
static inline void ph2l_hash_ref(const ph_key*K,const uint8_t*m,size_t n,uint8_t out[32]){
    size_t mmax=ph_nblocks(n)+1; ph_el*b=(ph_el*)malloc(mmax*sizeof(ph_el));
    size_t mb=ph2l_blocks(K,m,n,b); ph_el V=ph2l_outer(b,mb,(uint64_t)n,&K->Y,&K->z); free(b);
    ph_store(ph_finish(V,K),out); }
#endif
