/* xcheck.c -- randomized cross-check of ChainHash-256 against the normative reference ph2l_hash_ref:
 * many random RAW keys (tau / c_i limbs forced to 0, 1, ~0 at random: long carry chains in the twist), every runtime
 * backend, lengths 0..320 and around the block/region edges, random misalignment, one-shot and streaming (random
 * split points).  Uses the header's internal names (test only); -DCH256_HW_CLMUL speeds up the reference. */
#include <stdio.h>
#include <stdlib.h>
#include "chainhash256.h"
#ifndef NKEYS
#define NKEYS 600
#endif
static uint64_t rs=0x1234567; static uint64_t rnd(void){ return ch_splitmix(&rs); }
static uint64_t edge(void){ uint64_t v=rnd(); switch(rnd()%8){ case 0: return 0; case 1: return ~0ULL; case 2: return 1; case 3: return ~0ULL-1; case 4: return 1ULL<<63; default: return v; } }
static void rraw(ph2_raw*R){ ph_el*e=&R->s; for(int i=0;i<9;i++) for(int l=0;l<4;l++) e[i].w[l]=rnd();
    for(int l=0;l<4;l++){ R->t.w[l]=edge(); for(int j=0;j<5;j++) R->c[j].w[l]=(rnd()%4)? rnd() : edge(); }
    if(!(R->y.w[0]|R->y.w[1]|R->y.w[2]|R->y.w[3])) R->y.w[0]=1; }
int main(void){
    size_t maxn=2*PH_REGION+PH_BLOCK+300; uint8_t*buf=malloc(maxn+64); for(size_t i=0;i<maxn+64;i++) buf[i]=(uint8_t)rnd();
    int best=ph256_backend(), list[3], nb=0; list[nb++]=PH256_PORTABLE;
    if(best==PH256_AVX512){ list[nb++]=PH256_PCLMUL; list[nb++]=PH256_AVX512; } else if(best!=PH256_PORTABLE) list[nb++]=best;
    static ph256_key H; static ph_key R; static ph2_raw W; long fails=0, n_one=0, n_str=0;
    size_t edges[]={PH_BLOCK-1,PH_BLOCK,PH_BLOCK+1,2*PH_BLOCK-1,2*PH_BLOCK+64,PH_REGION-1,PH_REGION,PH_REGION+1,PH_REGION+63,PH_REGION+64,PH_REGION+65,PH_REGION+PH_BLOCK+7,2*PH_REGION+300};
    for(int bi=1;bi<nb;bi++){ int be=list[bi];
        for(int kk=0;kk<NKEYS;kk++){ rraw(&W); ph2_derive(&W,&R); ph256_init_raw2(&H,&W,be);
            for(int s=0;s<24;s++){ size_t n = s<16? (size_t)(rnd()%321) : s<20? (size_t)(64+rnd()%4100) : edges[rnd()%(sizeof edges/sizeof*edges)];
                if(kk%8==0 && s==0) n=(size_t)(kk/8)%321;                     /* every short length with several keys */
                size_t off=rnd()%64; uint8_t r[32],o[32]; ph2l_hash_ref(&R,buf+off,n,r);
                ph256(&H,buf+off,n,o); n_one++; if(memcmp(r,o,32)){ if(fails++<8) printf("  one-shot mismatch %s n=%zu\n",PH256_BACKEND_NAME[be],n); }
                ph256_stream S; ph256_stream_init(&S,&H); size_t i=0; while(i<n){ size_t c= rnd()%3? 1+rnd()%97 : 1+rnd()%(PH_REGION+5); if(c>n-i) c=n-i; ph256_update(&S,buf+off+i,c); i+=c; }
                ph256_final(&S,o); n_str++; if(memcmp(r,o,32)){ if(fails++<8) printf("  stream mismatch %s n=%zu\n",PH256_BACKEND_NAME[be],n); } } }
        printf("[xcheck] backend %-18s random raw keys %d: %s\n",PH256_BACKEND_NAME[be],NKEYS,fails?"FAIL":"PASS"); }
    printf("[xcheck] %ld one-shot + %ld streaming vs ph2l_hash_ref: %s\n",n_one,n_str,fails?"FAIL":"PASS"); return fails!=0; }
