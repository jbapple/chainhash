/* ChainHash-128: every reachable evaluation schedule gives the reference digest.
 * For each random input (random key, length, misalignment and contents) the
 * independent bit-serial evaluator chainhash128_portable is compared with
 * chainhash128_with_prefetch for every backend (and an unavailable id) x prefetch
 * hint x hint step, each with a random distance, and with streaming for every
 * backend x product method x stride 1..8 x eager/lazy, each with random splits.
 * Usage: schedule_knobs [cases=20000] [seed=1]. */
#include <stdio.h>
#include <stdlib.h>
#include "chainhash128.h"
static uint64_t rs;
static uint64_t rnd(void) { uint64_t z=(rs+=UINT64_C(0x9e3779b97f4a7c15)); z=(z^(z>>30))*UINT64_C(0xbf58476d1ce4e5b9); z=(z^(z>>27))*UINT64_C(0x94d049bb133111eb); return z^(z>>31); }
static size_t pick_len(void) {
    unsigned r=(unsigned)(rnd()%100); size_t base, d[5]={0,1,16,17,4095};
    if(r<30) return (size_t)(rnd()%600);
    if(r<70) return (size_t)(rnd()%40000);
    if(r<95) { base=4096*(size_t)(1+rnd()%12); return rnd()&1 ? base+d[rnd()%5] : base-d[rnd()%5]; }
    return (size_t)(rnd()%(300*1024));
}
static uint32_t pick_dist(void) { uint32_t d[8]={0,64,512,4096,8192,16384,65536,1u<<20}; return rnd()&1 ? d[rnd()%8] : (uint32_t)(rnd()%(1u<<18)); }
int main(int argc,char **argv) {
    long cases=argc>1?atol(argv[1]):20000, c, checks=0, bad=0;
    int backs[8], nb=0, b, hint, step;
    uint8_t *mem=(uint8_t *)malloc(300*1024+128);
    rs=argc>2?(uint64_t)atoll(argv[2]):1;
    for(b=0;b<=4;b++) if(chainhash128_has_backend(b)) backs[nb++]=b;
    backs[nb++]=77;
    printf("backends:"); for(b=0;b<nb;b++) printf(" %d",backs[b]); printf("\n");
    for(c=0;c<cases;c++) {
        size_t len=pick_len(), off=(size_t)(rnd()%64), i; uint8_t *p=mem+off; chainhash128_key k; ch128_word ref;
        uint8_t kb[128]; unsigned kind=(unsigned)(rnd()%100);
        for(i=0;i<128;i++) kb[i]=(uint8_t)rnd();
        if(kind==0) memset(kb+16,0,16); else if(kind==1) { memset(kb+16,0,16); kb[16]=1; }   /* y = 0, y = 1 */
        k=chainhash128_key_from_bytes(kb);
        for(i=0;i<len;i++) p[i]=(uint8_t)rnd();
        if(rnd()%8==0) memset(p,(int)(rnd()&0xff),len);
        ref=chainhash128_portable(&k,p,len);
        if(!ch128_equal(chainhash128(&k,p,len),ref)) { printf("FAIL chainhash128 len=%zu\n",len); bad++; }
        for(b=0;b<nb;b++) for(hint=0;hint<=4;hint++) for(step=64;step<=128;step+=64) {
            uint32_t d=pick_dist(); ch128_word h;
            if(backs[b]==0 && (len>8192 || hint || step==128)) continue;   /* portable ignores prefetch; keep it cheap */
            h=chainhash128_with_prefetch(&k,p,len,backs[b],hint,(unsigned)step,d); checks++;
            if(!ch128_equal(h,ref)) { if(bad<20) printf("FAIL len=%zu off=%zu b=%d hint=%d step=%d d=%u\n",len,off,backs[b],hint,step,d); bad++; }
        }
        { unsigned stride; int lazy, bi, school;
          for(bi=0;bi<nb-1;bi++) for(stride=1;stride<=8;stride++) for(lazy=0;lazy<2;lazy++) for(school=0;school<2;school++) {
            chainhash128_stream st; size_t done=0; int sb=backs[bi];
            if(sb==0 && (len>20000 || stride!=4 || !lazy)) continue;
            if(len>65536 && (stride&3)) continue;
            chainhash128_init(&st,&k,stride,lazy,sb,school);
            while(done<len) { size_t t=rnd()%4?(size_t)(rnd()%3000):(size_t)(rnd()%70000); if(t>len-done) t=len-done; chainhash128_update(&st,p+done,t); done+=t; if(rnd()%5==0) chainhash128_update(&st,p+done,0); }
            if(!ch128_equal(chainhash128_final(&st),ref)) { if(bad<20) printf("FAIL stream len=%zu stride=%u lazy=%d b=%d school=%d\n",len,stride,lazy,sb,school); bad++; } checks++; } }
    }
    printf("%s: %ld inputs, %ld evaluations, %ld mismatches\n",bad?"FAIL":"PASS",cases,checks,bad);
    free(mem); return bad!=0;
}
