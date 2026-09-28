/* Every reachable evaluation schedule gives the reference digest. For each random
 * input (random key, length, misalignment and contents) the independent serial
 * evaluator chainhash_portable is compared with chainhash_with_prefetch for every
 * backend (and an unavailable id) x prefetch hint x hint step, each with a random
 * distance; with chainhash_with_schedule through a random schedule; and with
 * streaming for every backend x stride 1..8 x eager/lazy, each with random splits.
 * Every 1000 inputs a calibrated and the table schedule are checked as well.
 * Usage: schedule_knobs [cases=20000] [seed=1]. */
#include <stdio.h>
#include <stdlib.h>
#include "chainhash.h"
#include "chainhash_calibrate.h"
static uint64_t rs;
static uint64_t rnd(void) { uint64_t z=(rs+=UINT64_C(0x9e3779b97f4a7c15)); z=(z^(z>>30))*UINT64_C(0xbf58476d1ce4e5b9); z=(z^(z>>27))*UINT64_C(0x94d049bb133111eb); return z^(z>>31); }
static size_t pick_len(void) {
    unsigned r=(unsigned)(rnd()%100); size_t base, d[5]={0,1,8,17,1023};
    if(r<35) return (size_t)(rnd()%3000);
    if(r<75) return (size_t)(rnd()%20000);
    if(r<95) { base=1024*(size_t)(1+rnd()%24); return rnd()&1 ? base+d[rnd()%5] : base-d[rnd()%5]; }
    return (size_t)(rnd()%(300*1024));
}
static uint32_t pick_dist(void) { uint32_t d[8]={0,64,1000,1024,4096,8192,65536,1u<<20}; return rnd()&1 ? d[rnd()%8] : (uint32_t)(rnd()%(1u<<18)); }
int main(int argc,char **argv) {
    long cases=argc>1?atol(argv[1]):20000, c, checks=0, bad=0;
    int backs[8], nb=0, b, hint, step;
    uint8_t *mem=(uint8_t *)malloc(300*1024+128);
    chainhash_schedule cal, tab;
    rs=argc>2?(uint64_t)atoll(argv[2]):1;
    for(b=0;b<=4;b++) if(chainhash_has_backend(b)) backs[nb++]=b;
    backs[nb++]=77;   /* an unavailable id falls back, it does not trap */
    printf("backends:"); for(b=0;b<nb;b++) printf(" %d",backs[b]); printf("\n");
    chainhash_calibrate(&cal,3000); chainhash_schedule_default(&tab);
    if(!chainhash_schedule_valid(&cal,64) || chainhash_schedule_valid(&cal,128) || sizeof(chainhash_schedule)!=64) { printf("FAIL schedule tag\n"); bad++; }
    for(c=0;c<cases;c++) {
        size_t len=pick_len(), off=(size_t)(rnd()%64), i; uint8_t *p=mem+off; chainhash_key k; uint64_t ref;
        uint64_t w[39]; unsigned kind=(unsigned)(rnd()%100);
        for(i=0;i<39;i++) w[i]=rnd();
        if(kind==0) w[32]=0; else if(kind==1) w[32]=1;   /* degenerate chain keys */
        k=chainhash_key_from_words(w);
        for(i=0;i<len;i++) p[i]=(uint8_t)rnd();
        if(rnd()%8==0) memset(p,(int)(rnd()&0xff),len);
        ref=chainhash_portable(&k,p,len);
        if(chainhash(&k,p,len)!=ref) { printf("FAIL chainhash len=%zu\n",len); bad++; }
        for(b=0;b<nb;b++) for(hint=0;hint<=4;hint++) for(step=64;step<=128;step+=64) {
            uint32_t d=pick_dist(); uint64_t h;
            if(backs[b]==0 && len>65536 && (hint|(step-64))) continue;   /* portable ignores prefetch; keep one */
            h=chainhash_with_prefetch(&k,p,len,backs[b],hint,(unsigned)step,d); checks++;
            if(h!=ref) { if(bad<20) printf("FAIL len=%zu off=%zu b=%d hint=%d step=%d d=%u\n",len,off,backs[b],hint,step,d); bad++; }
        }
        { chainhash_schedule s; unsigned j; memset(&s,0,sizeof(s));
          for(j=0;j<3;j++) s.limit[j]=(j?s.limit[j-1]:0)+rnd()%(64*1024);
          for(j=0;j<4;j++) { s.k[j].backend=(uint8_t)(rnd()%4==0?CHAINHASH_AUTO:backs[rnd()%nb]); s.k[j].hint=(uint8_t)(rnd()%5); s.k[j].step=(uint8_t)(rnd()&1?128:64); s.k[j].dist=pick_dist(); }
          if(chainhash_with_schedule(&k,p,len,&s)!=ref) { printf("FAIL schedule len=%zu\n",len); bad++; } checks++;
          if(c%1000==0) { checks+=2;
              if(chainhash_with_schedule(&k,p,len,&cal)!=ref || chainhash_with_schedule(&k,p,len,&tab)!=ref) { printf("FAIL calibrated/table schedule len=%zu\n",len); bad++; } } }
        { unsigned stride; int lazy, bi;
          for(bi=0;bi<nb-1;bi++) for(stride=1;stride<=8;stride++) for(lazy=0;lazy<2;lazy++) {
            chainhash_stream st; size_t done=0; int sb=backs[bi];
            if(sb==0 && len>20000 && (stride!=4 || !lazy)) continue;
            chainhash_init(&st,&k,stride,lazy,sb);
            while(done<len) { size_t t=rnd()%4?(size_t)(rnd()%3000):(size_t)(rnd()%70000); if(t>len-done) t=len-done; chainhash_update(&st,p+done,t); done+=t; if(rnd()%5==0) chainhash_update(&st,p+done,0); }
            if(chainhash_final(&st)!=ref) { if(bad<20) printf("FAIL stream len=%zu stride=%u lazy=%d b=%d\n",len,stride,lazy,sb); bad++; } checks++; } }
    }
    printf("%s: %ld inputs, %ld evaluations, %ld mismatches\n",bad?"FAIL":"PASS",cases,checks,bad);
    free(mem); return bad!=0;
}
