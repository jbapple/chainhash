/* chainhash_calibrate.h: opt-in, cacheable tuning of how ChainHash and
 * ChainHash-128 are evaluated. A schedule picks, per input-size class, the
 * backend and the software prefetch of the bulk kernels (hint, step, distance).
 * No schedule changes a digest: chainhash_with_schedule(k,m,n,&s) equals
 * chainhash(k,m,n) for every s, including a zeroed, stale or foreign one.
 *
 *   chainhash_schedule s;
 *   if(!load(&s) || !chainhash_schedule_valid(&s,64)) { chainhash_calibrate(&s,0); save(&s); }
 *   h=chainhash_with_schedule(&key,data,len,&s);
 *
 * chainhash_schedule_default (chainhash128_schedule_default) fills the built-in
 * table for this CPU family without measuring. chainhash_calibrate(&s,budget_us)
 * (chainhash128_calibrate) starts from the table and measures for about budget_us
 * of CPU time (0: 2000): the backend on a hot 64 KiB buffer, then, on x86, the
 * prefetch for inputs past L2 on a 2 MiB buffer evicted with CLFLUSHOPT before
 * every trial, which ranks the settings like a 256 MiB stream. Where eviction does
 * not reach DRAM (Apple: the lines stay in the SLC) the large classes keep the
 * table unless budget_us >= 50000, which probes a 64 MiB buffer directly. A
 * challenger must beat the table's setting by 3%. Classes below L2 always keep
 * the shipped kernels. The schedule is 64 bytes of plain data.
 * Include chainhash.h and/or chainhash128.h first (both if neither). C99. */
#ifndef CHAINHASH_CALIBRATE_H
#define CHAINHASH_CALIBRATE_H
#if !defined(CHAINHASH_H) && !defined(CHAINHASH128_H)
#include "chainhash.h"
#include "chainhash128.h"
#endif
#include <stdlib.h>
#include <string.h>
#include <time.h>
#if defined(__APPLE__)
#include <sys/types.h>
#include <sys/sysctl.h>
#elif defined(__linux__) && defined(__aarch64__)
#include <stdio.h>
#endif
#if !defined(CHAINHASH_PF_HINTS)
#error "include chainhash.h or chainhash128.h before chainhash_calibrate.h"
#endif
#if defined(CH_X86) || defined(CH128_X86)
#define CHC_X86 1
#include <cpuid.h>
#endif

#define CHAINHASH_AUTO 255                      /* knob backend: the detected one */
#define CHAINHASH_SCHEDULE_MAGIC(width) (0x43485300u|(uint32_t)(width))
typedef struct {
    uint8_t backend;    /* backend id (CH_XMM ... / CH128_XMM ...) or CHAINHASH_AUTO */
    uint8_t hint;       /* CH_PF_OFF, CH_PF_T0, CH_PF_T1, CH_PF_NTA */
    uint8_t step;       /* bytes between hints: 64 or 128 */
    uint8_t reserved;
    uint32_t dist;      /* bytes ahead of the region being hashed */
} chainhash_knobs;
typedef struct {
    uint32_t magic, cpu;          /* width tag; signature of the CPU it was tuned on */
    uint64_t limit[3];            /* class c holds len < limit[c]; class 3 the rest */
    chainhash_knobs k[4];
} chainhash_schedule;

/* ---- CPU signature and L2 size ---- */
typedef struct { uint32_t sig; int vendor; uint64_t l2; } chc_cpu;   /* vendor: 1 Intel, 2 AMD, 3 Apple, 4 other ARM */
static inline chc_cpu chc_cpuinfo(void) {
    chc_cpu c; c.sig=0; c.vendor=0; c.l2=0;
#ifdef CHC_X86
    { unsigned a,b,cc,d,i;
      __cpuid(0,a,b,cc,d);
      c.vendor = b==0x756e6547u ? 1 : b==0x68747541u ? 2 : 0;          /* "Genu", "Auth" */
      c.sig=(b*2654435761u)^cc;
      if(a>=1) { __cpuid(1,a,b,cc,d); c.sig^=a&0x0fff3ff0u; }         /* family and model, not stepping */
      if(c.vendor==1) for(i=0;i<8;i++) { unsigned w,bb,s,e; __cpuid_count(4,i,w,bb,s,e); if(!(w&31)) break;
          if(((w>>5)&7)==2) c.l2=(uint64_t)(((bb>>22)&0x3ff)+1)*(((bb>>12)&0x3ff)+1)*((bb&0xfff)+1)*((uint64_t)s+1); }
      else { __cpuid(0x80000000u,a,b,cc,d); if(a>=0x80000006u) { __cpuid(0x80000006u,a,b,cc,d); c.l2=(uint64_t)(cc>>16)*1024; } } }
#elif defined(__APPLE__)
    { uint64_t v=0,per=1; int32_t f=0; size_t n; c.vendor=3;
      n=sizeof(f); if(!sysctlbyname("hw.cpufamily",&f,&n,0,0)) c.sig=(uint32_t)f;
      n=sizeof(v); if(!sysctlbyname("hw.perflevel0.l2cachesize",&v,&n,0,0)) {
          n=sizeof(per); if(!sysctlbyname("hw.perflevel0.cpusperl2",&per,&n,0,0) && per>1) v/=per; c.l2=v; } }
#elif defined(__linux__) && defined(__aarch64__)
    { FILE *f=fopen("/sys/devices/system/cpu/cpu0/regs/identification/midr_el1","r"); unsigned long long m=0; unsigned long v; c.vendor=4;
      if(f) { if(fscanf(f,"%llx",&m)!=1) m=0; fclose(f); } c.sig=(uint32_t)(m&0xff0ffff0u);
      f=fopen("/sys/devices/system/cpu/cpu0/cache/index2/size","r"); if(f) { if(fscanf(f,"%luK",&v)==1) c.l2=(uint64_t)v*1024; fclose(f); } }
#endif
    if(c.l2<(256u<<10)) c.l2=(uint64_t)1<<20;
    c.sig|=1; return c;
}
static inline uint64_t chc_ticks(void) {
#ifdef CHC_X86
    unsigned lo,hi; __asm__ volatile("lfence\n\trdtsc\n\tlfence":"=a"(lo),"=d"(hi)::"memory"); return ((uint64_t)hi<<32)|lo;
#elif defined(__aarch64__) && (defined(__GNUC__) || defined(__clang__))
    uint64_t v; __asm__ volatile("isb\n\tmrs %0, cntvct_el0":"=r"(v)::"memory"); return v;
#else
    return (uint64_t)clock();
#endif
}
static inline double chc_cpu_us(void) { return (double)clock()*1e6/CLOCKS_PER_SEC; }
#ifdef CHC_X86
__attribute__((target("clflushopt"))) static inline void chc_flushopt(const uint8_t *p,size_t n) { size_t o; for(o=0;o<n;o+=64) _mm_clflushopt((void *)(p+o)); _mm_sfence(); }
__attribute__((target("sse2"))) static inline void chc_flush1(const uint8_t *p,size_t n) { size_t o; for(o=0;o<n;o+=64) _mm_clflush(p+o); _mm_mfence(); }
static inline int chc_flush(const uint8_t *p,size_t n) {   /* 1 if the range was evicted */
    unsigned a,b,c,d; __cpuid(0,a,b,c,d);
    if(a>=7) { __cpuid_count(7,0,a,b,c,d); if(b&(1u<<23)) { chc_flushopt(p,n); return 1; } }
    __cpuid(1,a,b,c,d); if(d&(1u<<19)) { chc_flush1(p,n); return 1; }
    return 0;
}
#else
static inline int chc_flush(const uint8_t *p,size_t n) { (void)p; (void)n; return 0; }
#endif

/* ---- schedules ---- */
static inline void chc_set(chainhash_knobs *k,int hint,int step,uint32_t dist) { k->hint=(uint8_t)hint; k->step=(uint8_t)step; k->dist=dist; }
/* The shipped evaluation: detected backend, no prefetch, in every class. */
static inline void chainhash_schedule_fixed(chainhash_schedule *s,int width) {
    unsigned c; memset(s,0,sizeof(*s)); s->magic=CHAINHASH_SCHEDULE_MAGIC(width);
    for(c=0;c<3;c++) s->limit[c]=UINT64_MAX;
    for(c=0;c<4;c++) { s->k[c].backend=CHAINHASH_AUTO; s->k[c].step=64; }
}
/* Classes: < 64 KiB | < L2 per core | < 64 x L2 | the rest. Measured on a Xeon
 * Platinum 8375C and an Apple M2 Pro; other CPUs keep the shipped evaluation.
 * Returns 1 if the table knows this CPU family. */
static inline int chc_default(chainhash_schedule *s,int width,int zmm) {
    chc_cpu c=chc_cpuinfo(); int known=0;
    chainhash_schedule_fixed(s,width); s->cpu=c.sig;
    s->limit[0]=64u<<10; s->limit[1]=c.l2; s->limit[2]=64*c.l2;
    if(c.vendor==1 && zmm) { known=1; chc_set(&s->k[2],CH_PF_T0,64,width==128?1024:4096); s->k[3]=s->k[2]; }
    else if(c.vendor==3) { known=1; if(width==64) { chc_set(&s->k[2],CH_PF_T1,128,65536); s->k[3]=s->k[2]; } }
    return known;
}
/* A stored schedule is usable if it has this width's tag and was tuned on this CPU. */
static inline int chainhash_schedule_valid(const chainhash_schedule *s,int width) {
    return s->magic==CHAINHASH_SCHEDULE_MAGIC(width) && s->cpu==chc_cpuinfo().sig;
}
static inline const chainhash_knobs *chainhash_schedule_class(const chainhash_schedule *s,uint64_t len) {
    unsigned c=0; while(c<3 && len>=s->limit[c]) c++; return &s->k[c];
}

/* ---- measurement ----
 * A probe hashes (p,n) under the caller's key with knobs k; any hash whose bulk
 * kernels take a prefetch schedule (ChainHash, ChainHash-128, CH-128/P) can be
 * calibrated with chainhash_calibrate_with. */
typedef uint64_t (*chainhash_probe)(const void *key,const uint8_t *p,size_t n,const chainhash_knobs *k);
/* Shortest of several timings per candidate, candidates interleaved; a round is
 * started only if one more (of the length measured so far) fits the deadline. */
static inline void chc_race(chainhash_probe fn,const void *key,uint8_t *p,size_t n,int reps,int cold,const chainhash_knobs *c,int nc,
                            uint64_t *best,double deadline,int minrounds,int *bad) {
    int r,i,j; uint64_t dig=0; double start=chc_cpu_us(),last=0;
    for(i=0;i<nc;i++) best[i]=UINT64_MAX;
    for(r=0;r<64;r++) {
        double t=chc_cpu_us(); if(r) last=(t-start)/r;
        if(r>=minrounds && t+last>deadline) break;
        for(i=0;i<nc;i++) { uint64_t t0,dt,h=0;
            if(cold) chc_flush(p,n);
            t0=chc_ticks(); for(j=0;j<reps;j++) h^=fn(key,p,n,&c[i]); dt=chc_ticks()-t0;
            if(dt<best[i]) best[i]=dt;
            if(reps%2==0) continue;
            if(!r && !i) dig=h; else if(h!=dig) *bad=1; }
    }
}
static inline int chc_pick(const uint64_t *best,int nc) {   /* index 0 is the incumbent */
    int i,b=0; for(i=1;i<nc;i++) if(best[i]<best[b]) b=i;
    return b && (double)best[b]*1.03<=(double)best[0] ? b : 0;
}
/* Refine the table-filled schedule s. has(b): backend b may be probed; prefetch:
 * the probe's kernels take prefetch hints on this target. Returns CPU us spent. */
static inline double chainhash_calibrate_with(chainhash_schedule *s,unsigned budget_us,chainhash_probe fn,const void *key,int (*has)(int),int prefetch) {
    double t0=chc_cpu_us(),deadline=t0+(budget_us?budget_us:2000);
    chainhash_knobs cand[8],table[4],base; uint64_t best[8],seed=UINT64_C(0x9e3779b97f4a7c15);
    int nc=0,i,bad=0,pick,big=budget_us>=50000;
    size_t hot=64u<<10,cold=big?(size_t)64<<20:(size_t)2<<20,n; uint8_t *buf;
    memcpy(table,s->k,sizeof(table));
    buf=(uint8_t *)malloc(cold); if(!buf) return chc_cpu_us()-t0;
    memset(buf+hot,0x5a,cold-hot);   /* the timing does not depend on the data */
    for(n=0;n<hot;n+=8) { seed^=seed<<13; seed^=seed>>7; seed^=seed<<17; memcpy(buf+n,&seed,8); }
    /* 1. backend on a hot buffer; index 0 is the shipped dispatch */
    base=table[0]; base.backend=CHAINHASH_AUTO; chc_set(&base,CH_PF_OFF,64,0); cand[nc++]=base;
    for(i=1;i<=4;i++) if(has(i)) { cand[nc]=base; cand[nc++].backend=(uint8_t)i; }
    chc_race(fn,key,buf,hot,5,0,cand,nc,best,t0+0.1*(deadline-t0),3,&bad);
    base=cand[chc_pick(best,nc)];
    for(i=0;i<4;i++) s->k[i].backend=base.backend;
    /* 2. prefetch past L2 on a cold buffer (x86 eviction, or a 64 MiB probe) */
    if(prefetch && (big || chc_flush(buf,4096))) {
#ifdef CHC_X86
        const int hint=CH_PF_T0,step=64; uint32_t d=1024;
#else
        const int hint=CH_PF_T1,step=128; uint32_t d=16384;
#endif
        nc=0; cand[nc++]=s->k[3];                             /* incumbent: the table's setting */
        if(s->k[3].hint) { cand[nc]=s->k[3]; chc_set(&cand[nc++],CH_PF_OFF,64,0); }
        for(i=0;i<3;i++,d*=4) if(!(s->k[3].hint==hint && s->k[3].step==step && s->k[3].dist==d)) { cand[nc]=s->k[3]; chc_set(&cand[nc++],hint,step,d); }
        if(big) fn(key,buf,cold,&cand[0]);                    /* fault it in; then it streams from DRAM */
        chc_race(fn,key,buf,cold,1,!big,cand,nc,best,deadline,big?2:1,&bad);
        pick=chc_pick(best,nc);
        s->k[3]=cand[pick];
        if(table[2].hint) s->k[2]=cand[pick];                 /* L2..64xL2 only where the table prefetches */
    }
    free(buf);
    if(bad) memcpy(s->k,table,sizeof(table));                 /* digests disagreed: keep the table */
    return chc_cpu_us()-t0;
}

#ifdef CHAINHASH_H
static inline uint64_t chainhash_with_schedule(const chainhash_key *k,const void *data,size_t len,const chainhash_schedule *s) {
    const chainhash_knobs *kn=chainhash_schedule_class(s,len);
    if(kn->backend==CHAINHASH_AUTO && !kn->hint) return chainhash(k,data,len);
    return chainhash_with_prefetch(k,data,len,kn->backend==CHAINHASH_AUTO?chainhash_backend():kn->backend,kn->hint,kn->step,kn->dist);
}
static inline int chainhash_schedule_default(chainhash_schedule *s) {
#ifdef CH_X86
    return chc_default(s,64,chainhash_backend()==CH_ZMM);
#else
    return chc_default(s,64,0);
#endif
}
static inline uint64_t chc_probe64(const void *key,const uint8_t *p,size_t n,const chainhash_knobs *k) {
    const chainhash_key *K=(const chainhash_key *)key;
    return chainhash_with_prefetch(K,p,n,k->backend==CHAINHASH_AUTO?chainhash_backend():k->backend,k->hint,k->step,k->dist);
}
static inline int chc_has64(int b) { return chainhash_has_backend(b); }
static inline double chainhash_calibrate(chainhash_schedule *s,unsigned budget_us) {
    chainhash_key k=chainhash_key_from_seed(7); chainhash_schedule_default(s);
#if defined(CH_X86) || defined(CH_ARM)
    return chainhash_calibrate_with(s,budget_us,chc_probe64,&k,chc_has64,1);
#else
    return chainhash_calibrate_with(s,budget_us,chc_probe64,&k,chc_has64,0);
#endif
}
#endif

#ifdef CHAINHASH128_H
static inline ch128_word chainhash128_with_schedule(const chainhash128_key *k,const void *data,size_t len,const chainhash_schedule *s) {
    const chainhash_knobs *kn=chainhash_schedule_class(s,len);
    if(kn->backend==CHAINHASH_AUTO && !kn->hint) return chainhash128(k,data,len);
    return chainhash128_with_prefetch(k,data,len,kn->backend==CHAINHASH_AUTO?chainhash128_backend():kn->backend,kn->hint,kn->step,kn->dist);
}
static inline int chainhash128_schedule_default(chainhash_schedule *s) {
#ifdef CH128_X86
    return chc_default(s,128,chainhash128_backend()==CH128_ZMM);
#else
    return chc_default(s,128,0);
#endif
}
static inline uint64_t chc_probe128(const void *key,const uint8_t *p,size_t n,const chainhash_knobs *k) {
    const chainhash128_key *K=(const chainhash128_key *)key; ch128_word h;
    h=chainhash128_with_prefetch(K,p,n,k->backend==CHAINHASH_AUTO?chainhash128_backend():k->backend,k->hint,k->step,k->dist);
    return h.lo^h.hi;
}
static inline int chc_has128(int b) { return chainhash128_has_backend(b); }
static inline double chainhash128_calibrate(chainhash_schedule *s,unsigned budget_us) {
    chainhash128_key k=chainhash128_key_from_seed(7); chainhash128_schedule_default(s);
#ifdef CH128_X86
    return chainhash_calibrate_with(s,budget_us,chc_probe128,&k,chc_has128,1);
#else
    return chainhash_calibrate_with(s,budget_us,chc_probe128,&k,chc_has128,0);
#endif
}
#endif
#endif
