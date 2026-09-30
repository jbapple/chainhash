/* cert_pk.c -- PH-256 v2 POWER KEY, level 1, certified on the AS-BUILT reference (ph_ref.h: ph2_derive + ph_block):
 *  the block value as a function of the power key s, c(s) = sum_i (x_i + s^(2i-1)) (y_i + s^(2i)) (pairs i = 1..m), is
 *  interpolated as a polynomial in s over GF(2^256) (4m nodes, degree <= 4m-1; checked at 8 more nodes) and must equal
 *      c(s) = sum_i x_i y_i  +  sum_i ( x_i s^(2i) + y_i s^(2i-1) )  +  sum_i s^(4i-1)        (exponent classes)
 *  so every data word has its own exponent in 1..2m and all data x data terms sit at s^0.  Adversarial differences
 *  (single x_i / y_i at several positions, a full pair, several pairs, dense, x-only with y = 0) must give a nonzero
 *  difference polynomial of degree <= 2m with the predicted coefficient; the 1-word case (pair 1, x only, y = 0) must be
 *  the monomial dx s^2 (one root: s = 0, squaring is bijective) -> d(1) = 1.
 *  Negative control (-DNEG): a derivation with a colliding exponent (l_j = s^(2j+1) = k_j) must be rejected.
 *  PH_M small by default (8) so the interpolation is quick; -DPH_M=64 runs the frozen size.
 *  -DDUMP: print 300 field products a b c (hex, limb 3 first) for ph_oracle.py's independent check of ph_mul. */
#include <stdio.h>
#include "ph_ref.h"
static uint64_t rs=77; static uint64_t rnd(void){ return ch_splitmix(&rs); }
static ph_el rel(void){ ph_el e; for(int i=0;i<4;i++) e.w[i]=rnd(); return e; }
static int eq(ph_el a,ph_el b){ return !memcmp(&a,&b,sizeof a); }
static int iszero(ph_el a){ return !(a.w[0]|a.w[1]|a.w[2]|a.w[3]); }
static ph_el inv(ph_el a){ ph_el r={{1,0,0,0}}, b=a; for(int i=1;i<256;i++){ b=ph_mul(b,b); r=ph_mul(r,b); } return r; }
#define NMAX (4*PH_M+8)
static void interp(const ph_el*x,const ph_el*v,int n,ph_el*c){
    static ph_el d[NMAX]; for(int i=0;i<n;i++) d[i]=v[i];
    for(int j=1;j<n;j++) for(int i=n-1;i>=j;i--) d[i]=ph_mul(ph_xor(d[i],d[i-1]),inv(ph_xor(x[i],x[i-j])));
    for(int i=0;i<n;i++) c[i]=(ph_el){{0,0,0,0}};
    for(int j=n-1;j>=0;j--){ static ph_el t[NMAX]; for(int i=0;i<n;i++) t[i]=(ph_el){{0,0,0,0}};
        for(int i=0;i<n-1;i++){ t[i+1]=ph_xor(t[i+1],c[i]); t[i]=ph_xor(t[i],ph_mul(c[i],x[j])); }
        t[0]=ph_xor(t[0],d[j]); for(int i=0;i<n;i++) c[i]=t[i]; } }
static ph_el peval(const ph_el*c,int n,ph_el x){ ph_el r={{0,0,0,0}}; for(int i=n-1;i>=0;i--) r=ph_xor(ph_mul(r,x),c[i]); return r; }
static int fails=0, checks=0;
#define CHECK(cond,...) do{ checks++; if(!(cond)){ fails++; if(fails<8){ printf("  FAIL: " __VA_ARGS__); printf("\n"); } } }while(0)
static ph_el Rtmp_mul(ph_el a,ph_el b){ return ph_mul(a,b); }
/* the as-built block value at key s */
static ph_el block_at(const ph_el*x,const ph_el*y,ph_el s){
    static ph_key K; ph2_raw R; memset(&R,0,sizeof R); R.s=s; R.y=(ph_el){{1,0,0,0}};
#ifdef NEG
    ph_el p=s; for(int j=0;j<PH_M;j++){ K.k[j]=p; K.l[j]=p; p=ph_mul(ph_mul(p,s),s); }   /* broken: l_j = k_j */
#else
    ph2_derive_with(&R,&K,Rtmp_mul);
#endif
    return ph_block(x,y,&K); }
static int N=4*PH_M;
static void poly(const ph_el*x,const ph_el*y,ph_el*c){
    static ph_el xs[NMAX],vs[NMAX]; for(int i=0;i<N+8;i++){ xs[i]=rel(); vs[i]=block_at(x,y,xs[i]); }
    interp(xs,vs,N,c); for(int i=N;i<N+8;i++) CHECK(eq(peval(c,N,xs[i]),vs[i]),"degree > %d",N-1); }
int main(void){
#ifdef DUMP
    for(int t=0;t<300;t++){ ph_el a=rel(), b=rel(); if(t<4){ memset(&a,0,sizeof a); a.w[3]=1ULL<<63; } ph_el c=ph_mul(a,b);
        for(int i=0;i<4;i++) printf("%016llx",(unsigned long long)a.w[3-i]); printf(" ");
        for(int i=0;i<4;i++) printf("%016llx",(unsigned long long)b.w[3-i]); printf(" ");
        for(int i=0;i<4;i++) printf("%016llx",(unsigned long long)c.w[3-i]); printf("\n"); }
    return 0;
#endif
    printf("[PH-256 v2] power-key level-1 certificate on the as-built reference (ph2_derive + ph_block), m=%d pairs\n",PH_M);
    static ph_el x[PH_M],y[PH_M],x2[PH_M],y2[PH_M],c[NMAX],c2[NMAX];
    int m=PH_M;
    /* exact expansion */
    for(int t=0;t<2;t++){ for(int i=0;i<m;i++){ x[i]=rel(); y[i]=rel(); } poly(x,y,c);
        static ph_el want[NMAX]; for(int e=0;e<N;e++) want[e]=(ph_el){{0,0,0,0}};
        for(int i=1;i<=m;i++){ want[0]=ph_xor(want[0],ph_mul(x[i-1],y[i-1])); want[2*i]=ph_xor(want[2*i],x[i-1]); want[2*i-1]=ph_xor(want[2*i-1],y[i-1]); want[4*i-1].w[0]^=1; }
        for(int e=0;e<N;e++) CHECK(eq(c[e],want[e]),"coefficient of s^%d differs from the exponent-class expansion",e); }
    printf("  exact expansion c(s) = sum x_i y_i + sum (x_i s^2i + y_i s^(2i-1)) + sum s^(4i-1): %s\n",fails?"FAIL":"PASS");
    /* adversarial differences */
    int pos[]={0,1,m/2,m-1}; int base=checks;
    for(int pi=0;pi<4;pi++){ int j=pos[pi];
        for(int kind=0;kind<3;kind++){ for(int i=0;i<m;i++){ x[i]=x2[i]=rel(); y[i]=y2[i]=rel(); }
            if(kind!=1) x2[j]=ph_xor(x2[j],rel()); if(kind!=0) y2[j]=ph_xor(y2[j],rel());
            poly(x,y,c); poly(x2,y2,c2); int nz=0; for(int e=0;e<N;e++) if(!eq(c[e],c2[e])) nz=1;
            CHECK(nz,"pair %d kind %d: difference polynomial is zero",j,kind);
            for(int e=2*m+1;e<N;e++) CHECK(eq(c[e],c2[e]),"pair %d: difference above s^%d",j,2*m);
            if(kind!=1) CHECK(eq(ph_xor(c[2*(j+1)],c2[2*(j+1)]),ph_xor(x[j],x2[j])),"pair %d: coefficient of s^%d is not dx",j,2*(j+1));
            if(kind!=0) CHECK(eq(ph_xor(c[2*(j+1)-1],c2[2*(j+1)-1]),ph_xor(y[j],y2[j])),"pair %d: coefficient of s^%d is not dy",j,2*(j+1)-1); } }
    { for(int i=0;i<m;i++){ x[i]=x2[i]=rel(); y[i]=y2[i]=rel(); } for(int i=0;i<m;i+=3) x2[i]=ph_xor(x2[i],rel());
      poly(x,y,c); poly(x2,y2,c2); int nz=0; for(int e=0;e<N;e++) if(!eq(c[e],c2[e])) nz=1; CHECK(nz,"several pairs: zero difference"); }
    { for(int i=0;i<m;i++){ x[i]=rel(); y[i]=rel(); x2[i]=rel(); y2[i]=rel(); } poly(x,y,c); poly(x2,y2,c2);
      int nz=0; for(int e=0;e<N;e++) if(!eq(c[e],c2[e])) nz=1; CHECK(nz,"dense: zero difference"); }
    printf("  adversarial differences (single x/y at pairs 1, 2, m/2+1, m; full pair; several; dense): %s [%d checks]\n",fails?"FAIL":"PASS",checks-base);
    /* d(1) = 1: pair 1 x only, y = 0 -> dx s^2 */
    { for(int i=0;i<m;i++){ x[i]=x2[i]=(ph_el){{0,0,0,0}}; y[i]=y2[i]=(ph_el){{0,0,0,0}}; } x2[0].w[0]=rnd()|1;
      poly(x,y,c); poly(x2,y2,c2); int only=1; for(int e=0;e<N;e++){ ph_el d=ph_xor(c[e],c2[e]); if(e==2) CHECK(eq(d,x2[0]),"1-word: s^2 coefficient"); else if(!iszero(d)) only=0; }
      CHECK(only,"1-word difference is not the monomial dx s^2");
      printf("  1-word difference (pair 1, x limb 0, y = 0) = dx * s^2 exactly -> one root (s = 0): d(1) = 1: %s\n",only?"PASS":"FAIL"); }
    printf("checks %d, failures %d\n%s\n",checks,fails,fails? "CERT PK FAIL" : "CERT PK PASS: power-key exponent classes are distinct (x_i -> s^2i, y_i -> s^(2i-1)), data x data at s^0; d(1) = 1");
    return fails!=0; }
