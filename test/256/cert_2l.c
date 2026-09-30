/* cert_2l.c -- the two-level (2L) outer stage of PH-256 v1, certified on the AS-BUILT reference (ph_ref.h):
 *  (A) region polynomial: for every q = 1..8 and random region values a, the as-built ph2l_region(a, q, y) is
 *      interpolated as a polynomial in y (4f+1 nodes, checked at 8 more) and must equal Lemma A's expansion exactly:
 *        c = D(a) + sum_{i<=f} (a_i y^(2i) + a_{h+i} y^(2i-1)) + kappa_q(y),  D = sum a_i a_{h+i} + [q odd] a_h,
 *        kappa_q = sum_{i<=f} y^(4i-1)
 *      (so distinct values sit at distinct linear exponents 1..2f and all data x data terms at y^0);
 *      adversarial differences: single value j (all j), only the bare a_h (odd q), pair with a zero partner,
 *      both members of a pair, dense -- the difference polynomial must be nonzero with the predicted top term;
 *  (B) outer polynomial: for several block counts m (m' = 1 and m' >= 2, odd last q), V(n, b; y, z) interpolated in z
 *      must be n z^{m'} + sum_rho c_rho z^{m'-rho} with c_rho = the as-built region values;
 *  (C) the <= 1-block path: for n <= BLOCK the full as-built hash equals finish(n z + b_1) (b_1 = the block value),
 *      and n = 0 gives finish(0).
 * Exit 0 iff all pass. */
#include <stdio.h>
#include "ph_ref.h"
#ifdef NEG   /* negative control: a broken region function (both operands masked with y^(2i): exponents collide) */
static inline ph_el ph2l_region_bad(const ph_el*a,int q,const ph_el*Y){
    int h=(q+1)/2, f=q/2; ph_el c={{0,0,0,0}}, yp[9]; yp[0]=(ph_el){{1,0,0,0}}; for(int e=1;e<=8;e++) yp[e]=ph_mul(yp[e-1],*Y);
    for(int i=1;i<=f;i++) c=ph_xor(c,ph_mul(ph_xor(a[i-1],yp[2*i]),ph_xor(a[h+i-1],yp[2*i])));
    if(q&1) c=ph_xor(c,a[h-1]); return c; }
#define ph2l_region ph2l_region_bad
#endif
static uint64_t rs=2025; static uint64_t rnd(void){ return ch_splitmix(&rs); }
static ph_el rel(void){ ph_el e; for(int i=0;i<4;i++) e.w[i]=rnd(); return e; }
static int eq(ph_el a,ph_el b){ return !memcmp(&a,&b,sizeof a); }
static int iszero(ph_el a){ return !(a.w[0]|a.w[1]|a.w[2]|a.w[3]); }
static ph_el inv(ph_el a){ ph_el r={{1,0,0,0}}, b=a; for(int i=1;i<256;i++){ b=ph_mul(b,b); r=ph_mul(r,b); } return r; }   /* a^(2^256-2) */
/* coefficients of the polynomial of degree < n through (x_i, v_i): Newton divided differences -> monomial form */
static void interp(const ph_el*x,const ph_el*v,int n,ph_el*c){
    ph_el d[40]; for(int i=0;i<n;i++) d[i]=v[i];
    for(int j=1;j<n;j++) for(int i=n-1;i>=j;i--) d[i]=ph_mul(ph_xor(d[i],d[i-1]),inv(ph_xor(x[i],x[i-j])));
    for(int i=0;i<n;i++) c[i]=(ph_el){{0,0,0,0}};
    for(int j=n-1;j>=0;j--){ /* c <- c * (X - x_j) + d_j */
        ph_el t[40]; for(int i=0;i<n;i++) t[i]=(ph_el){{0,0,0,0}};
        for(int i=0;i<n-1;i++){ t[i+1]=ph_xor(t[i+1],c[i]); t[i]=ph_xor(t[i],ph_mul(c[i],x[j])); }
        t[0]=ph_xor(t[0],d[j]); for(int i=0;i<n;i++) c[i]=t[i]; } }
static ph_el peval(const ph_el*c,int n,ph_el x){ ph_el r={{0,0,0,0}}; for(int i=n-1;i>=0;i--) r=ph_xor(ph_mul(r,x),c[i]); return r; }
static int fails=0, checks=0;
#define CHECK(cond,...) do{ checks++; if(!(cond)){ fails++; printf("  FAIL: " __VA_ARGS__); printf("\n"); } }while(0)
/* interpolate y -> ph2l_region(a,q,y) (degree <= 4f-1 < 16) */
static void region_poly(const ph_el*a,int q,ph_el*c,int n){
    ph_el xs[24],vs[24]; for(int i=0;i<n+8;i++){ xs[i]=rel(); vs[i]=ph2l_region(a,q,&xs[i]); }
    interp(xs,vs,n,c);
    for(int i=n;i<n+8;i++) CHECK(eq(peval(c,n,xs[i]),vs[i]),"q=%d region polynomial degree > %d",q,n-1); }
int main(void){
    printf("[PH-256 2L] outer-stage certificate on the as-built reference (ph2l_region / ph2l_outer / ph2l_hash_ref)\n");
    /* (A) expansion, every q */
    for(int q=1;q<=8;q++){ int h=(q+1)/2, f=q/2, n=(f? 4*f : 1); int base=checks;
        for(int t=0;t<3;t++){ ph_el a[8]; for(int i=0;i<q;i++) a[i]=rel(); ph_el c[24]; region_poly(a,q,c,n);
            ph_el D={{0,0,0,0}}; for(int i=1;i<=f;i++) D=ph_xor(D,ph_mul(a[i-1],a[h+i-1])); if(q&1) D=ph_xor(D,a[h-1]);
            ph_el want[24]; for(int i=0;i<n;i++) want[i]=(ph_el){{0,0,0,0}}; want[0]=D;
            for(int i=1;i<=f;i++){ want[2*i]=ph_xor(want[2*i],a[i-1]); want[2*i-1]=ph_xor(want[2*i-1],a[h+i-1]); want[4*i-1].w[0]^=1; }
            for(int e=0;e<n;e++) CHECK(eq(c[e],want[e]),"q=%d coefficient of y^%d differs from Lemma A",q,e); }
        /* adversarial differences */
        ph_el a[8],b[8],ca[24],cb[24];
        for(int j=0;j<q;j++){ for(int i=0;i<q;i++) a[i]=b[i]=rel(); b[j]=ph_xor(b[j],rel());
            region_poly(a,q,ca,n); region_poly(b,q,cb,n); int nz=0; for(int e=0;e<n;e++) if(!eq(ca[e],cb[e])) nz=1;
            CHECK(nz,"q=%d single value %d: difference polynomial is zero",q,j);
            int e_j = (q&1 && j==h-1)? 0 : (j<h? 2*(j+1) : 2*(j-h+1)-1);          /* exponent of a_j's linear term (0: bare) */
            CHECK(eq(ph_xor(ca[e_j],cb[e_j]),ph_xor(a[j],b[j])),"q=%d value %d: coefficient at y^%d",q,j,e_j);
            for(int e=2*f+1;e<n;e++) CHECK(eq(ca[e],cb[e]),"q=%d: difference has a term above y^%d",q,2*f); }
        if(f){ for(int i=0;i<q;i++) a[i]=b[i]=rel(); a[h]=b[h]=(ph_el){{0,0,0,0}}; b[0]=ph_xor(b[0],rel());          /* zero partner */
            region_poly(a,q,ca,n); region_poly(b,q,cb,n); CHECK(!eq(ca[2],cb[2]),"q=%d pair 1 with a zero partner: y^2 term cancels",q);
            for(int i=0;i<q;i++) a[i]=b[i]=rel(); b[0]=ph_xor(b[0],rel()); b[h]=ph_xor(b[h],rel());                     /* both members */
            region_poly(a,q,ca,n); region_poly(b,q,cb,n); CHECK(!eq(ca[2],cb[2]) && !eq(ca[1],cb[1]),"q=%d both members of pair 1",q); }
        for(int i=0;i<q;i++){ a[i]=rel(); b[i]=rel(); } region_poly(a,q,ca,n); region_poly(b,q,cb,n);
        { int nz=0; for(int e=0;e<n;e++) if(!eq(ca[e],cb[e])) nz=1; CHECK(nz,"q=%d dense difference is zero",q); }
        printf("  (A) q=%d: region polynomial = Lemma A expansion; adversarial differences nonzero with the predicted terms  [%d checks]\n",q,checks-base); }
    /* (B) outer polynomial in z */
    { size_t ms[]={1,2,3,7,8,9,13,16,17,24,31}; ph_el Y=rel();
      for(size_t t=0;t<sizeof ms/sizeof ms[0];t++){ size_t m=ms[t], mp=(m+7)/8; ph_el b[32]; for(size_t i=0;i<m;i++) b[i]=rel(); uint64_t n=rnd()>>20; int base=checks;
        ph_el xs[48],vs[48],c[48]; int N=(int)mp+1; for(int i=0;i<N+6;i++){ xs[i]=rel(); vs[i]=ph2l_outer(b,m,n,&Y,&xs[i]); }
        interp(xs,vs,N,c); for(int i=N;i<N+6;i++) CHECK(eq(peval(c,N,xs[i]),vs[i]),"m=%zu outer degree > m'",m);
        CHECK(eq(c[mp],(ph_el){{n,0,0,0}}),"m=%zu: z^{m'} coefficient is not the length",m);
        for(size_t r=0;r<mp;r++){ int q=(int)((m-8*r)<8? m-8*r : 8); CHECK(eq(c[mp-1-r],ph2l_region(b+8*r,q,&Y)),"m=%zu region %zu coefficient",m,r); }
        printf("  (B) m=%-2zu (m'=%zu, last q=%zu): V = n z^%zu + sum c_rho z^(m'-rho) exactly  [%d checks]\n",m,mp,m-8*(mp-1),mp,checks-base); } }
    /* (C) the <= 1-block path through the full hash */
    { static ph_key K; ph_key_from_seed(&K,31337); static uint8_t msg[PH_BLOCK+8]; for(size_t i=0;i<sizeof msg;i++) msg[i]=(uint8_t)rnd(); int base=checks;
      size_t ns[]={0,1,7,8,63,64,255,256,257,PH_BLOCK/2,PH_BLOCK-1,PH_BLOCK};
      for(size_t t=0;t<sizeof ns/sizeof ns[0];t++){ size_t n=ns[t]; uint8_t o[32],w[32]; ph2l_hash_ref(&K,msg,n,o);
          ph_el V={{0,0,0,0}};
          if(n){ uint8_t blk[PH_BLOCK]; memcpy(blk,msg,n); memset(blk+n,0,PH_BLOCK-n); ph_el x[PH_M],y[PH_M]; ph_tail_block(blk,x,y);
                 V=ph_xor(ph_mul((ph_el){{(uint64_t)n,0,0,0}},K.z),ph_block(x,y,&K)); }
          ph_store(ph_finish(V,&K),w); CHECK(!memcmp(o,w,32),"n=%zu: hash != finish(n z + b_1)",n); }
      printf("  (C) <=1-block path: hash(n) = finish(n z + b_1) for n in {0..%d} samples (n=0: V=0)  [%d checks]\n",PH_BLOCK,checks-base); }
    printf("checks %d, failures %d\n%s\n",checks,fails,fails? "CERT 2L FAIL" : "CERT 2L PASS: region polynomial = Lemma A, outer = length-leading Horner in z, <=1-block path, odd counts");
    return fails!=0; }
