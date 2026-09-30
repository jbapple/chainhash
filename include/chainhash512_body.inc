/* chainhash512_body.inc -- per-backend driver for chainhash512.h (included once per backend).
 * Parameters: PB(name) suffixes names, PBT is the target attribute, PB_POINTS(k,p,nch,P) is the
 * full-chunk kernel, and the PV* macros are the 128-bit vector primitives of the ISA.
 * All GF(2^512) arithmetic keeps elements as 4 x 128-bit vectors (limbs 2i, 2i+1).
 */
/* Karatsuba^3 points of X[4] in vector order (per 4-limb group: X0, X1, X0^X1, horizontal pair, half) */
PBT static inline __attribute__((always_inline)) void PB(ph1_evalv)(const PHV *X, PHV E[3][5]){
    PHV G[3][2]={{X[0],X[1]},{X[2],X[3]},{PVX(X[0],X[2]),PVX(X[1],X[3])}};
    for(int g=0;g<3;g++){ PHV x0=G[g][0], x1=G[g][1], S=PVX(x0,x1);
        E[g][0]=x0; E[g][1]=x1; E[g][2]=S; E[g][3]=PVX(PVULO(x0,x1),PVUHI(x0,x1)); E[g][4]=PVX(S,PVHSW(S)); }
}
/* P[27] ^= Ea (.) Eb, canonical point order */
PBT static inline __attribute__((always_inline)) void PB(ph1_pmacc)(PHV Ea[3][5], PHV Eb[3][5], PHV *Pv){
    for(int g=0;g<3;g++){ PHV *P=Pv+9*g;
        P[0]=PVX(P[0],PVC00(Ea[g][0],Eb[g][0])); P[1]=PVX(P[1],PVC11(Ea[g][0],Eb[g][0]));
        P[3]=PVX(P[3],PVC00(Ea[g][1],Eb[g][1])); P[4]=PVX(P[4],PVC11(Ea[g][1],Eb[g][1]));
        P[6]=PVX(P[6],PVC00(Ea[g][2],Eb[g][2])); P[7]=PVX(P[7],PVC11(Ea[g][2],Eb[g][2]));
        P[2]=PVX(P[2],PVC00(Ea[g][3],Eb[g][3])); P[5]=PVX(P[5],PVC11(Ea[g][3],Eb[g][3]));
        P[8]=PVX(P[8],PVC00(Ea[g][4],Eb[g][4])); }
}
PBT static inline __attribute__((always_inline)) void PB(ph1_eload)(const uint64_t e[3][5][2], PHV E[3][5]){
    for(int g=0;g<3;g++) for(int i=0;i<5;i++) E[g][i]=PVLD(e[g][i]); }
/* reduce 1024 -> 512: R[0..3] low (R[4] scratch), Hh[0..3] high */
PBT static inline __attribute__((always_inline)) void PB(ph1_red)(PHV *R, PHV *Hh, PHV *O){
    R[4]=PVZ();
    for(int k=0;k<4;k++){ PHV t=Hh[k];
        PHV L=PVX(PVX(t,PVSL(t,2)),PVX(PVSL(t,5),PVSL(t,8)));
        PHV C=PVX(PVX(PVSR(t,62),PVSR(t,59)),PVSR(t,56));
        R[k]=PVX(PVX(R[k],L),PVUP(C)); R[k+1]=PVX(R[k+1],PVDN(C)); }
    { PHV t=R[4]; R[0]=PVX(R[0],PVX(PVX(t,PVSL(t,2)),PVX(PVSL(t,5),PVSL(t,8)))); }
    O[0]=R[0]; O[1]=R[1]; O[2]=R[2]; O[3]=R[3];
}
PBT static inline __attribute__((always_inline)) void PB(ph1_r2)(const PHV *P, PHV *r){
    PHV m=PVX(PVX(P[0],P[1]),P[2]); r[0]=PVX(P[0],PVUP(m)); r[1]=PVX(P[1],PVDN(m)); }
PBT static inline __attribute__((always_inline)) void PB(ph1_r4)(const PHV *P, PHV *r){
    PHV a[2],b[2],c[2]; PB(ph1_r2)(P,a); PB(ph1_r2)(P+3,b); PB(ph1_r2)(P+6,c);
    r[0]=a[0]; r[1]=PVX(a[1],PVX(PVX(a[0],b[0]),c[0])); r[2]=PVX(b[0],PVX(PVX(a[1],b[1]),c[1])); r[3]=b[1]; }
/* recombine 27 point sums and reduce: O = R(P) */
PBT static inline void PB(ph1_finish)(const PHV *P, PHV *O){
    PHV a[4],b[4],c[4],m[4]; PB(ph1_r4)(P,a); PB(ph1_r4)(P+9,b); PB(ph1_r4)(P+18,c);
    for(int i=0;i<4;i++) m[i]=PVX(PVX(a[i],b[i]),c[i]);
    PHV R[5]={a[0],a[1],PVX(a[2],m[0]),PVX(a[3],m[1]),PVZ()}, Hh[4]={PVX(b[0],m[2]),PVX(b[1],m[3]),b[2],b[3]};
    PB(ph1_red)(R,Hh,O);
}
/* squaring (char 2: 8 self-products) */
PBT static inline void PB(ph1_sqr)(const PHV *A, PHV *O){
    PHV R[5]={PVC00(A[0],A[0]),PVC11(A[0],A[0]),PVC00(A[1],A[1]),PVC11(A[1],A[1]),PVZ()};
    PHV Hh[4]={PVC00(A[2],A[2]),PVC11(A[2],A[2]),PVC00(A[3],A[3]),PVC11(A[3],A[3])};
    PB(ph1_red)(R,Hh,O);
}
/* multiply: Karatsuba over 128-bit lanes (3 levels, 27 clmul) */
PBT static inline __attribute__((always_inline)) void PB(ph1_k2)(PHV a, PHV b, PHV *lo, PHV *hi){
    PHV ll=PVM00(a,b), hh=PVM11(a,b), m=PVX(PVX(PVM00(PVX(a,PVSWAP(a)),PVX(b,PVSWAP(b))),ll),hh);
    *lo=PVX(ll,PVUP(m)); *hi=PVX(hh,PVDN(m)); }
PBT static inline __attribute__((always_inline)) void PB(ph1_k4)(PHV a0, PHV a1, PHV b0, PHV b1, PHV *r){
    PHV q0l,q0h,q1l,q1h,q2l,q2h; PB(ph1_k2)(a0,b0,&q0l,&q0h); PB(ph1_k2)(a1,b1,&q1l,&q1h); PB(ph1_k2)(PVX(a0,a1),PVX(b0,b1),&q2l,&q2h);
    PHV ml=PVX(PVX(q0l,q1l),q2l), mh=PVX(PVX(q0h,q1h),q2h); r[0]=q0l; r[1]=PVX(q0h,ml); r[2]=PVX(q1l,mh); r[3]=q1h; }
PBT __attribute__((noinline)) static void PB(ph1_mul)(const PHV *A, const PHV *B, PHV *O){
    PHV q0[4],q1[4],q2[4];
    PB(ph1_k4)(A[0],A[1],B[0],B[1],q0); PB(ph1_k4)(A[2],A[3],B[2],B[3],q1); PB(ph1_k4)(PVX(A[0],A[2]),PVX(A[1],A[3]),PVX(B[0],B[2]),PVX(B[1],B[3]),q2);
    PHV m[4]; for(int i=0;i<4;i++) m[i]=PVX(PVX(q0[i],q1[i]),q2[i]);
    PHV R[5]={q0[0],q0[1],PVX(q0[2],m[0]),PVX(q0[3],m[1]),PVZ()}, Hh[4]={PVX(q1[0],m[2]),PVX(q1[1],m[3]),q1[2],q1[3]};
    PB(ph1_red)(R,Hh,O);
}
/* unreduced point sums of block t (full chunks by the kernel, tail pairs in registers) */
PBT static inline void PB(ph1_block_points)(const ph512v1_key *k, const uint8_t *p, size_t n, size_t t, PHV *Pv){
    size_t F=n/PH1_CHUNK, c0=t*PH1_CPB, c1=c0+PH1_CPB; if(c1>F) c1=F;
    if(c1>c0){ ph1_u128 P[28]; PB_POINTS(k,p+c0*PH1_CHUNK,(int)(c1-c0),P); for(int i=0;i<27;i++) Pv[i]=PVLD(&P[i]); }
    else for(int i=0;i<27;i++) Pv[i]=PVZ();
    if(c1==F && F<c0+PH1_CPB && F*PH1_CHUNK<n){
        size_t rem=n-F*PH1_CHUNK, q0=8*(F-c0); const uint8_t *b=p+F*PH1_CHUNK;
        for(size_t i=0;128*i<rem;i++){ size_t m=rem-128*i; const uint8_t *src=b+128*i; uint8_t buf[128];
            if(m<128){ memset(buf,0,128); memcpy(buf,src,m); src=buf; }
            PHV U[4],W[4],Eu[3][5],Ew[3][5];
            for(int j=0;j<4;j++){ U[j]=PVX(PVLD(src+16*j),PVLD(k->kup[q0+i]+2*j)); W[j]=PVX(PVLD(src+64+16*j),PVLD(k->kvp[q0+i]+2*j)); }
            PB(ph1_evalv)(U,Eu); PB(ph1_evalv)(W,Ew); PB(ph1_pmacc)(Eu,Ew,Pv); }
    }
}
/* region step: V <- V z + sum_i (b_i + y^(2i-1))(b_{h+i} + y^(2i)) [+ bare a_h], one reduction.
 * Pbare: unreduced points of the bare block (odd q) or NULL (then b[h-1] is added reduced).      */
PBT static inline void PB(ph1_region)(const ph512v1_key *k, PHV *V, PHV (*b)[4], int q, const PHV *Pbare){
    PHV P[27], E1[3][5], E2[3][5]; int h=(q+1)/2, f=q/2;
    if(Pbare) for(int i=0;i<27;i++) P[i]=Pbare[i]; else for(int i=0;i<27;i++) P[i]=PVZ();
    PB(ph1_evalv)(V,E1); PB(ph1_eload)(k->ezv,E2); PB(ph1_pmacc)(E1,E2,P);
    for(int i=1;i<=f;i++){ PHV X[4],W[4];
        for(int j=0;j<4;j++){ X[j]=PVX(b[i-1][j],PVLD(k->yp[2*i-1]+2*j)); W[j]=PVX(b[h+i-1][j],PVLD(k->yp[2*i]+2*j)); }
        PB(ph1_evalv)(X,E1); PB(ph1_evalv)(W,E2); PB(ph1_pmacc)(E1,E2,P); }
    PB(ph1_finish)(P,V);
    if((q&1) && !Pbare) for(int j=0;j<4;j++) V[j]=PVX(V[j],b[h-1][j]);
}
PBT __attribute__((noinline)) static void PB(ph1_final)(const ph512v1_key *k, const PHV *V, uint8_t *out){
    uint64_t w[8]; for(int i=0;i<4;i++) PVST(w+2*i,V[i]); ph1_add512(w,k->tau);
    PHV v[4],q[4],t1[4],t2[4],r[4],h[4]; for(int i=0;i<4;i++) v[i]=PVMK(w[2*i],w[2*i+1]);
    PB(ph1_sqr)(v,q);
    for(int i=0;i<4;i++){ t1[i]=PVX(q[i],PVLD(k->c[0]+2*i)); t2[i]=PVX(PVX(v[i],q[i]),PVLD(k->c[1]+2*i)); } PB(ph1_mul)(t1,t2,r);
    for(int i=0;i<4;i++){ t1[i]=PVX(v[i],PVLD(k->c[2]+2*i)); t2[i]=PVX(r[i],PVLD(k->c[3]+2*i)); } PB(ph1_mul)(t1,t2,h);
    for(int i=0;i<4;i++) PVST(out+16*i,PVX(h[i],PVLD(k->c[4]+2*i)));
}
/* V (the value before the finalizer) of the whole message */
PBT static void PB(ph1_hashV)(const ph512v1_key *k, const uint8_t *p, size_t n, PHV *V){
    size_t m=ph1_nblocks(n), mr=(m+PH1_R-1)/PH1_R;
    V[0]=PVMK((uint64_t)n,0); V[1]=V[2]=V[3]=PVZ();
    if(m==1){                                   /* <= 1 block: V = l z + b_1, one reduction (q = 1, bare) */
        PHV P[27], E1[3][5], E2[3][5]; PB(ph1_block_points)(k,p,n,0,P);
        PB(ph1_evalv)(V,E1); PB(ph1_eload)(k->ezv,E2); PB(ph1_pmacc)(E1,E2,P); PB(ph1_finish)(P,V);
        return; }
    for(size_t r=0;r<mr;r++){ int q = r+1<mr ? PH1_R : (int)(m-PH1_R*(mr-1)); int h=(q+1)/2;
        PHV b[PH1_R][4], Pb[27], Pt[27];
        for(int j=0;j<q;j++){ PB(ph1_block_points)(k,p,n,PH1_R*r+j,Pt);
            if((q&1) && j==h-1) for(int i=0;i<27;i++) Pb[i]=Pt[i]; else PB(ph1_finish)(Pt,b[j]); }
        PB(ph1_region)(k,V,b,q,(q&1)?Pb:0); }
}
PBT __attribute__((unused)) static void PB(ph1_hash)(const ph512v1_key *k, const uint8_t *p, size_t n, uint8_t *out){
    PHV V[4]; PB(ph1_hashV)(k,p,n,V); PB(ph1_final)(k,V,out); }
/* u64 interfaces (streaming) */
PBT static void PB(ph1_blockval)(const ph512v1_key *k, const uint8_t *p, size_t n, uint64_t *out){
    PHV Pt[27], O[4]; PB(ph1_block_points)(k,p,n,0,Pt); PB(ph1_finish)(Pt,O); for(int i=0;i<4;i++) PVST(out+2*i,O[i]); }
PBT static void PB(ph1_regionval)(const ph512v1_key *k, uint64_t *W, uint64_t bb[][8], int q){
    PHV V[4], b[PH1_R][4]; for(int i=0;i<4;i++) V[i]=PVLD(W+2*i);
    for(int j=0;j<q;j++) for(int i=0;i<4;i++) b[j][i]=PVLD(bb[j]+2*i);
    PB(ph1_region)(k,V,b,q,0); for(int i=0;i<4;i++) PVST(W+2*i,V[i]); }
PBT __attribute__((unused)) static void PB(ph1_finalu)(const ph512v1_key *k, const uint64_t *W, uint8_t *out){
    PHV V[4]; for(int i=0;i<4;i++) V[i]=PVLD(W+2*i); PB(ph1_final)(k,V,out); }
PBT static void PB(ph1_mulu)(const uint64_t *a, const uint64_t *b, uint64_t *o){
    PHV A[4],B[4],O[4]; for(int i=0;i<4;i++){ A[i]=PVLD(a+2*i); B[i]=PVLD(b+2*i); } PB(ph1_mul)(A,B,O); for(int i=0;i<4;i++) PVST(o+2*i,O[i]); }
