/* make calibrate: measure this CPU's ChainHash and ChainHash-128 schedules
 * (chainhash_calibrate.h) and print them as C initializers, to be compiled in
 * and passed to chainhash_with_schedule / chainhash128_with_schedule.
 * Usage: calibrate [budget_us=0, i.e. about 2000]. */
#include <stdio.h>
#include <stdlib.h>
#include "chainhash_calibrate.h"
static void print(const char *name,const chainhash_schedule *s,double us) {
    int c;
    printf("/* calibrated in %.0f us of CPU time */\n",us);
    printf("static const chainhash_schedule %s = { 0x%08xu, 0x%08xu,\n    { %lluu, %lluu, %lluu },\n    {",name,s->magic,s->cpu,
           (unsigned long long)s->limit[0],(unsigned long long)s->limit[1],(unsigned long long)s->limit[2]);
    for(c=0;c<4;c++) printf(" { %u, %u, %u, 0, %uu }%s",s->k[c].backend,s->k[c].hint,s->k[c].step,s->k[c].dist,c<3?",":" } };\n");
}
int main(int argc,char **argv) {
    unsigned budget=argc>1?(unsigned)atoi(argv[1]):0; chainhash_schedule s; double us;
    us=chainhash_calibrate(&s,budget); print("chainhash_tuned",&s,us);
    us=chainhash128_calibrate(&s,budget); print("chainhash128_tuned",&s,us);
    return 0;
}
