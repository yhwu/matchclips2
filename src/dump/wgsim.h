#ifndef WGSIM_H
#define WGSIM_H
#include <inttypes.h> 

void wgsim_core(FILE *fpout1, FILE *fpout2, const char *fn, int is_hap, uint64_t N, int dist, int std_dev, int size_l, int size_r, int POS1, int POS2, double er, double mr, double ir, double ie, double maxnr);

#endif
