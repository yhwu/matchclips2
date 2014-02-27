#ifndef _CHECK_RD_INSERT_H
#define _CHECK_RD_INSERT_H

/**** system headers ****/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <complex>
#include <algorithm>
#include <string>
#include <vector>
#include <unistd.h>
//#include <sys/types.h>
//#include <cstdlib>
using namespace std;

/* samtools header */
#include "samfunctions.h"

void check_rd_insert(samfile_t *fp_in,   bam_iter_t iter);
void optimize_resolution(samfile_t* fp_in, bam_index_t *bamidx, vector<cnv_st>& cnvlist);
void cnv_stat(samfile_t* fp_in, bam_index_t *bamidx, vector<cnv_st>& cnvlist);

    
#endif
