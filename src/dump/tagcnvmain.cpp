/*************************************************************************
BEGIN OF LICENSE
Copyright (c) Yinghua Wu and Hongzhe Li (rsicnv project).

This program is free software; you can redistribute and/or modify
the codes written by the authors under the terms of the GNU General 
Public License as published by the Free Software Foundation 
(www.fsf.org); either version 2 of the License, or (at your option) 
any later version. 

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; See the GNU General Public License for 
more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses
END OF LICENSE

CONTACT: wu_yinghua@hotmail.com; hongzhe@upenn.edu
*************************************************************************/
/**** system headers ****/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <complex>
#include <algorithm>
#include <string>
#include <cstring>
#include <vector>
#include <unistd.h>
#include <time.h>
using namespace std;

#include "functions.h"
#include "readref.h"
#include "matchreads.h"
#include "tagcnv.h"

int usage_main(int argc, char* argv[]) {
  cerr << "tagcnv is a combination of a few tools to check some properties of CNV regions.\n"
       << "\nUsage:\n"
       << "  " << argv[0] << " command options\n"
       << "\nCommands:\n"
       << "  pair   : check pairs that covers CNV\n"
       << "  rd     : check read depths outside and inside of CNV region\n"
       << "  het    : check zygosity of sites in CNV region\n"
       << "  hetsite: get heterozygous sites from a bam file\n"
       << "  noseq  : check if CNV overlap with no seq(N) region\n"
       << "  q0     : check proportion of reads with mapq==0\n"
       << "  sub    : subset CNVs within regions\n"
       << "  tab    : make a table of CNVs from multiple samples\n"
       << "  vcf    : modify vcf file according to CNVs\n"
       << "  fa     : convert CNV MERGE sequence to fastq file for mapping\n"
       << "\nNote: samtools and bcftools must be installed and in $PATH\n"
       << endl;
  return 0;
}


int main(int argc, char* argv[])
{
  time_t begin_T,end_T;
  int pid = getpid();
  
  time(&begin_T);
  
  if ( argc<2 ) exit( usage_main(argc,argv) );
  string func=argv[1];
  
  if ( func=="noseq" ) check_noseq_regions(argc, argv);
  else if ( func=="pair" ) check_pairs(argc, argv);
  else if ( func=="rd" ) check_readdepth(argc, argv);
  else if ( func=="q0" ) check_mapq0(argc, argv);
  else if ( func=="het" ) check_het(argc, argv);
  else if ( func=="hetsite" ) get_het_site(argc, argv);
  else if ( func=="sub" ) check_sample_subset(argc, argv);
  else if ( func=="tab" ) check_sample_table(argc, argv);
  else if ( func=="vcf" ) patch_vcf(argc, argv);
  else if ( func=="fa" ) variation_to_fastq(argc, argv);
  else usage_main(argc,argv);
  
  cerr << procpidstatus(pid,"VmPeak") ;
  time(&end_T);
  cerr << "Time elapsed: " << difftime(end_T,begin_T) << " seconds\n";
  exit(0);
} 

