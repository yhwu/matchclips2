#ifndef _TAGCNV_H
#define _TAGCNV_H

struct cnv_st {
  string RNAME;
  char T;
  int P1;
  int P2;
  int un;
  int64_t GP1;
  int64_t GP2;
  int isp5;       //0,1,9  0:notmatched, 1:matched
  int isp3;       //0,1,9  9:match continuous single locus
  string ID;
  cnv_st():RNAME("CHROM"),
           T('V'),
	   P1(0),
	   P2(0),
	   GP1(0),
	   GP2(0),
	   isp5(0),
	   isp3(0),
	   ID("ABCDEFG") {};
} ;

int check_noseq_regions(int argc, char* argv[]);
int check_readdepth(int argc, char* argv[]);
int check_pairs(int argc, char* argv[]);
int check_mapq0(int argc, char* argv[]);
int check_het(int argc, char* argv[]);
int get_het_site(int argc, char* argv[]);
int check_sample_table(int argc, char* argv[]);
int check_sample_subset(int argc, char* argv[]);
int patch_vcf(int argc, char* argv[]);
int variation_to_fastq(int argc, char* argv[]);

#endif
