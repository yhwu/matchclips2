#ifndef _SAMFUNCTIONS_H
#define _SAMFUNCTIONS_H

/**** samtools headers ****/
using namespace std;
#include <string>
#include <vector>
#include <bam.h>
#include <sam.h>

#define POS_BAM_CINS -1
#define POS_BAM_CPAD -2

struct POSCIGAR_st {     
  int tid;
  int pos;              // 1-based position in chromosome
  int base;
  int qual;
  int anchor;              // 1-based position in chromosome
  int iclip;
  int l_qseq;              // length of qseq
  vector<int> op;
  vector<int> nop;
  vector<int> cop;
  vector<int> qop;
  POSCIGAR_st(): tid(-1),
		 pos(0),
		 base(-1),
		 qual(0),
		 anchor(0),
		 iclip(-1),
		 l_qseq(0),
		 op(0),
		 nop(0),
		 cop(0),
		 qop(0) {};
};

/*! 
  @abstract a wrapper for bam api
  
  @field  ref   reference id in target 0-based
  @field  beg   position on reference 0-based
  @field  end   position on reference 1-based
  @field  bp1   position of CNV, 1-based
  @field  bp2   position of CNV, 1-based, bp1<bp2
  @field  T     type of CNV
  @field  c0    number of reads with mapping quality <=1
  @field  c1    number of reads with mapping quality >1
  @field  pair  number of reads whose pairs cover bp1 and bp2
                for deletion, pairs cover [bp1,bp2]
		for duplication, pairs with negative insert size
  @field  in    bam input
  @field  bamidx bam index
  @field  buf   buffer
  @field  n     vector of depths
*/
typedef struct {  
  int ref, beg, end;  // for bam api, ref and bed are 0-based, end is 1-based
  int bp1, bp2;       // 1-based break point
  char T;
  size_t c0, c1;  
  size_t pair;  
  samfile_t *in;  
  bam_index_t *bamidx;
  bam_plbuf_t *buf;  
  std::vector<int> n;
} bam_pileup_api_wrapper_t;  

bam1_t *bam_shorten(const bam1_t *src);
bam1_t *bam_shorten2(const bam1_t *src);

void resolve_cigar_pos(const bam1_t *b,  POSCIGAR_st& m, int base);
void resolve_cigar_pos(const bam1_t *b,  POSCIGAR_st& m);
void resolve_cigar_pos(int POS, string& CIGAR, POSCIGAR_st& m);

int calibrate_resolved_cigar_pos(string& FASTA, string& SEQ, POSCIGAR_st& m);
int calibrate_resolved_cigar_pos(string& FASTA, const bam1_t *b, POSCIGAR_st& m);

int calibrate_cigar_pos(const string& FASTA, bam1_t *b);

void get_cigar(const bam1_t *b,  string& cigar);  
string get_cigar(const bam1_t *b);
string get_cigar(const POSCIGAR_st& m)  ;

void get_qseq(const bam1_t *b,  std::string& seq);
string get_qseq(const bam1_t *b);

void get_rname(const bam_header_t *header, const bam1_t *b,  std::string& rname);  

void expand_cigar(const POSCIGAR_st& m, vector<int>& e_cigar );
void expand_cigar(const bam1_t *b, vector<int>& e_cigar );

void expand_pos(const bam1_t *b, vector<int>& e_pos);
void expand_pos(const POSCIGAR_st& m, vector<int>& e_pos);

string ref_projected_onto_qseq(const bam1_t *b, const string& FASTA);

int get_pos_for_base(const POSCIGAR_st& m, int p);
int get_pos_for_base(const bam1_t *b , int p);

void bampileup(bam_pileup_api_wrapper_t& tmp, string RNAME, int p1, int p2);

#endif
