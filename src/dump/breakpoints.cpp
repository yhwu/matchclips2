#include <stdio.h>
#include <string.h>
#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <complex>
#include <algorithm>
#include <string>
#include <vector>
#include <unistd.h>
using namespace std;

/**** samtools headers ****/
#include <bam.h>
#include <sam.h>

/**** user headers ****/
#include "samfunctions.h"
#include "functions.h"
#include "matchreads.h"

#include "preprocess.h"

/*! readMS and readSM overlap, find the break points
  ! 1. for F2 read, b1, get the reference position for p1
  ! 2. for R1 read, b2, get the reference position for 0
  | if ( cigar_MS == M && cigar_SM == M ) bp_F2 bp_R1 found
*/
void_get_break_points(const string& FASTA, const bam1_t *bF2, const bam1_t *bR1, int p1, int& F2, int& R1)
{
  F2=R1=-1;
  if ( bF2->core.n_cigar==0 || bF2->core.pos<0 ) return;
  if ( bR1->core.n_cigar==0 || bR1->core.pos<0 ) return;
  
  // 0 based positions on FASTA
  POSCIGAR_st bF2_m, bR1_m;
  resolve_cigar_pos(bF2, bF2_m, 0);
  resolve_cigar_pos(bR1, bR1_m, 0);
  
  int bpF2=get_pos_for_base(bF2_m, p1-1);
  int bpR1=get_pos_for_base(bR1_m, 0);
  if ( bpF2<1 || bpF2+bF2->core.l_qseq>FASTA.size() ) return;
  if ( bpR1<1 || bpR1+bR1->core.l_qseq>FASTA.size() ) return;
  
  if ( bR1_m.op.size()==1 &&  bR1_m.op[1]==BAM_CMATCH ) {
    F2=bpF2;
    R1=bpR1;
  }
  
  
  //! break points slide from (p1-1,0) to (l_qseq-1, )
  for(int q=p1-1, r=0; q<bF2->core.l_qseq; ++q,++r) {
    bpF2=get_pos_for_base(bmF2, q);
    bpR1=get_pos_for_base(bmR1, 0);
  }
  
  
  //! calculate break point from read position
  //! 5' end position of common string on reference
  REFPOS_com_beg=p1-MSREAD[i].Mrpos+MSREAD[i].pos;
  // 3' end position of common string on reference
  REFPOS_com_end=CL-1-SMREAD[k].Mrpos+SMREAD[k].pos;
  
  //! calculate break point from mismatched position
  //! 5' end position of common string on reference
  int REFPOS_com_beg1=MSREAD[i].sbeg-CL+MSREAD[i].S;
  //! 3' end position of common string on reference
  int REFPOS_com_end1=CL-1-SMREAD[k].S+SMREAD[k].pos;
  //! the difference is due to short indel in CIGAR
  //! an alternative is to calculate based on the CIGAR strings
  //! will be implemented later
  int dxindel_beg=(int)abs(REFPOS_com_beg1-REFPOS_com_beg);
  int dxindel_end=(int)abs(REFPOS_com_end1-REFPOS_com_end);
  /*
    if ( REFPOS_com_beg1 != REFPOS_com_beg ) 
    cerr << "REFPOS_com_beg1\t" 
    << REFPOS_com_beg1 << "\t" << REFPOS_com_beg << "\t"
    << REFPOS_com_beg1-REFPOS_com_beg <<  endl;
    if ( REFPOS_com_end != REFPOS_com_end ) 
    cerr << "REFPOS_com_end1\t" 
    << REFPOS_com_end1 << "\t" << REFPOS_com_end << "\t" 
    << REFPOS_com_end1 - REFPOS_com_end << endl;
  */
  REFPOS_com_beg=REFPOS_com_beg1;
  REFPOS_com_end=REFPOS_com_end1;
  
  //! check if _beg and _end positions are accurate
  string left10=MATCHED.substr(0,msc::minSNum);
  int ndiff=0;
  for(int i1=0;i1<(int)left10.size();++i1) 
    ndiff+= ( left10[i1] != FASTA[REFPOS_com_beg-1+i1] );
  if ( ndiff > 1 ) {
    dxindel_beg=max(10, dxindel_beg);
    int minndiff=ndiff;
    int mindx=0;
    for(int dx=-dxindel_beg;dx<=dxindel_beg;++dx) {
      int ndiffdx=0;
      for(int i1=0;i1<(int)left10.size();++i1) {
	ndiffdx+= ( left10[i1] != FASTA[REFPOS_com_beg+dx-1+i1] );
	if ( ndiffdx>minndiff ) break;
      }
      if ( ndiffdx<minndiff ) {
	minndiff=ndiffdx;
	mindx=dx;
      }
      if ( minndiff <=1 ) break;
    }
    if ( minndiff < ndiff ) {
      REFPOS_com_beg+=mindx;
      //cerr << "REFPOS_com_beg += " << mindx << "\t" 
      //     << ndiff << "\t" << minndiff << endl;
    }
  }
  
  string right10=MATCHED.substr(MATCHED.size()-msc::minSNum);
  ndiff=0;
  for(int i1=right10.size()-1, k1=REFPOS_com_end-1; i1>=0;--i1,--k1) 
    ndiff+= ( right10[i1] != FASTA[k1] );
  if ( ndiff > 1 ) {
    int minndiff=ndiff;
    int mindx=0;
    dxindel_end=max(10, dxindel_end);
    for(int dx=-dxindel_end;dx<=dxindel_end;++dx) {
      int ndiffdx=0;
      for(int i1=right10.size()-1, k1=REFPOS_com_end-1+dx; i1>=0;--i1,--k1) {
	ndiffdx+= ( right10[i1] != FASTA[k1] );
	if ( ndiffdx>minndiff ) break;
      }
      if ( ndiffdx<minndiff ) {
	minndiff=ndiffdx;
	mindx=dx;
      }
      if ( minndiff <=1 ) break;
    }
    if ( minndiff < ndiff ) {
      REFPOS_com_end+=mindx;
      //cerr << "REFPOS_com_end += " << mindx << "\t" 
      //     << ndiff << "\t" << minndiff << endl;
    }
  }
      
  //! find break points by minimizing edit distance between MATCHED and ref
  int mdis=CL;
  int mdx1=CL-p2-1;
  int uncertainty=0;
  int premdis=mdis;
  int premdx1=mdx1;
  int dx1,dx2,bp1,bp2,edistance;
  // break point 1, break point 2
  for(dx1=SMREAD[k].S;dx1<=CL-MSREAD[i].S;++dx1) {
    dx2=CL-dx1;
    bp1=REFPOS_com_beg+dx1-1;
    bp2=REFPOS_com_end-dx2+1;
    MATCHEDREF=FASTA.substr(REFPOS_com_beg-1, dx1)+FASTA.substr(bp2-1, dx2);
    if ( MATCHEDREF.length() != MATCHED.length() ) {
      cerr << "wrong MATCHEDREF length" << endl;
      exit(0);
    }
    edistance=edit_distance(MATCHEDREF, MATCHED);
    if (  edistance < mdis ) {
      premdis=edistance;
      premdx1=dx1;
    }
    if (  edistance <= mdis ) {
      mdis=edistance;
      mdx1=dx1;
    }
  }
  // right position
  //dx1=mdx1;
  //dx2=CL-mdx1;
  // left position
  dx1=premdx1;
  dx2=CL-premdx1;
  bp1=REFPOS_com_beg+dx1-1;
  bp2=REFPOS_com_end-dx2+1;
  //if ( premdis==mdis ) uncertainty=mdx1-premdx1;
  //if ( CL-SMREAD[k].S-MSREAD[i].S>0 ) uncertainty=mdx1-premdx1;
  uncertainty=CL-SMREAD[k].S-MSREAD[i].S;
  if ( uncertainty > abs(bp2-bp1) ) uncertainty=mdx1-premdx1;
  edistance=mdis;
  MATCHEDREF=FASTA.substr(REFPOS_com_beg-1, dx1)+FASTA.substr(bp2-1, dx2);
      
  //! according to the equation in paper; not so accurate
  // bp1=MSREAD[i].sbeg-CL+MSREAD[i].S+SMREAD[k].S+1;
  // bp2=SMREAD[k].pos;
  // uncertainty=CL-MSREAD[i].S-SMREAD[k].S;
      
  //! extend uncertainty according to repeats on reference
  //! does not seem to help much, but helps
  if ( bp1+1!=bp2 ) {
    dx1=0;
    int i5=bp1+uncertainty;
    int i3=bp2-1+uncertainty;
    if ( i5<(int)FASTA.size() && i3 <(int)FASTA.size() ) {
      while(FASTA[i5+dx1]==FASTA[i3+dx1] && FASTA[i5+dx1]!='N' ) {
	dx1++ ;
	if ( (i5+dx1)>(int)FASTA.size()-10 || (i3+dx1)>(int)FASTA.size()-10 ) break; 
	if ( dx1>CL ) break;
      }
    }
    //if ( dx1>0 ) cerr << "3' extended by " << dx1 << endl;
	
    dx2=0;
    i5=bp1-1;
    i3=bp2-2;
    if ( i5>0 && i3>0 ) {
      while(FASTA[i5+dx2]==FASTA[i3+dx2] && FASTA[i5+dx2]!='N' ) {
	dx2--;
	if ( (i5+dx2)<1 || (i3+dx2)<1 ) break; 
	if ( abs(dx2)>CL ) break;
      }
    }
    //if ( dx2<0 ) cerr << "5' extended by " << dx2 << endl;
    bp1+=dx1;
    bp2+=dx2;
    uncertainty+=dx1-dx2;
  }
      
}
