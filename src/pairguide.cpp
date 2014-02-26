#include <pthread.h>
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
#include "exhaustive.h"
#include "pairguide.h"

void check_read_pair_ends(const bam1_t *b )
{
  if ( !(b->core.flag & BAM_FPAIRED) ) return; 
  if ( b->core.mtid != b->core.tid )  return;
  
  int ISIZE;
  int F1=-1,F2=-1; // begin and end for Forward read
  int R1=-1,R2=-1; // begin and end for Reverse Complement read
  bool check=false;
  if ( bool(b->core.flag & BAM_FREVERSE)!=bool(b->core.flag & BAM_FMREVERSE) ){ 
    // either F or RC read in FR pair-end orientation
    // illumina pair-end has FR orientation
    // illumina mate pair has RF orientation, not considered yet
    // isize =   R2-F1+1  for F read
    // isize = -(R2-F1+1) for R read
    if ( (b->core.flag & BAM_FREVERSE) == 0 ) {
      // read is Forward 
      // all 4 positions are acurate 
      F1=b->core.pos;
      F2=bam_calend(&b->core, bam1_cigar(b));
      R1=b->core.mpos;
      R2=F1+b->core.isize-(b->core.isize>0 ? 1:-1);
      ISIZE = b->core.mpos+b->core.l_qseq-1-F1;
      ISIZE += (ISIZE>0 ? 1:-1);
    }
    else { 
      // read is Reverse Complement 
      // F2 is not accurate, approximate within the length of the read
      R1=b->core.pos;
      R2=bam_calend(&b->core, bam1_cigar(b));
      F1=b->core.mpos;
      F2=b->core.mpos+b->core.l_qseq-1;  //!!! approximate
      ISIZE=-(R2-F1+ (R2>F1 ? 1:-1) );
    }
  }
  else {
    // either read in FF or RR pair-end orientation
    // 454 mate pair has FF orientation
    // SOLID mat pair has RR orientation
    // don't know how isize is calculated yet, but certainly not F1 to R2
    if ( (b->core.flag & BAM_FREVERSE) == 0 ) {
      // read is Forward 
      F1=b->core.pos;
      F2=bam_calend(&b->core, bam1_cigar(b));
      R1=b->core.mpos;
      R2=b->core.mpos+b->core.l_qseq-1;
      //R2=F1+b->core.isize-(b->core.isize>0 ? 1:-1);
      ISIZE=R1-F1+(R1>F1 ? 1:-1);
    }
    else { 
      // read is Reverse Complement 
      R1=b->core.pos;
      R2=bam_calend(&b->core, bam1_cigar(b));
      F1=b->core.mpos;
      F2=b->core.mpos+b->core.l_qseq-1;  // approximate
      ISIZE=-(R1-F1+ (R1>F1 ? 1:-1) );
    }
    
  }
  
  if ( abs(ISIZE - b->core.isize) > 10 ) check=true;
  if ( check  ) {
    string cigar;
    get_cigar(b,  cigar);  
    
    cerr << "ends error "
	 << "FR"[bool(b->core.flag & BAM_FREVERSE)] 
	 << "FR"[bool(b->core.flag & BAM_FMREVERSE)] << "\t" 
	 << b->core.pos << "\t" 
	 << b->core.qual << "\t" 
	 << cigar << "\t" 
	 << b->core.mpos << "\t" 
	 << b->core.isize << "\t" 
	 << ISIZE << "\t" 
	 << ISIZE - b->core.isize << endl;
  }
  
}

//! report right forward and left reverse positions of a pair
//! if read is forward, both are accurate
//! if read is reverse complement, only R1 is accurate
//! 0 based.
void check_inner_pair_ends(const bam1_t *b, 
			   int& F2, bool& F2_acurate, int& R1, bool& R1_acurate)
{
  //*! BAM_FPROPER_PAIR may flag long/negative inserts as not properly paired
  F2=R1=-1;
  F2_acurate=R1_acurate=false;
  if ( !(b->core.flag & BAM_FPAIRED) ) return; 
  if ( b->core.mtid != b->core.tid )  return;
  if ( bool(b->core.flag & BAM_FREVERSE)==bool(b->core.flag & BAM_FMREVERSE) ) 
    return;
  
  if ( (b->core.flag & BAM_FREVERSE) ==0 ) {
    // read is Forward 
    F2= bam_calend(&b->core, bam1_cigar(b));
    F2_acurate=true;
    R1= b->core.mpos;
    R1_acurate=true;
  }
  else {
    // read is Reverse
    F2= b->core.mpos+b->core.l_qseq-1; // only approximate
    F2_acurate=false;
    R1= b->core.pos;
    R1_acurate=true;
  }
  
  //  check_read_pair_ends(b );
}

void check_outer_pair_ends(const bam1_t *b, 
			   int& F1, bool& F1_acurate, int& R2, bool& R2_acurate)
{
  //*! BAM_FPROPER_PAIR may flag long/negative inserts as not properly paired
  F1=R2=-1;
  F1_acurate=R2_acurate=false;
  if ( !(b->core.flag & BAM_FPAIRED) ) return; 
  if ( b->core.mtid != b->core.tid )  return;
  if ( bool(b->core.flag & BAM_FREVERSE)==bool(b->core.flag & BAM_FMREVERSE) ) 
    return;
  
  if ( (b->core.flag & BAM_FREVERSE) ==0 ) {
    // read is Forward 
    F1= b->core.pos;
    F1_acurate=true;
    R2=F1+b->core.isize-(b->core.isize>0 ? 1:-1);
    R2_acurate=true;
  }
  else {
    // read is Reverse
    F1=b->core.mpos;
    F1_acurate=true;
    R2=bam_calend(&b->core, bam1_cigar(b));
    R2_acurate=true;
  }
  return;
}

//! either F or RC read in FR pair-end orientation
//! illumina pair-end has FR orientation
//! illumina mate pair has RF orientation, not considered yet
//! isize =   R2-F1+1  for F read
//! isize = -(R2-F1+1) for R read
int get_isize(const bam1_t *b )
{
  if ( !(b->core.flag & BAM_FPAIRED) ) return -1; 
  if ( b->core.mtid != b->core.tid )  return -1;
  return  b->core.flag & BAM_FREVERSE ? -b->core.isize : b->core.isize;
}

int check_normalpairs_cross_pos(int ref, int end) 
{
  if ( ! msc::bam_is_paired ) return -1;

  bam1_t *b=NULL;   
  b = bam_init1();
  bam_iter_t iter=0;
  
  int dx = msc::bam_pe_insert+5*msc::bam_pe_insert_sd;
  
  //  int ref=msc::bam_ref;
  int beg=end-dx;
  if (beg<1) beg=1;
  
  iter = bam_iter_query(msc::bamidx, ref, beg, end);
  int d1=0;
  while( bam_iter_read(msc::fp_in->x.bam, iter, b)>0 ) {
    if ( b->core.tid!=ref || b->core.pos>end ) break;
    if ( ! is_read_count_for_pair(b) ) continue;
    // since we checked in [end-dx, end], we need FORWARD
    if (  b->core.flag & BAM_FREVERSE ) continue;
    if ( ! b->core.flag & BAM_FPROPER_PAIR ) continue;
    if ( abs(b->core.isize)>msc::bam_pe_insert+5*msc::bam_pe_insert_sd ||
	 abs(b->core.isize)<msc::bam_pe_insert-5*msc::bam_pe_insert_sd ) continue;
    int F1,R2; bool F1a,R2a;
    check_outer_pair_ends(b, F1, F1a, R2, R2a); 
    if ( F1<end && R2>end ) ++d1;
  }
  
  if ( d1<10 ) {
    int d2=0;
    iter = bam_iter_query(msc::bamidx, ref, end, end+dx);
    while( bam_iter_read(msc::fp_in->x.bam, iter, b)>0 ) {
      if ( b->core.tid!=ref || b->core.pos>end ) break;
      if ( ! is_read_count_for_pair(b) ) continue;
      // since we checked in [end, end+dx], we need REVERSE
      if ( ! ( b->core.flag & BAM_FREVERSE ) ) continue;
      if ( ! b->core.flag & BAM_FPROPER_PAIR ) continue;
      if ( abs(b->core.isize)>msc::bam_pe_insert+5*msc::bam_pe_insert_sd ||
	   abs(b->core.isize)<msc::bam_pe_insert-5*msc::bam_pe_insert_sd ) continue;
      int F1,R2; bool F1a,R2a;
      check_outer_pair_ends(b, F1, F1a, R2, R2a); 
      if ( F1<end && R2>end ) ++d2;
    }
    d1=max(d1, d2);
  }
  
  if ( b ) bam_destroy1(b);
  if ( iter ) bam_iter_destroy(iter);
  
  return d1;
}

//! note F2 and R1 are directional, reads are checked dx before F2 and after R1
//! regardless F2 and R1's order
//! pairs with BAM_FPROPER_PAIR bit are ignored
void check_normal_and_abnormalpairs_cross_region(int ref, int F2, int R1,
						 int& p_F2, int& p_R1, int& p_F2R1) 
{
  p_F2=p_R1=p_F2R1=0;
  if ( ! msc::bam_is_paired ) return;
  
  bam1_t *b=NULL;   
  b = bam_init1();
  bam_iter_t iter=0;
  int beg,end;
  
  int dx = msc::bam_pe_insert+7*msc::bam_pe_insert_sd;
  
  beg=max(1, F2-dx);
  end=F2;
  iter = bam_iter_query(msc::bamidx, ref, beg, end);
  while( bam_iter_read(msc::fp_in->x.bam, iter, b)>0 ) {
    if ( b->core.tid!=ref ) break;
    if ( !is_read_count_for_pair(b) ) continue; 
    
    int F1,R2; bool F1a,R2a;
    check_outer_pair_ends(b, F1, F1a, R2, R2a); 
    int isize=R2-F1;
    if ( !(b->core.flag & BAM_FREVERSE) &&  // since we checked -dx, we need FORWARD
	 isize > msc::bam_pe_insert-5*msc::bam_pe_insert_sd && 
	 isize < msc::bam_pe_insert+5*msc::bam_pe_insert_sd ) {
      // pair is normal
      p_F2 += ( F1<F2 && R2>F2 );
      continue;
    }
    
    if ( F1<=F2 && R2>=R1 ) {
      int bisize=F2-F1+R2-R1;
      if ( bisize>0 &&
	   bisize>msc::bam_pe_insert-7*msc::bam_pe_insert_sd &&
	   bisize<msc::bam_pe_insert+7*msc::bam_pe_insert_sd ) ++p_F2R1;
    }
  }
  
  beg=R1;
  end=R1+dx;
  iter = bam_iter_query(msc::bamidx, ref, beg, end);
  while( bam_iter_read(msc::fp_in->x.bam, iter, b)>0 ) {
    if ( b->core.tid!=ref ) break;
    if ( !is_read_count_for_pair(b) ) continue; 
    
    int F1,R2; bool F1a,R2a;
    check_outer_pair_ends(b, F1, F1a, R2, R2a); 
    int isize=R2-F1;
    if ( b->core.flag & BAM_FREVERSE && // since we checked +dx, we need REVERSE
	 isize > msc::bam_pe_insert-5*msc::bam_pe_insert_sd && 
	 isize < msc::bam_pe_insert+5*msc::bam_pe_insert_sd ) {
      // pair is normal
      p_R1 += ( F1<R1 && R2>R1 );
      continue;
    }
    
    if ( F1<=F2 && R2>=R1 ) {
      int bisize=F2-F1+R2-R1;
      if ( bisize>0 &&
	   bisize>msc::bam_pe_insert-7*msc::bam_pe_insert_sd &&
	   bisize<msc::bam_pe_insert+7*msc::bam_pe_insert_sd ) ++p_F2R1;
    }
    
  }
  // since we checked both -dx and +dx, both segments were counted
  p_F2R1/=2;
  
  if ( p_F2 < 10 ) p_F2=check_normalpairs_cross_pos(ref, F2); 
  if ( p_R1 < 10 ) p_R1=check_normalpairs_cross_pos(ref, R1); 
  
  bam_destroy1(b);
  bam_iter_destroy(iter);
  
  return;
}

//! check read depth at pair ends
//! check normal pairs across pair ends
void stat_and_filter_pair_group(vector<intpair_st>& bpi, vector<pairinfo_st>& bpinfo) 
{
  bpinfo.clear();
  if ( bpi.size()<1 ) return;
  
  vector<pairinfo_st> bp(0);
  for(size_t i=0; i<bpi.size(); ++i) {
    if ( bpi[i].pair_count<=5 ) continue;
    pairinfo_st ibpinfo;
    ibpinfo.F2=bpi[i].F2;
    ibpinfo.F2_acurate=bpi[i].F2_acurate;
    ibpinfo.R1=bpi[i].R1;
    ibpinfo.R1_acurate=bpi[i].R1_acurate;
    ibpinfo.pair_count=bpi[i].pair_count;
    
    bp.push_back(ibpinfo);
  }
  
  for(size_t i=0; i<bp.size(); ++i) {
    bp[i].F2_rp=check_normalpairs_cross_pos(msc::bam_ref, bp[i].F2); 
    bp[i].R1_rp=check_normalpairs_cross_pos(msc::bam_ref, bp[i].R1); 
  }
  
  // roughly approximate read pair, read depth, min read pairs needed to continue
  vector<int> t1(0);
  for(size_t i=0; i<bp.size(); ++i) {
    if ( bp[i].F2_rp>9 ) t1.push_back(bp[i].F2_rp);
    if ( bp[i].R1_rp>9 ) t1.push_back(bp[i].R1_rp);
  }
  std::nth_element(t1.begin(), t1.begin() + t1.size()/2, t1.end());  
  int medRP= t1.size()>0 ? t1[t1.size()/2] : 0 ;
  //int medRD=medRP*msc::bam_l_qseq*2/msc::bam_pe_insert;
  int minRP=medRP/16;
  if ( minRP<5 ) minRP=5;
  // 1. /2 for one chromosome, 
  // 2. /2 for variation, 
  // 3. /2 for mapping error
  // 4. /2 for mapq=0 reads
  
  for(size_t i=0; i<bp.size(); ++i) {
    
    if ( bp[i].pair_count<=minRP ) continue;
    
    int dx = msc::bam_l_qseq*5;
    if ( dx/3 > abs(bp[i].F2-bp[i].R1) ) dx=abs(bp[i].F2-bp[i].R1)*3;
    if ( dx < msc::bam_l_qseq ) dx = msc::bam_l_qseq;
    
    check_cnv_readdepth(msc::bam_ref, bp[i].F2, bp[i].R1, dx, 
			bp[i].F2_rd, bp[i].R1_rd, bp[i].rd);
    
    assess_rd_rp_sr_infomation(bp[i]);

  }
  
  
  for (size_t i=0; i<bp.size(); ++i) 
    if ( bp[i].rpscore>0 || bp[i].rdscore>0 ) bpinfo.push_back(bp[i]);
  
  return;
}

//! cluster the pairs according to F2 position and pair length
//! get the most possible brreak points
//! for F2, get the right most position
//! for R1, get the left most position
void kmean_pairgroup(vector<intpair_st>& pairs, vector<intpair_st>& bp)
{
  bp.clear();
  if ( pairs.size()<1 ) return;
  
  int pair_gap=msc::bam_pe_insert+msc::bam_pe_insert_sd*5;
  
  vector<int> len(pairs.size());
  sort(pairs.begin(), pairs.end(), sort_pair_len);
  for(size_t i=0; i<pairs.size(); ++i) len[i]=pairs[i].R1-pairs[i].F2;
  //  sort(len.begin(),len.end());
  
  // cluster according to pair length and remove outliers
  vector<int> g(0);       // cluster value
  vector<int> gc(0);      // cluster count
  vector<int> gi(len.size(),0); // cluster assignment
  g.push_back(len[0]);
  for(size_t i=0; i<len.size(); ++i) {
    if ( len[i]-g.back() <= pair_gap ) gi[i]=g.size()-1;
    else {
      g.push_back( len[i] );
      gi[i]=g.size()-1;
    }
  }
  
  // a cluster should have at least six pairs msc::minClusterSize=6
  // roughly 3 templates ( 6 segments ) 
  gc.resize(g.size());
  for(size_t i=0; i<len.size(); ++i) ++gc[ gi[i] ];
  if ( msc::verbose>1 ) 
    cerr << "density cluster\t" << g.size() << endl;
  for(int j=0; j<(int)g.size(); ++j) {
    bool del=false;
    if ( (double)gc[j]/(double)len.size() < 0.1 
	 || gc[j]<msc::minClusterSize ) del=true;
    if ( gc[j] > msc::minClusterSize ) del=false;
    if ( msc::verbose>1 ) { 
      cerr << "cluster " << j << "\t" 
	   << g[j] << "\t" 
	   << gc[j] << "\t" 
	   << (float)gc[j]/len.size();
      if ( del ) cerr << "\tx";
      cerr << endl;
    }
    if ( del ) {
      for(size_t i=0;i<gi.size();++i) if ( gi[i]==j ) gi[i]=-1;
    }
  }
  
  int good_pair_count=0;
  for(size_t i=0; i<pairs.size(); ++i) if ( gi[i]!=-1 ) ++good_pair_count;
  if ( msc::verbose>1 )     
    cerr << "cluster purified\t" << len.size() << "\t" << good_pair_count << endl;
  
  // get break points from the clusters
  for(int j=0; j<(int)g.size(); ++j) {
    vector<intpair_st> tmp(0);
    for(size_t i=0; i<pairs.size(); ++i) 
      if ( gi[i]==j ) tmp.push_back( pairs[i] );
    if ( tmp.size()==0 ) continue;
    intpair_st ibp=tmp[0];
    for(size_t i=0; i<tmp.size(); ++i) {
      // infer F2
      if ( tmp[i].F2_acurate && (!ibp.F2_acurate) ) {
	ibp.F2=tmp[i].F2;
	ibp.F2_acurate=tmp[i].F2_acurate;
      }
      if ( tmp[i].F2_acurate && tmp[i].F2>=ibp.F2 ) {
	ibp.F2=tmp[i].F2;
	ibp.F2_acurate=tmp[i].F2_acurate;
      }
      // note R1 must be always acurate
      if ( tmp[i].R1_acurate && tmp[i].R1<=ibp.R1 ) {
	ibp.R1=tmp[i].R1;
	ibp.R1_acurate=tmp[i].R1_acurate;
      }
    }
    ibp.pair_count=tmp.size();
    bp.push_back(ibp);
  }
  
  // make sure all pairs counted
  for(int i=0; i<(int)bp.size(); ++i) {
    good_pair_count-=bp[i].pair_count;
    if ( msc::verbose>1 )     
      cerr << "cluster identified " << i << "\t"
	   << bp[i].F2_acurate << "\t" << bp[i].F2 << "\t"
	   << bp[i].R1_acurate << "\t" << bp[i].R1 << "\t"
	   << bp[i].pair_count << "\n" << endl;
  }
  if ( good_pair_count != 0 ) 
    cerr << "not all pairs counted from the clusters" << endl;
  
}

//! sort pairs according to F2
//! roughly break them to clusters according to gap between pairs
//! call kmean to do fine clustering
void check_pair_group(vector<intpair_st>& pairs, vector<pairinfo_st>& bpinfo)
{
  bpinfo.clear();
  vector<intpair_st> bp(0);
  if ( pairs.size()<1 ) return;
  if ( msc::verbose>0 ) cerr << "total pairs:" << pairs.size() << endl;
  
  int pair_gap=msc::bam_pe_insert+msc::bam_pe_insert_sd*5-msc::bam_l_qseq/2;
  
  vector<intpair_st> pairgroup(0);    
  sort(pairs.begin(), pairs.end(), sort_pair);
  
  vector<intpair_st> bp0(0);
  pairgroup.push_back( pairs[0] );
  for(size_t i=1; i<pairs.size(); ++i) {
    if ( pairs[i].F2 - pairs[i-1].F2 > pair_gap ) {
      // cerr << "------------\t" << pairgroup.size() << endl;
      if ( (int)pairgroup.size()<msc::minClusterSize ) pairgroup.clear();
      kmean_pairgroup(pairgroup, bp0);
      pairgroup.clear();
      if ( bp0.size()>0 ) bp.insert(bp.end(), bp0.begin(), bp0.end());
    }
    
    if ( msc::verbose>2 ) {
      cerr << pairs[i].F2_acurate << "\t" << pairs[i].F2 << "\t" 
	   << pairs[i].R1_acurate << "\t" << pairs[i].R1 << "\t" 
	   << pairs[i].F2-pairs[i].R1 << endl;
    }
    pairgroup.push_back( pairs[i] );
  }
  if ( pairgroup.size()<3 ) pairgroup.clear();
  kmean_pairgroup(pairgroup, bp0);
  if ( bp0.size()>0 ) bp.insert(bp.end(), bp0.begin(), bp0.end());
  pairgroup.clear();
  
  if ( msc::verbose>0 ) {
    cerr << "clusters found --- " << bp.size() << endl;
  }
  
  // the above code counts reads
  // convert pair reads to pair templates
  for(size_t i=0; i<bp.size(); ++i) bp[i].pair_count/=2;

  // check read depth and filter
  for(size_t i=0; i<bp.size(); ++i) {
    
    // if ( bp[i].pair_count<=minRP ) continue;
    if ( bp[i].pair_count<=4 ) continue;
    
    pairinfo_st ibp;
    ibp.F2=bp[i].F2;
    ibp.F2_acurate=bp[i].F2_acurate;
    ibp.R1=bp[i].R1;
    ibp.R1_acurate=bp[i].R1_acurate;
    ibp.pair_count=bp[i].pair_count;
    
    ibp.F2_rp=check_normalpairs_cross_pos(msc::bam_ref, ibp.F2); 
    ibp.R1_rp=check_normalpairs_cross_pos(msc::bam_ref, ibp.R1); 
    
    if ( ibp.pair_count<=min(ibp.F2_rp, ibp.R1_rp)/16 ) continue;
    
    int dx = msc::bam_l_qseq*5;
    if ( dx/3 > abs(bp[i].F2-bp[i].R1) ) dx=abs(bp[i].F2-bp[i].R1)*3;
    if ( dx < msc::bam_l_qseq ) dx = msc::bam_l_qseq;
    
    check_cnv_readdepth(msc::bam_ref, ibp.F2, ibp.R1, dx, 
			ibp.F2_rd, ibp.R1_rd, ibp.rd);
    
    assess_rd_rp_sr_infomation(ibp);
    
    if ( ibp.rpscore>0 || ibp.rdscore>0 ) bpinfo.push_back(ibp);
  }
  
  
  // stat_and_filter_pair_group(bp, bpinfo);
  return;
}  


//! find possible match between reads around F2 and R1
//! get the best possible break points
//! reads are not filtered based on CIGAR score
//! reads around F2 are collected
//! reads around R1 are collected
//! check if two reads from each group match
//! get the break points
void match_reads_for_pairs(pairinfo_st& ipairbp, string& FASTA, int dx, bool pointmode)
{
  // these are return values
  ipairbp.MS_F2=-1;
  ipairbp.MS_R1=-1;
  ipairbp.MS_ED=1000;
  ipairbp.MS_ED_count=-1;
  
  bool with_matching_reads=false;

  bam1_t *b=NULL; b = bam_init1();
  bam_iter_t iter=NULL;
  
  vector<uint8_t> bdata(0); bdata.reserve(3000000);
  vector<bam1_t> b_F2(0);
  vector<bam1_t> b_R1(0);
  bam1_t ibam;
  
  if ( dx<1 ) dx=msc::bam_l_qseq/4;
  int beg,end, MS_F2_rd=0, MS_R1_rd=0;
  
  // get reads around F2
  beg=max(1, ipairbp.F2-dx);
  end=ipairbp.F2+dx;
  iter = bam_iter_query(msc::bamidx, msc::bam_ref, beg, end);
  while( bam_iter_read(msc::fp_in->x.bam, iter, b)>0 ) {
    if ( b->core.tid != msc::bam_ref ) break;
    if ( b->core.pos > end ) break;
    if ( ! is_read_count_for_depth(b) ) continue;
    //if ( pointmode && 
    //	 b->core.pos+msc::bam_l_qseq/5 > ipairbp.F2 ) break;
    POSCIGAR_st bm; resolve_cigar_pos(b, bm, 0);  
    if ( bm.pos<=0 || bm.op.size()<1 ) continue;
    bool cover= bm.cop[0] <= ipairbp.F2 && 
      bm.cop.back()+bm.nop.back()*(bm.op.back()!=BAM_CHARD_CLIP) >= ipairbp.F2; 
    MS_F2_rd+=cover;
    
    if ( pointmode && (! cover)  ) continue;
    
    _save_read_in_vector(b, ibam, bdata);
    b_F2.push_back(ibam);
  }
  ipairbp.MS_F2_rd=MS_F2_rd;
  
  // get reads around R1
  beg=max(1, ipairbp.R1-dx);
  end=ipairbp.R1+dx;
  iter = bam_iter_query(msc::bamidx, msc::bam_ref, beg, end);
  while( bam_iter_read(msc::fp_in->x.bam, iter, b)>0 ) {
    if ( b->core.tid != msc::bam_ref ) break;
    if ( b->core.pos > end ) break;
    if ( ! is_read_count_for_depth(b) ) continue;
    //if ( pointmode && 
    //	 b->core.pos+msc::bam_l_qseq/5*4 < ipairbp.R1 ) continue;   
    POSCIGAR_st bm; resolve_cigar_pos(b, bm, 0);  
    if ( bm.pos<=0 || bm.op.size()<1 ) continue;
    bool cover= bm.cop[0] <= ipairbp.R1 && 
      bm.cop.back()+bm.nop.back()*(bm.op.back()!=BAM_CHARD_CLIP) >= ipairbp.R1; 
    MS_R1_rd+=cover;
    
    if ( pointmode && ( ! cover ) ) continue;
    
    _save_read_in_vector(b, ibam, bdata);
    b_R1.push_back(ibam);
  }
  ipairbp.MS_R1_rd=MS_R1_rd;
  
  if ( b ) bam_destroy1(b);
  if ( iter ) bam_iter_destroy(iter);
  
  for(size_t i=0; i<b_F2.size(); ++i) b_F2[i].data = &bdata[(size_t)b_F2[i].data];
  for(size_t i=0; i<b_R1.size(); ++i) b_R1[i].data = &bdata[(size_t)b_R1[i].data];
  
  if ( b_F2.size()<1 || b_R1.size()<1 ) with_matching_reads=false;
  
  vector<ED_st> bp(0); 
  multithreads_read_matching(b_F2, b_R1, FASTA, -1, bp);
  
  int min_ED=msc::bam_l_qseq*2;
  for(size_t i=0; i<bp.size(); ++i) {
    if ( abs(bp[i].F2-bp[i].R1)<10 ) {
      bp[i].F2=-1;
      bp[i].R1=-1;
      bp[i].ED=msc::bam_l_qseq*2+1;
    }
    if ( bp[i].ED<min_ED ) min_ED=bp[i].ED;
  }
  if ( min_ED < msc::bam_l_qseq/2 ) with_matching_reads=true;
  
  vector<ED_st> best_ED(0); 
  reduce_matched_break_points(bp, best_ED);
  
  int it=0;
  if ( best_ED.size()>0 ) {
    for(int i=0; i<(int)best_ED.size(); ++i) {
      if ( best_ED[i].ED > msc::bam_l_qseq/2 ) continue;
      if ( best_ED[i].count > best_ED[it].count ) it=i;
      /*
      cerr << best_ED[i].F2 << "\t" 
	   << best_ED[i].R1 << "\t" 
	   << best_ED[i].ED << "\t"
	   << best_ED[i].count
	   << endl;
      */
    }
    int it_max_count=it;
    for(int i=0; i<(int)best_ED.size(); ++i) {
      if ( best_ED[i].ED > msc::bam_l_qseq/2 ) continue;
      if ( abs(best_ED[i].F2-ipairbp.F2) + abs(best_ED[i].R1-ipairbp.R1) <
	   abs(best_ED[it].F2-ipairbp.F2) + abs(best_ED[it].R1-ipairbp.R1) ) it=i;
    }
    int it_closest=it;
    
    int dis_closest=abs(best_ED[it_closest].F2-ipairbp.F2) 
      + abs(best_ED[it_closest].R1-ipairbp.R1);
    int dis_mac_count=abs(best_ED[it_max_count].F2-ipairbp.F2) 
      + abs(best_ED[it_max_count].R1-ipairbp.R1);
    
    if ( dis_closest*2 < dis_mac_count ) {
      if ( best_ED[it_closest].count*2 > best_ED[it_max_count].count &&
	   best_ED[it_closest].ED/2 <= best_ED[it_max_count].count ) it=it_closest;
    }
    else { it=it_max_count; }
    // cerr << "Final "
    //	 << best_ED[it].F2 << "\t" 
    //	 << best_ED[it].R1 << "\t" 
    //	 << best_ED[it].ED << "\t"
    //	 << best_ED[it].count
    //	 << endl;
  }
  
  vector<int> b_MS_m(b_F2.size(), 0);
  vector<int> b_SM_m(b_R1.size(), 0);
  if ( best_ED.size()>0 ) {
    for(size_t i=0; i<bp.size(); ++i) {
      if ( abs(bp[i].F2-best_ED[it].F2)<10 &&
	   abs(bp[i].R1-best_ED[it].R1)<10 &&
	   bp[i].ED <= best_ED[it].ED+msc::bam_l_qseq/15 &&
	   bp[i].ED <= msc::bam_l_qseq/8 ) {
	b_MS_m[ bp[i].iL ]=1;
	b_SM_m[ bp[i].iR ]=1;
      }
    }
  }
  ipairbp.F2_sr=0;
  for(size_t k=0; k<b_MS_m.size(); ++k) ipairbp.F2_sr += (b_MS_m[k]>0) ;
  ipairbp.R1_sr=0;
  for(size_t k=0; k<b_SM_m.size(); ++k) ipairbp.R1_sr += (b_SM_m[k]>0) ;
  // count for skipping reads
  ipairbp.F2_sr=(double)ipairbp.F2_sr*ED_st::linc;
  ipairbp.R1_sr=(double)ipairbp.R1_sr*ED_st::rinc;
  
  if ( ipairbp.F2_sr==0 || ipairbp.R1_sr==0 ) with_matching_reads=false;
  
  if ( with_matching_reads && best_ED.size()>0 ) {
    ipairbp.MS_F2=best_ED[it].F2;
    ipairbp.MS_R1=best_ED[it].R1;
    ipairbp.MS_ED=best_ED[it].ED;
    ipairbp.MS_ED_count=best_ED[it].count;
  }
  
  if ( msc::verbose ) {
    cerr << "good rdpr\t"
	 << ipairbp.F2_acurate << " " << ipairbp.F2 << "\t"
	 << ipairbp.R1_acurate << " " << ipairbp.R1 << "\t"
	 << ipairbp.R1-ipairbp.F2 << "\t"
	 << "d " << ipairbp.F2_rd << " " << ipairbp.R1_rd << " " << ipairbp.rd << " " << ipairbp.rdscore << "\t"
	 << "p " << ipairbp.F2_rp << " " << ipairbp.R1_rp << " " << ipairbp.pair_count << " " << ipairbp.rpscore << "\t"
	 << "m " << ipairbp.MS_F2 << " " << ipairbp.MS_R1 << " "
	 << ipairbp.F2_sr << " " << ipairbp.R1_sr << " " << ipairbp.MS_ED << " " << ipairbp.MS_ED_count
	 << endl;
  }
  
  //  vector<uint8_t> (0).swap(msc::bpdata);
  vector<uint8_t> (0).swap(bdata);
  vector<bam1_t> (0).swap(b_F2);
  vector<bam1_t> (0).swap(b_R1);
  return;
}

void match_reads_for_pairs(vector<pairinfo_st>& pairbp, string& FASTA, int dx, bool pointmode)
{
  if ( dx<1 ) dx=msc::bam_l_qseq/4;
  int pair_supported=0;
  for ( int i=0; i<(int)pairbp.size(); ++i) {
    match_reads_for_pairs(pairbp[i], FASTA, dx, pointmode);
    pair_supported += pairbp[i].sr_count>0 ;
  }
  cerr << "match found " << pair_supported << " out of " << pairbp.size() << endl;
  return;
}

void get_abnormal_pairs(int ref, int beg, int end, int min_pair_length,
			vector<intpair_st>& pairs)
{
  pairs.clear();
  if ( ! msc::bam_is_paired ) return;
  
  bam1_t *b=NULL; b = bam_init1();
  bam_iter_t iter=0;
  
  //! collect read pairs
  cerr << "processing pairs in:\t" 
       << msc::fp_in->header->target_name[ ref ]  << "\t" 
       << beg << "\t" << end << endl
       << "min_pair_length:\t" << min_pair_length 
       << endl;
  
  intpair_st ipair;
  
  iter = bam_iter_query(msc::bamidx, ref, beg, end);
  size_t count=0;
  while( bam_iter_read(msc::fp_in->x.bam, iter, b)>0 ) {
    count++;
    if ( count%1000000==0 ) {
      cerr << "#processed " << commify(count) << " reads at pos " 
	   << string(msc::fp_in->header->target_name[ref]) 
	   << "@" << commify(b->core.pos) 
	   << endl;
    }
    if ( abs(b->core.isize) < min_pair_length ) continue;
    if ( ! is_read_count_for_pair(b) ) continue;
    check_inner_pair_ends(b, ipair.F2, ipair.F2_acurate, ipair.R1, ipair.R1_acurate);
    if ( ipair.F2>0 && ipair.R1>0 ) pairs.push_back(ipair);
  }
  
  bam_destroy1(b);
  bam_iter_destroy(iter);
  return;
}

void pair_guided_search(vector<intpair_st>& pairs, string& FASTA, 
			vector<pairinfo_st>& pairbp) 
			
{
  check_pair_group(pairs, pairbp); 
  vector<intpair_st> (0).swap(pairs);
  
  for(int i=0; i<(int) pairbp.size(); ++i) pairbp[i].tid=msc::bam_ref;

  pairinfo_st ibp;
  
  size_t count=0;
  for (size_t i=0; i<pairbp.size(); ++i) {
    string cnv=cnv_format1(pairbp[i]);
    for(size_t t=0; t<cnv.size(); ++t) if ( cnv[t]=='\t' ) cnv[t]=' ';
    cerr << "checking: " << cnv << endl;
    
    int old_minOverlap = msc::minOverlap;
    msc::minOverlap += msc::minOverlapPlus;  
    match_reads_for_pairs(pairbp[i], FASTA, 0, false);
    msc::minOverlap = old_minOverlap;  
    
    assess_rd_rp_sr_infomation( pairbp[i] );
    cerr << "matched:  " << mr_format1(pairbp[i]) << endl;
    // update pos and check rd and rp info
    //if ( pairbp[i].srscore>0 && 
    //	 abs(pairbp[i].F2-pairbp[i].MS_F2)<msc::bam_l_qseq/2 && 
    //	 abs(pairbp[i].R1-pairbp[i].MS_R1)<msc::bam_l_qseq/2 ) {
    if ( pairbp[i].srscore>0 ) {
      ibp=pairbp[i];
      ibp.F2=ibp.MS_F2;
      ibp.R1=ibp.MS_R1;
      stat_region( ibp, FASTA, msc::bam_l_qseq*4 );
      assess_rd_rp_sr_infomation( ibp );
      cnv=cnv_format1(ibp);
      for(size_t t=0; t<cnv.size(); ++t) if (cnv[t]=='\t') cnv[t]=' ';
      cerr << "matched:  " << cnv << endl;
      if ( ibp.rdscore>=2 || 
	   ibp.rpscore>=2 ||
	   ibp.rpscore>=pairbp[i].rpscore ||
	   ibp.srscore>=3  ) pairbp[i]=ibp;
    }
    count += ( pairbp[i].srscore>0 )  ;
    
    cnv=cnv_format1(pairbp[i]);
    for(size_t t=0; t<cnv.size(); ++t) if (cnv[t]=='\t') cnv[t]=' ';
    cerr << "updated:  " << cnv << "\n" << endl;
  }
  cerr << "matching support " << count << " out of " << pairbp.size() << "\n"
       << "pair end mode done\n"
       << endl;;
  
  // now information is complete
  return;
}

void pair_guided_search(int ref, int beg, int end, int min_pair_length, string& FASTA,
			vector<pairinfo_st>& pairbp) 
{
  pairbp.clear();
  if ( ! msc::bam_is_paired ) return;
  
  vector<intpair_st> pairs(0);
  get_abnormal_pairs(ref, beg, end, min_pair_length, pairs);

  pair_guided_search(pairs, FASTA, pairbp);
  return;
}

