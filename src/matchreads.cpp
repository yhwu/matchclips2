#include <pthread.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <complex>
#include <map>
#include <algorithm>
#include <string>
#include <vector>
#include <unistd.h>
using namespace std;

/**** user samtools headers ****/
#include "samfunctions.h"

/**** user headers ****/
#include "functions.h"
#include "readref.h"
#include "matchreads.h"
#include "preprocess.h"
#include "exhaustive.h"
#include "pairguide.h"
//#include "statcnv.h"


pthread_mutex_t nout;

double ED_st::linc=1;
double ED_st::rinc=1;

int msc::pid=getpid();
long msc::seed=137;
bool msc::debug=false;
int msc::start=0;
int msc::end=0;
int msc::dx=1000;
ofstream msc::fout;
string msc::mycommand="";
string msc::execinfo="";
string msc::bamFile="";
vector<string> msc::bamRegion(0);
string msc::refFile="";
string msc::cnvFile="";
string msc::outFile="STDOUT";
string msc::logFile="";
string msc::function="";
int msc::verbose=0;
int msc::numThreads=1;
int msc::maxMR=4000;
int msc::errMatch=2;
int msc::minSNum=11;
int msc::minOverlap=25;
int msc::minOverlapPlus=10;
int msc::minMAPQ=10;
int msc::minBASEQ=0;
int msc::minClusterSize=6;
bool msc::search_length_set_by_user=false;
int msc::maxDistance=(int)1E6;
bool msc::noSecondary=true;
bool msc::dumpBam=false;
int msc::bam_ref=-1;
int msc::bam_l_qseq=0;
bool msc::bam_is_paired=false;
bool msc::bam_pe_disabled=false;
bool msc::bam_pe_set_by_user=false;
int msc::bam_pe_insert=500;
int msc::bam_pe_insert_sd=50;
int msc::bam_rd=0;
int msc::bam_rd_sd=0;
int msc::bam_tid=-1;        //! if msc::tid>0 only process msc::tid
string msc::chr="chr";
string msc::FASTA="";
vector<string> msc::bam_target_name(0);
vector<uint8_t> msc::bdata(0); 
vector<uint8_t> msc::bpdata(0); 
vector<uint32_t> msc::rd(0); 
samfile_t *msc::fp_in = NULL;
bam_index_t *msc::bamidx=NULL;
samfile_t *msc::fp_out = NULL;

bool sort_pair(const intpair_st& p1, const intpair_st& p2)
{
  if ( p1.F2 != p2.F2 ) return p1.F2<p2.F2;
  return p1.R1<p2.R1;
}
bool sort_pair_len(const intpair_st& p1, const intpair_st& p2)
{
  return (p1.R1-p1.F2) < (p2.R1-p2.F2) ;
}
bool sort_pair_info(const pairinfo_st& p1, const pairinfo_st& p2)
{
  if ( p1.F2 != p2.F2 ) return p1.F2<p2.F2;
  return p1.R1<p2.R1;
}
bool sort_pair_info_output(const pairinfo_st& p1, const pairinfo_st& p2)
{
  int p1F=p1.F2;
  int p1R=p1.R1;
  int p2F=p2.F2;
  int p2R=p2.R1;
  
  if ( p1F>p1R ) swap(p1F, p1R);
  if ( p2F>p2R ) swap(p2F, p2R);
  
  if ( p1F != p2F ) return p1F<p2F;
  return p1R<p2R;
}

string cnv_format1()
{
  std::stringstream ss;
  ss << "#CHR" << "\t" 
     << "BEGIN" << "\t" 
     << "END" << "\t" 
     << "TYPE" << "\t"
     << "LENGTH" << "\t"
     << "REF_REPEAT:" << "#" << "\t"
     << "READDEPTH:" 
     << "LSIDE" << ";"
     << "RSIDE" << ";"
     << "BETWEEN" << ":"
     << "RDSCORE" << "\t"
     << "READPAIR:"
     << "CROSSL" << ";"
     << "CROSSR" << ";"
     << "ENVELOPEBOTH" << ":"
     << "RPSCORE" << "\t"
     << "MATCHINGREADS:"
     << "LSIDE" << ";"
     << "RSIDE" << ";"
     << "EDITDISTANCE" << ":"
     << "MRSCORE" << "\t"
     << "SPLITREAD:#;#\t"
     << "Q0:q<=0;q<=10";
  
  return ss.str();
}

string cnv_format1(pairinfo_st &bp1)
{
  std::stringstream ss;
  pairinfo_st bp=bp1;
  
  string RNAME="chr";
  if ( bp.tid>=0 && bp.tid<(int)msc::bam_target_name.size() )
    RNAME=msc::bam_target_name[bp.tid];
  string TYPE= bp.F2 < bp.R1 ? "DEL" : "DUP";
  int len=bp.R1-bp.F2;
  if ( bp.F2>bp.R1 ) { 
    swap(bp.F2, bp.R1);
    swap(bp.F2_rd, bp.R1_rd);
    swap(bp.F2_rp, bp.R1_rp);
    swap(bp.F2_sr, bp.R1_sr);
  }
  
  ss << RNAME << "\t" 
     << bp.F2+1 << "\t" 
     << bp.R1+1 << "\t" 
     << TYPE << "\t"
     << len << "\t"
     << "UN:" << bp.un << "\t"
     << "RD:" 
     << bp.F2_rd << ";"
     << bp.R1_rd << ";"
     << bp.rd << ":"
     << bp.rdscore << "\t"
     << "RP:"
     << bp.F2_rp << ";"
     << bp.R1_rp << ";"
     << bp.FRrp << ":"
     << bp.rpscore << "\t"
     << "MR:"
     << bp.F2_sr << ";"
     << bp.R1_sr << ";"
     << bp.MS_ED << ":"
     << bp.srscore << "\t"
     << "SR:"
     << bp.sr_ed << ";"
     << bp.sr_count;
  // << bp.sr_count << "\t"
  // << BPMARKER ;
  
  return ss.str() ;
}

string cnv_format_all(pairinfo_st &bp1)
{
  std::stringstream ss;
  
  int F2=bp1.F2;
  int R1=bp1.R1;
  if ( F2> R1 ) swap(F2, R1);
  
  double q0, q10;
  check_map_quality(bp1.tid, F2, R1, q0, q10);  
  
  int i0=q0*100+0.5;
  int i10=q10*100+0.5;
  
  ss << cnv_format1(bp1) << "\tQ0:" << i0 << ";" << i10;
  
  return ss.str() ;
}

string mr_format1(pairinfo_st &ibp)
{
  std::stringstream ss;
  pairinfo_st bp=ibp;
  
  string RNAME="chr";
  if ( bp.tid>=0 && bp.tid<(int)msc::bam_target_name.size() )
    RNAME=msc::bam_target_name[bp.tid];
  string TYPE= bp.MS_F2 < bp.MS_R1 ? "DEL" : "DUP";
  int len=bp.MS_R1-bp.MS_F2;
  if ( bp.MS_F2>bp.MS_R1 ) {
    swap (bp.MS_F2, bp.MS_R1);
    swap (bp.MS_F2_rd, bp.MS_R1_rd);
    swap (bp.F2_sr, bp.R1_sr);
  }
  
  ss << RNAME << " " 
     << bp.MS_F2+1 << " " 
     << bp.MS_R1+1 << " " 
     << TYPE << " "
     << len << " "
     << "UN:" << bp.un << " "
     << "MD:"
     << bp.MS_F2_rd << ";"
     << bp.MS_R1_rd << " "
     << "MR:"
     << bp.F2_sr << ";"
     << bp.R1_sr << ";"
     << bp.MS_ED << ":"
     << bp.srscore << " "
     << "SR:"
     << bp.sr_ed << ";"
     << bp.sr_count;
  
  return ss.str() ;
}

void write_cnv_to_file(vector<pairinfo_st>& bp)
{
  for(size_t i=0;i<bp.size();++i) cout << cnv_format_all(bp[i]) << endl;
  return;
}

void write_cnv_to_file(vector<pairinfo_st>& bp, string fn)
{
  static vector<string> files(0);
  if ( fn=="STDOUT" ) {
    for(size_t i=0;i<bp.size();++i) cout << cnv_format_all(bp[i]) << endl;
    return;
  }
  
  bool is_new=false;
  is_new = find(files.begin(), files.end(), fn) == files.end() ? true : false ;
  
  ofstream FOUT;
  if ( is_new ) {
    FOUT.open(fn.c_str());
    files.push_back(fn);
    FOUT << cnv_format1() << endl;
  }    
  else FOUT.open(fn.c_str(), std::ofstream::app);
  
  for(size_t i=0;i<bp.size();++i) FOUT << cnv_format_all(bp[i]) << endl;
  FOUT.close();
  return;
}

string tmpfile(int thread_id) {
  string tmps="matchclipstmpdata."+
    to_string(msc::pid)+"."+to_string(msc::numThreads)+"."+to_string(thread_id);
  return tmps;
}

void write_ED_st_to_file(vector<ED_st>& ed, string fn)
{
  static vector<string> files(0);
  
  bool is_new=false;
  is_new = find(files.begin(), files.end(), fn) == files.end() ? true : false ;
  
  ofstream FOUT;
  if ( is_new ) {
    FOUT.open(fn.c_str());
    files.push_back(fn);
  }    
  else FOUT.open(fn.c_str(), std::ofstream::app);
  
  for(size_t i=0;i<ed.size();++i) 
    FOUT << ed[i].iL << "\t" << ed[i].F2 << "\t" << ed[i].iR << "\t" << ed[i].R1 
	 << ed[i].ED << "\t" << ed[i].count << "\n";
  
  FOUT.flush();
  return;
}

void remove_N_regions(string& fasta, vector<pairinfo_st>& bp)
{
  vector<int> N_beg(0), N_end(0);
  get_N_regions(fasta, N_beg, N_end);
  for(int i=bp.size()-1; i>=0; --i) {
    int F2=bp[i].F2;
    int R1=bp[i].R1;
    if ( F2>R1 ) swap(F2, R1);
    
    bool is_N=false;
    for(size_t k=0; k<N_beg.size(); ++k) {
      int r1=max(N_beg[k], F2);
      int r2=min(N_end[k], R1);
      if ( r2-r1 > (R1-F2)/3 ) { is_N=true; break; }
    }
    
    if ( is_N && bp[i].rdscore>=1 ) bp.erase( bp.begin()+i );
  }
  
  return;
}


void finalize_output(vector<pairinfo_st>& bp, 
		     vector<pairinfo_st>& strong, 
		     vector<pairinfo_st>& weak) 
{
  strong.clear();
  weak.clear();
  if ( bp.size() < 1 ) return;
  
  for(size_t i=1; i<bp.size(); ++i) cerr << cnv_format1(bp[i]) << endl;
  
  // remove duplicated 
  for(size_t i=1; i<bp.size(); ++i) {
    if ( bp[i].F2==bp[i-1].F2 && bp[i].R1==bp[i-1].R1 ) {
      int ikeep=i;
      
      ikeep = bp[i-1].sr_count > bp[i].sr_count ? i-1 : i;
      bp[i].un=bp[ikeep].un;
      bp[i].sr_ed=bp[ikeep].sr_ed;
      bp[i].sr_count=bp[ikeep].sr_count;
      
      if ( bp[i-1].rdscore != bp[i].rdscore )
	ikeep = bp[i-1].rdscore > bp[i].rdscore ? i-1 : i;
      bp[i].F2_rd=bp[ikeep].F2_rd;
      bp[i].R1_rd=bp[ikeep].R1_rd;
      bp[i].rd=bp[ikeep].rd;
      bp[i].rdscore=bp[ikeep].rdscore;
      
      if ( bp[i-1].rpscore != bp[i].rpscore )
	ikeep = bp[i-1].rpscore > bp[i].rpscore ? i-1 : i;
      bp[i].F2_rp=bp[ikeep].F2_rp;
      bp[i].R1_rp=bp[ikeep].R1_rp;
      bp[i].FRrp=bp[ikeep].FRrp;
      bp[i].rpscore=bp[ikeep].rpscore;
      
      if ( bp[i-1].srscore != bp[i].srscore )
	ikeep = bp[i-1].srscore > bp[i].srscore ? i-1 : i;
      bp[i].F2_sr=bp[ikeep].F2_sr;
      bp[i].R1_sr=bp[ikeep].R1_sr;
      bp[i].MS_ED=bp[ikeep].MS_ED;
      bp[i].MS_ED_count=bp[ikeep].MS_ED_count;
      bp[i].srscore=bp[ikeep].srscore;
      
      bp[i-1].sr_count=-9;
      bp[i-1].rdscore=bp[i-1].rpscore=bp[i-1].srscore=-1;
    }
  }
  
  // decide signal strength
  for(size_t i=0; i<bp.size(); ++i) {
    if ( bp[i].sr_count==-9 ) continue;
    
    // to discard
    if ( bp[i].rdscore<=0 && 
	 bp[i].rpscore<=0 &&
	 bp[i].srscore<=0 &&
	 bp[i].sr_count<=3 ) continue;
    
    if ( bp[i].rdscore<=0 && 
	 bp[i].rpscore<=0 &&
	 bp[i].srscore<=0 &&
	 bp[i].un>=msc::minOverlap &&
	 bp[i].un>=abs(bp[i].F2-bp[i].R1) ) continue;
    
    // RD and RP signal
    // if ( bp[i].rdscore>0 || bp[i].rpscore>0 ) {
    // strong.push_back( bp[i] );
    // continue;
    //}
    // combined score >=2
    if ( max(0, bp[i].rpscore)+
	 max(0, bp[i].rdscore)+
	 max(0, bp[i].srscore) >=2 ) {
      strong.push_back( bp[i] );
      continue;
    }
    
    // short CNV, strong SR signal
    if ( bp[i].rpscore<0 &&
	 bp[i].srscore>2 &&
	 bp[i].un< (float)msc::minOverlap*1.5  ) {
      strong.push_back( bp[i] );
      continue;
    }
    if ( abs(bp[i].F2-bp[i].R1)<msc::bam_pe_insert_sd*5 &&
	 bp[i].srscore>2 &&
	 bp[i].un< (float)msc::minOverlap*1.5  ) {
      strong.push_back( bp[i] );
      continue;
    }
    
    
    // strong softclips signal
    int medRD=(bp[i].F2_rd +  bp[i].R1_rd)/2;
    double expected_pairs= (double)medRD * 0.5 * 0.5 *
      (1.0 -(double)msc::minOverlap/(double)msc::bam_l_qseq ) *
      (1.0 - 2.0*(double)msc::minSNum/(double)msc::bam_l_qseq );
    if ( expected_pairs<20 ) expected_pairs=20;
    if ( bp[i].un< (float)msc::minOverlap*1.5 && 
	 bp[i].sr_count>=expected_pairs ) {
      strong.push_back( bp[i] );
      continue;
    }
    
    if ( bp[i].sr_count<3 && abs(bp[i].F2-bp[i].R1)<=bp[i].un ) continue;
    
    weak.push_back( bp[i] );
  }
  
  for(int i=strong.size()-1; i>=0; --i) {
    bool is_weak=false;
    
    if ( strong[i].rdscore<=0 && strong[i].rpscore==0 &&
	 abs(strong[i].F2-strong[i].R1) > msc::bam_pe_insert_sd*7 ) is_weak=true;
    
    if ( strong[i].rdscore<=0 && strong[i].rpscore==0 &&
	 strong[i].un*2>msc::bam_l_qseq ) is_weak=true;
    
    if ( is_weak ) {
      weak.push_back(strong[i]);
      strong.erase(strong.begin()+i);
    }
  }
  
  return;
}


void get_pairend_info(int ref, int beg, int end)
{
  bam1_t *b=NULL; b = bam_init1();
  bam_iter_t iter;
  
  size_t count=0;
  int l_qseq=100; double l_qseq_sum=0.0;
  double isize=0.0, isize2=0.0, isize_c=0, isize_sd=0.0;
  
  iter = bam_iter_query(msc::bamidx, ref, beg, end);
  
  vector<double> pe(0); pe.reserve(15000);
  while(  bam_iter_read(msc::fp_in->x.bam, iter, b)  > 0 ) {
    ++count;
    l_qseq_sum+=b->core.l_qseq;
    if ( (int)b->core.tid < 0 ) continue;
    if ( b->core.mtid != b->core.tid && b->core.mtid>0 ) continue;
    if ( bool(b->core.flag&BAM_FREVERSE) ==
	 bool(b->core.flag&BAM_FMREVERSE) ) continue;
    if ( b->core.flag & BAM_DEF_MASK ) continue;
    if ( !(b->core.flag & BAM_FPROPER_PAIR) ) continue;
    if ( b->core.qual < msc::minMAPQ ) continue;
    
    isize   += abs(b->core.isize);
    isize2  += (double)b->core.isize * (double)b->core.isize;
    isize_c += 1;
    pe.push_back( abs(b->core.isize) );
    if ( isize_c > 10000 ) break;
    if ( count > 200000 ) break;
  }
  cerr << endl;
  
  if (count>0) l_qseq=l_qseq_sum/count+0.5;
  else {
    cerr << "cannot find any reads in region: "
	 << msc::fp_in->header->target_name[ ref ] << "\t" 
	 << beg << "\t" << end 
	 << endl; 
    return;
  }
  
  //update main control
  msc::bam_ref=ref;
  msc::bam_l_qseq=l_qseq;
  if ( isize_c>2 ) {
    isize/=isize_c;
    isize_sd =sqrt( (isize2-isize*isize*isize_c)/isize_c );
    if ( isize_sd<30 ) {
      cerr << "isize sd is too small " << isize_sd << " changed to 30" << endl;  
      isize_sd=30;
    }
    
    sort( pe.begin(), pe.end() );
    // for(int i=0; i<pe.size(); ++i) cerr << pe[i] << "\t"; cerr << endl;
    isize=pe[ pe.size()/2 ];
    isize_sd= ( pe[pe.size()/4*3] - pe[pe.size()/4] ) / 1.35;
    
    msc::bam_is_paired=true;
    msc::bam_pe_insert=(int)isize;
    msc::bam_pe_insert_sd=(int)isize_sd;
  }
  cerr << "sampled from " << isize_c << " reads" << endl
       << "bam_target\t" << msc::fp_in->header->target_name[ ref ] << "\n"
       << "bam_l_qseq\t" << msc::bam_l_qseq << "\n"
       << "bam_is_paired\t" << msc::bam_is_paired << "\n"
       << "bam_pe_insert\t" << msc::bam_pe_insert << "\n"
       << "bam_pe_insert_sd\t" << msc::bam_pe_insert_sd << "\n"
       << endl;
  
  
  if ( isize>1000 || isize_sd>2*isize || isize_sd<0 ) {
    msc::bam_pe_insert_sd=msc::bam_pe_insert/2;
    cerr << "unusual behavior\n"
	 << "continue with default values " << msc::bam_pe_insert << " "
	 << msc::bam_pe_insert_sd << endl;
  }

  if ( msc::minOverlap*4<msc::bam_l_qseq ) {
    msc::minOverlap=msc::bam_l_qseq/4;
    cerr << "length of minimum overlap changed to: " << msc::minOverlap << endl;
  }

  if ( b ) bam_destroy1(b);
  if ( iter) bam_iter_destroy(iter);
  return;
}


void get_bam_info()
{
  
  bam1_t *b=NULL;   
  b = bam_init1();
  bam_iter_t iter=0;
  int ref=0, beg=0, end=0x7fffffff;
  
  for(int i=0;i<msc::fp_in->header->n_targets;++i) 
    msc::bam_target_name.push_back( string(msc::fp_in->header->target_name[i]) );
  
  if ( msc::bam_target_name.size() == 0 ) 
    cerr << "#no header in bam file? " << msc::bamFile << endl;
  
  if ( msc::bamRegion.size() == 0 ) msc::bamRegion=msc::bam_target_name;
  bool have_reads=true;
  for( int i=0; i<(int)msc::bamRegion.size(); ++i) {
    have_reads=true;
    if ( msc::bamRegion[i].find(":") == string::npos ) msc::bamRegion[i]+=":";
    ref=-1;
    int is_solved=bam_parse_region(msc::fp_in->header, msc::bamRegion[i].c_str(), &ref, &beg, &end); 
    if ( is_solved<0 || ref<0 || ref>=msc::fp_in->header->n_targets ) {
      cerr << "Cannt resolve region " << msc::bamRegion[i] << endl;
      have_reads=false;
    }
    else {
      iter = bam_iter_query(msc::bamidx, ref, beg, end);
      if ( bam_iter_read(msc::fp_in->x.bam, iter, b)>0 ) {
	string RNAME=msc::fp_in->header->target_name[b->core.tid];
	if ( msc::bamRegion[i].find(RNAME) == 0 ) 
	  cerr << msc::bamFile << " has reads on " << msc::bamRegion[i] << endl;
	else have_reads=false;
      }
      else have_reads=false;;
    }
    
    if ( have_reads==false ) {
      // cerr << msc::bamFile << " does not have reads on " << msc::bamRegion[i] << endl;
      msc::bamRegion[i]="NA";
    }
  }
  
  if ( b ) bam_destroy1(b);
  if ( iter ) bam_iter_destroy(iter);
  
  return;
}

int usage_match_MS_SM_reads(int argc, char* argv[]) {
  string app=string(argv[0]);
  if ( app.rfind('/') != string::npos ) app=app.substr(1+app.rfind('/'));
  cerr << "Usage:\n" 
       << "  " << app << " <options> -f REFFILE -b BAMFILE [REGION]\n"
       << "\nOptions:\n"
       << "  -t  INT  number of threads, INT=1 \n"
       << "  -e  INT  max allowed mismatches when matching strings, INT=2 \n"
       << "  -l  INT  minimum length of overlap, INT=25 \n"
       << "  -s  INT  minimum number of soft clipped bases, INT=10 \n"
       << "  -q  INT  minimum mapping score, INT=10 \n"
       << "  -Q  INT  minimum base read quality, INT=0 \n"
       << "  -L  INT  check reads before and after INT bases, INT=auto \n"
       << "  -L  0    pair end mode only \n"
       << "  -se      single end mode, do not use pair end distances\n"
       << "  -pe INT INT provide insert and s.d. of insert, otherwise calculate them\n"
       << "  -o  STR  outputfile, STR=STDOUT \n"
       << "   REGION  if given should be in samtools's region format \n"
       << "\nExamples:\n"
       << "  " << app << "      -f hg19.fasta -b A.bam -o A.txt\n"
       << "  " << app << " -t 4 -f hg19.fasta -b A.bam chr1 -o A.txt\n"
       << "\nDiscussion:\n"
       << "  1. max allowed mismatches when matching reads\n"
       << "     This number is also adjusted according to the length of\n"
       << "     overlap. If this numer multipled by 12 is greated than the length\n"
       << "     false is returned.\n"
       << "  2. minimum length of overlap\n"
       << "     This number is also adjusted so that it is no less than\n"
       << "     1/4 of the reads.\n"
       << "  3. mapping quality\n"
       << "     Mapping quality is checked for reads matching and pair end distances\n"
       << "     but not for read depth. It is so to be consistent with samtools\n"
       << "     mpileup. The BAM_DEF_MASK filter is used for all reads. Reads that are\n"
       << "     BAM_FUNMAP or BAM_FSECONDARY or BAM_FQCFAIL or BAM_FDUP are ignored.\n"
       << "  4. pair-end insert and s.d.\n"
       << "     Pair-end insert and s.d. are calculated from the bam file if not given.\n"
       << "     They can be overwritten by the -pe INT INT option. Inserts longer than\n"
       << "     insert+6s.d. are used to find structure variations. Reads matching is\n"
       << "     performed for reads whose POS's are with in insert+8s.d. So basically,\n"
       << "     pair end distances and read depths are used to find longer variations,\n"
       << "     and split reads are used to find shorter variations. This allows for\n"
       << "     fast processing of high coverage whole genome data. The range of reads\n"
       << "     matching can be overwritten by the -L INT option.\n"
       << "\nOutput:\n"
       << "  Please check https://github.com/yhwu/matchclips2\n"
       << "\nReference: doi: 10.3389/fgene.2013.00157"
       << endl;
  
  return(0);
}

int get_parameters(int argc, char* argv[])
{
#define _next3 ARGV[i]=""; ARGV[i+1]=""; ARGV[i+2]=""; continue;
#define _next2 ARGV[i]=""; ARGV[i+1]=""; continue;
#define _next1 ARGV[i]=""; continue;
  
  size_t i;
  vector<string> ARGV(0);

  for(i=0; (int)i<argc; ++i) ARGV.push_back( string(argv[i]) );
  
  msc::mycommand=ARGV[0];
  for(i=1;i<ARGV.size();++i) msc::mycommand+=" "+ARGV[i];
  
  for(i=1;i<ARGV.size();++i) {
    if ( ARGV[i]=="-b" ) {    // input bam file
      msc::bamFile=ARGV[i+1]; // input regions
      for( int k=i+2; k<(int)ARGV.size(); ++k ) {
	if ( ARGV[k][0] == '-' ) break;
	msc::bamRegion.push_back( ARGV[k] );
	ARGV[k]="";
      }
      _next2;
    }
    if ( ARGV[i]=="-f" ) { msc::refFile=ARGV[i+1]; _next2; }
    if ( ARGV[i]=="-o" ) { msc::outFile=ARGV[i+1]; _next2; }
    if ( ARGV[i]=="-t" ) { msc::numThreads=atoi(ARGV[i+1].c_str()); _next2; }
    if ( ARGV[i]=="-e" ) { msc::errMatch=atoi(ARGV[i+1].c_str()); _next2; }
    if ( ARGV[i]=="-l" ) { msc::minOverlap=atoi(ARGV[i+1].c_str()); _next2; }
    if ( ARGV[i]=="-s" ) { msc::minSNum=atoi(ARGV[i+1].c_str()); _next2; }
    if ( ARGV[i]=="-q" ) { msc::minMAPQ=atoi(ARGV[i+1].c_str()); _next2; }
    if ( ARGV[i]=="-Q" ) { msc::minBASEQ=atoi(ARGV[i+1].c_str()); _next2; }
    if ( ARGV[i]=="-2" ) { msc::noSecondary=false; _next1; }
    if ( ARGV[i]=="-L" ) { 
      msc::search_length_set_by_user=true;
      msc::maxDistance=atoi(ARGV[i+1].c_str()); _next2; 
    }
    if ( ARGV[i]=="-v" ) { msc::verbose=1; _next1; }
    if ( ARGV[i]=="-vv" ) { msc::verbose=2; _next1; }
    if ( ARGV[i]=="-vvv" ) { msc::verbose=3; _next1; }
    if ( ARGV[i]=="-dump" ) { msc::dumpBam=true; _next1; }
    if ( ARGV[i]=="-cnv" ) { msc::cnvFile=ARGV[i+1]; _next2; }
    if ( ARGV[i]=="-d" ) { msc::dx=atoi(ARGV[i+1].c_str()); _next2; }
    if ( ARGV[i]=="-se" ) { msc::bam_pe_disabled=true; _next1; }
    if ( ARGV[i]=="-pe" ) { 
      msc::bam_pe_set_by_user=true;
      msc::bam_pe_insert=atoi(ARGV[i+1].c_str());
      msc::bam_pe_insert_sd=atoi(ARGV[i+2].c_str());
      if ( to_string(abs(msc::bam_pe_insert)) == ARGV[i+1] && 
	   to_string(abs(msc::bam_pe_insert_sd)) == ARGV[i+2] ) _next3;
      cerr << "pair end parameters are given by -pe isize isize_std" << endl;
      exit(0);
    }
  }
  
  if ( msc::bamFile=="" ) {
    cerr << "Need bam file\n";
    exit( usage_match_MS_SM_reads(argc, argv) );
  }
  if ( msc::refFile=="" ) {
    cerr << "Need reference file\n";
    exit( usage_match_MS_SM_reads(argc, argv) );
  }
  bool is_unknown_parameter=false;
  for(i=1;i<ARGV.size();++i) {
    if ( ARGV[i]!="" ) {
      cerr << "unknown argument:\t" << ARGV[i] << endl;
      is_unknown_parameter=true;
    }
  }
  if ( is_unknown_parameter ) exit( usage_match_MS_SM_reads(argc, argv) );
  
  string tmps="";
  for(i=0;i<msc::bamRegion.size();++i) {
    tmps+=msc::bamRegion[i]+" ";
  }
  msc::execinfo="#Command Line     : " + msc::mycommand + "\n" +
    "#Input bamfile    : " + msc::bamFile + "\n" +
    "#Bamfile region   : " + tmps + "\n" +
    "#FASTAfile        : " + msc::refFile + "\n" +
    "#Num of S bases   : " + to_string(msc::minSNum) + "\n" +
    "#Minumum overlap  : " + to_string(msc::minOverlap) + "\n" +
    "#Minimum mapq     : " + to_string(msc::minMAPQ) + "\n" +
    "#Allowed mismatch : " + to_string(msc::errMatch) + "\n" +
    "#Maximum distance : " + to_string(msc::maxDistance) + "\n" +
    "#Output           : " + msc::outFile + "\n"+
    "#Output           : " + msc::outFile + ".weak\n";
  
  cerr << msc::execinfo << endl;
  
  return 1;
}

void match_MS_SM_reads(int argc, char* argv[])
{
  if ( argc<3 )  exit( usage_match_MS_SM_reads(argc, argv) );
  get_parameters(argc, argv);
  
  string FASTA; string fastaname="";
  
  int ref=0, beg=0, end=0x7fffffff;
  
  // load BAM and index; check regions
  msc::fp_in=samopen(msc::bamFile.c_str(), "rb", 0);
  if ( ! msc::fp_in ) {
    cerr << msc::bamFile << " not found!" << endl;
    exit(0);
  }
  msc::bamidx=bam_index_load(msc::bamFile.c_str()); 
  if ( ! msc::bamidx ) {
    cerr << msc::bamFile << " idx not found!" << endl;
    exit(0);
  }
  get_bam_info();
  
  if ( msc::dumpBam ) {
    msc::fp_out = msc::outFile=="STDOUT" ? 
      samopen("-", "w", msc::fp_in->header) :
      samopen(msc::outFile.c_str(), "wb", msc::fp_in->header) ;
  }
  
  for(int ichr=0; ichr<(int)msc::bamRegion.size(); ++ichr ) {
    if ( msc::bamRegion[ichr]=="NA" ) continue;
    cerr << "processing region:\t" << msc::bamRegion[ichr] << endl;
    
    ref=-1; beg=0; end=0x7fffffff;
    int is_solved=bam_parse_region(msc::fp_in->header, msc::bamRegion[ichr].c_str(), &ref, &beg, &end); 
    if ( is_solved<0 || ref<0 || ref>=(int)msc::bam_target_name.size() ) continue;
    msc::bam_ref=ref;
    
    if ( !msc::bam_pe_set_by_user ) get_pairend_info(ref, beg, end);
    
    //! load reference sequence
    if ( msc::bam_target_name[msc::bam_ref] != fastaname ) {
      fastaname=msc::bam_target_name[ref];
      load_reference(msc::refFile, fastaname, FASTA);
      if ( FASTA.size() != msc::fp_in->header->target_len[msc::bam_ref] )
	cerr << "not exactly the same reference, expected length " 
	     << msc::fp_in->header->target_len[msc::bam_ref]  
	     << " loaded " << FASTA.size() 
	     << endl;
    }
    
    vector<intpair_st> pairs(0);
    vector<bam1_t> b_MS(0);
    vector<bam1_t> b_SM(0);
    
    int min_pair_length=msc::bam_pe_insert+msc::bam_pe_insert_sd*6;
    //if ( min_pair_length<1000 ) min_pair_length=1000;
    prepare_pairend_matchclip_data(ref, beg, end, min_pair_length, FASTA,
				   pairs, b_MS, b_SM);
    
    vector<pairinfo_st> pairbp_pe(0);
    if (! msc::bam_pe_disabled ) {
      pair_guided_search(pairs, FASTA, pairbp_pe) ;
      // pair_guided_search(ref, beg, end, min_pair_length, FASTA, pairbp_pe);
      vector<intpair_st> (0).swap(pairs);
    }
    
    vector<pairinfo_st> pairbp_mc(0);
    int search_length=msc::bam_pe_insert+msc::bam_pe_insert_sd*8;
    if ( msc::bam_pe_disabled || msc::search_length_set_by_user ) 
      search_length=msc::maxDistance;
    exhaustive_search(b_MS, b_SM, search_length, FASTA, pairbp_mc);
    // exhaustive_search(ref, beg, end, search_length, FASTA, pairbp_mc);
    
    pairbp_mc.insert(pairbp_mc.end(), pairbp_pe.begin(), pairbp_pe.end() ); 
    sort(pairbp_mc.begin(), pairbp_mc.end(), sort_pair_info);
    
    vector<pairinfo_st> strong, weak;
    remove_N_regions(FASTA, pairbp_mc);
    finalize_output(pairbp_mc, strong, weak);    
    sort(strong.begin(), strong.end(), sort_pair_info_output);
    sort(weak.begin(), weak.end(), sort_pair_info_output);
    write_cnv_to_file(strong, msc::outFile);
    write_cnv_to_file(weak, string(msc::outFile+".weak"));    
    
  } // done
  
  if ( msc::fp_in ) samclose(msc::fp_in);
  if ( msc::bamidx )bam_index_destroy(msc::bamidx);
  if ( msc::fp_out ) samclose(msc::fp_out);
  if ( msc::dumpBam && msc::outFile!="" ) bam_index_build(msc::outFile.c_str());
  
  cerr << msc::execinfo << endl;
  return;
}

