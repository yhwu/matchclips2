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
#include <vector>
#include <unistd.h>
#include "color.h"
using namespace std;

/**** user samtools headers ****/
#include "samfunctions.h"

/**** user headers ****/
#include "functions.h"
#include "readref.h"
#include "read_CIGAR.h"
#include "filters.h"
#include "matchreads.h"

vector<string> BREAKSEEK_st::TARGET;
map<string, int> BREAKSEEK_st::RNAMEINT;
int BREAKSEEK_st::ctid;

int MINIMUMS=11;     // minimum S number
int MINIMUMOVERLAP=25;     // minimum length of common string
int MINIMUMMAPQ=10;  // minumum mapping quality
int MINIMUMBASEQ=0; // minumum base quality
int ERR_EDGE=1;      // ignore error in first EDGE and last EDGE bases in a read
int ERR_MATCH=3;     // allowed mismatch when matching ms sm reads
bool SHOWMAPPED=false; // only print mappped part; otherwise print all reads

int maxMismatch=2;
long int BUFFERSIZE=(int)600E6; // buffer size to hold reads, two buffers used
string SEQBUFFER="";

//! commonString : common string with a tolorance of ERR_MATCH with ERR_EDGE ignored
//! p1 : position of the beginning of match in readMS, 0-based
//! p2 : position of the beginning of match in readSM, 0-based; 0 is matched 
bool is_MS_SM_overlap_fuzzy_lite(const string& readMS, const string& readSM, string& commonString, int& p1, int& p2) 
{
  p1=p2=-1;
  commonString="";
  
  int i,j,k;
  int maxerror=ERR_MATCH;
  
  bool match=false;
  int ndiff=0;
  
  int minOver=20;
  for(i=readMS.length()-minOver+1-1;i>=0;--i) {
    ndiff=0;
    for(k=i+ERR_EDGE,j=ERR_EDGE;
	k<(int)readMS.length()-ERR_EDGE && j<(int)readSM.length();
	++k,++j) {
      ndiff += ( readMS[k]!=readSM[j] );
      if (ndiff>maxerror) break;
    }
    if (ndiff<=maxerror) { 
      match=true;
      p1=i;
      p2=0;
      //commonString=readMS.substr(p1);
      string comMS=readMS.substr(p1);
      string comSM=readSM.substr(0,comMS.length());
      int L_MS=comMS.length()/2;
      commonString=comMS.substr(0,L_MS)+comSM.substr(L_MS);
      //commonString=comSM.substr(0,L_MS)+comMS.substr(L_MS);
      return(match);
    }
  }
  
  return(match);
  
}

void readsMSandSMmatchfuzzy(const string& readMS, int Sms, 
			    const string& readSM, int Ssm, 
			    string& commonString, int& p1, int& p2) 
{
  p1=p2=-1;
  commonString="";
  
  int i,j,k;
  int maxerror=ERR_MATCH;
  
  bool match=false;
  int ndiff=0;
  
  //! look for deletion and tandem duplication
  //! for those, overlap is equal or longer than Sms + Ssm
  //int minOver=Sms+Ssm;
  int minOver=20;
  for(i=readMS.length()-minOver+1-1;i>=0;--i) {
    ndiff=0;
    for(k=i+ERR_EDGE,j=ERR_EDGE;
	k<(int)readMS.length()-ERR_EDGE && j<(int)readSM.length();
	++k,++j) {
      ndiff += ( readMS[k]!=readSM[j] );
      if (ndiff>maxerror) break;
    }
    if (ndiff<=maxerror) { 
      match=true;
      p1=i;
      p2=0;
      //commonString=readMS.substr(p1);
      string comMS=readMS.substr(p1);
      string comSM=readSM.substr(0,comMS.length());
      int L_MS=comMS.length()/2;
      commonString=comMS.substr(0,L_MS)+comSM.substr(L_MS);
      //commonString=comSM.substr(0,L_MS)+comMS.substr(L_MS);
      break;
    }
  }
  
  return;
  
}


//! print out BP
string BP_format(BREAKSEEK_st &BP, bool h)
{
  std::stringstream ss;

  if ( h ) { // header?
    ss << "#RNAME\t" 
       << "START\t"
       << "END\t"
       << "TYPE\t"
       << "UNCERTAINTY\t"
       << "LENGTH\t"
       << "INSERT\t"
       << "MATCHINFO\t"
       << "MATCHCOUNT\t"
       << "VALID\t" 
       << "SVSEQ\t"
       << "MERGE\t"
       << "BP";
  }
  else {
    string BPMARKER="BP_"+to_string(BP.T)+"_"+
      BP.RNAME()+":"+
      to_string(BP.P1)+"-"+to_string(BP.P2); 
    ss << BP.RNAME() << "\t" 
       << BP.P1 << "\t" 
       << BP.P2 << "\t" 
       << BP.T << "\t" 
       << BP.UN << "\t" 
       << BP.len << "\t"
       << BP.INSEQ << "\t"
       << "CL" << BP.CL << ","
       << "ED" << BP.ED << ","
       << "L" << BP.L1 << ","
       << "S" << BP.S1 << ","
       << "N" << BP.n1 << ","
       << "Q" << BP.q1 << ","
       << "L" << BP.L2 << ","
       << "S" << BP.S2 << ","
       << "N" << BP.n2 << ","
       << "Q" << BP.q2 << ";"
       << BP.TAG << "\t"
       << BP.count << "\t"
       << BP.isit << "\t"
       << BP.SVSEQ << "\t"
       << BP.MERGE << "\t"
       << BPMARKER;
    //	 << BPSEEK[i].ED5p << "\t"
    //	 << BPSEEK[i].ED3p << "\t"
  }
  
  return ss.str() ;
}

void load_cnv_from_file(vector<BREAKSEEK_st>& BPSEEK, string cnvFile)
{
  size_t i;
  
  BPSEEK.clear();
  BREAKSEEK_st ABP;
  ifstream FIN(cnvFile.c_str());
  if ( !FIN ) { cerr << "Can't open file " << cnvFile << endl; exit(0); }
  
  while ( !FIN.eof() ) {
    string chr,tmps,comment;
    getline(FIN,tmps);
    if ( tmps.size()<3 ) continue;
    if ( tmps[0]=='#' ) continue;
    
    istringstream iss(tmps);
    
    ABP.P1=ABP.P2=ABP.UN=0;
    iss >> chr
	>> ABP.P1
	>> ABP.P2
	>> ABP.T
	>> ABP.UN
	>> ABP.len
	>> ABP.INSEQ
	>> comment
	>> ABP.count
	>> ABP.isit
	>> ABP.SVSEQ
	>> ABP.MERGE;
    
    if ( ABP.RNAMEINT.count(chr)==0 ) 
      cerr << "#" << chr << " not found in bam file" << endl;
    ABP.tid=ABP.RNAMEINT[chr];
    
    for(i=0;i<comment.size();++i) {
      if ( comment[i]==',' ) comment[i]='\t';
      if ( comment[i]=='C' ) comment[i]=' ';
      if ( comment[i]=='L' ) comment[i]=' ';
      if ( comment[i]=='E' ) comment[i]=' ';
      if ( comment[i]=='D' ) comment[i]=' ';
      if ( comment[i]=='S' ) comment[i]=' ';
      if ( comment[i]=='N' ) comment[i]=' ';
      if ( comment[i]=='Q' ) comment[i]=' ';
    }
    iss.clear();
    iss.str(comment);
    
    iss >> ABP.CL
	>> ABP.ED
	>> ABP.L1
	>> ABP.S1
	>> ABP.n1
	>> ABP.q1
	>> ABP.L2
	>> ABP.S2
	>> ABP.n2
	>> ABP.q2;
    
    BPSEEK.push_back(ABP);
  }
  
  // cerr << "laoded " << BPSEEK.size() << " cnv" << endl;
  FIN.close();
  return;
}


void push_alignment_to_string(POSCIGAR_st &bm, string &BUF, vector<RSAI_st>& READ)
{
  RSAI_st iread;
  
  iread.tid=bm.tid;
  iread.pos=bm.pos;
  iread.q1=bm.qual;
  iread.len=bm.l_qseq;
  iread.len_cigar=bm.cigar_s.size();
  iread.M=bm.nop[bm.anchor];
  iread.Mrpos=bm.qop[bm.anchor];
  iread.S=iread.Spos=iread.sbeg=iread.send=0;
  if ( bm.iclip>=0 ) {
    iread.S=bm.nop[bm.iclip];
    iread.Spos=bm.cop[bm.iclip];
    iread.sbeg=bm.cop[bm.iclip];
    iread.send=bm.cop[bm.iclip]+bm.nop[bm.iclip]-1;
  }
  iread.pos_beg=bm.cop[0];
  iread.pos_end=bm.cop.back()+bm.nop.back()-1;
  
  iread.p1 = BUF.size();
  BUF += bm.seq + bm.cigar_s + "\t";
  READ.push_back(iread); 
  
  return;
}

int usage_plot_reads(int argc, char* argv[]) {
  cerr << "Usage:\n" 
       << "  " << argv[0] << " <options> -f REFFILE -b BAMFILE -cnv CNVFILE\n"
       << "\nOptions:\n"
       << "  -M       only print matched parts of reads; default is to print all parts\n"
       << "  -o  STR  outputfile, STR=STDOUT \n"
       << "  -l  INT  length of reads to print out, default is full length \n"
       << "  -ee INT  maximum allowed mismatches at 3' edges, INT=1 \n"
       << "  -me INT  max allowed mismatches when matching strings, INT=1 \n"
       << "\nExamples:\n"
       << "  " << argv[0] << " -f human_g1k_v37.fasta -b bwap40X.bam -cnv bwap40X.bam.txt -o tmp.txt\n"
       << endl;
  
  return(0);
}
int main(int argc, char* argv[])
{
  if ( argc<3 )  exit( usage_plot_reads(argc, argv) );
  
  string mycommand="";
  int i, k;
  int readlength=-1;
  string QNAME,RNAME,POS,MAPQ,CIGAR,MRNM,MPOS,ISIZE,SEQ,QUAL,OPT;
  string QSEQ;
  string CIGAR1,SEQ1,QUAL1,OPT1;
  string REFSEQ=string(200,'C');
  string REFCLIP=string(200,'C');
  string CLIPPEDSEQ=string(200,'C');
  string CLIPPEDQUAL=string(200,'C');
  string CLIPPEDCIGAR=string(200,'C');
  vector<size_t>mismatchCount(10000,0);
  vector<size_t>mapqCount(10000,0);
  string oper;
  //FILE *fp;
  bool NOSECONDARY=true;
  string bamFile="",bamRegion="",fastaFile="",outputFile="STDOUT";
  string cnvFile="";
  vector<string> inputArgv;
  vector<BREAKSEEK_st> BPSEEK(0);
  
  for(i=0;i<argc;++i) inputArgv.push_back(string(argv[i]));
  mycommand=inputArgv[0];
  for(i=1;i<(int)inputArgv.size();++i) mycommand+=" "+inputArgv[i];
  for(i=1;i<(int)inputArgv.size();++i) {
    if ( inputArgv[i]=="-b" ) {  // input bam file
      bamFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-f" ) {  // input reference file
      fastaFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-cnv" ) {   // load cnv file
      cnvFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    } 
    if ( inputArgv[i]=="-o" ) {  // output file
      outputFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-l" ) {  // output file
      stringstream(inputArgv[i+1]) >> readlength;
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-M" ) {  // show all parts of reads or just mapped
      SHOWMAPPED=true;
      inputArgv[i]="";
    }
    if ( inputArgv[i]=="-ee" ) {   // number of bases at edge that can be ignored
      ERR_EDGE=atoi(inputArgv[i+1].c_str());
      inputArgv[i]="";
      inputArgv[i+1]="";
    } 
    if ( inputArgv[i]=="-me" ) {   // err_rate when fuzzy matching
      ERR_MATCH=atoi(inputArgv[i+1].c_str());
      inputArgv[i]="";
      inputArgv[i+1]="";
    } 
    if ( inputArgv[i]=="-s" ) {   // number of softclipped bases
      MINIMUMS=atoi(inputArgv[i+1].c_str());
      inputArgv[i]="";
      inputArgv[i+1]="";
    } 
    if ( inputArgv[i]=="-q" ) {   // min mapq score
      MINIMUMMAPQ=atoi(inputArgv[i+1].c_str());
      inputArgv[i]="";
      inputArgv[i+1]="";
    } 
    if ( inputArgv[i]=="-Q" ) {   // min mapq score
      MINIMUMBASEQ=atoi(inputArgv[i+1].c_str());
      inputArgv[i]="";
      inputArgv[i+1]="";
    } 
    if ( inputArgv[i]=="-2" ) {   // min mapq score
      NOSECONDARY=false;
      inputArgv[i]="";
    } 
  }
  if ( bamFile=="" ) {
    cerr << "Need bam file\n------------------\n";
    exit( usage_plot_reads(argc, argv) );
  }
  if ( cnvFile=="" ) {
    cerr << "Need cnv file\n------------------\n";
    exit( usage_plot_reads(argc, argv) );
  }
  if ( fastaFile=="" ) {
    cerr << "Need reference file\n------------------\n";
    exit( usage_plot_reads(argc, argv) );
  }
  for(i=1;i<(int)inputArgv.size();++i) 
    if ( inputArgv[i]!="" ) cerr << "unknown argument:\t" << inputArgv[i] << endl;
  for(i=1;i<(int)inputArgv.size();++i) 
    if ( inputArgv[i]!="" ) exit( usage_plot_reads(argc, argv) );
  
  streambuf* sbuf = cout.rdbuf();
  ofstream FOUT;
  if ( outputFile != "STDOUT" ) {
    FOUT.open(outputFile.c_str());
    cout.rdbuf(FOUT.rdbuf());
  }
  
  string comment="#Command Line     : " + mycommand + "\n" +
    "#Input bamfile    : " + bamFile + "\n" +
    "#cnvFile          : " + cnvFile + "\n" +
    "#FASTAfile        : " + fastaFile + "\n" +
    //    "#Num of S bases   : " + to_string(MINIMUMS) + "\n" +
    //    "#Minumum overlap  : " + to_string(MINIMUMOVERLAP) + "\n" +
    //    "#Minimum mapq     : " + to_string(MINIMUMMAPQ) + "\n" +
    //    "#3' bases ignored : " + to_string(ERR_EDGE) + "\n" +
    //    "#Allowed mismatch : " + to_string(ERR_MATCH) + "\n" +
    //    "#Maximum distance : " + to_string(MAXDISTANCE) + "\n" +
    "#Output           : " + outputFile + "\n";
  cerr << comment << endl;
  
  samfile_t *fp_in = NULL;
  fp_in = samopen(bamFile.c_str(), "rb", 0);
  bam1_t *b=NULL; b=bam_init1();
  bam_index_t *bamidx=NULL;
  bamidx = bam_index_load(bamFile.c_str()); // load BAM index
  bam_iter_t iter=0;
  int ref=0;
  
  for(int i=0;i<fp_in->header->n_targets;++i) {
    BREAKSEEK_st::TARGET.push_back( string(fp_in->header->target_name[i]) );
    BREAKSEEK_st::RNAMEINT[ string(fp_in->header->target_name[i]) ]=i;
  }
  if ( BREAKSEEK_st::TARGET.size() == 0 ) 
    cerr << "#no header in bam file? " << bamFile << endl;
  
  vector<RSAI_st> MSREAD(0);         // info used to find the read
  vector<RSAI_st> SMREAD(0);
  RSAI_st iread;
  string FASTA;
  string fastaname="";
  
  SEQBUFFER.reserve(BUFFERSIZE+1); // buffer for reads with MS type of CIGAR
  
  load_cnv_from_file(BPSEEK, cnvFile);
  RD_PAIR_filters(BPSEEK, bamFile);
  //exit(0);

  /*
  for(size_t icnv=0; icnv<BPSEEK.size(); ++icnv) {
    int P5,P3;
    ref=BPSEEK[icnv].tid;
    P5=BPSEEK[icnv].P1;
    P3=BPSEEK[icnv].P2;
    if ( BPSEEK[icnv].T == 'A' ) swap(P5,P3);
    cerr << ref << "\t" << BPSEEK[icnv].RNAME() << "\t"  << P5 << "\t" << P3 << "\t" << BPSEEK[icnv].T << endl;
  }
  */
  
  for(size_t icnv=0; icnv<BPSEEK.size(); ++icnv) {
    
    int P5,P3;
    ref=BPSEEK[icnv].tid;
    P5=BPSEEK[icnv].P1;
    P3=BPSEEK[icnv].P2;
    
    if ( BPSEEK[icnv].T == 'A' ) swap(P5,P3);
    
    cerr << ref << "\t" << BPSEEK[icnv].RNAME() << "\t"  << P5 << "\t" << P3 << "\t" << BPSEEK[icnv].T << endl;
    
    /* load reference if necessary and process data for previous chromosome*/
    RNAME=BPSEEK[icnv].RNAME();
    if ( RNAME != fastaname ) {
      read_fasta(fastaFile,RNAME,FASTA);
      fastaname=RNAME;
      cerr << "#" << fastaFile << "\t" 
	   << fastaname << "\t" 
	   << FASTA.size() << " loaded" << endl;
      if ( FASTA.size()<1 ) {
	cerr << "something loading reference" << endl;
	exit(0);
      }
    }
    
    MSREAD.clear();
    SMREAD.clear();
    SEQBUFFER="";
    
    iter = bam_iter_query(bamidx, ref, P5-100, P5+100);
    while ( bam_iter_read(fp_in->x.bam, iter, b) > 0 ) {
      POSCIGAR_st bm;
      resolve_cigar_pos(b, bm);  
      if ( bm.pos==0 ) continue;
      if ( (int)bm.cop[0] > P5 ) continue;
      if ( (int)bm.cop.back() + (int)bm.nop.back()-1 < P5 ) continue;
      push_alignment_to_string(bm, SEQBUFFER, MSREAD);
    }
    
    iter = bam_iter_query(bamidx, ref, P3-100, P3+100);
    while ( bam_iter_read(fp_in->x.bam, iter, b) > 0 ) {
      POSCIGAR_st bm;
      resolve_cigar_pos(b, bm);  
      if ( bm.pos==0 ) continue;
      if ( (int)bm.cop[0] > P3 ) continue;
      if ( (int)bm.cop.back()+(int)bm.nop.back()-1 < P3 ) continue;
      push_alignment_to_string(bm, SEQBUFFER, SMREAD);
    }
    
    string readMS, MSCIGAR, readSM, SMCIGAR, MATCHED;
    int p1, p2;
    POSCIGAR_st bm1;
    vector<bool> is_MS_matched(MSREAD.size(), false);
    vector<bool> is_SM_matched(SMREAD.size(), false);
    for(i=0;i<(int)MSREAD.size();++i) {
      readMS=SEQBUFFER.substr(MSREAD[i].p1, MSREAD[i].len);
      MSCIGAR=SEQBUFFER.substr(MSREAD[i].p1+MSREAD[i].len, MSREAD[i].len_cigar);
      
      for(k=0;k<(int)SMREAD.size();++k) {
	readSM=SEQBUFFER.substr(SMREAD[k].p1, SMREAD[k].len);
	SMCIGAR=SEQBUFFER.substr(SMREAD[k].p1+SMREAD[k].len, SMREAD[k].len_cigar);
	
	bool is_overlap = false;
	is_overlap = is_MS_SM_overlap_fuzzy_lite(readMS, readSM, MATCHED, p1, p2); 
	if ( ! is_overlap ) 
	  is_overlap = is_MS_SM_overlap_fuzzy_lite(readSM, readMS, MATCHED, p1, p2); 
	is_MS_matched[i] = is_MS_matched[i] || is_overlap;
	is_SM_matched[k] = is_SM_matched[k] || is_overlap;
	
      }
    }
    int is_MS_matched_count=0;
    int is_SM_matched_count=0;
    for(i=0;i<(int)MSREAD.size();++i) is_MS_matched_count+=is_MS_matched[i];
    for(i=0;i<(int)SMREAD.size();++i) is_SM_matched_count+=is_SM_matched[i];
    cerr << icnv << "\t" << ref << "\t" << P5 << "\t" << P3 << "\n"
	 << "P5 reads " << MSREAD.size() << "\t" << is_MS_matched_count << "\n"
	 << "P3 reads " << SMREAD.size() << "\t" << is_SM_matched_count << "\n"
	 << endl;
    
    int reflen=MSREAD[0].len;
    for(i=0;i<(int)MSREAD.size();++i) if ( MSREAD[i].len>reflen ) reflen=MSREAD[i].len;
    
    string left100=FASTA.substr(P5-reflen,reflen); //start=P5-reflen+1
    string right100=FASTA.substr(P3-1,reflen);  //start=P3
    string divider=string(10,'-');
    string refheader=left100+divider+right100;
    string leftheader=BPSEEK[icnv].RNAME() + " "  + to_string(P5);
    string rightheader=to_string(P3)+ " " + to_string(BPSEEK[icnv].T);
    string header="";
    
    header = string(reflen-leftheader.size() + to_string(P5).size()-1,' ')+leftheader 
      + string(divider.size()-to_string(P5).size()+1, ' ') + rightheader;
    
    if ( header.size() < refheader.size() ) header += string(refheader.size()-header.size() ,' ');
    //cout << header << endl;
    //cout << refheader << endl;
    
    vector<string> MS_left(MSREAD.size()), MS_right(MSREAD.size());
    //goto PRINTRIGHT;
    for(i=0;i<(int)MSREAD.size();++i) {
      readMS=SEQBUFFER.substr(MSREAD[i].p1, MSREAD[i].len);
      MSCIGAR=SEQBUFFER.substr(MSREAD[i].p1+MSREAD[i].len, MSREAD[i].len_cigar);
      resolve_cigar_string(MSREAD[i].pos, MSCIGAR, bm1);
      
      string leftside="",rightside="";
      int pad=bm1.cop[0]-P5+reflen-1;
      
      int start;
      if ( pad>=0 ) { 
	//cout << string(pad,' '); 
	leftside+=string(pad,' '); 
	start=0; 
      }
      else start=-pad;
      string o_string,l_string,r_string;
      for(int k=0;k<(int)bm1.op.size();++k) {
	o_string=l_string=r_string="";
	if ( start >= (int)bm1.qop[k]+(int)bm1.nop[k] ) continue;
	
	o_string=readMS.substr(bm1.qop[k], bm1.nop[k]);
	
	if ( start>(int)bm1.qop[k] ) o_string=o_string.substr(start-bm1.qop[k]);
	
	if ( (int)bm1.cop[k] > P5 ) { 
	  if (  bm1.op[k] != BAM_CDEL &&
		bm1.op[k] != BAM_CREF_SKIP ) r_string = o_string; 
	  if ( ( bm1.op[k] == BAM_CINS ) && ( !is_MS_matched[i] ) ) r_string="";
	  rightside += r_string;
	  continue; 
	}
	
	if ( (int)bm1.cop[k]+(int)bm1.nop[k]>P5 ) {
	  l_string = o_string.substr(0, o_string.size()-bm1.cop[k]-bm1.nop[k]+P5+1);
	  r_string = o_string.substr(o_string.size()-bm1.cop[k]-bm1.nop[k]+P5+1);
	}
	else {
	  l_string = o_string;
	}
	
	if (  bm1.op[k] == BAM_CDEL ||
	      bm1.op[k] == BAM_CREF_SKIP ) {
	  l_string=std::string( l_string.size(), '_');
	  r_string="";
	}
	if ( bm1.op[k] == BAM_CINS ) l_string="";
	//cout << l_string ;
	leftside += l_string ;
	rightside += r_string;
      }
      //cout << "\t" << MSCIGAR;
      //cout << "\n";
      MS_left[i]=leftside;
      if ( rightside.size() < right100.size() ) rightside+=string(right100.size()-rightside.size(), ' ');
      MS_right[i]=rightside;
      //cout << leftside << divider << rightside << endl;
    }
    
    //cout << left100 << divider << right100 << "\n";
    vector<string> SM_left(SMREAD.size()), SM_right(SMREAD.size());
    for(i=0;i<(int)SMREAD.size();++i) {
      
      //cout << string(left100.size(), ' ');
      //cout << string(divider.size(), ' ');
      
      readSM=SEQBUFFER.substr(SMREAD[i].p1, SMREAD[i].len);
      SMCIGAR=SEQBUFFER.substr(SMREAD[i].p1+SMREAD[i].len, SMREAD[i].len_cigar);
      resolve_cigar_string(SMREAD[i].pos, SMCIGAR, bm1);
      if ( bm1.pos==0 ) continue;
      
      string leftside="",rightside="";
      string o_string,l_string,r_string;
      for(int k=0;k<(int)bm1.op.size();++k) {
	o_string=readSM.substr(bm1.qop[k], bm1.nop[k]);
	l_string=r_string="";
	
	if ( (int)bm1.cop[k]+(int)bm1.nop[k] <= P3 ) { 
	  if (  bm1.op[k] != BAM_CDEL &&
		bm1.op[k] != BAM_CREF_SKIP ) l_string = o_string; 
	  if ( bm1.op[k] == BAM_CINS && (!is_SM_matched[i]) ) l_string="";
	  leftside += l_string;
	  continue;
	}
	
	if ( (int)bm1.cop[k]>=P3 ) { 
	  l_string = "";
	  r_string = o_string;
	}
	else {
	  l_string = o_string.substr(0, o_string.size()-bm1.cop[k]-bm1.nop[k]+P3);
	  r_string = o_string.substr(o_string.size()-bm1.cop[k]-bm1.nop[k]+P3);
	}
	
	if (  bm1.op[k] == BAM_CDEL ||
	      bm1.op[k] == BAM_CREF_SKIP ) {
	  l_string="";
	  r_string=std::string( l_string.size(), '_');
	}
	
	if ( bm1.op[k] == BAM_CINS ) r_string="";
	//cout << r_string ;
	leftside += l_string;
	rightside += r_string ;
      }
      //cout << "\t" << SMCIGAR;
      //cout << "\n";
      
      if ( leftside.size() < left100.size() ) leftside=string(left100.size()-leftside.size(), ' ')+leftside;
      if ( rightside.size() < right100.size() ) rightside+=string(right100.size()-rightside.size(), ' ');
      //cout << leftside << string(divider.size(), ' ') << rightside << endl;
      
      SM_left[i]=leftside;
      SM_right[i]=rightside;
    }
    
    //! cut to fit screen size
    if ( readlength>0 && readlength<reflen ) {
      left100=left100.substr(left100.size()-readlength);
      right100=right100.substr(0,readlength);
      refheader=left100+divider+right100;
      header = header.substr(reflen-readlength);
      header = header.substr(0,header.size()-reflen+readlength);
      for(i=0;i<(int)MSREAD.size();++i) { 
	MS_left[i]=MS_left[i].substr(MS_left[i].size()-readlength);
	MS_right[i]=MS_right[i].substr(0,readlength);
      }
      for(i=0;i<(int)SMREAD.size();++i) { 
	SM_left[i]=SM_left[i].substr(SM_left[i].size()-readlength);
	SM_right[i]=SM_right[i].substr(0,readlength);
      }
    }
    
    //! begin print cnv
    string info=string( (header.size()-BPSEEK[icnv].TAG.size())/2, ' ');
    info+=BPSEEK[icnv].TAG;
    info+=string(header.size()-info.size(), ' ');
    cout << "-----" << string(refheader.size(),'-') << "-----" << endl;
    cout << "POS |" << header << "|    "  << endl;
    cout << "STAT|" << info << "|    "  << endl; 
    cout << "REF |" << refheader << "|    "  << endl;
    cout << "----|" << string(refheader.size(),'-') << "|----" << endl;
    //! print the full matched reads with breaks at the break points
    i=0;k=SMREAD.size()-1;
    while( i<(int)MSREAD.size() || k>=0 ) {
      
      string qual5="  ",qual3="  ";
      if ( i<(int)MSREAD.size() ) {
	while ( ! is_MS_matched[i] ) { 
	  ++i;
	  if ( i==(int)MSREAD.size() ) break;
	}
      }
      if ( i<(int)MSREAD.size() ) {
	if ( ! is_MS_matched[i] ) { cerr << "loop error\n"; exit(0); }
	if ( i<(int)MSREAD.size() ) {
	  if ( MSREAD[i].q1>99 ) MSREAD[i].q1=99;
	  qual5=to_string(MSREAD[i].q1);
	  if ( qual5.length()<2 ) qual5 = string(2-qual5.length(), ' ') +qual5;
	}
	cout << qual5 << "@5|" << MS_left[i] << divider  <<  MS_right[i] 
	     << "|    " << endl;
	//cout << "@5' |" << MS_left[i] << divider  <<  MS_right[i] 
	//   << "|    " << endl;
	++i;
      }
      
      if ( k>=0 ) {
	while ( ! is_SM_matched[k] ) { 
	  --k;
	  if ( k<0 ) break;
	}
      }
      if ( k>=0 ) {
	if ( ! is_SM_matched[k] ) { cerr << "loop error\n"; exit(0); }
	if ( k<(int)SMREAD.size() ) {
	  if ( SMREAD[k].q1>99 ) SMREAD[k].q1=99;
	  qual3=to_string(SMREAD[k].q1);
	  if ( qual3.length()<2 ) qual3 = qual3 + string(2-qual3.length(), ' ') ;
	}
	cout << "    |" << SM_left[k] << divider  <<  SM_right[k] 
	     << "|3@" << qual3 << endl;
	//   << "| @3'" << endl;
	--k;
      }
      
    }
    
    //! print other reads
    if ( SHOWMAPPED ) {
      cout << "....|" << string(header.size(), '.') << "|...." << endl;
      i=0;k=0;
      while( i<(int)MSREAD.size() || k<(int)SMREAD.size() ) {
	string leftside="",rightside="";
	
	if ( i<(int)MSREAD.size() ) {
	  while ( is_MS_matched[i] ) { 
	    ++i;
	    if ( i==(int)MSREAD.size() ) break;
	  }
	}
	if ( i<(int)MSREAD.size() ) leftside=MS_left[i];
	else leftside=string(reflen,' ');
	++i;
	
	if ( k<(int)SMREAD.size() ) {
	  while ( is_SM_matched[k] ) { 
	    ++k;
	    if ( k==(int)SMREAD.size() ) break;
	  }
	}
	if ( k<(int)SMREAD.size() ) rightside=SM_right[k];
	else rightside=string(reflen,' ');
	++k;
	
	if ( leftside==string(reflen,' ') && rightside==string(reflen,' ') ) continue;
	
	if ( readlength>0 && readlength<reflen ) {
	  leftside=leftside.substr(leftside.size()-readlength);
	  rightside=rightside.substr(0,readlength);
	}
	
	string qual5="  ",qual3="  ";
	if ( i<(int)MSREAD.size() ) {
	  if ( MSREAD[i].q1>99 ) MSREAD[i].q1=99;
	  qual5=to_string(MSREAD[i].q1);
	  if ( qual5.length()<2 ) qual5 = string(2-qual5.length(), ' ') +qual5;
	}
	if ( k<(int)SMREAD.size() ) {
	  if ( SMREAD[k].q1>99 ) SMREAD[k].q1=99;
	  qual3=to_string(SMREAD[k].q1);
	  if ( qual3.length()<2 ) qual3 = qual3 + string(2-qual3.length(), ' ') ;
	}
	cout << qual5 << "@5|" << leftside << string(divider.size(), ' ')  <<  rightside 
	     << "|3@" << qual3 << endl;
	//cout << "@5' |" << leftside << string(divider.size(), ' ')  <<  rightside 
	//   << "| @3'" << endl;
      }
    }
    else {
      
      leftheader=BPSEEK[icnv].RNAME() + " "  + to_string(P5);
      rightheader=to_string(P5+1);
      header = string(reflen-leftheader.size() + to_string(P5).size()-1,' ')+leftheader 
	+ string(divider.size()-to_string(P5).size()+1, ' ') 
	+ rightheader+string(reflen-rightheader.size() ,' ');
      if ( readlength>0 && readlength<reflen ) {
	header = header.substr(reflen-readlength);
	header = header.substr(0,header.size()-reflen+readlength);
      }
      for(i=0;i<(int)header.size();++i) if ( header[i]==' ' ) header[i]='.';
      cout << "....|" << header << "|...." << endl;
      string qual5="  ",qual3="  ";
      for(i=0;i<(int)MSREAD.size();++i) {
	if ( is_MS_matched[i] ) continue;
	if ( MSREAD[i].q1>99 ) MSREAD[i].q1=99;
	qual5=to_string(MSREAD[i].q1);
	if ( qual5.length()<2 ) qual5 = string(2-qual5.length(), ' ') +qual5;
	cout << qual5 << "@5|" << MS_left[i] << string(divider.size(), ' ')  <<  MS_right[i] 
	     << "|    " << endl;
	//cout << "@5' |" << MS_left[i] << string(divider.size(), ' ')  <<  MS_right[i] 
	//   << "|    " << endl;
      }
      
      leftheader=BPSEEK[icnv].RNAME() + " "  + to_string(P3-1);
      rightheader=to_string(P3);
      header = string(reflen-leftheader.size() + to_string(P3-1).size()-1,' ')+leftheader 
	+ string(divider.size()-to_string(P3-1).size()+1, ' ') 
	+ rightheader+string(reflen-rightheader.size() ,' ');
      if ( readlength>0 && readlength<reflen ) {
	header = header.substr(reflen-readlength);
	header = header.substr(0,header.size()-reflen+readlength);
      }
      for(i=0;i<(int)header.size();++i) if ( header[i]==' ' ) header[i]='.';
      cout << "....|" << header << "|...." << endl;
      for(i=SMREAD.size()-1;i>=0;--i) {
	if ( is_SM_matched[i] ) continue;
	if ( SMREAD[i].q1>99 ) SMREAD[i].q1=99;
	qual3=to_string(SMREAD[i].q1);
	if ( qual3.length()<2 ) qual3 = qual3 + string(2-qual3.length(), ' ') ;
	cout << "    |" << SM_left[i] << string(divider.size(), ' ')  <<  SM_right[i] 
	     << "|3@" << qual3 << endl;
	//     << "| @3'" << endl;
      }
    }
    
    cout << "-----" << string(refheader.size(),'-') << "-----" << endl;
    cout << endl << endl;
  }
  
  samclose(fp_in);
  bam_destroy1(b);
  bam_iter_destroy(iter);
  bam_index_destroy(bamidx);
  
  FOUT.close();
  cout.rdbuf(sbuf);
  
  exit(0);
}

