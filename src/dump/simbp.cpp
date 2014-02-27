#include <stdio.h>
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
#include <time.h>
using namespace std;

/**** user headers ****/
#include "readref.h"
template <class T>
inline std::string to_string (const T& t)
{
  std::stringstream ss;
  ss << t;
  return ss.str();
}

extern "C" {
#include "wgsim.h"
}


/* Simple normal random number generator, copied from genran.c;
   Simple normal random number generator, copied from wgsim.c */
double unifrand() { return rand() / double(RAND_MAX); }
double ran_normal()
{ 
  static int iset = 0; 
  static double gset; 
  double fac, rsq, v1, v2; 
  if (iset == 0) {
    do { 
      v1 = 2.0 * unifrand() - 1.0;
      v2 = 2.0 * unifrand() - 1.0; 
      rsq = v1 * v1 + v2 * v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0 * log(rsq) / rsq); 
    gset = v1 * fac; 
    iset = 1;
    return v2 * fac;
  } else {
    iset = 0;
    return gset;
  }
}


int usage_simbp(int argc, char* argv[]) {
  string me=string(argv[0]);
  me=me.substr(me.find("/")+1);
  cerr << me << ":    simulate reads around CNVs within certain range. Half of the reads are selected from the normal chromosome and the other half is from the chromosome with CNV. \n"
       << "\nUsage:\n" 
       << "  " << me << " <options> -f ref.fasta \n"
       << "\nOptions:\n"
       << "  -cnv    STR  read CNV from file STR \n"
       << "  -p5     INT  5 prime position, INT=RANDOM \n"
       << "  -p3     INT  3 prime position, INT=p5+lcnv \n"
       << "  -lcnv   INT  length of CNV, can be negative, INT=1000 \n"
       << "  -lrange INT  range in which reads are simulated, INT=4000 \n"
       << "  -lrd    INT  length of reads, INT=100 \n"
       << "  -X      INT  read depth, INT=10 \n"
       << "  -n      INT  number of reads, if not given, determined by -X \n"
       << "  -I      INT  insert size, INT=500 \n"
       << "  -m      FLT  mutation rate, FLT=0.001 \n"
       << "  -Idis   INT  insert size standard deviation, INT=50 \n"
       << "  -s      INT  seed for rand, INT=time(NULL) \n"
       << "  -o      STR  outputfile, STR=STDOUT \n"
       << "  -a           append output, default=false \n"
       << "\nExamples:\n"
       << "  " << me << " -f human_g1k_v37.fasta -cnv NA12878hg19.txt -o tmp \n"
       << "\nNote:\n"
       << "  The reads around the CNV are randomly simulated from the string\n"
       << "     RNAME(p5-lrange/2+1 ~ p5) + RNAME(p3 ~ p3+lrange/2-1)\n"
       << endl;
  return 0;
}

struct CNV_ST {
  string RNAME;
  int P5;
  int P3;
  string TYPE;
  CNV_ST():RNAME(""),
	    P5(-1),
	    P3(-1),
	    TYPE(""){};
};

int simbp(int argc, char* argv[])
{
  
  string mycommand="";
  size_t i, seed=0;
  string QNAME,FLAG,RNAME="",POS,MAPQ,CIGAR,MRNM,MPOS,ISIZE,SEQ,QUAL,OPT;
  
  string fastaFile="",outputFile="STDOUT",logFile="CIGARPOS.log";
  string cnvfile="";
  vector<string> inputArgv;
  
  size_t nread=0;
  int lcnv=1000, lrange=5000, lread=100, readX=10;
  int insert=500, idis=50;
  bool appendmode=false;
  int POS1=-1, POS2=-1;
  string wgsim_command="";
  double er=0.02;
  double mr=0.01;
  double ir=0.15;
  double ie=0.3;
  double maxnr=0.05;

  RNAME="19";
  
  for(i=0;i<(size_t) argc;++i) inputArgv.push_back(string(argv[i]));
  mycommand=inputArgv[0];
  for(i=1;i<inputArgv.size();++i) mycommand+=" "+inputArgv[i];
  for(i=1;i<inputArgv.size();++i) {
    if ( inputArgv[i]=="-f" ) {  // input reference file
      fastaFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-r" ) {  // RNAME
      RNAME=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-p5" ) {  // number of reads needed
      istringstream( inputArgv[i+1] ) >> POS1;
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-p3" ) {  // number of reads needed
      istringstream( inputArgv[i+1] ) >> POS2;
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-cnv" ) {  // distance between break point
      cnvfile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-lcnv" ) {  // distance between break point
      lcnv=atoi(inputArgv[i+1].c_str());
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-lrange" ) {  // length of CNV
      lrange=atoi(inputArgv[i+1].c_str());
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-lrd" ) {  // length of read
      lread=atoi(inputArgv[i+1].c_str());
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-X" ) {  // read depth
      readX=atoi(inputArgv[i+1].c_str());
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-m" ) {  // read depth
      istringstream( inputArgv[i+1] ) >> mr;
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-e" ) {  // read depth
      istringstream( inputArgv[i+1] ) >> er;
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-n" ) {  // number of reads needed
      istringstream( inputArgv[i+1] ) >> nread;
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-I" ) {  // insert size
      istringstream( inputArgv[i+1] ) >> insert;
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-Idis" ) {  // insert size distribution
      istringstream( inputArgv[i+1] ) >> idis;
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-s" ) {  // seed
      istringstream( inputArgv[i+1] ) >> seed;
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-o" ) {  // output file
      outputFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-a" ) {  // append mode
      appendmode=true;
      inputArgv[i]="";
    }
  }
  if ( fastaFile=="" ) exit( usage_simbp(argc, argv) );
  if ( seed==0 ) seed=time(NULL);
  logFile=cnvfile+".log";
  srand(seed);
  srand48(seed);
  
  cerr << "#Command Line     : " << mycommand << "\n"
       << "#fastafile        : " << fastaFile << "\n"
       << "#length of cnv    : " << lcnv << "\n"
       << "#length of read   : " << lread << "\n"
       << "#Read depth       : " << readX << "\n"
       << "#Number of reads  : " << nread << "\n"
       << "#wgsim command    : " << wgsim_command << "\n"
       << "#seed             : " << seed << "\n"
       << endl;
  
  vector<CNV_ST> CNV(0);
  CNV_ST icnv;
  
  if ( cnvfile!="" ) {
    ifstream FIN(cnvfile.c_str());
    while ( !FIN.eof() ) {
      string tmps;
      getline(FIN,tmps);
      if ( tmps.length()<3 ) continue;
      if ( tmps[0]=='#' ) continue;
      
      
      istringstream iss(tmps);
      
      iss >> icnv.RNAME
	  >> icnv.P5 
	  >> icnv.P3
	  >> icnv.TYPE;
      CNV.push_back(icnv);
    }
    FIN.close();
  }
  else {
    icnv.RNAME=RNAME;
    icnv.P5=POS1;
    icnv.P3=POS2;
    icnv.TYPE="DEL";
    CNV.push_back(icnv);
  }
  
  string FASTA;
  string fastaname="";
  
  int nDeleted=1;
  while ( nDeleted>0 ) {
    vector<bool> iskept(CNV.size(), true);
    for(size_t k=0;k<CNV.size();++k) {
      bool pass=true;
      if ( k>0 ) {
	int com_beg=max(CNV[k].P5, CNV[k-1].P5);
	int com_end=min(CNV[k].P3, CNV[k-1].P3);
	if ( com_beg<=com_end ) { 
	  cerr << "everlap\t" 
	       << CNV[k-1].RNAME << "\t" 
	       << CNV[k-1].P5 << "\t" << CNV[k-1].P3 << "\t" << CNV[k-1].TYPE << "\t"
	       << CNV[k].P5 << "\t" << CNV[k].P3 << "\t" << CNV[k].TYPE << "\t"
	     << endl;
	  pass=false;
	}
	else if ( com_beg<=com_end+lrange/2 ) {
	  cerr << "close\t" 
	       << CNV[k-1].RNAME << "\t" 
	       << CNV[k-1].P5 << "\t" << CNV[k-1].P3 << "\t" << CNV[k-1].TYPE << "\t"
	       << CNV[k].P5 << "\t" << CNV[k].P3 << "\t" << CNV[k].TYPE << "\t"
	     << endl;
	  pass=false;
	}
      }
      if ( pass==false ) {
	if ( abs(CNV[k].P5-CNV[k].P3) < abs(CNV[k-1].P5-CNV[k-1].P3) ) {
	  iskept[k]=true;
	  iskept[k-1]=false;
	}
	else {
	  iskept[k]=false;
	  iskept[k-1]=true;
	}
      }
      /*
	WRITETOFILE:
	if ( !pass ) continue;
	cout << CNV[k].RNAME << "\t" 
	<< CNV[k].P5 << "\t" 
	<< CNV[k].P3 << "\t" 
	   << CNV[k].TYPE << endl;
      */
    }
    vector<CNV_ST> tmpCNV(0);
    nDeleted=0;
    for(size_t k=0;k<CNV.size();++k) {
      if ( iskept[k] ) tmpCNV.push_back(CNV[k]);
      else ++nDeleted;
    }
    CNV=tmpCNV;
    cerr << CNV.size() << endl;
  }
  //  exit(0);
  
  ofstream FLOG(logFile.c_str());
  
  for(size_t k=0;k<CNV.size();++k) {
    
    RNAME=CNV[k].RNAME;
    POS1=CNV[k].P5;
    POS2=CNV[k].P3;
    if ( CNV[k].TYPE=="DUP" ) swap(POS1, POS2);
    
    cerr << k+1 << "\t" << RNAME << "\t" << POS1 << "\t" << POS2 << endl;
    
    if ( RNAME != fastaname ) {
      string RNAMEtmp=RNAME;
      if ( RNAME.find("chr") != string::npos ) RNAMEtmp=RNAME.substr(3); 
      read_fasta(fastaFile,RNAMEtmp,FASTA);
      fastaname=RNAME;
      cerr << "Load ref  "
	   << fastaFile << "\t"
	   << fastaname << "\t"
	   << FASTA.size() << endl;
      if ( FASTA=="" ) {
	cerr << "reference " << RNAME << " not read sucessfully, exit\n"
	     << endl;
	exit(0);
      }
    }
    
    // fasta around cnv 
    int lrange2=lrange/2;
    int fasta_beg=min(POS1,POS2)-lrange2;
    int fasta_end=max(POS1,POS2)+lrange2;
    string FASTA1,FASTA0;
    if ( POS1>0 && POS2>0 ) {
      if ( POS2+lrange > (int)FASTA.size() || POS2-lrange < 0 ||
	   POS1+lrange > (int)FASTA.size() || POS1-lrange < 0 ) {
	cerr << "bad positions\t" << POS1 << "\t" << POS2 << endl;
	exit(0);
      }
      
      if ( fasta_beg < 1 ) fasta_beg=1; 
      if ( fasta_end >= (int)FASTA.size() ) fasta_end=FASTA.size()-1; 
      FASTA0=FASTA.substr(fasta_beg-1, fasta_end-fasta_beg+1);
      FASTA1=FASTA.substr(fasta_beg-1, POS1-fasta_beg+1)+
	FASTA.substr(POS2-1, fasta_end-POS2+1);
    }
    else {
      do {
	POS1=unifrand()*(FASTA.size()-lcnv-2)+lrange/2;
	POS2=POS1+lcnv-1;
	while ( FASTA[POS1]==FASTA[POS2-1] ) { POS1++; POS2++; }
	FASTA1 = FASTA.substr(POS1-lrange/2+1-1, lrange/2)
	  +FASTA.substr(POS2-1, lrange/2);
      } 
      while ( FASTA1.find("N") != string::npos || 
	      POS1-lrange<0 || 
	      POS2+lrange>(int)FASTA.size() )  ;
    }
    
    int Ncount=0;
    for(int j=0;j<(int)FASTA0.size();++j) 
      if (FASTA0[j]=='N') Ncount++;
    if (Ncount>2) continue;
    Ncount=0;
    //for(int j=0;j<(int)FASTA1.size();++j) 
    //  if (FASTA1[j]=='N') Ncount++;
    //if (Ncount>2) continue;
    
    FLOG << CNV[k].RNAME << "\t" 
	 << CNV[k].P5 << "\t"
	 << CNV[k].P3 << "\t"
	 << CNV[k].TYPE << endl;
    
    cerr << "#break points: " 
	 << POS1 << " " << FASTA[POS1-1] << "\t"
	 << POS2 << " " << FASTA[POS2-1] << "\n"
	 << "#Ends of two segments: " 
	 << fasta_beg << " " << FASTA[fasta_beg-1] << "\t"
	 << POS1 << " " << FASTA[POS1-1] << "\t"
	 << POS2 << " " << FASTA[POS2-1] << "\t"
	 << fasta_end << " " << FASTA[fasta_end-1] 
	 << endl;
    
    
    string tmpcnvfasta=outputFile+".cnv.fasta";
    ofstream FOUT(tmpcnvfasta.c_str());
    string header="";
    header=">"+RNAME+" dna:" + 
      RNAME+":"+to_string(fasta_beg)+"-"+to_string(POS1) + "+"
      +
      RNAME+":"+to_string(POS2)+"-"+to_string(fasta_end) +
      " chromosome:SIM:1:1:"+to_string(FASTA1.length()) + ":1";    
    FOUT << header << "\n";
    for(i=0;i<FASTA1.length();++i) {
      FOUT << FASTA1[i];
      if ( (i+1)%60 == 0 ) FOUT << "\n";
    }
    if ( i%60 != 0 ) FOUT << "\n";
    FOUT.close();
    cerr << "#temp cnv fasta written to " << tmpcnvfasta << endl;
    
    string tmpfasta=outputFile+".fasta";
    FOUT.open(tmpfasta.c_str());
    header=">"+RNAME+" dna:" + 
      RNAME+":"+to_string(fasta_beg)+"-"+to_string(fasta_end)
      +
      " chromosome:SIM:1:1:"+to_string(FASTA0.length()) + ":1";
    FOUT << header << "\n";
    for(i=0;i<FASTA0.length();++i) {
      FOUT << FASTA0[i];
      if ( (i+1)%60 == 0 ) FOUT << "\n";
    }
    if ( i%60 != 0 ) FOUT << "\n";
    FOUT.close();
    cerr << "#temp fasta written to " << tmpfasta << endl;
    
    
    
    FILE* fpout1, *fpout2;
    const char* mode= appendmode ? "a":"w";
    if ( k>0 ) mode="a";
    fpout1 = fopen((outputFile+"_1.fq").c_str(), mode);
    fpout2 = fopen((outputFile+"_2.fq").c_str(), mode);
    
    /* simulate reads on the cnv chromosome */
    nread=FASTA1.length()*readX/lread/4+1;
    wgsim_core(fpout1, fpout2, 
	       tmpcnvfasta.c_str(), 0, nread, insert, idis, lread, lread, POS1, POS2, er,mr,ir,ie,maxnr );  
    
    
    fclose(fpout1); fclose(fpout2);
    
    /* simulate reads on the normal chromosome */
    mode="a";
    fpout1 = fopen((outputFile+"_1.fq").c_str(), mode);
    fpout2 = fopen((outputFile+"_2.fq").c_str(), mode);
    POS1=fasta_beg;
    POS2=0;
    nread=FASTA0.length()*readX/lread/4+1;
    wgsim_core(fpout1, fpout2, 
	       tmpfasta.c_str(), 0, nread, insert, idis, lread, lread, POS1, POS2, er,mr,ir,ie,maxnr);  
    fclose(fpout1); fclose(fpout2);
  }  
  
  return 0;
}

int main(int argc, char* argv[]) 
{
  return simbp(argc, argv);
}
