/*********************************************************************
6. CIGAR: CIGAR string. 
The CIGAR operations are given in the following table (set ‘*’ if unavailable):
    Op    BAM     Description
    M      0     alignment match (can be a sequence match or mismatch)
    I      1     insertion to the reference
    D      2     deletion from the reference
    N      3     skipped region from the reference
    S      4     soft clipping (clipped sequences present in SEQ)
    H      5     hard clipping (clipped sequences NOT present in SEQ)
    P      6     padding (silent deletion from padded reference)
    =      7     sequence match
    X      8     sequence mismatch

   H can only be present as the first and/or last operation.
   S may only have H operations between them and the ends of the CIGAR string.
   For mRNA-to-genome alignment, an N operation represents an intron. For other types of alignments, the interpretation of N is not defined.
   Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ.
from sam1.pdf
*********************************************************************/
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
  cerr << "\nUsage:\n" 
       << "  " << me << " <options> -f ref.fasta -c CHR -o PREFIX \n"
       << "\nOptions:\n"
       << "  -i      INT  insert size, INT=250 \n"
       << "  -l      INT  length of reads, INT=100 \n"
       << "  -d      INT  increament, INT=11 \n"
       << "\nExamples:\n"
       << "  " << me << " -f human_g1k_v37.fasta -c 20 -o tmp \n"
       << "\nNote:\n"
       << "        This subroutine samples reads on chromosome CHR every 11 bases\n" 
       << endl;
  return 0;
}

int simwg(int argc, char* argv[])
{
  
  string mycommand="";
  size_t i, seed=0;
  string QNAME,FLAG,RNAME="",POS,MAPQ,CIGAR,MRNM,MPOS,ISIZE,SEQ,QUAL,OPT;
  
  string fastaFile="",outputFile="STDOUT";
  string cnvfile="";
  vector<string> inputArgv;
  
  size_t nread=0;
  int lcnv=1000, lrange=5000, lread=100, dx=11;
  int insert=300, idis=50;
  bool appendmode=false;
  int POS1=-1, POS2=-1;
  double errormultiplier=0.001;
  string wgsim_command="";
  
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
    if ( inputArgv[i]=="-c" ) {  // RNAME
      RNAME=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-l" ) {  // length of read
      lread=atoi(inputArgv[i+1].c_str());
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-i" ) {  // insert size
      istringstream( inputArgv[i+1] ) >> insert;
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-d" ) {  // number of reads needed
      istringstream( inputArgv[i+1] ) >> dx;
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-o" ) {  // output file
      outputFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
  }
  if ( fastaFile=="" ) exit( usage_simbp(argc, argv) );
  
  cerr << "#Command Line     : " << mycommand << "\n"
       << "#fastafile        : " << fastaFile << "\n"
       << "#length of read   : " << lread << "\n"
       << "#Insert           : " << insert << "\n"
       << "#Dx               : " << dx << "\n"
       << "#Output           : " << outputFile << "\n"
       << endl;
  
  string FASTA;
  string fastaname="";
  
  read_fasta(fastaFile,RNAME,FASTA);
  fastaname=RNAME;
  
  string fread1=outputFile+"_1.fq";
  string fread2=outputFile+"_2.fq";
  ofstream FOUT1(fread1.c_str());
  ofstream FOUT2(fread2.c_str());
  
  string read1,read2,readq;
  readq=string(lread,'I');
  
  size_t nRead=0;
  for(int i=1;i<FASTA.size()-insert-1;i+=dx) {
    string frag=FASTA.substr(i,insert);
    read1=frag.substr(0,lread);
    read2=frag.substr(insert-lread,lread);
    std::reverse(read2.begin(), read2.end());    
    for(int k=0;k<read2.length();++k) 
      switch ( read2[k] ) {
      case 'A': read2[k]='T'; break;
      case 'T': read2[k]='A'; break;
      case 'C': read2[k]='G'; break;
      case 'G': read2[k]='C'; break;
      }
    
    int Ncount=0;
    for(int k=0;k<read1.length();++k) if ( read1[k]=='N' ) Ncount++;
    if ( Ncount>=2 ) continue;
    Ncount=0;
    for(int k=0;k<read2.length();++k) if ( read2[k]=='N' ) Ncount++;
    if ( Ncount>=2 ) continue;
    
    FOUT1 << "@" << RNAME << "_" << i+1 << "_" << i+insert-lread+1 << "/1\n" 
	  << read1 << "\n"
	  << "+\n"
	  << readq << endl;
    
    FOUT2 << "@" << RNAME << "_" << i+1 << "_" << i+insert-lread+1 << "/2\n" 
	  << read2 << "\n"
	  << "+\n"
	  << readq << endl;

    nRead+=2;
    if ( nRead%1000000==0 ) cerr << nRead << " reads simulated" << endl;
  }
  
  cerr << "#Total reads : " << nRead << endl;
  
  return 0;
}

int main(int argc, char* argv[]) 
{
  return simwg(argc, argv);
}
