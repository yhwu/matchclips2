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

int main(int argc, char* argv[])
{
  
  string mycommand="";
  size_t i, k, seed=0;
  string QNAME,FLAG,RNAME="",POS,MAPQ,CIGAR,MRNM,MPOS,ISIZE,SEQ,QUAL,OPT;
  
  string fastaFile="",outputFile="STDOUT",logFile="CIGARPOS.log";
  string cnvfile="";
  vector<string> inputArgv;
  
  size_t nread=0;
  int lcnv=1000, lrange=5000, lread=100, readX=10;
  int insert=500, idis=50;
  bool appendmode=false;
  int POS1=-1, POS2=-1;
  double errormultiplier=0.001;
  string wgsim_command="";
  
  for(i=0;i<(size_t) argc;++i) inputArgv.push_back(string(argv[i]));
  mycommand=inputArgv[0];
  for(i=1;i<inputArgv.size();++i) mycommand+=" "+inputArgv[i];
  SEQ=inputArgv[1];
  
  int maxRepeat=6;
  int maxLength=0;
  string rep="";
  string maxrep="";
  for(i=0;i<SEQ.size();++i) {
    for(k=1;k<=maxRepeat;++k) {
      if ( i+k>=SEQ.size() ) continue;
      //cerr << i << "\t" << k << endl;
      rep=SEQ.substr(i,k);
      
      int ndiff=0;
      int lRepeat=0;
      int idx=i+k;
      while(ndiff<=1 && idx+k<SEQ.size() ) {
	for(int i1=idx, i2=0; i1<idx+k; ++i1,++i2)
	  if ( rep[i2]!=SEQ[i1] ) ++ndiff;
	if ( ndiff<=1 ) lRepeat+=k;
	idx+=k;
      }
      
      if ( lRepeat > maxLength ) {
	maxLength=lRepeat;
	maxrep=rep;
      }
    }
  }
  
  cerr << "#Command Line     : " << mycommand << "\n"
       << "#SEQ              : " << SEQ << "\n"
       << maxLength << "\t" << maxrep << endl
       << endl;
  
  exit(0);
}

