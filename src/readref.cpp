#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <vector>
#include <cstdlib>
#include <unistd.h>
using namespace std;

bool read_fasta(string fastaFile, string chr, string& ref)
{
  if ( chr=="*" ) { 
    cerr << "Warning\nWarning\nWarning\n"
	 << chr << " not a valid RNAME\n"; 
    ref="";
    return false;
  } 
  
  string faiFile=fastaFile+".fai";
  ifstream FAI;
  FAI.close();
  FAI.clear();
  FAI.open(faiFile.c_str(), ifstream::in);
  int failcount=0;
  while ( !FAI && failcount<5 ) {
    usleep(3000000);
    failcount++;
    FAI.close();
    FAI.clear();
    FAI.open(faiFile.c_str(), ifstream::in);
    cerr << "[read_fasta] Failed to open " << faiFile << " " 
	 << failcount << " time" << endl; 
  }
  if ( !FAI ) {
    cerr << "[read_fasta] Index file " << faiFile << " not found\n"; 
    exit(0); 
  }
  
  long i,k;
  string ichr;
  long len, offset, nbases, lwidth;
  bool found=false;
  while ( !FAI.eof() ) {
    string tmps;
    getline(FAI,tmps);
    istringstream iss(tmps);
    iss >> ichr >> len >> offset >> nbases >> lwidth;
    if ( ichr==chr ) { found=true; break; }
  }
  //  FAI.close();
  if ( !found ) { 
    ref="";
    return false;
  } 
  
  ifstream FIN(fastaFile.c_str());
  
  long flen;
  flen=len+(len/nbases)*(lwidth-nbases);
  ref.reserve(len+1);
  ref.resize(len);
  char *buffer = new char [flen+1];
  FIN.seekg(offset, std::ios::beg);

  //cerr << "seeking to offset " << offset << endl;
  FIN.read(buffer,flen);
  FIN.close();
  buffer[flen]=0;
  for(k=0,i=0; k<flen ;++k) 
    if ( buffer[k] != '\n' ) { 
      ref[i]=buffer[k];
      ++i;
    }
  delete[] buffer;
  if ( i!=len ) {
    cerr << "Error reading the reference fasta\n"
	 << "read " << i << " bases\n"
	 << "expecting " << len << " bases"
	 << endl;
    exit(0);
  }
  
  return true;
}

void get_N_regions(string& fasta, vector<int>& beg, vector<int>& end)
{
  
  beg.clear(); 
  end.clear();
  
  for(int i=0; i<(int)fasta.size();++i) {
    if( fasta[i]!='N' ) continue;
    if ( end.size()==0 ) {
      beg.push_back(i);
      end.push_back(i);
      continue;
    }
    
    if ( i==(end.back()+1) ) {
      end.back()=i;
      continue;
    }
    
    beg.push_back(i);
    end.push_back(i);
  }
  
  return;
}

void load_reference(string fastaFile, string chr, string& ref)
{
  
  ref.clear();
  if ( read_fasta(fastaFile, chr, ref) ) {
    cerr << "#loaded " << chr << " from " << fastaFile << "\tlength=" 
	 << ref.size() << endl;
    return;
  }
  
  string tmps=chr;
  if ( chr.find("chr")==0 || chr.find("CHR")==0 ) chr=chr.substr(3); 
  else  chr="chr"+chr;
  if ( read_fasta(fastaFile, chr, ref) ) {
    cerr << "#loading " << tmps << " failed\n"
	 << "#loaded " << chr << " from " << fastaFile << "\tlength=" 
	 << ref.size() 
	 << endl;
    return;
  }
  else {
    cerr << "loading " << tmps << " failed" << endl
	 << "loading " << chr << " failed" << endl;
  }
  
  return;
}
