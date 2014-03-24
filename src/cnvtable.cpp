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
using namespace std;

#include "functions.h"

struct cnv_st {
  string RNAME;
  string type;
  int P1;
  int P2;
  int un;
  cnv_st():RNAME("CHROM"),
	   type("var"),
	   P1(0),
	   P2(0) {};
} ;

string cnv_marker(const cnv_st& icnv) { 
  return icnv.RNAME+":"+to_string(icnv.P1)+"-"+to_string(icnv.P2)+":"+icnv.type;
}

bool comp_cnv_st(const cnv_st& cnv1, const cnv_st& cnv2)
{
  if ( cnv1.RNAME != cnv2.RNAME ) {
    string R1=cnv1.RNAME;
    string R2=cnv2.RNAME;
    if ( R1.find("chr")==0 ) R1=R1.substr(3);
    if ( R1=="X" ) R1="23";
    if ( R1=="Y" ) R1="24";
    if ( R2.find("chr")==0 ) R2=R2.substr(3);
    if ( R2=="X" ) R2="23";
    if ( R2=="Y" ) R2="24";
    
    if ( (R1[0]>='0' && R1[0]<='9') && (R2[0]>='0' && R2[0]<='9') ) {
      return ( atoi(R1.c_str()) < atoi(R2.c_str()) );
    }
    else { return ( R1 < R2 ); }
    
  }
  else {
    if ( cnv1.P1 == cnv2.P1 ) return cnv1.P2 < cnv2.P2 ;
    return cnv1.P1 < cnv2.P1 ;
  }
}

bool comp_cnv_rname(const string& chr1, const string& chr2)
{
  if ( chr1==chr2 ) return false;
  
  string R1=chr1;
  string R2=chr2;
  if ( R1.find("chr")==0 ) R1=R1.substr(3);
  if ( R1=="X" ) R1="23";
  if ( R1=="Y" ) R1="24";
  if ( R2.find("chr")==0 ) R2=R2.substr(3);
  if ( R2=="X" ) R2="23";
  if ( R2=="Y" ) R2="24";
  
  return (R1[0]>='0' && R1[0]<='9') && (R2[0]>='0' && R2[0]<='9') ?  
    atoi(R1.c_str()) < atoi(R2.c_str()) : R1 < R2 ;
}

bool is_cnv_overlap(cnv_st& cnv1, cnv_st& cnv2)
{
  if ( cnv1.RNAME!=cnv2.RNAME ) return false;
  else {
    int cnv1_p1=cnv1.P1;
    int cnv1_p2=cnv1.P2;
    if ( cnv1_p1 > cnv1_p2 ) swap(cnv1_p1, cnv1_p2);
    int cnv2_p1=cnv2.P1;
    int cnv2_p2=cnv2.P2;
    if ( cnv2_p1 > cnv2_p2 ) swap(cnv2_p1, cnv2_p2);
    
    return max(cnv1_p1,cnv2_p1) <= min(cnv1_p2, cnv2_p2 );
  }
}

int cnv_overlap_length(cnv_st& cnv1, cnv_st& cnv2)
{
  if ( cnv1.RNAME!=cnv2.RNAME ) return 0;
  else {
    int cnv1_p1=cnv1.P1;
    int cnv1_p2=cnv1.P2;
    if ( cnv1_p1 > cnv1_p2 ) swap(cnv1_p1, cnv1_p2);
    int cnv2_p1=cnv2.P1;
    int cnv2_p2=cnv2.P2;
    if ( cnv2_p1 > cnv2_p2 ) swap(cnv2_p1, cnv2_p2);
    
    if ( max(cnv1_p1,cnv2_p1) > min(cnv1_p2, cnv2_p2 ) ) return 0;
    else return min(cnv1_p2, cnv2_p2 ) - max(cnv1_p1,cnv2_p1)  + 1;
  }
}

void read_samplecnvfile(string& inpfile, 
			vector<string>& samplecnvfile, vector<string>& sampleid)
{
  ifstream FIN(inpfile.c_str()) ;
  if ( !FIN ) {
    cerr << "cannot find file " << inpfile << endl;
    exit(0);
  }
  
  while ( !FIN.eof() ) {
    string tmps;
    getline(FIN,tmps);
    if ( tmps[0] == '#' || tmps.length() <2 ) continue;
    istringstream iss(tmps);
    string filename, samplename;
    if ( iss >> filename ) samplecnvfile.push_back(filename);
    if ( iss >> samplename ) sampleid.push_back(samplename);
  }
  FIN.close();
  
  return;
}

int usage_check_sample_table(int argc, char* argv[]) {
  cerr << "Function:\n  This program tabulates CNV overlaps from multiple CNV files\n"
       << "\nUsage:\n" 
       << "  cnvtable <options>\n"
       << "\nOptions:\n"
       << "  -cnv  STR  file names separated by white spaces\n"
       << "  -i    STR  IDs for the files separated by white spaces \n"
       << "  -cnvf STR  a file with list of file names and(or) ids\n"
       << "  -l    INT  minimum length to include, INT=3 \n"
       << "  -L    INT  maximum length to include, INT=10000000 \n"
       << "  -t    STR  only process STR type of CNVs\n"
       << "  -O  FLOAT  minimum reciprocal overlap ratio [0.0, 1.0], FLOAT=0.5\n"
       << "  -o    STR  outputfile, STR=STDOUT \n"
       << "\nExamples :\n"
       << "  cnvtable -cnv 1.txt 2.txt 3.txt -L 10000 -O 0.5 -o 123.txt\n"
       << "\nNote :\n"
       << "  If the list of filenames is long, it is necessary to put them in a file,\n"
       << "  and use the -cnvf option. Each line should contain 1(filename) string \nor 2(filename idname) strings.\n"
       << endl;
  
  return(0);
}


int check_sample_table(int argc, char* argv[])
{
  
  if ( argc<3 ) exit( usage_check_sample_table(argc, argv) );
  
  string mycommand="";
  string QNAME,FLAG,RNAME,POS,MAPQ,CIGAR,MRNM,MPOS,ISIZE,SEQ,QUAL,OPT;
  string POS1S,POS2S;
  string CNVTYPE="";
  int minLength=3,maxLength=10000000;
  double rOverlap=0.5;
  string cnvFileList="", cnvFile="", outputFile="STDOUT";
  vector<string> inputArgv;
  vector<string> samplecnvfile(0);
  vector<string> sampleid(0);
  vector<int> samplegenotype(0);

#define _next2 inputArgv[i]=""; inputArgv[i+1]=""; continue;
#define _next1 inputArgv[i]=""; continue;
  
  for(size_t i=0;i<(size_t) argc;++i) inputArgv.push_back(string(argv[i]));
  mycommand=inputArgv[0];
  for(size_t i=1;i<inputArgv.size();++i) mycommand+=" "+inputArgv[i];
  for(size_t i=1;i<inputArgv.size();++i) {
    if ( inputArgv[i]=="-cnv" ) {  // input variation files
      inputArgv[i]="";
      size_t k=i+1;
      while( inputArgv[k][0] != '-' ) {
	samplecnvfile.push_back(inputArgv[k]);
	inputArgv[k]="";
	++k; 
	if ( k>=inputArgv.size() ) break;
      }
      continue;
    }
    if ( inputArgv[i]=="-i" ) {  // sample ids
      inputArgv[i]="";
      size_t k=i+1;
      while( inputArgv[k][0] != '-' ) {
	sampleid.push_back(inputArgv[k]);
	inputArgv[k]="";
	++k; 
	if ( k==inputArgv.size() ) break;
      }
      continue;
    }
    if ( inputArgv[i]=="-cnvf" ) {  // input variation files
      cnvFileList=inputArgv[i+1];
      _next2;
    }
    if ( inputArgv[i]=="-l" ) {  // minimum length
      minLength=atoi(inputArgv[i+1].c_str());
      _next2;
    }
    if ( inputArgv[i]=="-L" ) {  // max length
      maxLength=atoi(inputArgv[i+1].c_str());
      _next2;
    }
    if ( inputArgv[i]=="-t" ) {  // cnv type
      CNVTYPE=inputArgv[i+1];
      _next2;
    }
    if ( inputArgv[i]=="-O" ) {  // output file
      rOverlap=atof(inputArgv[i+1].c_str());
      _next2;
    }
    if ( inputArgv[i]=="-o" ) {  // output file
      outputFile=inputArgv[i+1];
      _next2;
    }
  }
  
  if ( cnvFileList!="" ) read_samplecnvfile(cnvFileList, samplecnvfile, sampleid);
  if ( rOverlap>1.0 || rOverlap<0.0 ) {
    cerr << "Overlap ration must be between 0.0 and 1.0" << endl;
    exit( usage_check_sample_table(argc, argv) );
  }
  if ( samplecnvfile.size()==0 ) exit( usage_check_sample_table(argc, argv) );
  
  cerr << "#CommandL : " << mycommand << "\n";
  
  // infer sample id from filenames
  if ( sampleid.size() != samplecnvfile.size() ) {
    sampleid.empty();
    for(size_t i=0;i<samplecnvfile.size();++i) {
      size_t found=samplecnvfile[i].rfind('/');
      string id=samplecnvfile[i].substr(found+1);
      found=id.find('.');
      if (found==string::npos) found=id.length();
      id=id.substr(0,found);
      sampleid.push_back(id);
    }
  }
  
  streambuf* sbuf = cout.rdbuf();
  ofstream FOUT;
  if ( outputFile != "STDOUT" ) {
    FOUT.open(outputFile.c_str());
    cout.rdbuf(FOUT.rdbuf());
  }
  
  vector<cnv_st> allcnvlist(0);
  cnv_st icnv;
  
  for(size_t i=0;i<samplecnvfile.size();++i) {
    
    if ( ! file_exist(samplecnvfile[i]) ) { 
      cerr << samplecnvfile[i] << " not found\n"; 
      exit(0); 
    }
    
    ifstream FIN(samplecnvfile[i].c_str()) ;
    vector<cnv_st> cnvlist(0);
    while ( !FIN.eof() ) {
      string tmps,tmps1;
      getline(FIN,tmps);
      if ( tmps[0] == '#' || tmps.length() <2 ) continue;
      istringstream iss(tmps);
      iss >> icnv.RNAME >> icnv.P1 >> icnv.P2 >> icnv.type;
      if ( CNVTYPE!="" && CNVTYPE!=icnv.type ) continue;
      if ( icnv.P1 > icnv.P2 ) swap(icnv.P1, icnv.P2);
      if ( icnv.P2 - icnv.P1 < minLength ) continue;
      if ( icnv.P2 - icnv.P1 > maxLength ) continue;
      
      cnvlist.push_back(icnv);
    }
    FIN.close();
    
    sort(cnvlist.begin(), cnvlist.end(), comp_cnv_st);
    cerr << "#" << sampleid[i] << "\t" << cnvlist.size() << endl;
    allcnvlist.insert(allcnvlist.end(), cnvlist.begin(),cnvlist.end() );
  }
  cerr << "#total cnv " << allcnvlist.size() << endl;
  
  sort(allcnvlist.begin(), allcnvlist.end(), comp_cnv_st);
  cerr << "#sorted cnv " << allcnvlist.size() << endl;
  
  // compact the list, remove duplicate and mostly overlapped
  // /*
  vector<cnv_st> allcnvlist_tmp(0);
  allcnvlist_tmp.push_back(allcnvlist[0]);
  for(size_t k=1;k<allcnvlist.size();++k) {
    
    if ( ! is_cnv_overlap(allcnvlist_tmp.back(), allcnvlist[k] )  ) {
      allcnvlist_tmp.push_back(allcnvlist[k]);
      continue;
    }
    
    //double overLen=cnv_overlap_length(allcnvlist[k], allcnvlist_tmp.back());
    // if ( overLen > ( allcnvlist_tmp.back().P2-allcnvlist_tmp.back().P1 )*0.9 &&
    //	 overLen > ( allcnvlist[k].P2-allcnvlist[k].P1 )*0.9 ) {
    if ( allcnvlist_tmp.back().RNAME==allcnvlist[k].RNAME && 
	 allcnvlist_tmp.back().P1==allcnvlist[k].P1 &&
	 allcnvlist_tmp.back().P2==allcnvlist[k].P2 ) {
      continue;
      icnv=allcnvlist_tmp.back();
      icnv.P1=min( allcnvlist[k].P1, allcnvlist_tmp.back().P1 );
      icnv.P2=max( allcnvlist[k].P2, allcnvlist_tmp.back().P2 );
      allcnvlist_tmp.back()=icnv;
    }
    else { allcnvlist_tmp.push_back(allcnvlist[k]); }
  }
  allcnvlist=allcnvlist_tmp;
  cerr << "#total after overlaps combined : " << allcnvlist.size() << endl;
  //*/
  
  // read cnv again and check overlap
  vector< vector<int> > cnvoverlap( allcnvlist.size(), 
				    vector<int>(samplecnvfile.size(), 0) );
  
  for(size_t i=0;i<samplecnvfile.size();++i) {
    
    cerr << i << "\t" << samplecnvfile[i] << endl;
    ifstream FIN(samplecnvfile[i].c_str()) ;
    vector<cnv_st> cnvlist(0);
    while ( !FIN.eof() ) {
      string tmps,tmps1;
      getline(FIN,tmps);
      if ( tmps[0] == '#' || tmps.length() <2 ) continue;
      istringstream iss(tmps);
      if (!( iss >> icnv.RNAME >> icnv.P1 >> icnv.P2 >> icnv.type) ) continue;
      
      if ( icnv.P1 > icnv.P2 ) swap(icnv.P1, icnv.P2);
      if ( icnv.P2 - icnv.P1 < minLength ) continue;
      if ( icnv.P2 - icnv.P1 > maxLength ) continue;
      if ( CNVTYPE!=""  && icnv.type != CNVTYPE) continue;
      
      cnvlist.push_back(icnv);
    }
    sort(cnvlist.begin(), cnvlist.end(), comp_cnv_st);
    
    size_t istart=0;
    size_t count=0;
    for(size_t k=0; k<cnvlist.size(); ++k) {
      
      for(size_t j=istart; j<allcnvlist.size(); ++j) {
	
	int overlap_len=cnv_overlap_length(allcnvlist[j] , cnvlist[k]);
	if ( comp_cnv_st(allcnvlist[j], cnvlist[k]) && overlap_len==0 ) {
	  if ( allcnvlist[j].RNAME==cnvlist[k].RNAME &&
	       cnvlist[k].P1>allcnvlist[j].P2+10000000 ) istart=j;
	  continue;
	}
	
	if ( overlap_len>0 ) { 
	  if ( overlap_len > (double)max(abs(allcnvlist[j].P2-allcnvlist[j].P1), abs(cnvlist[k].P2-cnvlist[k].P1))*rOverlap ) {
	    count+=1;
	    cnvoverlap[j][i]=overlap_len;
	  }
	  continue;
	}
	
	if ( (!comp_cnv_st(allcnvlist[j], cnvlist[k])) && overlap_len==0) break;
      }
      
    }
    cerr << i << "\t" << cnvlist.size()  << "\t" << count << endl;
    
  }
  
  vector<size_t> sumocnv(samplecnvfile.size(),0);
  vector<size_t> sumoid(allcnvlist.size(),0);
  for(size_t k=0;k<allcnvlist.size();++k) 
    for(size_t i=0;i<samplecnvfile.size();++i) 
      if ( cnvoverlap[k][i]>0 ) { 
	sumoid[k]++; 
	sumocnv[i]++;
      }
  
  cout << "##sampleid_and_inputfile\n";
  for(size_t i=0;i<samplecnvfile.size();++i) 
    cout << "##S" << i+1 << "\t" 
	 << sampleid[i] << "\t" 
	 << samplecnvfile[i] << "\n";
  cout << "##length of cnv included: [" << minLength << ", " << maxLength << "]" << endl;
  cout << "##portion of reciprocal overlap required: " << rOverlap << endl;
  
  cout << "#RNAME\t" << "POS1\t" << "POS2\t" << "TYPE\t" << "ROWSUM";
  for(size_t i=0;i<samplecnvfile.size();++i) cout << "\tS" << i+1;
  cout << "\n";
  cout << "##CSUM\t" << "POS1\t" << "POS2\t" << "TYPE\t" << "ROWSUM";
  for(size_t i=0;i<samplecnvfile.size();++i) cout << "\t" << sumocnv[i];
  cout << endl;
  
  for(size_t k=0;k<allcnvlist.size();++k) {
    cout << allcnvlist[k].RNAME << "\t"
	 << allcnvlist[k].P1 << "\t"
	 << allcnvlist[k].P2 << "\t"
	 << allcnvlist[k].type << "\t"
	 << sumoid[k];
    for(size_t i=0;i<samplecnvfile.size();++i) cout << "\t" << cnvoverlap[k][i];
    cout << endl;
  }
  
  cout.rdbuf(sbuf);
  return(0);
}

int main(int argc, char* argv[]) { return check_sample_table(argc, argv); }
  
