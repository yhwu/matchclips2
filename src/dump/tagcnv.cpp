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
#include <map>
#include <unistd.h>
using namespace std;

/**** user headers ****/
#include "functions.h"
#include "readref.h"
#include "read_CIGAR.h"
#include "matchreads.h"
#include "tagcnv.h"




int chrom2int(string RNAME)
{
  static map<string, int> RNAMEINT;
  if ( RNAME=="X" ) return 23;
  if ( RNAME=="Y" ) return 24;
  if ( RNAME.length()<3 ) return atoi(RNAME.c_str());
  RNAME=to_upper(RNAME);
  if ( RNAME.find("CHROM")!=string::npos ) RNAME=RNAME.substr(5);
  if ( RNAME.find("CHR")!=string::npos ) RNAME=RNAME.substr(3);
  if ( RNAME=="X" ) return 23;
  if ( RNAME=="Y" ) return 24;
  int ichr=atoi(RNAME.c_str());
  if (  ichr>0 && ichr<=22 ) return ichr;
  
  map<string,int>::const_iterator it = RNAMEINT.find(RNAME);
  if ( it!=RNAMEINT.end() ) return RNAMEINT[RNAME]+24;
  RNAMEINT[RNAME]=RNAMEINT.size()+1;
  return RNAMEINT[RNAME]+24;
}

int64_t chrompos2int64(string RNAME, int POS) { return chrom2int(RNAME)*(int64_t)1E10+POS; }


/* extract fields from genotype 
 * if OFORMAT=="" output is same as input
 */
string getvcfgenotypefield(string& igenotype, string& FORMAT, string& OFORMAT)
{
  string ogenotype="NF";
  if ( OFORMAT=="" ) return igenotype;
  vector<int> ofields;
  vector<string> ifmt,ofmt, ogt;
  stringsplit(FORMAT, ifmt, ':');
  stringsplit(OFORMAT, ofmt, ':');
  stringsplit(igenotype, ogt, ':');
  
  size_t i;
  for(i=0; i<ofmt.size(); ++i) {
    int iof= find(ifmt.begin(),ifmt.end(), ofmt[i]) - ifmt.begin();
    if ( iof!=(int)ifmt.size() ) ofields.push_back(iof);
  }
  
  if ( ofields.size()<1 ) return ogenotype;
  
  ogenotype = ofields[0]<(int)ogt.size() ? ogt[ ofields[0] ] : "";
  for(i=1;i<ofields.size();++i) {
    ogenotype+=":";
    if ( ofields[i] < (int)ogt.size() ) ogenotype+=ogt[ ofields[i] ];
  }
  
  return ogenotype;
}

int usage_check_noseq_regions(int argc, char* argv[]) {
  cerr << "This subroutine checks if the CNV regions cover the non-coded part in the reference fasta.\n"
       << "PASS or FAIL is attached at the end of a line\n" 
       << "\nUsage:\n" 
       << "  " << argv[0] << " " << argv[1] << " <options> -f REFERENCE -v CNVFILE\n"
       << "\nOptions:\n"
       << "  -o  STR  outputfile, STR=STDOUT \n"
       << "\nExamples:\n"
       << "  " << argv[0] <<  " " << argv[1] << " -f human_g1k_v37.fasta -b bwap40X.bam -v bwap40X.bam.bp\n"
       << endl;
  
  return(0);
}

int check_noseq_regions(int argc, char* argv[])
{
  
  bool NOSEQ=false;
  string mycommand="";
  size_t i, k;
  string QNAME,FLAG,RNAME,POS,MAPQ,CIGAR,MRNM,MPOS,ISIZE,SEQ,QUAL,OPT;
  int POS1,POS2;
  string POS1S,POS2S;
  string CIGAR1,SEQ1,QUAL1,OPT1;
  
  string cnvFile="",fastaFile="",outputFile="STDOUT";
  string FASTA="",fastaname="";
  vector<string> inputArgv;
  
  for(i=0;i<(size_t) argc;++i) inputArgv.push_back(string(argv[i]));
  mycommand=inputArgv[0];
  for(i=1;i<inputArgv.size();++i) mycommand+=" "+inputArgv[i];
  for(i=2;i<inputArgv.size();++i) {
    if ( inputArgv[i]=="-v" ) {  // input variation file
      cnvFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-f" ) {  // input reference file
      fastaFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-o" ) {  // output file
      outputFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
  }

  cerr << "#CommandL : " << mycommand << "\n";
  if ( fastaFile=="" || cnvFile=="" ) exit( usage_check_noseq_regions(argc, argv) );
  
  if ( ! file_exist(fastaFile) ) { cerr << fastaFile << " not found\n"; exit(0); }
  if ( ! file_exist(cnvFile) ) { cerr << cnvFile << " not found\n"; exit(0); }
  
  streambuf* sbuf = cout.rdbuf();
  ofstream FOUT;
  if ( outputFile != "STDOUT" ) {
    FOUT.open(outputFile.c_str());
    cout.rdbuf(FOUT.rdbuf());
  }
  
 
 
  vector<string> FILTER(2);
  FILTER[0]="PASS";
  FILTER[1]="FAIL";
  
  ifstream FIN(cnvFile.c_str());
  if ( !FIN ) { cerr << "Can't open file " << outputFile << endl; exit(0); }
  i=0;
  while ( !FIN.eof() ) {
    string tmps,tmps1;
    getline(FIN,tmps);
    istringstream iss(tmps);
    if ( tmps[0] == '#' ) {
      cout << tmps << "\tNOSEQ" << endl;
      continue;
    }
    if ( tmps.length() <2 ) continue;
    
    iss >> RNAME 
	>> POS1S 
	>> POS2S;
    
    POS1=atoi(POS1S.c_str());
    if ( to_string(POS1) != POS1S ) {
      cerr << "Integer overflow " << POS1S << "\t" << POS1 << endl;
      exit(0);
    }
    POS2=atoi(POS2S.c_str());
    if ( to_string(POS2) != POS2S ) {
      cerr << "Integer overflow " << POS2S << "\tas\t" << POS2 << endl;
      exit(0);
    }
    
    if ( RNAME != fastaname ) {
      string RNAMEtmp=RNAME;
      if ( ci_find(RNAME,"chr") != string::npos ) RNAMEtmp=RNAME.substr(3); 
      read_fasta(fastaFile,RNAMEtmp,FASTA);
      fastaname=RNAME;
      cerr << fastaFile << "\t"
	   << fastaname << "\t"
	   << FASTA.size() << endl;
    }
    NOSEQ=false;
    size_t p1=POS1;
    size_t p2=POS2;
    if ( p1>p2 ) swap(p1,p2); 
    if (p1==0) p1=1;
    if ( p1>=FASTA.length() || p2>=FASTA.length() ) NOSEQ=true;
    else {
      for(k=p1-1;k<p2;++k) if ( FASTA[k]=='N' ) { NOSEQ=true; break; }
    }
    
    cout << tmps << "\t" << FILTER[ NOSEQ ] << endl; 
  }
  FIN.close();

  FOUT.close();
  cout.rdbuf(sbuf);
  return(0);
}

void get_read_pos(string& read, int& pos1, int& pos2)
{
  string QNAME,FLAG,RNAME,POS,MAPQ,CIGAR,MRNM,MPOS,ISIZE,SEQ,QUAL,OPT;
  int misbeg;

  istringstream iss(read);
  iss >> QNAME
      >> FLAG
      >> RNAME 
      >> POS 
      >> MAPQ 
      >> CIGAR 
      >> MRNM 
      >> MPOS 
      >> ISIZE 
      >> SEQ 
      >> QUAL 
      >> OPT;

  pos1=atoi(POS.c_str());;
  pos2=pos1+SEQ.length()-1;
  
  vector<char> op(100);
  vector<int> opnum(100);
  vector<int> opposseq(100);
  vector<int> oppos(100);
  size_t iclip,anchor;
  
  read_POS_CIGAR(pos1,CIGAR,
		 anchor, iclip, 
		 op, opnum,
		 opposseq, oppos );
  
  if ( op[0]=='*' ) return;
  
  misbeg=opposseq.back()+opnum.back();
  if ( misbeg !=(int)SEQ.size() ) {
    cerr << "something wrong reading CIGAR while " 
	 << "calculate the position of the operators in the string SEQ\n"
	 << read << "\n"
	 << CIGAR << "\t" << op.size() << "\n"
	 << misbeg << "\t" << SEQ.size()
	 << endl;
    exit(0);
  }
  
  pos1=oppos[0];
  pos2=oppos[0]+SEQ.length()-1;

  return;
}


int usage_check_readdepth(int argc, char* argv[]) {
  cerr << "This subroutine checks the read depths outside, inside, and on the end points of CNV regions. Read depth info is appended at the end of a line in the format\n"
       << "  RDO=average_RD_of_outter_region;\n  RDI=average_RD_inside_CNV_region;\n  RD1=RD_of_the_lower_end;\n  RD2=RD_of_the_upper_end .\n"
       << "\nUsage:\n" 
       << "  " << argv[0] << " " << argv[1] << " <options> -b BAMFILE -v CNVFILE \n"
       << "Options  :\n"
       << "  -d  INT  distance before and after a CNV to be checked, INT=1000 \n"
       << "  -o  STR  outputfile, STR=STDOUT \n"
       << "\nExamples :\n"
       << "  " << argv[0] << " " << argv[1] <<  " -f human_g1k_v37.fasta -b bwap40X.bam -q 10 -v bwap40X.bam.bp \n"
       << endl;
  
  return(0);
}
int check_readdepth_bak(int argc, char* argv[])
{
  
  string mycommand="";
  size_t i;
  string QNAME,FLAG,RNAME,POS,MAPQ,CIGAR,MRNM,MPOS,ISIZE,SEQ,QUAL,OPT;
  int POS1,POS2,DIS=1000;
  string POS1S,POS2S;
  string CIGAR1,SEQ1,QUAL1,OPT1;
  
  string cnvFile="",outputFile="STDOUT";
  string bamFile="";
  string fastaFile="",FASTA="",fastaname="";
  string samtoolsOpt1="", samtoolsOpt2="", samtoolsCommand="";
  vector<string> inputArgv;
  
  for(i=0;i<(size_t) argc;++i) inputArgv.push_back(string(argv[i]));
  mycommand=inputArgv[0];
  for(i=1;i<inputArgv.size();++i) mycommand+=" "+inputArgv[i];
  for(i=2;i<inputArgv.size();++i) {
    if ( inputArgv[i]=="-v" ) {  // input variation file
      cnvFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-b" ) {  // input variation file
      bamFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-f" ) {  // input reference file
      fastaFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-o" ) {  // output file
      outputFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-d" ) {  // output file
      DIS=atoi(inputArgv[i+1].c_str());
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
  }
  for(i=2;i<inputArgv.size();++i) 
    if ( inputArgv[i]!="" ) samtoolsOpt1+=inputArgv[i]+" ";
  samtoolsCommand="samtools view " + samtoolsOpt1 + bamFile;
  
  if ( bamFile=="" || cnvFile=="" ) exit( usage_check_readdepth(argc, argv) );
  if ( ! file_exist(bamFile) ) { cerr << bamFile << " not found\n"; exit(0); }
  if ( ! file_exist(cnvFile) ) { cerr << cnvFile << " not found\n"; exit(0); }
  
  cerr << "#CommandL  : " << mycommand << "\n"
       << "#CNV file  : " << cnvFile << "\n"
       << "#BAM file  : " << bamFile << "\n"
       << "#Distance  : " << DIS << "\n" 
       << "#SAMT opts : " << samtoolsOpt1 << "\n"
       << "#SAMT exec : " << samtoolsCommand << "\n"
       << "#OUTPUT    : " << outputFile << "\n"
       << endl;
  
  streambuf* sbuf = cout.rdbuf();
  ofstream FOUT;
  if ( outputFile != "STDOUT" ) {
    FOUT.open(outputFile.c_str());
    cout.rdbuf(FOUT.rdbuf());
  }
  
  ifstream FIN(cnvFile.c_str());
  if ( !FIN ) { cerr << "Can't open file " << outputFile << endl; exit(0); }
  i=0;
  while ( !FIN.eof() ) {
    string tmps,tmps1;
    getline(FIN,tmps);
    istringstream iss(tmps);
    if ( tmps[0] == '#' ) {
      cout << tmps << "\tRD[DIS=+-" << DIS << "]" << endl;
      continue;
    }
    if ( tmps.length() <2 ) continue;
    
    iss >> RNAME 
	>> POS1S 
	>> POS2S;
    
    POS1=atoi(POS1S.c_str());
    if ( to_string(POS1) != POS1S ) {
      cerr << "Integer overflow " << POS1S << "\t" << POS1 << endl;
      exit(0);
    }
    POS2=atoi(POS2S.c_str());
    if ( to_string(POS2) != POS2S ) {
      cerr << "Integer overflow " << POS2S << "\tas\t" << POS2 << endl;
      exit(0);
    }
    
    int p1=POS1, p2=POS2;
    if ( p1>p2 ) swap(p1,p2);
    int p1e=p1-DIS, p2e=p2+DIS;
    if ( p1e<1 ) p1e=1;
    
    string region=RNAME+":"+to_string(p1)+"-"+to_string(p2);
    string regionouter=RNAME+":"+to_string(p1e)+"-"+to_string(p2e);
    //    string SAMEXEC=samtoolsCommand+" "+region;
    string SAMEXEC=samtoolsCommand+" "+regionouter;

    FILE *fp;
    fp=popen( SAMEXEC.c_str(),"r");
    if(!fp) {
      cerr << "can't open pipe for input \n"
	   << SAMEXEC << endl;
      exit(0);
    }

    size_t C_inner=0,C_outer=0,C_p1=0,C_p2=0;
    vector<size_t> RC(p2e-p1e+2,0);
    while ( !feof(fp) ) {
      string  read="";
      fgetline(fp, read);
      //cout << endl;
      //cout << read << endl;
      //cout << SAMEXEC << endl;
      if ( read.length() < 3 ) continue;
      int r1,r2;
      get_read_pos(read, r1, r2);
      //get_read_pos_bak(read, r1, r2);
      for(int i=r1;i<=r2;++i) if ( i>=p1e && i<=p2e ) RC[i-p1e]++;
    }
    fclose(fp);
    
    for(int i=p1e;i<p1;++i) C_outer+=RC[i-p1e];
    for(int i=p2+1;i<=p2e;++i) C_outer+=RC[i-p1e];
    C_outer/=(p1-p1e+p2e-p2);
    for(int i=p1;i<=p2;++i) C_inner+=RC[i-p1e];
    C_inner/=(p2-p1+1);
    C_p1=RC[p1-p1e];
    C_p2=RC[p2-p1e];
    
    string tag="RDO="+to_string(C_outer)+";"+
      "RDI="+to_string(C_inner)+";"+
      "RD1="+to_string(C_p1)+";"+
      "RD2="+to_string(C_p2);
    
    /*
    cout << region << "\t" 
	 << p2-p1+1 << "\t"
	 << C_outer << "\t"
	 << C_inner << "\t"
	 << C_p1 << "\t"
	 << C_p2 << "\t"
	 << tag
	 << endl;
    */
    cout << tmps << "\t" << tag << endl;
  }
  FIN.close();

  FOUT.close();
  cout.rdbuf(sbuf);
  return(0);
}

int check_readdepth(int argc, char* argv[])
{
  
  string mycommand="";
  size_t i;
  string QNAME,FLAG,RNAME,POS,MAPQ,CIGAR,MRNM,MPOS,ISIZE,SEQ,QUAL,OPT;
  int POS1,POS2,DIS=500;
  string POS1S,POS2S;
  string CIGAR1,SEQ1,QUAL1,OPT1;
  
  string cnvFile="",outputFile="STDOUT";
  string bamFile="";
  string fastaFile="",FASTA="",fastaname="";
  string samtoolsOpt1="", samtoolsOpt2="", samtoolsCommand="";
  vector<string> inputArgv;
  
  for(i=0;i<(size_t) argc;++i) inputArgv.push_back(string(argv[i]));
  mycommand=inputArgv[0];
  for(i=1;i<inputArgv.size();++i) mycommand+=" "+inputArgv[i];
  for(i=2;i<inputArgv.size();++i) {
    if ( inputArgv[i]=="-v" ) {  // input variation file
      cnvFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-b" ) {  // input variation file
      bamFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-f" ) {  // input reference file
      fastaFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-o" ) {  // output file
      outputFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-d" ) {  // output file
      DIS=atoi(inputArgv[i+1].c_str());
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
  }
  for(i=2;i<inputArgv.size();++i) 
    if ( inputArgv[i]!="" ) samtoolsOpt1+=inputArgv[i]+" ";
  samtoolsCommand="samtools view " + samtoolsOpt1 + bamFile;
  
  if ( bamFile=="" || cnvFile=="" ) exit( usage_check_readdepth(argc, argv) );
  if ( ! file_exist(bamFile) ) { cerr << bamFile << " not found\n"; exit(0); }
  if ( ! file_exist(cnvFile) ) { cerr << cnvFile << " not found\n"; exit(0); }
  
  cerr << "#CommandL  : " << mycommand << "\n"
       << "#CNV file  : " << cnvFile << "\n"
       << "#BAM file  : " << bamFile << "\n"
       << "#Distance  : " << DIS << "\n" 
       << "#SAMT opts : " << samtoolsOpt1 << "\n"
       << "#SAMT exec : " << samtoolsCommand << "\n"
       << "#OUTPUT    : " << outputFile << "\n"
       << endl;
  
  streambuf* sbuf = cout.rdbuf();
  ofstream FOUT;
  if ( outputFile != "STDOUT" ) {
    FOUT.open(outputFile.c_str());
    cout.rdbuf(FOUT.rdbuf());
  }
  
  vector<char> op(100);    // CIGAR operator
  vector<int> opnum(100);  // CIGAR operator bases
  vector<int> oppos(100);  // positions of operators in REF
  vector<int> opposseq(100);  // positions of operators in SEQ
  size_t anchor,iclip;
  string  read = "";
  
  
  ifstream FIN(cnvFile.c_str());
  if ( !FIN ) { cerr << "Can't open file " << outputFile << endl; exit(0); }
  i=0;
  while ( !FIN.eof() ) {
    string tmps,tmps1;
    getline(FIN,tmps);
    istringstream iss(tmps);
    if ( tmps[0] == '#' ) {
      cout << tmps << "\tRD[DIS=+-" << DIS << "]" << endl;
      continue;
    }
    if ( tmps.length() <2 ) continue;
    
    iss >> RNAME 
	>> POS1S 
	>> POS2S;
    
    POS1=atoi(POS1S.c_str());
    if ( to_string(POS1) != POS1S ) {
      cerr << "Integer overflow " << POS1S << "\t" << POS1 << endl;
      exit(0);
    }
    POS2=atoi(POS2S.c_str());
    if ( to_string(POS2) != POS2S ) {
      cerr << "Integer overflow " << POS2S << "\tas\t" << POS2 << endl;
      exit(0);
    }
    
    int p1=POS1, p2=POS2;
    if ( p1>p2 ) swap(p1,p2);
    int p1e=p1-DIS, p2e=p2+DIS;
    if ( p1e<1 ) p1e=1;
    
    string regionouter=RNAME+":"+to_string(p1e)+"-"+to_string(p2e);
    string SAMEXEC=samtoolsCommand+" "+regionouter;
    
    FILE *fp;
    fp=popen( SAMEXEC.c_str(),"r");
    if(!fp) {
      cerr << "can't open pipe for input \n"
	   << SAMEXEC << endl;
      exit(0);
    }
    
    size_t dout=0,din=0,dp1=0,dp2=0;
    while ( !feof(fp) ) {
      fgetline(fp, read);
      
      istringstream iss(read);
      iss >> QNAME
	  >> FLAG
	  >> RNAME 
	  >> POS 
	  >> MAPQ 
	//  >> mapq 
	  >> CIGAR 
	  >> MRNM 
	  >> MPOS 
	  >> ISIZE 
	  >> SEQ 
	  >> QUAL 
	  >> OPT;
      
      if ( QNAME[0]=='@' || QNAME[0]=='#' ) continue;
      if ( CIGAR=="*" ) continue;
      if ( RNAME=="*" ) continue;
      
      POS1=atoi(POS.c_str());
      read_POS_CIGAR(POS1,CIGAR,
		     anchor, iclip, 
		     op, opnum,
		     opposseq, oppos );
      if ( op[0]=='*' ) continue;
      
      vector<int> rpos(SEQ.size(),0);
      int idx=0;
      for(size_t k=0;k<op.size();++k) {
	if ( op[k]=='D' || op[k]=='N' || op[k]=='P' ) continue;
	for(int j=0;j<opnum[k];++j) {
	  if ( idx==(int)SEQ.size() ) break;
	  if ( op[k]!='I' ) rpos[idx]=oppos[k]+j;
	  else rpos[idx]=0;
	  // else rpos[idx]=oppos[k];
	  ++idx;
	}
      }
      if ( idx!=(int)SEQ.size() ) {
	cerr << "Error understanding CIGAR and POS " 
	     << POS << "\t" << CIGAR << "\t" << idx << "\t" << SEQ.size() << endl;
	for(size_t k=0;k<op.size();++k) 
	  cerr << op[k] << "\t" << oppos[k] << "t" << opnum[k] << endl;
      }
      
      for(size_t k=0;k<SEQ.size();++k) {
	if ( rpos[k]>p1-SEQ.size() && rpos[k]<p1 ) ++dout;
	if ( rpos[k]>p2 && rpos[k]<p2+SEQ.size() ) ++dout;
	if ( rpos[k]==p1 ) ++dp1;
	if ( rpos[k]==p2 ) ++dp2;
	if ( rpos[k]>p1 && rpos[k]<p2 ) ++din;
      }
      
    }
    fclose(fp);

    int dx=p2-p1-1;
    if ( dx<=0 ) dx=1;
    
    dout /= (2*SEQ.size()-2);
    din  /= dx;
    
    string tag="RDO="+to_string(dout)+";"+
      "RDI="+to_string(din)+";"+
      "RD1="+to_string(dp1)+";"+
      "RD2="+to_string(dp2);
    
    cout << tmps << "\t" << tag << endl;
    
    
  }
  
  FIN.close();
  
  FOUT.close();
  cout.rdbuf(sbuf);
  return(0);
}

int usage_check_mapq0(int argc, char* argv[]) {
  cerr << "This subroutine checks the proportion of ambiguous reads(mapq=0) within a CNV region. A number between 0-1 is attached at the end of a line.\n";
  cerr << "\nUsage:\n" 
       << "   " << argv[0] << " " << argv[1] << " <options> -v CNVFILE -b BAMFILE \n"
       << "\nOptions:\n"
       << "  -q  INT  minimum mapq, INT=0 \n"
       << "  -o  STR  outputfile, STR=STDOUT \n"
       << "\nExamples:\n"
       << "  " << argv[0] << " " << argv[1] << " -b bwap40X.bam -v bwap40X.bam.bp\n"
       << "  " << argv[0] << " " << argv[1] << " -q 5 -b bwap40X.bam -v bwap40X.bam.bp\n"
       << "\nNote:\n"
       << "  Different mapping software may have different minimum mapping quality score. \n"
       << "  Use -q to change to desired value.\n"
       << endl;
  
  return(0);
}

int check_mapq0(int argc, char* argv[])
{
  int MINIMUMMAPQ=0;
  
  string mycommand="";
  size_t i;
  string QNAME,FLAG,RNAME,POS,MAPQ,CIGAR,MRNM,MPOS,ISIZE,SEQ,QUAL,OPT;
  int POS1,POS2;
  string POS1S,POS2S;
  string CIGAR1,SEQ1,QUAL1,OPT1;
  
  string cnvFile="",outputFile="STDOUT";
  string bamFile="";
  string fastaFile="",FASTA="",fastaname="";
  string samtoolsOpt1="", samtoolsOpt2="", samtoolsCommand="";
  vector<string> inputArgv;
  
  for(i=0;i<(size_t) argc;++i) inputArgv.push_back(string(argv[i]));
  mycommand=inputArgv[0];
  for(i=1;i<inputArgv.size();++i) mycommand+=" "+inputArgv[i];
  for(i=2;i<inputArgv.size();++i) {
    if ( inputArgv[i]=="-q" ) {  // mapq0
      MINIMUMMAPQ=atoi(inputArgv[i+1].c_str());
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-v" ) {  // input variation file
      cnvFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-b" ) {  // input variation file
      bamFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-o" ) {  // output file
      outputFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
  }
  for(i=2;i<inputArgv.size();++i) 
    if ( inputArgv[i]!="" ) samtoolsOpt1+=inputArgv[i]+" ";
  samtoolsCommand="samtools view " + bamFile;
  
  if ( bamFile=="" || cnvFile=="" ) exit( usage_check_mapq0(argc, argv) );
  
  if ( ! file_exist(bamFile) ) { cerr << bamFile << " not found\n"; exit(0); }
  if ( ! file_exist(cnvFile) ) { cerr << cnvFile << " not found\n"; exit(0); }
  
  cerr << "#CommandL : " << mycommand << "\n"
       << "#CNV file : " << cnvFile << "\n"
       << "#BAM file : " << bamFile << "\n"
       << "SAMT exec : " << samtoolsCommand << "\n"
       << endl;
  
  streambuf* sbuf = cout.rdbuf();
  ofstream FOUT;
  if ( outputFile != "STDOUT" ) {
    FOUT.open(outputFile.c_str());
    cout.rdbuf(FOUT.rdbuf());
  }
  
  ifstream FIN(cnvFile.c_str());
  if ( !FIN ) { cerr << "Can't open file " << outputFile << endl; exit(0); }
  i=0;
  while ( !FIN.eof() ) {
    string tmps,tmps1;
    getline(FIN,tmps);
    istringstream iss(tmps);
    if ( tmps[0] == '#' ) {
      cout << tmps << "\tMAPQ0" << endl;
      continue;
    }
    if ( tmps.length() <2 ) continue;
    
    iss >> RNAME 
	>> POS1S 
	>> POS2S;
    
    POS1=atoi(POS1S.c_str());
    if ( to_string(POS1) != POS1S ) {
      cerr << "Integer overflow " << POS1S << "\t" << POS1 << endl;
      exit(0);
    }
    POS2=atoi(POS2S.c_str());
    if ( to_string(POS2) != POS2S ) {
      cerr << "Integer overflow " << POS2S << "\tas\t" << POS2 << endl;
      exit(0);
    }
    
    int p1=POS1, p2=POS2;
    if ( p1>p2 ) swap(p1,p2);
    
    string region=RNAME+":"+to_string(p1)+"-"+to_string(p2);
    string SAMEXEC=samtoolsCommand+" "+region;
    
    FILE *fp;
    fp=popen( SAMEXEC.c_str(),"r");
    if(!fp) {
      cerr << "can't open pipe for input \n"
	   << SAMEXEC << endl;
      exit(0);
    }
    
    size_t C_q0=0,C_all=0;
    while ( !feof(fp) ) {
      string  read="";
      fgetline(fp, read);
      if ( read.length() < 3 ) continue;
      istringstream iss(read);
      iss >> QNAME >> FLAG >> RNAME >> POS >> MAPQ;
      int q=atoi(MAPQ.c_str());
      C_all++;
      if (q<=MINIMUMMAPQ) C_q0++;
    }
    fclose(fp);
    
    if ( C_all==0 ) C_all=1;
    float r=(float)C_q0/(float)C_all;
    cout << tmps << "\t" << setprecision(2) << setw(4) << r << endl;
  }
  FIN.close();
  
  FOUT.close();
  cout.rdbuf(sbuf);
  return(0);
}

void subsetbam(string bamFile, string regions, string tmpbam)
{
  
  string samtoolsCommand="samtools view -h " + bamFile + " "+ regions;
  cerr << "#" << samtoolsCommand << endl;
  
  string tmpsam=tmpbam+".sam";
  if ( ci_equal(tmpbam.substr(tmpbam.length()-4),".bam") )
    tmpbam=tmpbam.substr(0,tmpbam.length()-4);
  
  FILE *fp;
  fp=popen( samtoolsCommand.c_str(),"r");
  if(!fp) {
    cerr << "can't open pipe for input \n"
	 << samtoolsCommand << endl;
    exit(0);
  }
  
  ofstream FOUT(tmpsam.c_str());
  if ( !FOUT ) {
    cerr << "Cann't write to " << tmpbam << endl;
    exit(0);
  }
  
  size_t C_all=0;
  while ( !feof(fp) ) {
    string  read="";
    fgetline(fp, read);
    if ( read.length() < 3 ) continue;
    C_all++;
    FOUT << read << endl;
  }
  fclose(fp);
  FOUT.close();
  
  samtoolsCommand="samtools view -bS " + tmpsam + " | samtools sort - " + tmpbam; 
  
  cerr << "Converting to bam format " << endl;
  cerr << "executing " << samtoolsCommand << endl;
  system(samtoolsCommand.c_str());
  samtoolsCommand="samtools index " + tmpbam+".bam";
  cerr << "executing " << samtoolsCommand << endl;
  system(samtoolsCommand.c_str());
  string command="rm " + tmpsam;
  remove(tmpsam.c_str());
  
  return;
}

void subsetbam_long(string bamFile, string regions, string tmpbam)
{
  
  // this is too large for system and pipe commands 
  size_t arg_max = sysconf(_SC_ARG_MAX);
  arg_max = 10000;
  size_t cuts= int(arg_max*0.8);
  
  string samtoolsCommand="samtools view -h " + bamFile + " "+ regions;
  
  string tmpsam=tmpbam+".sam";
  // tmpbam prefix
  if ( ci_equal(tmpbam.substr(tmpbam.length()-4),".bam") )
    tmpbam=tmpbam.substr(0,tmpbam.length()-4);
  
  size_t i=0;
  bool end=false;
  vector<string> REG(0);
  string tmps="";
  while ( !end ) {
    size_t breakup=regions.find(" ",cuts-1);
    if ( breakup==string::npos ) { 
      end=true; 
      tmps=regions; 
    }
    else { 
      tmps=regions.substr(0,breakup); 
      regions=regions.substr(breakup);
    }
    if ( tmps.size()>2 ) REG.push_back(tmps);
  }
  
  if ( REG.size() > 1 )
    cerr << "Regions broken up into " << REG.size() << " parts" << endl;
  
  // subset into tmp sam file
  for( i=0; i<REG.size(); ++i ) {
    if ( i==0 ) 
      samtoolsCommand="samtools view -h " 
	+ bamFile + " "+ REG[i] + " > "; 
    else  
      samtoolsCommand="samtools view " 
	+ bamFile + " "+ REG[i] + " >> ";
    samtoolsCommand += tmpsam;
    cerr << samtoolsCommand << endl;
    system(samtoolsCommand.c_str());
  }
  
  // convert to sorted bam file
  cerr << "Converting to bam format " << endl;
  samtoolsCommand="samtools view -uS " + tmpsam + " | samtools sort - " + tmpbam; 
  cerr << samtoolsCommand << endl;
  system(samtoolsCommand.c_str());
  samtoolsCommand="samtools index " + tmpbam + ".bam";
  cerr << "executing " << samtoolsCommand << endl;
  system(samtoolsCommand.c_str());
  string command="rm -f " + tmpsam;
  system(command.c_str());
  
  return;
}


int usage_check_het(int argc, char* argv[]) {
  cerr  << "\nUsage:\n" 
	<< "  " << argv[0] << " " << argv[1] << " <options> -f REFERENCE -b BAMFILE -v CNVFILE \n"
	<< "\nOptions:\n"
	<< "  -t  FLT  mutton rate, samtools default FLT=0.001 \n"
	<< "  -d  INT  check INT before and after for read depth \n"
       << "   -q  INT  minimum mapq, INT=0 \n"
	<< "  -o  STR  outputfile, STR=STDOUT \n"
	<< "   *   *   all others are passed to samtools \n"
	<< "\nExamples:\n"
	<< "  " << argv[0] << " " << argv[1] << " -f human_g1k_v37.fasta -b bwap40X.bam -q 10 -v bwap40X.bam.bp \n"
	<< "\nNote:\n"
	<< "This subroutine checks the number of heterozyguous sites within a CNV region. If a region is homozyguous, HOM is appended. Otherwise, heteozyguous sites are given in the format\n"
	<< "  HET=T(MSRQ):POS(Q)[,POS(Q),...]\n"
	<< "where T=total number of heterozyguous sites,\n"
	<< "      MSRQ=mean square root of quality scores of zygosity as called by samtools/bcftools\n"
	<< "      POS=position, Q=quality score for each site.\n"
	<< "It really does not matter for duplications. \n"
	<< endl;   
  return(0);
}

int check_het(int argc, char* argv[])
{
  /* samtools mpileup -Q 13 -q 1 -uIf ../data/human_g1k_v37.fasta -l tmpcnvNA19240.19.txt  tmpreads.NA19240.19.bam | bcftools view -t 0.5 -g -
     FQ	int Consensus quality. Positive: sample genotypes different; negative: otherwise
     http://samtools.sourceforge.net/samtools.shtml
  */
  
  
  string mycommand="";
  size_t i;
  string QNAME,FLAG,RNAME,POS,MAPQ,CIGAR,MRNM,MPOS,ISIZE,SEQ,QUAL,OPT;
  int POS1,POS2;
  int dx=200;
  int MINIMUMDX=200;
  int MINIMUMMAPQ=0;
  string POS1S,POS2S;
  string CIGAR1,SEQ1,QUAL1,OPT1;
  double muttion_rate=0.001;
  string cnvFile="",outputFile="STDOUT";
  string bamFile="",hetFile="";
  string fastaFile="",FASTA="",fastaname="";
  string samtoolsOpt1="", samtoolsOpt2="", samtoolsCommand="", bcftoolsCommand="";
  string hetCommand="";
  string rdCommand="";
  string q0Command="";
  vector<string> inputArgv;
  char hostname[1024];
  hostname[1023] = '\0';
  gethostname(hostname, 1023);
  string host=hostname;
  i=host.find(".");
  if ( i!=string::npos ) host=host.substr(0,i);
  
  //cerr << "Running on "<< "\t" << host << endl;
  
  for(i=0;i<(size_t) argc;++i) inputArgv.push_back(string(argv[i]));
  mycommand=inputArgv[0];
  for(i=1;i<inputArgv.size();++i) mycommand+=" "+inputArgv[i];
  for(i=2;i<inputArgv.size();++i) {
    if ( inputArgv[i]=="-v" ) {  // input variation file
      cnvFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-b" ) {  // input bam file
      bamFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-f" ) {  // reference file
      fastaFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-z" ) {  // reference file
      hetFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-o" ) {  // output file
      outputFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-t" ) {  // mutation rate 
      muttion_rate=atof(inputArgv[i+1].c_str());
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-d" ) {  // outer region for RD file
      MINIMUMDX=atoi(inputArgv[i+1].c_str());
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-q" ) {  // minimum mapq
      MINIMUMMAPQ=atoi(inputArgv[i+1].c_str());
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
  }
  for(i=2;i<inputArgv.size();++i) 
    if ( inputArgv[i]!="" ) samtoolsOpt1+=inputArgv[i]+" ";
  
  cerr << "#CommandL   : " << mycommand << "\n"
       << "#CNV file   : " << cnvFile << "\n"
       << "#BAM file   : " << bamFile << "\n"
       << "#HET file   : " << hetFile << "\n"
       << "#FASTA      : " << fastaFile << "\n"
    //<< "#SAMT exec  : " << samtoolsCommand << "\n"
       << "#OUTPUT     : " << outputFile <<"\n"
       << endl;
  
  if ( bamFile=="" || cnvFile=="" || fastaFile=="" ) 
    exit( usage_check_het(argc, argv) );
  if ( ! file_exist(fastaFile) ) { cerr << fastaFile << " not found\n"; exit(0); }
  if ( ! file_exist(bamFile) ) { cerr << bamFile << " not found\n"; exit(0); }
  if ( ! file_exist(cnvFile) ) { cerr << cnvFile << " not found\n"; exit(0); }

  streambuf* sbuf = cout.rdbuf();
  ofstream FOUT;
  if ( outputFile != "STDOUT" ) {
    FOUT.open(outputFile.c_str());
    cout.rdbuf(FOUT.rdbuf());
  }
  
  ifstream FIN(cnvFile.c_str());
  if ( !FIN ) { cerr << "Can't open file " << outputFile << endl; exit(0); }
  
  vector<int> vbeg(0);
  vector<int> vend(0);
  string region="";
  i=0;
  while ( !FIN.eof() ) {
    string tmps,tmps1;
    getline(FIN,tmps);
    if ( tmps[0] == '#' ) {
      cout << tmps 
	   << "\tNC"
	   << "\tq0"
	   << "\tRD[" << dx << "]"
	   << "\tZYGOSITY" << endl;
      continue;
    }
    if ( tmps.length() <2 ) continue;
    istringstream iss(tmps);
    
    iss >> RNAME 
	>> POS1 
	>> POS2;
    if ( POS1>POS2 ) swap(POS1,POS2);
    
    region=RNAME+":"+to_string(POS1)+"-"+to_string(POS2);
    bcftoolsCommand= " bcftools view -t " + to_string(muttion_rate) + " -g - ";
    samtoolsCommand = "samtools view " + samtoolsOpt1 + " -u "+ bamFile
      + " " + region + 
      " | samtools mpileup -uIf " + fastaFile + " - | " +
      bcftoolsCommand;
    
    if ( hetFile!="" ) 
      hetCommand="awk '$1==\""+RNAME+"\""
	+" && $2>="+to_string(POS1)
	+" && $2<="+to_string(POS2)+"' "
	+ hetFile;
    else hetCommand=samtoolsCommand;
    
    dx=MINIMUMDX;
    if ( POS2-POS1>MINIMUMDX ) dx=POS2-POS1;
    int p1=POS1-dx;
    int p2=POS2+dx;
    if ( p1<1 ) p1=1;
    region=RNAME+":"+to_string(p1)+"-"+to_string(p2);
    rdCommand="samtools mpileup -A -BQ0 " + bamFile + " -r " + region +
      " | cut -f1-4 " ;
    
    region=RNAME+":"+to_string(POS1)+"-"+to_string(POS2);
    q0Command = "samtools view " + bamFile + " " + region;
    
    /*
    cerr << hetCommand << endl
	 << rdCommand << endl
	 << q0Command 
	 << endl;
	 exit(0);
    */
    
    int Ncount=0;
    vector<string> hetrname(0);
    vector<int> hetpos(0);
    vector<float> hetscore(0);
    FILE *fp;
    fp=popen( hetCommand.c_str(),"r");
    while ( !feof(fp) ) {
      string  read="";
      fgetline(fp, read);
      if ( read.length() < 3 ) continue;
      int pos,gq;
      char ref,alt,c1,c2;
      string info;
      istringstream iss(read);
      iss >> RNAME 
	  >> pos
	  >> c1
	  >> ref
	  >> alt
	  >> gq
	  >> c2
	  >> info;
      if ( pos<p1 ) continue;
      if ( pos>p2 ) break;
      
      size_t FQPOS=read.find("FQ=");
      if (FQPOS==string::npos) continue;
      string numchar="";
      for(i=FQPOS+3;i<read.size();++i) 
	if ( read[i]==';'  || 
	     read[i]==' '  || 
	     read[i]=='\t' || 
	     read[i]=='\n' || 
	     read[i]=='\r' ) break;
	else { numchar+=read[i]; }
      float FQscore=atof(numchar.c_str());
      if ( FQscore>0 && pos>POS1 && pos<POS2 ) {
	hetrname.push_back(RNAME);
	hetpos.push_back(pos);
	hetscore.push_back(FQscore);
	cerr << read << "\n" << FQscore << endl;
      }
      
    }
    fclose(fp);
    int C_het=hetpos.size();
    double FQ2=0.0;
    string taghet="";
    for(i=0;i<hetpos.size();++i) FQ2+=hetscore[i]*hetscore[i];
    if ( C_het==0 ) { taghet="HOM"; }
    else {
      FQ2=sqrt(FQ2/(double)C_het);
      taghet="HET="+to_string(C_het)+"("+to_string(FQ2)+")";
      if ( hetpos.size()<=10 ) {
	for(i=0;i<hetpos.size();++i) 
	  taghet+=","+to_string(hetpos[i])+"("+to_string(hetscore[i])+")";
      }
    }
    
    int dout=0,din=0,dp1=0,dp2=0;
    fp=popen( rdCommand.c_str(),"r");
    while ( !feof(fp) ) {
      string  read="";
      fgetline(fp, read);
      if ( read.length() < 3 ) continue;
      int pos,gq,rd;
      char ref,alt,c1,c2;
      string info;
      istringstream iss(read);
      iss >> RNAME 
	  >> pos
	  >> ref
	  >> rd;
      
      if ( pos<p1 ) continue;
      if ( pos>p2 ) break;
      
      int DP=rd;
      if ( pos < POS1 ) dout+=DP;
      if ( pos==POS1 ) dp1=DP;
      if ( pos==POS2 ) dp2=DP;
      if ( pos > POS2 ) dout+=DP;
      if ( pos>POS1 && pos<POS2 ) din+=DP;
    }
    fclose(fp);
    int Lin=POS2-POS1-1;
    if ( Lin<=1 ) Lin=1;
    dout /= (dx*2);
    din /= Lin;
    string tagrd="RDO="+to_string(dout)+";"+
      "RDI="+to_string(din)+";"+
      "RD1="+to_string(dp1)+";"+
      "RD2="+to_string(dp2);
    
    fp=popen( q0Command.c_str(),"r");
    int rCount=0;
    int q0Count=0;
    while ( !feof(fp) ) {
      string  read="";
      fgetline(fp, read);
      if ( read.length() < 3 ) continue;
      istringstream iss(read);
      int mapq;
      iss >> QNAME
	  >> FLAG
	  >> RNAME 
	  >> POS 
	  >> mapq ;
      
      ++rCount;
      if ( mapq<=MINIMUMMAPQ ) ++q0Count;
    }
    fclose(fp);
    float q0P=(float)q0Count/(float)rCount;
    string tagq0=to_string(q0P);
    
    if ( RNAME != fastaname ) {
      string RNAMEtmp=RNAME;
      if ( ci_find(RNAME,"chr") != string::npos ) RNAMEtmp=RNAME.substr(3); 
      read_fasta(fastaFile,RNAMEtmp,FASTA);
      get_N_regions(FASTA, vbeg, vend);
      
      fastaname=RNAME;
      cerr << fastaFile << "\t"
	   << fastaname << "\t"
	   << FASTA.size() << endl;
      
      cerr << "N regions:\n";
      for(int k=0;k<(int)vbeg.size();++k)
	cerr << fastaname << "\t" << vbeg[k] << "\t" << vend[k] << endl;
    }
    bool NOSEQ=false;
    if ( POS1>=FASTA.length() || POS2>=FASTA.length() ) NOSEQ=true;
    else {
      //for(int k=POS1;k<=POS2;++k) if ( FASTA[k]=='N' ) { NOSEQ=true; break; }
      for(int k=0;k<(int)vbeg.size();++k) 
	if ( POS1<=(vend[k]+1) && POS2>=(vbeg[k]+1) ) { NOSEQ=true; break; }
    }
    
    string tagN= NOSEQ ? "N":"P";
    
    cout << tmps << "\t"
	 << tagN << "\t"
	 << tagq0 << "\t" 
	 << tagrd << "\t" 
	 << taghet 
	 << endl;
    
  }
  FIN.close();
  FOUT.close();
  cout.rdbuf(sbuf);
  return(0);
}

int usage_get_het_site(int argc, char* argv[]) {
  cerr  << "\nUsage:\n" 
	<< "  " << argv[0] << " " << argv[1] << " -f REFERENCE -b BAMFILE -o OUTPUT \n"
	<< "\nOptions:\n"
	<< "  -t  FLT  mutton rate, samtools default FLT=0.001 \n"
	<< "   *   *   all others are passed to samtools \n"
	<< "\nExamples:\n"
	<< "  " << argv[0] << " " << argv[1] << " -f human_g1k_v37.fasta -b bwap40X.bam -o bwap40X.bam.het \n"
	<< "\nNote:\n"
	<< "This subroutine runs a samtools/bcftools script to call zygosity to get heterozygous sites only in the output file.\n"
	<< endl;   
  return(0);
}
int get_het_site(int argc, char* argv[])
{
  /* samtools mpileup -Q 13 -q 1 -uIf ../data/human_g1k_v37.fasta -l tmpcnvNA19240.19.txt  tmpreads.NA19240.19.bam | bcftools view -t 0.5 -g -
     FQ	int Consensus quality. Positive: sample genotypes different; negative: otherwise
     http://samtools.sourceforge.net/samtools.shtml
  */
  
  string mycommand="";
  size_t i;
  string QNAME,FLAG,RNAME,POS,MAPQ,CIGAR,MRNM,MPOS,ISIZE,SEQ,QUAL,OPT;
  int MINIMUMMAPQ=0;
  string POS1S,POS2S;
  string CIGAR1,SEQ1,QUAL1,OPT1;
  double muttion_rate=0.001;
  string outputFile="STDOUT";
  string bamFile="", region="";
  string fastaFile="",FASTA="",fastaname="";
  string samtoolsOpt1="", samtoolsOpt2="", samtoolsCommand="", bcftoolsCommand="";
  vector<string> inputArgv;
  char hostname[1024];
  hostname[1023] = '\0';
  gethostname(hostname, 1023);
  string host=hostname;
  i=host.find(".");
  if ( i!=string::npos ) host=host.substr(0,i);
  
  for(i=0;i<(size_t) argc;++i) inputArgv.push_back(string(argv[i]));
  mycommand=inputArgv[0];
  for(i=1;i<inputArgv.size();++i) mycommand+=" "+inputArgv[i];
  for(i=2;i<inputArgv.size();++i) {
    if ( inputArgv[i]=="-b" ) {  // input bam file
      bamFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
      if ( i+2<inputArgv.size() ) {
	if ( inputArgv[i+2][0]!='-' ) {
	  region=inputArgv[i+2];
	  inputArgv[i+2]="";
	}
      }
    }
    if ( inputArgv[i]=="-f" ) {  // reference file
      fastaFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-t" ) {  // mutation rate file
      muttion_rate=atof(inputArgv[i+1].c_str());
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-q" ) {  // minimum mapq
      MINIMUMMAPQ=atoi(inputArgv[i+1].c_str());
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-o" ) {  // output file
      outputFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
  }
  for(i=2;i<inputArgv.size();++i) 
    if ( inputArgv[i]!="" ) samtoolsOpt1+=inputArgv[i]+" ";
  
  if ( bamFile=="" ) {
    cerr << "Need BAM file\n";
    exit( usage_get_het_site(argc, argv) );
  }
  if ( fastaFile=="" ) {
    cerr << "Need reference file\n";
    exit( usage_get_het_site(argc, argv) );
  }
  if ( ! file_exist(fastaFile) ) { cerr << fastaFile << " not found\n"; exit(0); }
  if ( ! file_exist(bamFile) ) { cerr << bamFile << " not found\n"; exit(0); }
  
  streambuf* sbuf = cout.rdbuf();
  ofstream FOUT;
  if ( outputFile != "STDOUT" ) {
    FOUT.open(outputFile.c_str());
    cout.rdbuf(FOUT.rdbuf());
  }
  
  samtoolsCommand = "samtools view " + samtoolsOpt1 + " -u "+ bamFile;
  bcftoolsCommand= " bcftools view -t " + to_string(muttion_rate) + " -g - ";
  samtoolsCommand+=" " + region + 
    " | samtools mpileup -uIf " + fastaFile + " - | " +
    bcftoolsCommand;
  
  cerr << "#CommandL   : " << mycommand << "\n"
       << "#BAM file   : " << bamFile << "\n"
       << "#FASTA      : " << fastaFile << "\n"
       << "#SAMT exec  : " << samtoolsCommand << "\n"
       << "#OUTPUT     : " << outputFile <<"\n"
       << endl;
  
  FILE *fp;
  fp=popen( samtoolsCommand.c_str(),"r");
  while ( !feof(fp) ) {
    string  read="";
    fgetline(fp, read);
    if ( read.length() < 3 ) continue;
    int pos,gq;
    char ref,alt,c1,c2;
    string info;
    istringstream iss(read);
    iss >> RNAME 
	>> pos
	>> c1
	>> ref
	>> alt
	>> gq
	>> c2
	>> info;
    
    size_t FQPOS=read.find("FQ=");
    if (FQPOS==string::npos) continue;
    string numchar="";
    for(i=FQPOS+3;i<read.size();++i) 
      if ( read[i]==';'  || 
	   read[i]==' '  || 
	   read[i]=='\t' || 
	   read[i]=='\n' || 
	   read[i]=='\r' ) break;
      else { numchar+=read[i]; }
    float FQscore=atof(numchar.c_str());
    if ( FQscore>0 ) cout << read << "\n";
  }
  fclose(fp);
  
  FOUT.close();
  cout.rdbuf(sbuf);
  return(0);
}




int usage_check_sample_table(int argc, char* argv[]) {
  cerr << "This subroutine tabulates CNV overlaps from multiple CNV files\n"
       << "0 for overlap, D for DEL, A for DUP, V for type unknown\n"
       << "\nUsage:\n" 
       << "  " << argv[0] << " " << argv[1] << " <options> -v CNVFILE(s)\n"
       << "\nOptions  :\n"
       << "  -i  STR  IDs for CNVFILE(S) \n"
       << "  -l  INT  minimum length, INT=3 \n"
       << "  -L  INT  minimum length, INT=100000000 \n"
       << "  -t  STR  only process STR type of CNVs, if given, D for DEL or A for DUP \n"
       << "  -o  STR  outputfile, STR=STDOUT \n"
       << "\nExamples :\n"
       << "  " << argv[0] << " " << argv[1] << " -v DSE1.txt DSE2.txt -o tmp1\n"
       << endl;
  
  return(0);
}

int check_sample_table(int argc, char* argv[])
{
  
  if ( argc<3 ) exit( usage_check_sample_table(argc, argv) );
  
  string mycommand="";
  size_t i, k;
  string QNAME,FLAG,RNAME,POS,MAPQ,CIGAR,MRNM,MPOS,ISIZE,SEQ,QUAL,OPT;
  string POS1S,POS2S;
  char CNVTYPE='V';
  int minLength=3,maxLength=100000000;
  string cnvFile="",outputFile="STDOUT";
  vector<string> inputArgv;
  vector<string> samplecnv(0);
  vector<string> sampleid(0);
  vector<int> samplegenotype(0);
  
  for(i=0;i<(size_t) argc;++i) inputArgv.push_back(string(argv[i]));
  mycommand=inputArgv[0];
  for(i=1;i<inputArgv.size();++i) mycommand+=" "+inputArgv[i];
  for(i=2;i<inputArgv.size();++i) {
    if ( inputArgv[i]=="-v" ) {  // input variation file
      inputArgv[i]="";
      k=i+1;
      while( inputArgv[k][0] != '-' ) {
	samplecnv.push_back(inputArgv[k]);
	inputArgv[k]="";
	++k; if ( k==inputArgv.size() ) break;
      }
    }
    if ( inputArgv[i]=="-i" ) {  // sample ids
      inputArgv[i]="";
      k=i+1;
      while( inputArgv[k][0] != '-' ) {
	sampleid.push_back(inputArgv[k]);
	inputArgv[k]="";
	++k; if ( k==inputArgv.size() ) break;
      }
    }
    if ( inputArgv[i]=="-l" ) {  // minimum length
      minLength=atoi(inputArgv[i+1].c_str());
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-L" ) {  // max length
      maxLength=atoi(inputArgv[i+1].c_str());
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-t" ) {  // cnv type
      CNVTYPE=inputArgv[i+1][0];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-o" ) {  // output file
      outputFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
  }
  
  if ( samplecnv.size()==0 ) exit( usage_check_sample_table(argc, argv) );
  
  cerr << "#CommandL : " << mycommand << "\n";
  
  // sample id
  if ( sampleid.size() != samplecnv.size() ) {
    sampleid.empty();
    for(i=0;i<samplecnv.size();++i) {
      size_t found=samplecnv[i].rfind('/');
      string id=samplecnv[i].substr(found+1);
      found=id.find('.');
      if (found==string::npos) found=id.length();
      id=id.substr(0,found);
      sampleid.push_back(id);
      //cerr << id << "\t" << samplecnv[i] << endl;
    }
  }
  
  streambuf* sbuf = cout.rdbuf();
  ofstream FOUT;
  if ( outputFile != "STDOUT" ) {
    FOUT.open(outputFile.c_str());
    cout.rdbuf(FOUT.rdbuf());
  }
  
  
  vector<cnv_st> cnvlist(0);
  vector<cnv_st> allcnvlist(0);
  cnv_st icnv;
  ifstream FIN;
  for(i=0;i<samplecnv.size();++i) {
    if ( ! file_exist(samplecnv[i]) ) { 
      cerr << samplecnv[i] << " not found\n"; 
      exit(0); 
    }
    
    FIN.open(samplecnv[i].c_str());
    cnvlist.clear();
    while ( !FIN.eof() ) {
      string tmps,tmps1;
      getline(FIN,tmps);
      if ( tmps[0] == '#' || tmps.length() <2 ) continue;
      istringstream iss(tmps);
      iss >> icnv.RNAME >> icnv.P1 >> icnv.P2 >> tmps1;
      if ( tmps1=="D" || tmps1=="A" || tmps1=="V" ) icnv.T=tmps1[0];
      if ( icnv.P1 > icnv.P2 ) swap(icnv.P1, icnv.P2);
      if ( icnv.P2 - icnv.P1 < minLength ) continue;
      if ( icnv.P2 - icnv.P1 > maxLength ) continue;
      if ( icnv.T != CNVTYPE && CNVTYPE!='V' ) continue;
      
      icnv.ID=sampleid[i];
      cnvlist.push_back(icnv);
    }
    FIN.close();
    
    allcnvlist.insert(allcnvlist.end(), cnvlist.begin(),cnvlist.end() );
  }
  cerr << "#total cnv " << cnvlist.size() << endl;

  // sort the cnv list according to P1
  vector<int64_t> GPOSID(allcnvlist.size());
  for(k=0;k<allcnvlist.size();++k) {
    allcnvlist[k].GP1=chrompos2int64(allcnvlist[k].RNAME,allcnvlist[k].P1);
    allcnvlist[k].GP2=chrompos2int64(allcnvlist[k].RNAME,allcnvlist[k].P2);
    GPOSID[k]=allcnvlist[k].GP1;
  }
  vector<int> idx(allcnvlist.size());
  arrayindex(GPOSID,idx,1);
  vector<cnv_st> allcnvlist_tmp(allcnvlist.size());
  for(k=0;k<allcnvlist.size();++k) allcnvlist_tmp[k]=allcnvlist[ idx[k] ];
  allcnvlist=allcnvlist_tmp;
  /*
  for(k=0;k<allcnvlist.size();++k) 
    cerr << allcnvlist[k].RNAME << "\t" 
	 << allcnvlist[k].P1 << "\t"
	 << allcnvlist[k].GP1 << endl;
  */
  cerr << "#total from all samples : " << allcnvlist.size() << endl;
  
  // compact the list, remove duplicate and mostly overlapped
  allcnvlist_tmp.clear();
  allcnvlist_tmp.push_back(allcnvlist[0]);
  for(k=1;k<allcnvlist.size();++k) {
    
    bool overlap=false;
    if ( max( allcnvlist[k].GP1, allcnvlist_tmp.back().GP1 ) <=
	 min( allcnvlist[k].GP2, allcnvlist_tmp.back().GP2 ) ) overlap=true;
    if ( !overlap ) {
      allcnvlist_tmp.push_back(allcnvlist[k]);
      continue;
    }
    
    double overLen=min( allcnvlist[k].GP2, allcnvlist_tmp.back().GP2 ) - 
      max( allcnvlist[k].GP1, allcnvlist_tmp.back().GP1 );
    if ( overLen > ( allcnvlist_tmp.back().GP2-allcnvlist_tmp.back().GP1 )*0.8 &&
	 overLen > ( allcnvlist[k].GP2-allcnvlist[k].GP1 )*0.8 ) {
      icnv=allcnvlist_tmp.back();
      icnv.P1=min( allcnvlist[k].P1, allcnvlist_tmp.back().P1 );
      icnv.GP1=min( allcnvlist[k].GP1, allcnvlist_tmp.back().GP1 );
      icnv.P2=max( allcnvlist[k].P2, allcnvlist_tmp.back().P2 );
      icnv.GP2=max( allcnvlist[k].GP2, allcnvlist_tmp.back().GP2 );
      allcnvlist_tmp.back()=icnv;
    }
    else { allcnvlist_tmp.push_back(allcnvlist[k]); }
  }
  allcnvlist=allcnvlist_tmp;
  /*
  for(k=0;k<allcnvlist.size();++k) 
    cerr << allcnvlist[k].RNAME << "\t" 
	 << allcnvlist[k].P1 << "\t" 
	 << allcnvlist[k].P2 << "\t" 
	 << allcnvlist[k].GP1 << "\t" 
	 << allcnvlist[k].GP2 << "\t" 
	 << endl;
  */
  cerr << "#total after overlaps combined : " << allcnvlist.size() << endl;
  // exit(0);
  
  // read cnv again and check overlap
  size_t ncnv=allcnvlist.size();
  size_t nids=samplecnv.size();
  vector< vector<char> > cnvtype(nids);
  for(i=0;i<nids;++i) {
    cnvtype[i].resize(ncnv);
    cnvtype[i]=vector<char>(ncnv,'0');
  }
  vector< vector<char> > cnvoverlap(ncnv, vector<char>(nids,'0'));
  
  for(i=0;i<samplecnv.size();++i) {
    
    FIN.open(samplecnv[i].c_str());
    while ( !FIN.eof() ) {
      string tmps,tmps1;
      getline(FIN,tmps);
      if ( tmps[0] == '#' || tmps.length() <2 ) continue;
      istringstream iss(tmps);
      
      iss >> icnv.RNAME >> icnv.P1 >> icnv.P2 >> tmps1;
      if ( tmps1=="D" || tmps1=="A" || tmps1=="V" ) icnv.T=tmps1[0];
      if ( icnv.P1 > icnv.P2 ) swap(icnv.P1, icnv.P2);
      if ( icnv.P2 - icnv.P1 < minLength ) continue;
      if ( icnv.P2 - icnv.P1 > maxLength ) continue;
      if ( icnv.T != CNVTYPE && CNVTYPE!='V' ) continue;
      icnv.GP1=chrompos2int64(icnv.RNAME,icnv.P1);
      icnv.GP2=chrompos2int64(icnv.RNAME,icnv.P2);
      
      bool overlap=false;
      for(k=0;k<allcnvlist.size();++k) {
	if ( max( allcnvlist[k].GP1, icnv.GP1 ) <=
	     min( allcnvlist[k].GP2, icnv.GP2 ) ) {
	  overlap=true;
	  cnvoverlap[k][i]=icnv.T;
	  cnvtype[i][k]=icnv.T;
	}
	if ( allcnvlist[k].GP1 > icnv.GP2 ) break;
      }

    }
    FIN.close();
    
  }
  
  vector<size_t> sumocnv(samplecnv.size(),0);
  vector<size_t> sumoid(allcnvlist.size(),0);
  for(k=0;k<allcnvlist.size();++k) 
    for(i=0;i<samplecnv.size();++i) 
      if ( cnvoverlap[k][i]!='0' ) { 
	sumoid[k]++; 
	sumocnv[i]++;
      }
  
  cout << "##sampleid_and_inputfile\n";
  for(i=0;i<samplecnv.size();++i) cout << "##" << sampleid[i] << "\t" << samplecnv[i] << "\n";
  
  cout << "#RNAME\t" << "POS1\t" << "POS2\t" << "TYPE\t" << "ROWSUM";
  for(i=0;i<samplecnv.size();++i) cout << "\t" << sampleid[i];
  cout << "\n";
  cout << "##CSUM\t" << "POS1\t" << "POS2\t" << "TYPE\t" << "ROWSUM";
  for(i=0;i<samplecnv.size();++i) cout << "\t" << sumocnv[i];
  cout << endl;
  
  for(k=0;k<allcnvlist.size();++k) {
    cout << allcnvlist[k].RNAME << "\t"
	 << allcnvlist[k].P1 << "\t"
	 << allcnvlist[k].P2 << "\t"
      // << allcnvlist[k].GP1 << "\t"
	 << allcnvlist[k].T << "\t"
	 << sumoid[k];
    for(i=0;i<samplecnv.size();++i) cout << "\t" << cnvoverlap[k][i];
    cout << endl;
  }
  
  FOUT.close();
  cout.rdbuf(sbuf);
  return(0);
}


int usage_check_sample_subset(int argc, char* argv[]) {
  cerr << "This subroutine subset CNVs from input files within given regions\n\n";
  cerr << "Usage:\n" 
       << "  " << argv[0] << " " << argv[1] << " <options> \n"
       << "\nOptions:\n"
       << "  -overlap true if overlap, default is true if inside \n"
       << "  -v  STR  input variation files \n"
       << "  -r  STR  regions given as chr:start-end or filename  \n"
       << "  -o  STR  outputfile, STR=STDOUT \n"
       << "  -d       verbose \n"
       << "\nExamples :\n"
       << "  " << argv[0] << " " << argv[1] << " -v DSE*.pm DSE2.pm -r chr19:123-456\n"
       << "  " << argv[0] << " " << argv[1] << " -v DSE*.pm DSE2.pm -r reg.bed\n"
       << endl;
  
  return(0);
}
int check_sample_subset(int argc, char* argv[])
{
  
  if ( argc<3 ) exit( usage_check_sample_subset(argc, argv) );
  
  string mycommand="";
  size_t i, k;
  string QNAME,FLAG,RNAME,POS,MAPQ,CIGAR,MRNM,MPOS,ISIZE,SEQ,QUAL,OPT;
  string POS1S,POS2S;
  string cnvFile="",outputFile="STDOUT";
  bool OVERLAP=false;
  bool VERBOSE=false;
  vector<string> inputArgv;
  vector<string> samplecnv(0);
  vector<string> sampleregion(0);
  
  for(i=0;i<(size_t) argc;++i) inputArgv.push_back(string(argv[i]));
  mycommand=inputArgv[0];
  for(i=1;i<inputArgv.size();++i) mycommand+=" "+inputArgv[i];
  for(i=2;i<inputArgv.size();++i) {
    if ( inputArgv[i]=="-v" ) {  // input variation file
      inputArgv[i]="";
      for(k=i+1;k<inputArgv.size();++k) {
	if ( inputArgv[k][0]=='-' ) break;
	samplecnv.push_back(inputArgv[k]);
      }
    }
    if ( inputArgv[i]=="-r" ) {  // region file
      inputArgv[i]="";
      for(k=i+1;k<inputArgv.size();++k) {
	if ( inputArgv[k][0]=='-' ) break;
	sampleregion.push_back(inputArgv[k]);
      }
    }
    if ( inputArgv[i]=="-o" ) {  // output file
      outputFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-overlap" ) {  // output file
      inputArgv[i]="";
      OVERLAP=true;
    }
    if ( inputArgv[i]=="-d" ) {  // output file
      inputArgv[i]="";
      VERBOSE=true;
    }
  }
  
  cerr << "#CommandL : " << mycommand << "\n";

  if ( samplecnv.size()==0 ) exit( usage_check_sample_subset(argc, argv) );
  if ( sampleregion.size()==0 ) exit( usage_check_sample_subset(argc, argv) );
  
  streambuf* sbuf = cout.rdbuf();
  ofstream FOUT;
  if ( outputFile != "STDOUT" ) {
    FOUT.open(outputFile.c_str());
    cout.rdbuf(FOUT.rdbuf());
  }
  
  /* load regions either from command line or file */  
  cnv_st icnv;
  ifstream FIN;
  vector<cnv_st> region(0);
  for(i=0;i<sampleregion.size();++i) {
    if ( sampleregion[i].find(":") != string::npos ) {
      size_t foundc=sampleregion[i].find(":");
      size_t foundt=sampleregion[i].find("-");
      icnv.RNAME=sampleregion[0].substr(0,foundc);
      icnv.P1=0;
      icnv.P2=255555555;
      if ( foundt!=string::npos ) {
	icnv.P1=atoi(sampleregion[0].substr(foundc+1,foundt-foundc-1).c_str());
	icnv.P2=atoi(sampleregion[0].substr(foundt+1).c_str());
      }
      if ( icnv.P1 > icnv.P2 ) swap(icnv.P1, icnv.P2);
      region.push_back(icnv);
    }
    else {
      if ( ! file_exist(sampleregion[i]) ) {
	cerr << "#File " << sampleregion[i] << " not found\n";
	exit(0);
      }
      FIN.open(sampleregion[i].c_str());
      
      while ( !FIN.eof() ) {
	string tmps,tmps1;
	getline(FIN,tmps);
	if ( tmps[0] == '#' || tmps.length() <2 ) continue;
	istringstream iss(tmps);
	icnv.P1=0;
	icnv.P2=255555555;
	iss >> icnv.RNAME >> icnv.P1 >> icnv.P2;
	if ( icnv.P1 > icnv.P2 ) swap(icnv.P1, icnv.P2);
	region.push_back(icnv);
      }
      FIN.close();
    }
  }
  
  if ( VERBOSE ) {
    cerr << "#REGIONS:\n" ;
    for(k=0;k<region.size();++k) {
      cerr << region[k].RNAME << "\t"
	   << region[k].P1 << "\t"
	   << region[k].P2 << endl;
    }
  }
  
  for(i=0;i<samplecnv.size();++i) {
    
    cout << "##SUBSETFROM " << samplecnv[i] << endl;
    
    if ( ! file_exist(samplecnv[i]) ) {
      cerr << "File " << samplecnv[0] << " not found\n";
      exit(0);
    }
    FIN.open(samplecnv[i].c_str());
    while ( !FIN.eof() ) {
      string tmps,tmps1;
      getline(FIN,tmps);
      if ( tmps.length() <2 ) continue;
      if ( tmps[0] == '#' ) {
	cout << tmps << endl; 
	continue;
      }
      istringstream iss(tmps);
      iss >> icnv.RNAME >> icnv.P1 >> icnv.P2;
      if ( icnv.P1 > icnv.P2 ) swap(icnv.P1, icnv.P2);
      
      bool accept=false;
      for(k=0;k<region.size();++k) {
	if ( region[k].RNAME != icnv.RNAME ) continue;
	// inside
	if ( icnv.P1>=region[k].P1 && icnv.P2<=region[k].P2 ) {
	  accept=true;
	  break;
	}
	// overlap
	if ( OVERLAP ) {
	  if ( max( region[k].P1, icnv.P1 ) <=
	       min( region[k].P2, icnv.P2 ) ) {
	    accept=true;
	    break;
	  }
	}
      }

      if ( accept ) cout << tmps << endl;
    }
    FIN.close();
    
  }
  
  
  
  FOUT.close();
  cout.rdbuf(sbuf);
  return(0);
}

/* sort cnv according to P1 */
void sort_cnv(vector<cnv_st>& cnvlist, vector<cnv_st>& cnv_sorted)
{
  
  vector<cnv_st> allcnvlist(0);
  cnv_st icnv;
  size_t i;
  
  // convert to numbers
  vector<int64_t> GPOSID(cnvlist.size());
  for(i=0;i<cnvlist.size();++i) {
    cnvlist[i].GP1=chrompos2int64(cnvlist[i].RNAME,cnvlist[i].P1);
    cnvlist[i].GP2=chrompos2int64(cnvlist[i].RNAME,cnvlist[i].P2);
    GPOSID[i]=cnvlist[i].GP1;
  }
  
  cnv_sorted.empty();
  cnv_sorted.resize(cnvlist.size());
  vector<int> idx(cnvlist.size());
  arrayindex(GPOSID,idx,1);
  for(i=0;i<cnvlist.size();++i) cnv_sorted[i]=cnvlist[ idx[i] ];
  
  return;
  
}

int recode(string& genotype, string FORMAT, int64_t gipos, 
	    string ID, vector<cnv_st>& cnvlist, int lastcnv)
{
  int i;
  string tmps;
  char cnvtype='N';
  char GT[3]={'.','.',0};
  vector<int> AD(2); 
  vector<float> PL(3); 
  vector<string> vgt(1);// vgt[0]=".";
  vector<string> vfmt; 
  int iAD=-1,iPL=-1;
  
  // find cnv type
  for(i=lastcnv;i<(int)cnvlist.size();++i) {
    if ( cnvlist[i].GP1>gipos ) break;
    if ( cnvlist[i].GP1<=gipos && cnvlist[i].GP2>=gipos && cnvlist[i].ID==ID ) 
      cnvtype=cnvlist[i].T;
  }
  
  //  if ( cnvtype=='N' ) goto recodedone;
  if ( cnvtype=='N' ) return 0;
  
  stringsplit(genotype, vgt, ':');
  //  if ( vgt[0].find(".") != string::npos ) goto recodedone;
  if ( vgt[0].find(".") != string::npos ) return 0;

  stringsplit(FORMAT, vfmt, ':');
  if ( vgt.size()!=vfmt.size() ) {
    if ( vgt[0].find(".") == string::npos ) {
      cerr << "Error read geotype and FORMAT\n"
	   << genotype << "\n" << FORMAT 
	   << endl;
      exit(0);
    }
  }
  
  GT[0]=vgt[0][0];
  GT[1]=vgt[0][2];

  if ( GT[0]==GT[1] ) {
    if ( cnvtype == 'A' ) vgt[0][1]=vgt[0][0];
    else if ( cnvtype == 'D' ) vgt[0]=vgt[0][0];
    goto recodedone;
  }
  
  iAD= find(vfmt.begin(),vfmt.end(),"AD") - vfmt.begin();
  if ( iAD != (int)vfmt.size() && (int)vgt.size()>iAD ) {
    stringsplit(vgt[iAD], AD, ',');
    if ( AD.size() != 2 ) {
      cerr << "does not look like AD:\t" << vgt[iAD] << "\n";
      exit(0);
    }
    if ( cnvtype=='D' ) {
      if ( AD[0]>AD[1] ) vgt[0]=GT[0];
      else if ( AD[1]>AD[0] ) vgt[0]=GT[1];
    }
    if ( cnvtype=='A' ) {
      if ( AD[0]>AD[1] ) vgt[0][1]=GT[0];
      else if ( AD[1]>AD[0] ) vgt[0][1]=GT[1];
    }
    goto recodedone;
  }
  
  iPL= find(vfmt.begin(),vfmt.end(),"PL") - vfmt.begin();
  if ( iPL != (int)vfmt.size() && (int)vgt.size()>iPL ) {
    stringsplit(vgt[iPL], PL, ',');
    if ( PL.size() != 3 ) {
      cerr << "does not look like PL:\t" << vgt[iPL] << "\n";
      exit(0);
    }
    if ( cnvtype=='D' ) {
      if ( PL[0]<PL[2] ) vgt[0]=GT[0];
      else if ( PL[0]>PL[2] ) vgt[0]=GT[1];
    }
    if ( cnvtype=='A' ) {
      if ( PL[0]<PL[2] ) vgt[0][1]=GT[0];
      else if ( PL[0]>PL[2] ) vgt[0][1]=GT[1];
    }
    goto recodedone;
  }
  
 recodedone:  
  string genotype1=vgt[0];
  for(i=1;i<(int)vgt.size();++i) genotype1+=":"+vgt[i];
  // cout << cnvtype << "\t" << GT << "\t" << vgt[0] << "\t" << genotype << "\t" 
  //    << AD[0] << "/" << AD[1] << "\t" << genotype1 << endl;
  
  genotype=genotype1;
  return 1;
}


int usage_patch_vcf(int argc, char* argv[]) {
  cerr << "This subroutine modifies a VCF file according to samples' CNV files\n"
       << "For DEL, 0/0->0, 1/1->1, 0/1->0(or 1), depending on AL or PL\n"
       << "For DUP, 0/0->000, 1/1->111, 0/1->00(or 1)1, depending on AL or PL\n"
       << "\nUsage:\n" 
       << "  " << argv[0] << " " << argv[1] 
       << " -vcf VCFFILE -i ID -v CNVFILE <options> \n"
       << "\nOptions:\n"
       << "  -i  STR  sample ID \n"
       << "  -v  STR  sample ID's cnv \n"
       << "  -L  INT  CNVs longer than INT are ignored, INT=100000000\n"
       << "  -s       subset modified; default is append at next line\n"
       << "  -f  STR  recode genotype in STR formart\n, default is same as input\n"
       << "  -o  STR  outputfile, STR=STDOUT \n"
       << "\nExamples:\n"
       << "  " << argv[0] << " " << argv[1] << " -vcf AMD.vcf -i DSE1 -v DSE1.pm -i DSE2 -v DSE2.pm -o tmp.vcf \n"
       << "\nNote :\n"
       << "  -v must immediately follow -i\n"
       << "  ID must be in the VCF file's header file\n"
       << "  multiple -i -v combinations can be given; order does not matter\n"
       << "  SNPIDs for cnv adjusted SNPs are changed to YW$SNPID\n"
       << endl;
  
  return(0);
}
int patch_vcf(int argc, char* argv[])
{
  if ( argc<3 ) exit( usage_patch_vcf(argc, argv) );
  
  string CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,OFORMAT="";
  string mycommand="";
  size_t i;
  string vcfFile="",cnvFile="",outputFile="STDOUT";
  bool SUBSET=false;
  string TAG="YW";
  int Lmax=100000000;
  int NSNP=0,CSNP=0;
  vector<string> inputArgv;
  vector<string> cnvid(0);
  vector<string> cnvfile(0);
  
  vector<string> vcfid(0);
  
  vector<string> COLUMNID(9);
  COLUMNID[0]="CHROM";
  COLUMNID[1]="POS";
  COLUMNID[2]="ID";     
  COLUMNID[3]="REF";
  COLUMNID[4]="ALT";   
  COLUMNID[5]="QUAL";   
  COLUMNID[6]="FILTER"; 
  COLUMNID[7]="INFO";
  COLUMNID[8]="FORMAT";
  
  for(i=0;i<(size_t) argc;++i) inputArgv.push_back(string(argv[i]));
  mycommand=inputArgv[0];
  for(i=1;i<inputArgv.size();++i) mycommand+=" "+inputArgv[i];
  for(i=2;i<inputArgv.size();++i) {
    if ( inputArgv[i]=="-vcf" ) {  // input vcf file
      vcfFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
      continue;
    }
    if ( inputArgv[i]=="-i" ) {  // sample id
      if ( inputArgv[i+2] != "-v" ) exit( usage_patch_vcf(argc, argv) );
      cnvid.push_back(inputArgv[i+1]);
      cnvfile.push_back(inputArgv[i+3]);
      inputArgv[i]="";
      inputArgv[i+1]="";
      inputArgv[i+2]="";
      inputArgv[i+3]="";
      continue;
    }
    if ( inputArgv[i]=="-L" ) {  // CNV length
      Lmax=atoi(inputArgv[i+1].c_str());
      inputArgv[i]="";
      inputArgv[i+1]="";
      continue;
    }
    if ( inputArgv[i]=="-s" ) {  // region file
      SUBSET=true;
      inputArgv[i]="";
      continue;
    }
    if ( inputArgv[i]=="-f" ) {  // output format
      OFORMAT=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
      continue;
    }
    if ( inputArgv[i]=="-o" ) {  // output file
      outputFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
      continue;
    }
  }
  
  cerr << "#CommandL : " << mycommand << "\n";
  
  if ( vcfFile=="" ) exit( usage_check_sample_subset(argc, argv) );
  if ( cnvfile.size()==0 ) exit( usage_check_sample_subset(argc, argv) );
  
  streambuf* sbuf = cout.rdbuf();
  ofstream FOUT;
  if ( outputFile != "STDOUT" ) {
    FOUT.open(outputFile.c_str());
    cout.rdbuf(FOUT.rdbuf());
  }
  
  /* load CNV regions from files */  
  cnv_st icnv;
  ifstream FIN;
  vector<cnv_st> cnvlist(0);
  for(i=0; i<cnvfile.size(); ++i) {
    if ( ! file_exist( cnvfile[i] ) ) {
      cerr << cnvfile[i] << " not found\n";
      exit(0);
    }
    FIN.open( cnvfile[i].c_str() );
    while ( !FIN.eof() ) {
      string tmps,tmps1;
      getline(FIN,tmps);
      if ( tmps[0] == '#' || tmps.length() <2 ) continue;
      istringstream iss(tmps);
      icnv.ID=cnvid[i];
      icnv.P1=0;
      icnv.P2=255555555;
      iss >> icnv.RNAME >> icnv.P1 >> icnv.P2 >> icnv.T >> icnv.un;
      if ( icnv.P1 > icnv.P2 ) swap(icnv.P1, icnv.P2);
      icnv.P1-=icnv.un;
      if ( icnv.P2 - icnv.P1 > Lmax ) continue;
      cnvlist.push_back(icnv);
    }
    FIN.close();
  }
  
  vector<cnv_st> cnvlist_sort;
  sort_cnv(cnvlist, cnvlist_sort);
  cnvlist=cnvlist_sort;
  
  vector<string> cell(0);
  vector<int> cols(0);
  int64_t pregpos=0;

  int lastcnv=0;
  if ( ! file_exist( vcfFile ) ) {
    cerr << vcfFile << " not found\n";
    exit(0);
  }
  FIN.open( vcfFile.c_str() );
  while ( !FIN.eof() ) {
    string tmps,tmps1;
    getline(FIN,tmps);
    if ( tmps.length() < 2 ) continue;
    if ( tmps.substr(0,2) == "##" ) {
      cout << tmps << endl;
      continue;
    }
    
    istringstream iss(tmps);
    
    int icol=0;
    if ( tmps[0] == '#' ) {
      cout << tmps << endl; 
      if ( ci_find(tmps, "CHROM") == string::npos ) continue;
      bool is_ID=false;
      while( iss >> tmps1 ) {
	if ( is_ID ) {
	  vcfid.push_back(tmps1);
	  ++icol;
	  COLUMNID.push_back(tmps1);
	  int ifind=find(cnvid.begin(),cnvid.end(),tmps1)-cnvid.begin();
	  if ( ifind != (int)cnvid.size() ) cols.push_back(icol); 
	}
	else if ( ci_equal(tmps1,"FORMAT") ) { 
	  is_ID=true; 
	  icol=8; 
	}
      }
      if ( cols.size() == 0 ) {
	cerr << "check IDs in vcf file and IDs in command line\n";
	exit(0);
      }
      cerr << cols.size() << " id found out of " << vcfid.size() << endl;
      //for(i=0;i<cols.size();++i) 
      //cerr << cols[i] << "\t" << COLUMNID[ cols[i] ] << endl;
      continue;
    }
    
    icol=0;
    cell.resize(0);
    while ( iss >> tmps1 ) cell.push_back(tmps1);
    
    CHROM=cell[0];
    POS=cell[1];
    ID=cell[2];     
    REF=cell[3];
    ALT=cell[4];   
    QUAL=cell[5];   
    FILTER=cell[6]; 
    INFO=cell[7];
    FORMAT=cell[8]; 
    
    int ichr=chrom2int(CHROM);
    int ipos=atoi(POS.c_str());
    int64_t gipos=ichr*(int)1E10+ipos;
    
    bool snpincnv=false;
    for(i=lastcnv; i<cnvlist.size(); ++i) {
      if ( gipos>cnvlist[i].GP1 && gipos<cnvlist[i].GP2 )  { 
	snpincnv=true; 
	lastcnv=i;
	break; 
      }
      if ( gipos < cnvlist[i].GP1 ) break;
      lastcnv=i;
    }
    if ( gipos<pregpos ) { // vcf not ordered
      for(i=0;i<=(size_t)lastcnv;++i) {
	if ( gipos>cnvlist[i].GP1 && gipos<cnvlist[i].GP2 )  { 
	  snpincnv=true; 
	  lastcnv=i;
	  break; 
	}
	if ( gipos < cnvlist[i].GP1 ) break;
	lastcnv=i;
      }
    }
    
    if ( !SUBSET ) {
      if ( OFORMAT=="" ) cout << tmps << "\n";
      else {
	for(i=0;i<8;++i) cout << cell[i] << "\t";
	cout << OFORMAT; //cout << cell[8];
	for(i=9;i<cell.size();++i)
	  cout << "\t" << getvcfgenotypefield(cell[i], FORMAT, OFORMAT);
	cout << endl;
      }
    }
    
    ++NSNP;
    if ( !snpincnv ) continue;
    ++CSNP;
    
    //cout << ichr << "\t" << ipos << "\t" << snpincnv;
    //for(i=0;i<cols.size();++i) cout << "\t" << cell[cols[i]]; 
    //cout << endl;
    int modifiedcount=0;
    for(i=0;i<cols.size();++i) {
      string genotype=cell[ cols[i] ]; 
      modifiedcount+=
	recode(genotype, FORMAT, gipos, COLUMNID[ cols[i] ], cnvlist, lastcnv);
      cell[ cols[i] ] = genotype;
    }
    //if (modifiedcount==0) cerr << "snp not modified\n" << endl; 
    //cout << ichr << "\t" << ipos << "\t" << snpincnv;
    //for(i=0;i<cols.size();++i) cout << "\t" << cell[cols[i]]; 
    //cout << endl;
    
    if (modifiedcount>0) { 
      cout << CHROM << "\t"
	   << POS << "\t"
	   << ID << TAG << "\t"
	   << REF << "\t"
	   << ALT << "\t"   
	   << QUAL << "\t"   
	   << FILTER << "\t"
	   << INFO << "\t";
      
      if ( OFORMAT=="" ) {
	cout << FORMAT; 
	for(i=9;i<cell.size();++i) cout << "\t" << cell[i];
      }
      else  {
	cout << OFORMAT; 
	for(i=9;i<cell.size();++i)
	  cout << "\t" << getvcfgenotypefield(cell[i], FORMAT, OFORMAT);
      }
      
      cout << endl;
    }
    
  }
  FIN.close();
  FOUT.close();
  cout.rdbuf(sbuf);
  
  //  for(i=0; i<COLUMNID.size(); ++i) cerr << i << "\t" << COLUMNID[i] << endl;  
  cerr << CSNP << " found out of " << NSNP << endl;
  cerr << cols.size() << " id found out of " << vcfid.size() << endl;
  return(0);
}


/*!-------------------------------------------------------------------!*/
int usage_check_pairs(int argc, char* argv[]) {
  cerr << "This subroutine checks if there are pairs that cover the cnv break points\n"
       << "For deletion, we look for templates whose 5' point is lower than CNV 5' point and whose 3' is higher than CNV's 3'\n"
       << "For duplication, we look for templates with negative(positive) TLENgth for forward(reverse compliment) reads.\n"
       << "The result is append at the end of a line in the format\n"
       << "   evident_reads,total_paired_reads,total_reads \n"
       << "\nUsage:\n" 
       << "  " << argv[0] << " " << argv[1] << " <options> -b BAMFILE -v CNVFILE \n"
       << "\nOptions:\n"
       << "  -d  INT  check reads INT before and after a CNV, INT=1000 \n"
       << "  -o  STR  outputfile, STR=STDOUT \n"
       << "   *       all other options are passed to samtools \n" 
       << "\nExamples :\n"
       << "  " << argv[0] << " " << argv[1] <<  " -b bwap40X.bam -v bwap40X.bam.bp -q 10\n"
       << "\nNote:\n"
       << "   This subroutine is useless for single ended mapping\n"
       << "   use  | awk '$NF~/^0,/' to find CNVs with no pairs support\n"
       << endl;
  
  return(0);
}
int check_pairs(int argc, char* argv[])
{
  
  string mycommand="";
  size_t i;
  string QNAME,RNAME,CIGAR,MRNM,SEQ,QUAL,OPT;
  int FLAG,POS,MAPQ,MPOS,ISIZE;
  char TYPE;
  int POS1,POS2,UN,DIS=1000;
  string POS1S,POS2S;
  string CIGAR1,SEQ1,QUAL1,OPT1;
  
  string cnvFile="",bamFile="",outputFile="STDOUT";
  string samtoolsCommand="",samtoolsOpt1="", samtoolsOpt2="";
  vector<string> inputArgv;
  
  for(i=0;i<(size_t) argc;++i) inputArgv.push_back(string(argv[i]));
  mycommand=inputArgv[0];
  for(i=1;i<inputArgv.size();++i) mycommand+=" "+inputArgv[i];
  for(i=2;i<inputArgv.size();++i) {
    if ( inputArgv[i]=="-v" ) {  // input variation file
      cnvFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
      continue;
    }
    if ( inputArgv[i]=="-b" ) {  // input bam file
      bamFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
      continue;
    }
    if ( inputArgv[i]=="-d" ) {  // input bam file
      DIS=atoi(inputArgv[i+1].c_str());
      inputArgv[i]="";
      inputArgv[i+1]="";
      continue;
    }
    if ( inputArgv[i]=="-o" ) {  // output file
      outputFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
      continue;
    }
  }
  for(i=2;i<inputArgv.size();++i) 
    if ( inputArgv[i]!="" ) samtoolsOpt1+=inputArgv[i]+" ";
  samtoolsCommand="samtools view " + samtoolsOpt1 + bamFile;
  
  if ( bamFile=="" || cnvFile=="" ) exit( usage_check_pairs(argc, argv) );
  if ( ! file_exist(bamFile) ) { cerr << bamFile << " not found\n"; exit(0); }
  if ( ! file_exist(cnvFile) ) { cerr << cnvFile << " not found\n"; exit(0); }
  
  cerr << "#CommandL  : " << mycommand << "\n"
       << "#CNV file  : " << cnvFile << "\n"
       << "#BAM file  : " << bamFile << "\n"
       << "#SAMT opts : " << samtoolsOpt1 << "\n"
       << "#SAMT exec : " << samtoolsCommand << "\n"
       << endl;
  
  streambuf* sbuf = cout.rdbuf();
  ofstream FOUT;
  if ( outputFile != "STDOUT" ) {
    FOUT.open(outputFile.c_str());
    cout.rdbuf(FOUT.rdbuf());
  }
  
  ifstream FIN(cnvFile.c_str());
  if ( !FIN ) { cerr << "Can't open file " << outputFile << endl; exit(0); }
  i=0;
  while ( !FIN.eof() ) {
    string tmps,tmps1;
    getline(FIN,tmps);
    istringstream iss(tmps);
    if ( tmps[0] == '#' ) {
      cout << tmps << "\tPAIRS" << endl;
      continue;
    }
    if ( tmps.length() <2 ) continue;
    
    iss >> RNAME >> POS1 >> POS2 >> TYPE >> UN;
    
    if ( POS1>POS2 ) swap(POS1,POS2);
    int p1e=POS1-DIS, p2e=POS2+DIS;
    if ( p1e<1 ) p1e=1;
    
    string region=RNAME+":"+to_string(p1e)+"-"+to_string(p2e);
    string SAMEXEC=samtoolsCommand+" "+region;
    
    FILE *fp;
    fp=popen( SAMEXEC.c_str(),"r");
    if(!fp) {
      cerr << "can't open pipe for input \n"
	   << SAMEXEC << endl;
      exit(0);
    }
    
    int readcount=0;
    int paircount=0;
    int supportcount=0;
    while ( !feof(fp) ) {
      string  read="";
      fgetline(fp, read);
      if ( read.length() < 3 ) continue;
      iss.clear();
      iss.str(read);
      iss >> QNAME
	  >> FLAG
	  >> RNAME 
	  >> POS 
	  >> MAPQ 
	  >> CIGAR 
	  >> MRNM 
	  >> MPOS 
	  >> ISIZE
	  >> SEQ;
      //>> QUAL 
      //>> OPT;
      ++readcount;
      
      if ( QNAME == "*" ) continue;
      if ( POS == 0 ) continue;
      if ( MPOS == 0 ) continue;
      if ( ISIZE == 0 ) continue;
      if ( (MRNM != "=") && (MRNM != QNAME) ) continue;
      ++paircount;
      
      int r1,r2;
      if ( (FLAG & 0x10)==0 ) {   // read is Forward 
	r1= POS;
	r2= MPOS+SEQ.size();
      }
      else {  // read is Reverse
	r1= MPOS; 
	r2=POS+SEQ.size();
	ISIZE=-ISIZE;
      }
      
      //r1=min(POS,MPOS);
      //r2=max(POS+SEQ.size(),MPOS+SEQ.size());
      if ( TYPE=='D' && r1<POS1 && r2>POS2 ) ++supportcount;
      if ( TYPE=='A' && ISIZE<0  ) ++supportcount;
    }
    fclose(fp);
    
    string tag=to_string(supportcount)+","
      +to_string(paircount)+","
      +to_string(readcount);
    
    cout << tmps << "\t" << tag << endl;
  }
  FIN.close();

  FOUT.close();
  cout.rdbuf(sbuf);
  return(0);
}

/*!-------------------------------------------------------------------!*/
int usage_variation_to_fastq(int argc, char* argv[]) {
  cerr << "This subroutine extract merged reads in a CNV file and convert them to fastq file for mapping to check whether the merged sequence belongs to a different position than suggested by the CNV.\n"
       << "\nUsage:\n" 
       << "  " << argv[0] << " " << argv[1] << " <options> -f BWAREF -v CNVFILE \n"
       << "\nOptions:\n"
       << "  -o  STR  outputfile, STR=BAMFILE.map{.fa .sam .log} \n"
       << "  bwa=STR  bwa executable \n"
       << "  -s  STR  aligned SAM file, if given, alignment not performed \n"
       << "  -0       do not use secondary alignments \n"
       << "  -d       output debug information \n"
       << "   *       all other options are passed to bwa \n" 
       << "\nExamples :\n"
       << "  " << argv[0] << " " << argv[1] <<  " bwa=bwa070 -f $REF -v -v bwap40X.bam.bp \n"
       << endl;
  
  return(0);
}
int variation_to_fastq(int argc, char* argv[])
{
  
  size_t i;
  string mycommand,BWA="";
  bool ISBWA=false;
  string QNAME,RNAME,CIGAR,MRNM,SEQ,OPT;
  int POS,FLAG,MAPQ,MPOS,ISIZE;
  string CNVNAME;
  string POS1S,POS2S;
  
  bool VERBOSE=false, SECONDARY=true;
  string COMMENT="notMatch";
  string HEADER="";
  string cnvFile="",bamFile="",samFile="",outputFile="STDOUT",refFile="";
  string bwaoutputFile="", fastqFile="", logFile="";
  int iField=-1;
  string bwaCommand="",bwaOpt1="";
  vector<string> inputArgv;
  
  for(i=0;i<(size_t) argc;++i) inputArgv.push_back(string(argv[i]));
  mycommand=inputArgv[0];
  for(i=1;i<inputArgv.size();++i) mycommand+=" "+inputArgv[i];
  for(i=2;i<inputArgv.size();++i) {
    //    if ( inputArgv[i].find("=") != string::npos ) {  // bwa command exec
    if ( ci_find(inputArgv[i], "bwa=") != string::npos ) {  // bwa command exec
      BWA=inputArgv[i].substr(inputArgv[i].find("=")+1);
      inputArgv[i]="";
      continue;
    }
    if ( inputArgv[i]=="-v" ) {  // input variation file
      cnvFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
      continue;
    }
    if ( inputArgv[i]=="-b" ) {  // input bam file
      bamFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
      continue;
    }
    if ( inputArgv[i]=="-s" ) {  // alignment file for fa
      samFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
      continue;
    }
    if ( inputArgv[i]=="-f" ) {  // reference file
      refFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
      continue;
    }
    if ( inputArgv[i]=="-i" ) {  // input field
      iField=atoi(inputArgv[i+1].c_str())-1;
      inputArgv[i]="";
      inputArgv[i+1]="";
      continue;
    }
    if ( inputArgv[i]=="-o" ) {  // output file
      outputFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
      continue;
    }
    if ( inputArgv[i]=="-d" ) {  // debug
      VERBOSE=true;
      inputArgv[i]="";
      continue;
    }
    if ( inputArgv[i]=="-0" ) {  // debug
      SECONDARY=false;
      inputArgv[i]="";
      continue;
    }
  }
  if ( outputFile == "STDOUT" ) outputFile=cnvFile+".map";
  fastqFile=outputFile+".fa";
  bwaoutputFile=outputFile+".sam";
  logFile=outputFile+".log";
  if ( BWA=="" ) BWA=to_string(std::getenv("BWA"));
  if ( BWA=="" ) BWA="bwa";
  bwaOpt1+=" -M "; // flag secondary alignment
  // bwaOpt1="";
  for(i=2;i<inputArgv.size();++i) 
    if ( inputArgv[i]!="" ) bwaOpt1+=inputArgv[i]+" ";
  bwaCommand=BWA+" bwasw "
    + bwaOpt1 + " -f " + bwaoutputFile
    + " " + refFile + " " + fastqFile;
  
  
  if ( cnvFile=="" ) exit( usage_variation_to_fastq(argc, argv) );
  //if ( ! file_exist(bamFile) ) { cerr << bamFile << " not found\n"; exit(0); }
  if ( ! file_exist(cnvFile) ) { cerr << cnvFile << " not found\n"; exit(0); }
  
  string testbwa=BWA+" bwasw > /dev/null 2>&1";
  int bwaexit=system(testbwa.c_str());
  if ( bwaexit == 0 ) ISBWA=true;
  else cerr << "#BWA bwasw not executable " << bwaexit << endl;
  
  cerr << "#CommandL  : " << mycommand << "\n"
       << "#CNV file  : " << cnvFile << "\n"
       << "#BWA exec  : " << BWA << "\t" << std::boolalpha << ISBWA << "\n"
       << "#BWA ref   : " << refFile << "\n"
       << "#BWA opt   : " << bwaOpt1 << "\n"
       << "#BWA exec  : " << bwaCommand << "\n"
       << "#Output    : " << outputFile << ", "  
       << fastqFile << ", " << bwaoutputFile << "\n"
       << endl;

  
  ofstream FLOG(logFile.c_str());
  
  ofstream FOUT(fastqFile.c_str());
  if ( !FOUT ) { cerr << "Can't write to file " << outputFile << endl; exit(0); }
  
  ifstream FIN(cnvFile.c_str());
  if ( !FIN ) { cerr << "Can't open file " << cnvFile << endl; exit(0); }
  
  // find MERGE field
  istringstream iss;
  vector<string> cell(0);
  if ( iField<0 ) 
    while ( !FIN.eof() ) {
      string tmps,tmps1;
      getline(FIN,tmps);
      //if ( tmps[0] != '#' ) continue;
      if ( tmps.find("MERGE") == string::npos ) continue;
      iss.clear();
      iss.str(tmps);
      cell.clear();
      while( iss>>tmps1 ) cell.push_back(tmps1);
      for(i=0;i<cell.size();++i) if ( cell[i]=="MERGE" ) iField=i;
      //if (iField>0) break;
    }
  if ( iField<0 ) {
    cerr << "MERGE field not found, use -i INT to specify" << endl;
    exit(0);
  } 
  FIN.clear();
  FIN.seekg(0,FIN.beg);
  
  i=0;
  cnv_st icnv;
  vector<cnv_st> CNV(0);
  while ( !FIN.eof() ) {
    string tmps,tmps1;
    getline(FIN,tmps);
    if ( tmps[0] == '#' ) { HEADER=tmps; continue; }
    if ( tmps.length() <2 ) continue;
    
    iss.clear();
    iss.str(tmps);
    cell.clear();
    while( iss>>tmps1 ) cell.push_back(tmps1);
    
    if ( (int)cell.size()<iField+1 ) {
      cerr << "Input line has fewer fields than " << iField+1 << "\n"
	   << tmps 
	   << endl;
      continue;
    }
    
    icnv.T=cell[3][0];
    icnv.RNAME=cell[0];
    icnv.P1=atoi(cell[1].c_str());
    icnv.P2=atoi(cell[2].c_str());
    icnv.un=atoi(cell[4].c_str());
    icnv.isp5=0;
    icnv.isp3=0;
    icnv.ID=tmps;
    SEQ=cell[iField];
    if ( icnv.T=='A' || icnv.T=='J' ) {
      if ( icnv.P1<icnv.P2 ) swap(icnv.P1,icnv.P2);
    }
    else {
      if ( icnv.P1>icnv.P2 ) swap(icnv.P1,icnv.P2);
    }
    CNV.push_back(icnv);
    
    string COMMENT=to_string(CNV.size()-1)+"|"
      +to_string(icnv.T)
      +"_"+icnv.RNAME
      +":"+to_string(icnv.P1)
      +"-"+to_string(icnv.P2);
    
    FOUT << "@" << COMMENT << "\n"
	 << SEQ << "\n"
	 << "+" << "\n"
	 << string(SEQ.size(),'I') << endl;
    
    SEQ=cell[iField].substr(0, int(cell[iField].length()*0.7));
    FOUT << "@" << COMMENT << "/left\n"
	 << SEQ << "\n"
	 << "+" << "\n"
	 << string(SEQ.size(),'I') << endl;
    
    SEQ=cell[iField].substr(int(cell[iField].length()*0.3));
    FOUT << "@" << COMMENT << "/right\n"
	 << SEQ << "\n"
	 << "+" << "\n"
	 << string(SEQ.size(),'I') << endl;
    
  }
  FIN.close();
  FOUT.close();
  
  cerr << CNV.size() << " sequence written to " << fastqFile << endl;
  
  if ( samFile=="" ) {
    if ( ISBWA ) bwaexit=system(bwaCommand.c_str());
    cerr << "bwa bwasw exit code " << bwaexit << endl;
    if ( bwaexit != 0 ) return(0);
    cerr << bwaexit << " : alignment written to " << bwaoutputFile << endl;
  }
  else bwaoutputFile=samFile;

  vector<char> op(100);    // CIGAR operator
  vector<int> opnum(100);  // CIGAR operator bases
  vector<int> oppos(100);  // positions of operators in REF
  vector<int> opposseq(100);  // positions of operators in SEQ
  size_t anchor,iclip;
  
  FIN.open(bwaoutputFile.c_str());
  while ( !FIN.eof() ) {
    string tmps,tmps1;
    getline(FIN,tmps);
    if ( tmps[0] == '@' ) continue;
    if ( tmps.length() <2 ) continue;
    
    iss.clear();
    iss.str(tmps);
    iss >> QNAME
	>> FLAG
	>> RNAME 
	>> POS 
	>> MAPQ 
	>> CIGAR
	>> MRNM 
	>> MPOS 
	>> ISIZE
	>> SEQ;
    //>> QUAL 
    //>> OPT;
    
    int cnvid=atoi( QNAME.substr(0,QNAME.find("|")).c_str() );
    
    read_POS_CIGAR(POS,CIGAR,
		   anchor, iclip, 
		   op, opnum,
		   opposseq, oppos );
    // determine read MS orientation, MS, SM, or M
    string readType= (anchor<=iclip) ? "MS":"SM";
    int nM=0,nS=0,nID=0;
    for(int i=0;i<(int)op.size();++i) {
      if ( op[i]=='M' ) nM+=opnum[i];
      if ( op[i]=='S' ) nS+=opnum[i];
      if ( op[i]=='I' ) nID+=0;
      if ( op[i]=='D' || op[i]=='N' ) nID+=opnum[i];
    }
    if ( nS<SEQ.length()*0.08 )  readType="M";
    
    // check which side, left or right, the read comes from
    char SIDE='.';
    if ( QNAME.find("left")!=string::npos ) SIDE='L';
    if ( QNAME.find("right")!=string::npos ) SIDE='R';
    
    COMMENT="";
    
    // is it secondary read
    int SEC=0;
    if ( (FLAG & 0x100)>0 ) SEC=1;
    if ( (!SECONDARY) && SEC==1 ) { COMMENT=":SEC"; goto Note; }
    
    // read is reversed
    if ( (FLAG & 0x10)>0 ) { COMMENT+=":revComp"; goto Note; }
    
    if ( RNAME != CNV[cnvid].RNAME ) { COMMENT+=":notMatch"; goto Note; }
    
    // if read is from 5' side, should be MS or M
    if ( SIDE=='L' ) {
      if ( readType=="SM" ) { COMMENT=":wrongSide"; goto Note; }
      if ( POS<CNV[cnvid].P1 && CNV[cnvid].P1<POS+nM+10 ) {
	COMMENT+=":match5'";
	if ( CNV[cnvid].isp5!=9 && CNV[cnvid].isp3!=9 ) CNV[cnvid].isp5=1;
	goto Note;
      }
      COMMENT+=":notMatch";
      goto Note;
    }
    
    // if read is from 3' side, should be SM or M
    if ( SIDE=='R' ) {
      if ( readType=="MS" ) { COMMENT+=":wrongSide"; goto Note; }
      if ( abs(CNV[cnvid].P2-POS) < 11 ) {
	COMMENT+=":match3'";
	if ( CNV[cnvid].isp5!=9 && CNV[cnvid].isp3!=9 ) CNV[cnvid].isp3=1;
	goto Note;
      }
      COMMENT+=":notM";
      goto Note;
    }
    
    // if whole read and MS, should match 5'
    if ( readType=="MS" ) { 
      if ( POS<CNV[cnvid].P1 && CNV[cnvid].P1<POS+nM+10 ) {
	COMMENT+=":match5'";
	if ( CNV[cnvid].isp5!=9 && CNV[cnvid].isp3!=9 ) CNV[cnvid].isp5=1;
	goto Note;
      }
      COMMENT+=":notMatch";
      goto Note;
    }
    
    // if whole read and SM, should match 3'
    if ( readType=="SM" ) {
      if ( abs(CNV[cnvid].P2-POS) < 11 ) {
	COMMENT+=":match3'";
	if ( CNV[cnvid].isp5!=9 && CNV[cnvid].isp3!=9 ) CNV[cnvid].isp3=1;
	goto Note;
      }
      COMMENT+=":notMatch";
      goto Note;
    }
    
    // if whole read and M, should cover both
    if ( readType=="M" ) {
      COMMENT+=":singleLoc"; 
      CNV[cnvid].isp5=9;
      CNV[cnvid].isp3=9;
      if ( POS<CNV[cnvid].P1 && CNV[cnvid].P1<POS+nM ) {
	CNV[cnvid].isp5=1;
	COMMENT+=":match5'";
      }
      if ( POS<CNV[cnvid].P2 && CNV[cnvid].P2<POS+nM+nID ) {
	CNV[cnvid].isp3=1;
	COMMENT+=":match3'";
      }
      if ( CNV[cnvid].isp5==1 && CNV[cnvid].isp3==1 ) COMMENT+=":singleIndel"; 
      if ( CNV[cnvid].isp5==9 && CNV[cnvid].isp3==9 ) COMMENT+=":notMatch"; 
    } 
    
  Note:
    if ( VERBOSE ) {
      cerr << cnvid << "\t" << CNV[cnvid].T << "_" 
	   << CNV[cnvid].RNAME << ":" 
	   << CNV[cnvid].P1 << "-" << CNV[cnvid].P2 << "\t"
	   << SIDE << SEC << "|" << RNAME << ":" 
	   << POS << "\t" << CIGAR << "\t:"
	   << COMMENT << endl;
      
    }

    FLOG << cnvid << "\t" << CNV[cnvid].T << "_" 
	 << CNV[cnvid].RNAME << ":" 
	 << CNV[cnvid].P1 << "-" << CNV[cnvid].P2 << "\t"
	 << SIDE << SEC << "|" << RNAME << ":" 
	 << POS << "\t" << CIGAR << "\t"
	 << COMMENT << endl;
    
  }
  
  FOUT.open(outputFile.c_str());
  if ( HEADER!="" ) FOUT << HEADER << "\tMAP" << endl;
  for(int i=0;i<(int)CNV.size();++i) 
    FOUT << CNV[i].ID << "\t" << CNV[i].isp5 << CNV[i].isp3 << endl;
  FOUT.close();
  FLOG.close();
  
  return(0);
}

