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

   • H can only be present as the first and/or last operation.
   • S may only have H operations between them and the ends of the CIGAR string.
   • For mRNA-to-genome alignment, an N operation represents an intron. For other types of alignments, the interpretation of N is not defined.
   • Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ.
from sam1.pdf
*********************************************************************/
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
using namespace std;

void read_POS_CIGAR(int POS, string& CIGAR,
		    size_t& anchor, size_t& iclip, 
		    vector<char>& op, vector<int>& opnum,
		    vector<int>& opposseq, vector<int>& oppos )
{
  string ops="";
  char numchar[]="999999999999999999999999999999";
  size_t i,k;
  int baseidx;
  
  anchor=0;
  iclip=string::npos;
  op.resize(0);
  opnum.resize(0);

  for(i=0, k=0;i<CIGAR.size();++i) {
    if ( CIGAR[i]>='0' && CIGAR[i]<='9'  ) {
      numchar[k]=CIGAR[i];
      ++k;
    }
    else {
      if ( CIGAR[i]=='*' ) { 
	op.resize(1);
	op[0]='*';
	opnum.resize(1);
	opnum[0]=0;
	opposseq.resize(1);
	opposseq[0]=0;
	oppos.resize(1);
	oppos[0]=POS;
	return;
	//break;
      }
      if ( k>9 ) { cerr << "CIGAR error" << CIGAR << endl; exit(0); }
      ops+=CIGAR[i];
      numchar[k]=0;
      if ( CIGAR[i]!='H' ) {              // just ignore H hard clip
	opnum.push_back(atoi(numchar));
	op.push_back(CIGAR[i]);
      }
      k=0;
    }
  }
  
  opposseq.resize(op.size());
  oppos.resize(op.size());
  
  anchor=string::npos;
  iclip=0;  
  size_t maxclipped=0;
  baseidx=0;
  for(i=0;i<op.size();++i) {
    opposseq[i]=baseidx;
    if ( op[i]!='D' && op[i]!='N' && op[i]!='P' && op[i]!='H') 
      baseidx+=opnum[i]; 
    if ( (op[i]=='S' || op[i]=='H') && opnum[i]>(int)maxclipped ) {
      maxclipped=opnum[i];
      iclip=i;
    }
    if ( op[i]=='M' && anchor==string::npos) anchor=i;
  }
  if ( op.size()>1 && op[0]=='H' && op.back()=='H' ) {
    cerr << "H can only be present as the first and/or last operation\t"
	 << CIGAR << endl;
    op.resize(1);
    op[0]='*';
    opnum.resize(1);
    opnum[0]=0;
    opposseq.resize(1);
    opposseq[0]=0;
    oppos.resize(1);
    oppos[0]=POS;
    return;
  }
  if ( iclip>0 && iclip<op.size()-1 ) {
    cerr << "[read_CIGAR.cpp] clipped part not at ends \n"
	 << CIGAR << "\t" << iclip << endl;
    op.resize(1);
    op[0]='*';
    opnum.resize(1);
    opnum[0]=0;
    opposseq.resize(1);
    opposseq[0]=0;
    oppos.resize(1);
    oppos[0]=POS;
    return;
  }
  if ( anchor==string::npos && CIGAR!="*" ) { 
    //cerr << "[read_CIGAR.cpp] Error CIGAR does not has M " << CIGAR << endl; 
    op.resize(1);
    op[0]='*';
    opnum.resize(1);
    opnum[0]=0;
    opposseq.resize(1);
    opposseq[0]=0;
    oppos.resize(1);
    oppos[0]=POS;
    return;
  }
  if ( op[iclip]!='S' && op[iclip]!='H' ) iclip=string::npos;
  
  /* find the position of the operators in the REF 
     I,P don't change positions in REF 
     find the first match then calculate positions before and after */
  oppos[anchor]=POS;
  for(i=anchor-1;i>=0;--i) {
    if ( anchor==0 ) break;      // be careful with unsigned and --
    if ( op[i]=='I' || op[i]=='P' ) oppos[i]=oppos[i+1];
    else { oppos[i]=oppos[i+1]-opnum[i]; }
    if ( i==0 ) break;
  }
  for(i=anchor+1;i<op.size();++i) 
    if (op[i-1]=='I' || op[i-1]=='P' ) oppos[i]=oppos[i-1];
    else { oppos[i]=oppos[i-1]+opnum[i-1]; }
  
  
  return;
}

/* note H clipped bases not present in read SEQ 
   no actual bases are produced in SEQ by H,D,P,N operators 
*/

//this subroutine only adjust first S and/or last S parts
int calibrate_cigar(string& SEQ, int& POS, string& FASTA,
		    size_t& anchor, size_t& iclip, 
		    vector<char>& op, vector<int>& opnum,
		    vector<int>& opposseq, vector<int>& oppos )
{
  anchor=0;
  iclip=string::npos;
  int nChanged=0;
  
  if ( op[0]=='H' ) {
    op.erase(op.begin());
    opnum.erase(opnum.begin());
    opposseq.erase(opposseq.begin());
    oppos.erase(oppos.begin());
  }
  if ( op.back()=='H' ) {
    op.pop_back();
    opnum.pop_back();
    opposseq.pop_back();
    oppos.pop_back();
  }
  if ( op.size()<=1 ) return nChanged;
  
  // first, adjust S at the 3' end, POS don't change, only S number
  if ( op.back()=='S' ) {
    string CLIPPEDSEQ=SEQ.substr(opposseq.back(),opnum.back());
    string REFCLIP=FASTA.substr(oppos.back()-1,opnum.back());
    
    int iclip=op.size()-1;
    int p1misMatch=REFCLIP.size();
    for(int i=0;i<(int)REFCLIP.size();++i) {
      if ( REFCLIP[i]==CLIPPEDSEQ[i] ) continue;
      p1misMatch=i;
      break;
    }
    int dx=p1misMatch;
    nChanged+=dx;
    
    opnum[iclip-1]+=dx;
    opnum[iclip]-=dx;
    oppos[iclip]+=dx;
    opposseq[iclip]+=dx;
  }
  
  // second, adjust S at the 5' end, POS change too 
  if ( op[0]=='S' ) {
    string CLIPPEDSEQ=SEQ.substr(opposseq[0],opnum[0]);
    string REFCLIP=FASTA.substr(oppos[0]-1,opnum[0]);
    
    int p2misMatch=0;
    for(int i=REFCLIP.size()-1;i>=0;--i) {
      if ( REFCLIP[i]==CLIPPEDSEQ[i] ) continue;
      p2misMatch=i+1;
      break;
    }
    int dx=REFCLIP.size()-p2misMatch;
    nChanged+=dx;
    
    opnum[0]-=dx;
    opnum[1]+=dx;
    oppos[1]-=dx;
    opposseq[1]-=dx;
    if ( op[1]=='M' ) POS=oppos[1];
  }
  
  iclip=string::npos;
  if ( op[0]=='S' ) iclip=0;
  if ( op.back()=='S' ) iclip=op.size()-1; 
  if ( iclip==string::npos ) {
    cerr << "error reading cigar \t";
    for(int i=0;i<op.size();++i) cerr << opnum[i] << op[i]; 
    cerr << endl;
  }
  if ( op.back()=='S' && opnum.back()>opnum[iclip] ) iclip=op.size()-1;
  if ( op[0]=='S' && opnum[0]>opnum[iclip] ) iclip=0;
  anchor=string::npos;
  for(int i=0;i<op.size();++i) if ( op[i]=='M' ) { anchor=i; break; }
  if ( anchor==string::npos ) {
    cerr << "error reading cigar \t";
    for(int i=0;i<op.size();++i) cerr << opnum[i] << op[i]; 
    cerr << endl;
  }
  
  
  return nChanged;
}

