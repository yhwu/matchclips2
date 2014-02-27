#include <iostream>
#include <string>
#include <vector>
#include <map>
using namespace std;

#include "functions.h"
#include "samfunctions.h"
#include "matchreads.h"
#include "filters.h"

void RD_PAIR_filters(vector<BREAKSEEK_st>& BP, string bamFile)
{
  
  bam_pileup_api_wrapper_t tmp;  
  tmp.in = NULL;
  tmp.bamidx = NULL;
  
  if ( is_binary(bamFile) ) {
    tmp.in = samopen(bamFile.c_str(), "rb", 0);
    tmp.bamidx = bam_index_load(bamFile.c_str()); // load BAM index
    if (tmp.bamidx == 0) {  
      cerr << "BAM indexing file for " << bamFile << " is not available." << endl ;  
      return ;  
    }  
  }
  else return;
  
  // cerr << "Filtering " << BP.size() << " cnvs" << endl;
  
  size_t i, k;
  string bamRegion="";
  
  for(i=0;i<BP.size();++i) {
    size_t dp1=0,dp2=0,din=0,dout=0;
    int dx=max(BP[i].L1,BP[i].L2)*2;
    if ( dx<(int)abs(BP[i].P2-BP[i].P1) ) dx=(int)abs(BP[i].P2-BP[i].P1);
    if ( dx>400) dx=400;
    int p1=min(BP[i].P1,BP[i].P2) -dx;
    int p2=max(BP[i].P1,BP[i].P2) +dx;
    if ( p1<1 ) p1=1;
    int bp1=BP[i].P1+BP[i].UN;
    int bp2=BP[i].P2;
    if ( BP[i].UN > BP[i].len/2 ) {
      bp1=BP[i].P1;
      bp2=BP[i].P2;
    }
    int Lout=bp1-p1-1 + p2-bp2-1;
    int Lin=bp2-bp1-1;
    if ( Lin<1 ) Lin=1;
    
    int klen=kmerlength(BP[i].MERGE);
    
    //if ( BP[i].count>10 && klen<40 ) goto FILTERDONE;
    //if ( BP[i].count==1 && BP[i].CL<30 && BP[i].len>100000 ) BP[i].isit=false;
    //if ( BP[i].isit==false ) goto FILTERDONE;
    
    tmp.T=BP[i].T;
    tmp.bp1=BP[i].P1;
    tmp.bp2=BP[i].P2;
    if ( tmp.bp1 > tmp.bp2 ) swap(tmp.bp1,tmp.bp2);
    //    bampileup(tmp, BP[i].RNAME(), p1, p2);
    bampileup(tmp, msc::bam_target_name[BP[i].tid], p1, p2);
    
    for(k=p1,dout=0; (int)k<bp1; ++k) dout+=tmp.n[k-p1];
    for(k=bp2+1; (int)k<=p2; ++k) dout+=tmp.n[k-p1];
    for(k=bp1+1, din=0; (int)k<BP[i].P2; ++k) din+=tmp.n[k-p1];
    dp1=tmp.n[bp1-p1];
    dp2=tmp.n[bp2-p1];
    
    dout /= Lout;
    din  /=  Lin;
    
    float q0=(float)tmp.c0/(float)(tmp.c0+tmp.c1);
    
    //if ( BP[i].count==1 && BP[i].CL<30 ) BP[i].isit=false;
    //if ( BP[i].count==2 && BP[i].CL<30 ) BP[i].isit=false;
    
    if ( BP[i].count==1 || klen>40 ) {
      if ( BP[i].T=='D' && din>=dout*3/4 ) BP[i].isit=false;
      if ( BP[i].T=='A' && din<=dout*5/4 ) BP[i].isit=false;
    }
    
    if ( BP[i].count==2 ) {
      if ( BP[i].T=='D' && din>=dout ) BP[i].isit=false;
      if ( BP[i].T=='A' && din<=dout ) BP[i].isit=false;
    }
    
    if (dout>30 || din>30)  {
      if ( BP[i].T=='D' && din>=dout ) BP[i].isit=false;
      if ( BP[i].T=='A' && din<=dout ) BP[i].isit=false;
    }
    if (dout>50 || din>50)  {
      if ( BP[i].T=='D' && din>=dout*3/4 ) BP[i].isit=false;
      if ( BP[i].T=='A' && din<=dout*5/4 ) BP[i].isit=false;
    }

    if ( klen > 33 && BP[i].count<5 ) BP[i].isit=false;
    
    //FILTERDONE:    
    BP[i].TAG="DO="+to_string(dout)+","
      +"DI="+to_string(din)+","
      +"D1="+to_string(dp1)+","
      +"D2="+to_string(dp2)+";"
      +"PR="+to_string(tmp.pair/2)+";"
      +"KM="+to_string(klen)+";"
      +"Q0="+to_string((int)(q0*100+0.499))+"%";
    
    //    cerr << BP[i].RNAME() << "\t" 
    cerr << msc::bam_target_name[ BP[i].tid ] << "\t" 
	 << BP[i].P1 << "\t"
	 << BP[i].P2 << "\t"
 	 << BP[i].T << "\t" 
 	 << BP[i].UN << "\t" 
 	 << BP[i].len << "\t" 
	 << BP[i].TAG << "\t"
	 << BP[i].count << "\t"
	 << klen << "\t" 
	 << BP[i].isit << endl; 
    
  }
  
  samclose(tmp.in);
  bam_index_destroy(tmp.bamidx);  
  
  return;
}

