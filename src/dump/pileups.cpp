#include <pthread.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <complex>
#include <algorithm>
#include <vector>
using namespace std;

/**** samtools headers ****/
#include <bam.h>
#include <sam.h>

/**** user headers ****/
#include "functions.h"
#include "samfunctions.h"
#include "matchreads.h"

#include "wu2.h"
#include "pileups.h"

//! check if readMS and readSM overlap with miinimum minOver charactors
//! maxErr mismatches are tolorated
//! readMS is always fully overlaped with or in front of readSM
//! p1 is the position on readMS where overlap begins
//! p1=-1 if not overlaped
//! p_err stores indices mismatched charactors in readMS
bool string_align(const string& FASTA, const string& read, 
		  int beg, int end,
		  const int maxErr,
		  //		  int& p1, vector<int>& p_err) 
		  int& p1, int& p1_err) 
{
  p1=-1;
  vector<int> p_err(0);
  p_err.clear();
  
  bool match=false;
  bool founddel=false;
  bool foundadd=false;
  vector<int> p_match(0);
  vector<int> p_match_err(0);
  
  int i,j,l,idel,iadd,ndiff;
  
  for(i=beg; i<end; ++i) {
    match=false;
    p_err.clear();
    
    // check base match
    ndiff=0;
    for(j=0; j<(int)read.length(); ++j) {
      if ( FASTA[i+j]==read[j] ) continue;
      ++ndiff;
      p_err.push_back(j);
      if (ndiff>maxErr) break;
    }
    if (ndiff<=maxErr) { 
      match=true;
      p1=i;
      p_match.push_back(i); 
      p_match_err.push_back(ndiff); 
      // break;
    } 
    // if ( match ) {cerr << "match\t" << p_match.back() << "\t" << p_match_err.back() << endl;}
    if ( match ) continue;
    // continue;
    
    // now base match failed; check base deletion
    bool front_mismatch=false;
    for(l=0; l<(int)p_err.size(); ++l) if ( p_err[l]<3 ) front_mismatch=true;
    if ( front_mismatch ) continue;

    founddel=false;
    for(l=0; l<(int)p_err.size(); ++l) {
      //if ( p_err[i]==0 ) continue;
      for(idel=1; idel<=maxErr; ++idel ) {
	ndiff=l+idel;
	for(j=p_err[l]; j<(int)read.length(); ++j) {
	  if ( FASTA[i+j+idel]==read[j] ) continue;
	  ndiff += 1;
	  if (ndiff>maxErr) break;
	}
	if (ndiff<maxErr) { 
	  p_match.push_back(i); 
	  p_match_err.push_back(ndiff); 
	  founddel=true;
	  break;
	} 
      }
      if ( founddel ) break;
    }
    if ( founddel ) {cerr << "deletion\t" << p_match.back() << "\t" << p_match_err.back() << endl;}
    if ( founddel ) continue;
    
    // now base match failed; check base insertion
    foundadd=false;
    for(l=0; l<(int)p_err.size(); ++l) {
      for(iadd=1; iadd<=maxErr; ++iadd ) {
	ndiff=l+1;
	for(j=p_err[l]; j<(int)read.length(); ++j) {
	  if ( FASTA[i+j]==read[j+iadd] ) continue;
	  ndiff += 1;
	  if (ndiff>maxErr) break;
	}
	if (ndiff<maxErr) { 
	  p_match.push_back(i); 
	  p_match_err.push_back(ndiff); 
	  foundadd=true;
	  break;
	} 
      }
      if ( founddel ) break;
    }
    if ( foundadd ) {cerr << "insertion\t" << p_match.back() << "\t" << p_match_err.back() << endl;}
    
    if (ndiff<=maxErr) { 
      match=true;
      p1=i;
      break;
    }
    
  }
  
  if ( p_match.size()>0 ) {
    match=true;
    int imin=0;
    for(i=0;i<(int)p_match.size();++i) 
      if ( p_match_err[i]<p_match_err[imin] ) imin=i;
    p1=p_match[imin];
    p1_err=p_match_err[imin];
  }
  
  return(match);
}

void get_continuous_soft_clipped_seq(const string& FASTA, vector<pairinfo_st>& pair)
{
  vector<unsigned short int> rd( FASTA.size(), 0 );
  vector<char> c( FASTA.size(), 0 );
  vector<unsigned short int> rd_ma( FASTA.size(), 0 ); // moving_average
  
  ifstream fin("sc_acgt1.txt"); 
  
  double s2=0.0, s=0.0, n=0.0;
  double s_mean=0, s_sd=0;
  while ( !fin.eof() ) {
    string tmps;
    getline(fin,tmps);
    if ( tmps.length()<5 ) continue;
    for(int i=tmps.size()-1; i>=0; --i ) 
      if ( tmps[i]==' ' ) { tmps[i]='\t'; break; }
    istringstream iss(tmps);
    int pos, c_rd=0;
    char ref_char, s_char;
    string tmps1;
    if ( ! ( iss >> tmps1 >> pos >> ref_char >> s_char >> c_rd ) ) {
      cerr << "input error\n" << tmps << "\n"
	   << pos << "\t" << ref_char << "\t" << s_char << "\t" << c_rd 
	   << endl;
    }
    rd[pos]=c_rd;
    c[pos]=s_char;
    
    if ( c_rd>0 && ref_char != 'N' ) {
      s+=c_rd;
      s2+=c_rd*c_rd;
      n+=1;
    }
  }
  fin.close();
  if ( n>10 ) {
    s_mean=s/n;
    s_sd=sqrt( s2/n - s/n*s/n );
    cerr << "mean and std\t" << s/n << "\t" << sqrt( s2/n - s/n*s/n ) << endl;
  }
  else { return; }
  
  /*
    double ms=0;
    for(int i=0; i<10; ++i) ms+=rd[i];
    for(int i=5; i<(int)rd.size()-6; ++i) {
    ms-=rd[i-5];
    ms+=rd[i+5];
    rd_ma[i]=ms/10;
    }
    runmed(&rd[0], &rd_ma[0], rd.size(), 11);
  */
  
  
  rd_ma=rd;
  s2=0.0; s=0.0; n=0.0;
  for(int i=0; i<(int)rd_ma.size(); ++i) {
    if ( rd_ma[i]>0 ) {
      s+=rd_ma[i];
      s2+=rd_ma[i]*rd_ma[i];
      n+=1;
    }
  }  
  if ( n>10 ) {
    s_mean=s/n;
    s_sd=sqrt( s2/n - s/n*s/n );
    cerr << "mean and std\t" << s/n << "\t" << sqrt( s2/n - s/n*s/n ) << endl;
  }
  
  vector<int> p1(0),p2(0);
  int S_beg=-1, S_end=-1;
  for(int i=0; i<(int)rd_ma.size(); ++i) {
    if ( rd_ma[i] > s_mean+s_sd*2 ) {
      if ( S_beg<0 ) { S_beg=i; S_end=i; }
      else { S_end=i; }
    }
    else {
      if (S_beg>0 && S_end>S_beg+11 ) {
	p1.push_back(S_beg);
	p2.push_back(S_end);
      }
      S_beg=-1;
      S_end=-1;
    }
  }
  
  int p_pos=-1;
  //  vector<int> p_err(0);
  int p_err=-1;
  
  for(int i=0; i<(int)p1.size(); ++i) {
    // for(int i=0; i<20; ++i) {
    string comS=string( &c[ p1[i] ], p2[i]-p1[i]+1 );
    int ED=edit_distance(comS, FASTA.substr(p1[i], p2[i]-p1[i]+1 ) );
    // if ( ED*12 < comS.length() ) continue;
    
    cerr << i << "\t" << p1[i] << "\t" << rd_ma[ p1[i] ] 
	 << "\t" << p2[i] << "\t" << rd_ma[ p2[i] ] << "\t"
	 << p2[i]-p1[i]
	 << endl;
    cerr << comS << "\n"
	 << FASTA.substr(p1[i], p2[i]-p1[i]+1 )  << "\t"
	 << p2[i]-p1[i] << "\t" << edit_distance(comS, FASTA.substr(p1[i], p2[i]-p1[i]+1 ) )
	 << endl;
    
    string_align(FASTA, comS, p1[i]-1000, p1[i]+1000, 3, p_pos, p_err);
    if ( p_pos>0 && p_pos!=p1[i] ) {
      //if ( ( p_pos>0 && p_pos!=p1[i] )  || 1 ) {
      cerr << "msoft\t" 
	   << p1[i] << "\t" << p2[i]-p1[i]+1 << "\t" << ED << "\t"
	   << p_pos << "\t" << p_err << "\t"
	   << p_pos - p1[i] << endl;
      if ( ( p_pos-p1[i] )>100 ) {
	pairinfo_st ipair;
	ipair.F2=p1[i];
	ipair.R1=p_pos;
	pair.push_back(ipair);
	ipair.F2=p_pos;
	ipair.R1=p1[i];
	pair.push_back(ipair);
      }
    }
    cerr << "-------------\n";
  }
  
}

void S_pileup(const int ref, const int beg, const int end, const string& FASTA, 
	      vector<pairinfo_st>& pair)
{  
  pair.clear();
  get_continuous_soft_clipped_seq(FASTA, pair);
   return;

  POSCIGAR_st bm;
  bam1_t *b=NULL; b = bam_init1();
  bam_iter_t iter=NULL;
  
  const string ACGTN="ACGTN";
  
  //  Mat<unsigned short int > PU(msc::fp_in->header->target_len[ref], 5, (unsigned short int)0 );  
  Mat<unsigned short int > PU(FASTA.size(), 5, (unsigned short int)0 );  
  
  iter = bam_iter_query(msc::bamidx, ref, beg, end);
  int iMillion=-1;
  while( bam_iter_read(msc::fp_in->x.bam, iter, b)>0 ) {
    
    string SEQ=get_qseq(b);
    resolve_cigar_pos(b, bm, 0);  
    if ( (b->core.pos - beg)/1000000>iMillion ) {
      iMillion=(b->core.pos - beg)/1000000;
      cerr << msc::fp_in->header->target_name[ref] << "@" << commify(b->core.pos) << '\xD';
    }
    
    if ( b->core.tid != ref ) break;
    if ( bm.pos<=0 ) continue;;
    if ( bm.cop[0]+bm.nop[0] < bm.cop[0] ) continue;
    if ( bm.cop.back()+bm.nop.back() >= FASTA.size() ) continue;
    if ( bm.cop[0]<1 ) continue;
    if ( bm.iclip<0 ) continue;
    
    if ( bm.nop[bm.iclip]<=3 ) continue;
    
    calibrate_resolved_cigar_pos((string&)FASTA, SEQ, bm);  
    vector<int> e_pos(0);
    expand_pos(bm, e_pos );
    
    for(int i=0; i<(int)bm.nop[bm.iclip]; ++i) {
      int q=bm.qop[bm.iclip]+i;
      size_t pos=ACGTN.find( SEQ[q] );
      if ( pos!=string::npos && e_pos[q]<(int)FASTA.size() ) {
	PU[e_pos[q]][pos]+=1;
	if ( PU[e_pos[q]][pos] > 65000 ) PU[e_pos[q]][pos]=65000;
      }
    }
  }
  cerr << endl;
  
  ofstream fout("sc.txt");
  
  double s2=0.0, s=0.0, n=0.0;
  for(int i=beg; i<end; ++i) {
    if ( i<0 ) continue;
    if ( i>=(int)FASTA.size() ) break;
    
    int rd=0;
    int imax=4;
    for(int j=0;j<5; ++j) {
      //cerr << i << "\t" << j << endl;
      if ( PU[i][j] > PU[i][imax] ) imax=j;
      rd+=PU[i][j];
    }
    int inext=4;
    for(int j=0;j<5; ++j) {
      if ( j==imax ) continue;
      if ( PU[i][j] > PU[i][inext] ) inext=j;
    }
    
    char var_char='N';
    
    if ( ACGTN[ imax ] != FASTA[i] ) var_char=ACGTN[ imax ];
    else {
      if ( PU[i][inext]*3 >= PU[i][imax] ) var_char=ACGTN[ inext ];
    }
    
    if ( msc::verbose>1 || 1 ) 
      fout << "SC\t" << i << "\t" << FASTA[i] << "\t"
	//cerr << "SC\t" << i << "\t" << FASTA[i] << "\t"
	   << var_char << " " << rd << "\t"
	   << PU[i][0] << " "
	   << PU[i][1] << " "
	   << PU[i][2] << " "
	   << PU[i][3] << " "
	   << PU[i][4] 
	   << endl;
    
    if ( PU[i][imax]>0 || FASTA[i]!='N' ) {
      s+=PU[i][imax];
      s2+=PU[i][imax]*PU[i][imax];
      n+=1;
    }
  }
  
  if ( n>0 ) {
    cerr << "mean and std\t" << s/n << "\t" << sqrt( s2/n - s/n*s/n ) << endl;
  }

  if ( b!=NULL ) bam_destroy1(b);
  if ( iter != NULL ) bam_iter_destroy(iter);

  return;
}

