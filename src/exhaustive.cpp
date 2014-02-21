#include <pthread.h>
#include <stdio.h>
#include <string.h>
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

/**** samtools headers ****/
#include <bam.h>
#include <sam.h>

/**** user headers ****/
#include "samfunctions.h"
#include "functions.h"
#include "matchreads.h"

#include "preprocess.h"
#include "pairguide.h"
#include "exhaustive.h"

extern pthread_mutex_t nout;

static bool sort_bp(const ED_st& p1, const ED_st& p2)
{
  if ( p1.F2 != p2.F2 ) return p1.F2<p2.F2;
  return p1.R1<p2.R1;
}
static bool sort_count(const ED_st& p1, const ED_st& p2)
{
  if ( p1.count != p2.count ) return p1.count>p2.count;
  return p1.ED<p2.ED;
}

void compact_pairgroup(vector<ED_st>& pairgroup, vector<ED_st>& bp0)
{
  bp0.clear();
  if ( pairgroup.size()<1 ) return;
  
  // merge same break points
  vector<ED_st> cluster(0);
  for(int i=0; i<(int)pairgroup.size(); ++i) {
    // cerr << pairgroup[i].F2 << "\t" << pairgroup[i].R1 << "\t" << pairgroup[i].ED << endl;  
    if ( pairgroup[i].F2+1==pairgroup[i].R1 ) continue;
    bool is_in_cluster=false;
    for(size_t k=0; k<cluster.size(); ++k) {
      if ( pairgroup[i].F2==cluster[k].F2 && 
	   pairgroup[i].R1==cluster[k].R1 ) {
	is_in_cluster=true;
	cluster[k].count+=pairgroup[i].count;
	cluster[k].ED+=pairgroup[i].ED;
	break;
      }
    }
    if ( ! is_in_cluster ) cluster.push_back(pairgroup[i]);
  }
  for(int i=0; i<(int)cluster.size(); ++i) {
    cluster[i].ED /= cluster[i].count;
    //cerr << "cluster " << cluster[i].F2 << "\t" << cluster[i].R1 << "\t" << cluster[i].ED << "\t" << cluster[i].count << endl;
  }
  //cerr << "-------" << endl;
  
  // merge simutaneous displacements
  bool merged=false;
  for(int i=0; i<(int)cluster.size()-1; ++i) {
    if ( cluster[i].count==0 ) continue;
    
    vector<int> idx(0);
    for(int k=i+1; k<(int)cluster.size(); ++k) {
      if ( cluster[k].count==0 ) continue;
      if ( cluster[i].F2-cluster[k].F2 == cluster[i].R1-cluster[k].R1 &&
	   abs(cluster[i].F2-cluster[k].F2) < abs(cluster[i].F2-cluster[i].R1) && 
	   abs(cluster[i].F2-cluster[k].F2) < 11 ) 
	idx.push_back(k);
    }
    if ( idx.size()==0 ) continue;
    
    idx.push_back(i);
    int isel=i;
    double new_count=0;
    double new_ED=0;
    for(int k=0; k<(int)idx.size(); ++k) {
      new_count+=cluster[ idx[k] ].count;
      new_ED+=(double)cluster[ idx[k] ].ED * (double)cluster[ idx[k] ].count;
      if ( cluster[isel].count < cluster[ idx[k] ].count ) isel=idx[k];
    }
    new_ED = new_ED/new_count;
    
    for(int k=0; k<(int)idx.size(); ++k) cluster[idx[k]].count=0;
    cluster[isel].count=(int)new_count;
    cluster[isel].ED=(int)(new_ED+0.5);
    merged=true;
  }
  
  // merge close break points
  // merged=false;
  for(int i=0; i<(int)cluster.size()-1; ++i) {
    if ( cluster[i].count==0 ) continue;
    
    vector<int> idx(0);
    for(int k=i+1; k<(int)cluster.size(); ++k) {
      if ( cluster[k].count==0 ) continue;
      if ( abs(cluster[i].F2-cluster[k].F2)<5 &&
	   abs(cluster[i].R1-cluster[k].R1)<5 ) 
	idx.push_back(k);
    }
    if ( idx.size()==0 ) continue;
    
    idx.push_back(i);
    int isel=i;
    double new_count=0;
    double new_ED=0;
    for(int k=0; k<(int)idx.size(); ++k) {
      new_count+=cluster[ idx[k] ].count;
      new_ED+=(double)cluster[ idx[k] ].ED * (double)cluster[ idx[k] ].count;
      if ( cluster[isel].count < cluster[ idx[k] ].count ) isel=idx[k];
    }
    new_ED = new_ED/new_count;
    
    for(int k=0; k<(int)idx.size(); ++k) cluster[idx[k]].count=0;
    cluster[isel].count=(int)new_count;
    cluster[isel].ED=(int)(new_ED+0.5);
    merged=true;
  }
  
  if ( msc::verbose>1 ) {
    for(int i=0; i<(int)cluster.size(); ++i) {
      if ( cluster[i].count==0 ) continue;
      cerr << "cluster " << i << "\t" << cluster[i].F2 << "\t" << cluster[i].R1 << "\t" << cluster[i].ED << "\t" << cluster[i].count << endl;
    }
    cerr << "-------" << endl;
  }
  
  sort(cluster.begin(), cluster.end(), sort_count);
  for(int i=0; i<(int)cluster.size(); ++i) {
    if ( i>0 ) { if ( cluster[i-1].count>cluster[i].count*2 ) break; }
    bp0.push_back( cluster[i] );
    // cerr << "cluster " << cluster[i].F2 << "\t" << cluster[i].R1 << "\t" << cluster[i].ED << "\t" << cluster[i].count << endl;
    if ( bp0.size()>1 ) break;
  }
  // cerr << "Total : " << bp0.size() << endl;
  // cerr << endl;
  
  return;
}

void reduce_matched_break_points(vector<ED_st>& bp, vector<ED_st>& reduced)
{
  reduced.clear();
  if (bp.size()<1 ) return;
  
  int pair_gap=max(5, msc::bam_l_qseq/5);
  
  vector<ED_st> pairgroup(0), bp0(0);
  
  sort(bp.begin(), bp.end(), sort_bp);
  
  size_t im;
  for(im=0; im<bp.size(); ++im) if ( bp[im].F2>0 ) break;
  if ( im>=bp.size() ) return;
  
  pairgroup.push_back( bp[im] );
  for(size_t i=im+1; i<bp.size(); ++i) {
    if ( bp[i].F2<0 || bp[i].R1<0 ) continue;
    if ( bp[i].F2 - bp[i-1].F2 > pair_gap ) {
      compact_pairgroup(pairgroup, bp0);
      pairgroup.clear();
      if ( bp0.size()>0 ) 
	reduced.insert(reduced.end(), bp0.begin(), bp0.end());
    }
    pairgroup.push_back( bp[i] );
  }
  compact_pairgroup(pairgroup, bp0);
  pairgroup.clear();
  if ( bp0.size()>0 ) reduced.insert(reduced.end(), bp0.begin(), bp0.end());
  
  return;
}

void reduce_matched_break_points(vector<ED_st>& bp, vector<pairinfo_st>& mcbp)
{
  
  mcbp.clear();
  if (bp.size()<1 ) return;
  
  vector<ED_st> reduced(0);
  reduce_matched_break_points(bp, reduced);
  
  pairinfo_st ipair;
  for(int i=0; i<(int)reduced.size();++i) {
    ipair.tid=msc::bam_ref;
    ipair.F2=reduced[i].F2;
    ipair.MS_F2=reduced[i].F2;
    ipair.R1=reduced[i].R1;
    ipair.MS_R1=reduced[i].R1;
    ipair.sr_ed=reduced[i].ED;
    ipair.sr_count=reduced[i].count;
    mcbp.push_back(ipair);
    
    if ( msc::verbose>1 ) {
      cerr << "MC\t" << i << "\t" 
	   << reduced[i].F2 << "\t" 
	   << reduced[i].R1 << "\t" 
	   << reduced[i].ED << "\t" 
	   << reduced[i].count
	   << endl;
    }
  }
  
  return;
}

void compact_reads(vector<bam1_t>& bset)
{
  /*
    for(size_t i=1; i<bset.size(); ++i) {
    bool saveaspre=( get_qseq(&bset[i]) == get_qseq(&bset[i-1]) );
    cerr << b_MS[i].core.pos << "\t" 
    << get_cigar(&b_MS[i]) << "\t" 
    << get_qseq(&b_MS[i]) << "\t"
	 << saveaspre
	 << endl;
  }
  */
  if ( bset.size()<100 ) return;
  
  vector<size_t> posidx(0);
  vector<int> poss(0);
  poss.push_back(bset[0].core.pos);
  for(size_t i=0; i<bset.size(); ++i) {
    int pos=bset[i].core.pos;
    if ( pos==poss.back() ) continue;
    poss.push_back(pos);
  }
  
  size_t jstart=0;
  for(size_t i=0; i<poss.size(); ++i) {
    posidx.clear();
    map <string, size_t> qseqmap; qseqmap.clear();
    for(size_t j=jstart; j<bset.size(); ++j) {
      if ( bset[j].core.pos>poss[i] ) break;
      jstart=j;
      if ( bset[j].core.pos<poss[i] ) continue;
      posidx.push_back(j);
      string qseq=get_qseq( &bset[j] );
      if ( qseqmap.find( qseq ) == qseqmap.end() ) qseqmap[qseq]=1; 
      else qseqmap[qseq]+=1;
    }
    
    cerr << poss[i] << "\t" << posidx.size() << "\t" << qseqmap.size() << endl;
    //for(map<string,size_t>::iterator it = qseqmap.begin(); it!=qseqmap.end(); ++it) 
    // cerr << "\t" << it->first << "\t" << it->second << endl;
    
  }
  

}

void sampleidx(size_t N, size_t NS, vector<size_t>& id)
{
  
  id.clear();
  if( N==0 || NS==0 ) return;
  
  if ( NS>=N ) {
    id.resize(N);
    for(size_t i=0; i<N; ++i) id[i]=i;
    return;
  }
  
  id.resize(NS);
  double inc=((double)N-1.0)/((double)NS-1.0);
  double eps=0.01/((double)NS-1.0);
  for(size_t i=0; i<NS; ++i) {
    double di=(double)i*inc+eps;
    id[i]=(size_t) di;
    // cerr << id[i] << " ";
  }
  // cerr << N << endl;
  return;
}

// check_length > 0, check matching within the range 
// check_length = 0, check matching within the default range 
// check_length < 0, check matching among all reads 
void match_reads_for_exhaustive_search(int thread_id,
				       int NUM_THREADS,
				       vector<bam1_t>& b_MS,
				       vector<bam1_t>& b_SM,
				       string& FASTA, 
				       int check_length,
				       vector<ED_st>& ibp)
//				       vector<ED_st>& bp)
{
  ibp.clear();
  // bp.clear();
  if ( b_MS.size()<1 || b_SM.size()<1 || FASTA.size()<2 ) return;
  if ( check_length==0 ) return;
  if ( NUM_THREADS<1 ) {
    cerr << "thread error NUM_THREADS=" << NUM_THREADS << endl;
    exit(0);
  }
  
  if ( thread_id==0 ) {
    cerr << "matching " 
	 << b_MS.size() << " X " << b_SM.size()
	 << " reads within range " << check_length << endl;
  }

  //compact_reads(b_MS); return; 
  // reduce repeated reads
  // not useful, only reduced a few
  
  string readMS,readSM;
  
  int minOver=msc::minOverlap;
  int maxErr=msc::errMatch;
  size_t m_count=0;

  size_t bpreserve=1000000;
  vector<ED_st> bp(0); bp.reserve(bpreserve);
  
  // get indice of reads to be matched; limit reads to msc::maxNR
  vector<size_t> ii,kk;
  if ( check_length<0 ) {     // only skip when in all match mode
    sampleidx(b_MS.size(), msc::maxMR, ii);
    sampleidx(b_SM.size(), msc::maxMR, kk);
    ED_st::linc=max( 1.0, (double)b_MS.size()/(double)msc::maxMR );
    ED_st::rinc=max( 1.0, (double)b_SM.size()/(double)msc::maxMR );
    if ( (ii.size()<b_MS.size() || kk.size()<b_SM.size() ) && thread_id==0 ) 
      cerr << "sampleing " << ii.size() << " x " << kk.size() << " reads" << endl;
  }
  else {
    sampleidx(b_MS.size(), b_MS.size(), ii);
    sampleidx(b_SM.size(), b_SM.size(), kk);
    ED_st::linc=1.0;
    ED_st::rinc=1.0;
  }
  
  size_t ndiv= ii.size()/NUM_THREADS + bool(ii.size()%NUM_THREADS);
  size_t istart=thread_id*ndiv;
  size_t iend=min(thread_id*ndiv+ndiv, ii.size());
  
  int imm=b_MS[istart].core.pos/1000000;
  size_t kstart=0;
  
  for(size_t si=istart; si<iend; ++si) {
    size_t i=ii[si];
    if ( b_MS[i].core.pos /1000000 > imm ) {
      imm=b_MS[i].core.pos /1000000 ;
      pthread_mutex_lock(&nout);
      cerr << "#thread " << thread_id << "\t" 
	   << string(msc::fp_in->header->target_name[ b_MS[i].core.tid ]) 
	   << "@" << commify( b_MS[i].core.pos ) 
	   << endl; 
      pthread_mutex_unlock(&nout);
    }
    
    readMS=get_qseq( &b_MS[i] );
    
    for(size_t sk=kstart; sk<kk.size(); ++sk ) {
      size_t k=kk[sk];
      if ( check_length > 0 ) {
	if ( b_SM[k].core.pos + check_length < b_MS[i].core.pos ) {
	  kstart=sk+1;
	  continue;
	}
	if ( b_SM[k].core.pos  > b_MS[i].core.pos + check_length ) break;
      }
      
      readSM=get_qseq( &b_SM[k] );
      
      int p1=-1; // 0 based position on F2 where strings begin overlap
      vector<int> p_err(0); // positions on F2 where mismatch happens
      bool match=false;
      match=string_overlap(readMS, readSM, minOver, maxErr, p1, p_err);
      if ( p1<0 || !match ) continue;
      int cl=readMS.length() > p1+readSM.length() ?   // overlap length
	readSM.length() : readMS.length() - p1 ;
      if ( (int)p_err.size()*12 > cl ) continue;
      
      int F2, R1, e_dis;
      ED_st ipair;
      get_break_points(FASTA, &b_MS[i], &b_SM[k], p1, p_err, F2, R1, e_dis);
      int ml=p1+readSM.length();                      // merged length
      if ( e_dis*15 > ml ) continue;
      if ( R1-F2==1 ) continue; 
      
      size_t pre_found=0;
      if ( check_length < 0 ) pre_found=0;
      else {
	if ( bp.size()>20 ) {
	  for(size_t l=bp.size()-1; l>=bp.size()-10; --l) {
	    if ( F2==bp[l].F2 && R1==bp[l].R1 && e_dis==bp[l].ED ) {
	      pre_found=l;
	      break;
	    }
	  }
	}
      }
      
      if ( pre_found>0 ) bp[pre_found].count+=1;
      else {
	ipair.iL=i;
	ipair.F2=F2;
	ipair.iR=k;
	ipair.R1=R1;
	ipair.ED=e_dis;
	ipair.count=1;
	bp.push_back(ipair);
      }
      
      if ( bp.size()>bpreserve/4*3 ) {
	sort(bp.begin(), bp.end(), sort_bp);
	pthread_mutex_lock(&nout);
	ibp.insert(ibp.end(), bp.begin(), bp.end() );
	pthread_mutex_unlock(&nout);
	m_count+=bp.size();
	bp.clear();
      }
      
      // in case there are too many pairs, save to disk; 
      // this probem is solved by skipping  with iinc, kinc  
      // if ( bp.size()>bpreserve/4*3 ) {
      // write_ED_st_to_file(bp, tmpfile(thread_id) );
      // bp.clear();
      // }
      
      /*
      pthread_mutex_lock(&nout);
      cerr << "MSSM\t" << i << "\t" << k << "\t" << F2 << "\t" << R1 << "\t" << e_dis << endl;
      pthread_mutex_unlock(&nout);
      */
    }
  }
  
  if ( bp.size()>0 ) {
    sort(bp.begin(), bp.end(), sort_bp);
    pthread_mutex_lock(&nout);
    ibp.insert(ibp.end(), bp.begin(), bp.end() );
    pthread_mutex_unlock(&nout);
    m_count+=bp.size();
    bp.clear();
  }
  
  // write_ED_st_to_file(bp, tmpfile(thread_id) );
  // bp.clear();
  
  pthread_mutex_lock(&nout);
  cerr << "thread " << thread_id << " returned " << m_count << endl;
  pthread_mutex_unlock(&nout);
  
  return;
}

//! pass pointers to threads
void* multithreads_search_wrapper(void* threadarg)
{
  struct exhaustive_search_thread_data_t *my_data = 
    (struct exhaustive_search_thread_data_t *) threadarg;
  
  int thread_id = my_data->thread_id;
  int NUM_THREADS = my_data->NUM_THREADS;
  vector<bam1_t>* b_MS = my_data->b_MS;
  vector<bam1_t>* b_SM = my_data->b_SM;
  string* FASTA = my_data->FASTA;
  int check_length = my_data->check_length;
  vector<ED_st>* bp = my_data->bp;
  
  match_reads_for_exhaustive_search(thread_id,
				    NUM_THREADS,
				    *b_MS,
				    *b_SM,
				    *FASTA, 
				    check_length,
				    *bp );
  
  pthread_exit((void*) 0);
}

void* multithreads_calibrate_wrapper(void* threadarg)
{
  struct exhaustive_search_thread_data_t *my_data = 
    (struct exhaustive_search_thread_data_t *) threadarg;
  
  int me = my_data->thread_id;
  int NT = my_data->NUM_THREADS;
  vector<bam1_t>* b_MS = my_data->b_MS;
  vector<bam1_t>* b_SM = my_data->b_SM;
  string* FASTA = my_data->FASTA;
  
  for(size_t i=me; i<(*b_MS).size(); i+=NT) calibrate_cigar_pos((*FASTA), &(*b_MS)[i]);
  for(size_t i=me; i<(*b_SM).size(); i+=NT) calibrate_cigar_pos((*FASTA), &(*b_SM)[i]);
  
  pthread_exit((void*) 0);
}

void multithreads_calibrate(vector<bam1_t>& b_MS, vector<bam1_t>& b_SM, string& FASTA) 
{
  if ( FASTA.size()<2 ) return;
  
  pthread_attr_t attr;
  vector<pthread_t> thread_ca(msc::numThreads);
  vector<exhaustive_search_thread_data_t> thread_ca_arg(msc::numThreads);
  
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  
  for(int i=0; i< msc::numThreads; ++i) {
    thread_ca_arg[i].thread_id=i;
    thread_ca_arg[i].NUM_THREADS=msc::numThreads;
    thread_ca_arg[i].b_MS = &b_MS;
    thread_ca_arg[i].b_SM = &b_SM;
    thread_ca_arg[i].FASTA = &FASTA;
    int rc = pthread_create(&thread_ca[i], 
			    &attr, 
			    multithreads_calibrate_wrapper, 
			    &thread_ca_arg[i] );
    if (rc) {
      cerr << "ERROR; return code from pthread_create() is " << rc << endl; 
      exit(-1);
    }
  }
  for (int i=0; i<msc::numThreads; i++) pthread_join(thread_ca[i], NULL);

  return;
}


void multithreads_read_matching(vector<bam1_t>& b_MS,
				vector<bam1_t>& b_SM,
				string& FASTA, 
				int check_length,
				vector<ED_st>& bp)
{
  if ( b_MS.size()<1 || b_SM.size()<1 || FASTA.size()<5 ) return;
  
  multithreads_calibrate(b_MS, b_SM, FASTA) ;
  
  pthread_attr_t attr;
  vector<pthread_t> thread_es(msc::numThreads);
  vector<exhaustive_search_thread_data_t> thread_es_arg(msc::numThreads);
  
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  
  bp.clear();
  bp.reserve(4000000);
  
  for(int i=0; i< msc::numThreads; ++i) {
    thread_es_arg[i].thread_id=i;
    thread_es_arg[i].NUM_THREADS=msc::numThreads;
    thread_es_arg[i].b_MS = &b_MS;
    thread_es_arg[i].b_SM = &b_SM;
    thread_es_arg[i].FASTA = &FASTA;
    thread_es_arg[i].check_length = check_length;
    thread_es_arg[i].bp = &bp;
    int rc = pthread_create(&thread_es[i], 
			    &attr, 
			    multithreads_search_wrapper, 
			    &thread_es_arg[i] );
    if (rc) {
      cerr << "ERROR; return code from pthread_create() is " << rc << endl; 
      exit(-1);
    }
  }
  for (int i=0; i<msc::numThreads; i++) pthread_join(thread_es[i], NULL);
  
  if ( msc::numThreads>1 ) cerr << "all threads returned: " << bp.size() << endl;
  
  return;
}

void get_softclip_reads(int ref, int beg, int end, string& FASTA,
			vector<bam1_t>& b_MS, vector<bam1_t>& b_SM) 
{
  bam1_t *b=NULL; b = bam_init1();
  bam_iter_t iter=0;
  bam1_t ibam;
  
  msc::bdata.reserve(300000000);
  b_MS.clear();
  b_SM.clear();
  
  iter = bam_iter_query(msc::bamidx, ref, beg, end);
  size_t count=0;
  while( bam_iter_read(msc::fp_in->x.bam, iter, b)>0 ) {
    if ( b->core.tid != msc::bam_ref ) break;
    if ( b->core.pos > end ) break;
    RSAI_st iread;
    
    count++;
    if ( count%1000000==0 ) 
      cerr << "#processed " << commify(count) << " reads at pos " 
	   << string(msc::fp_in->header->target_name[ref]) 
	   << "@" << commify(b->core.pos) 
	   << endl;
    
    if ( ! is_keep_read(b, FASTA, iread) ) continue;
    
    if ( msc::dumpBam ) {
      bam_aux_append(b, "ns", 'i', 4, (uint8_t*)&iread.S);
      bam_aux_append(b, "nm", 'i', 4, (uint8_t*)&iread.mmm);
      samwrite(msc::fp_out, b);
      continue;
    }
    
    _save_read_in_vector(b, ibam, msc::bdata);

    if ( iread.sbeg > iread.pos ) b_MS.push_back(ibam); // type M...S
    else b_SM.push_back(ibam);  // type S...M
    
  }
  bam_destroy1(b);
  bam_iter_destroy(iter);
  
  cerr << b_MS.size() << "\t" << b_SM.size() << "\t" << msc::bdata.size() << endl;
  
  return;
}

void exhaustive_search(vector<bam1_t>& b_MS, vector<bam1_t>& b_SM,
		       int min_pair_length, string& FASTA,
		       vector<pairinfo_st>& mcbp) 
{
  mcbp.clear();
  
  vector<ED_st> bp(0);
  
  multithreads_read_matching(b_MS, b_SM, FASTA, min_pair_length, bp);
  reduce_matched_break_points(bp, mcbp);
  cerr << "Done softclips matching\n" << endl;
  
  // release memory
  vector<ED_st>(0).swap(bp); 
  vector<uint8_t> (0).swap(msc::bdata);
  vector<bam1_t> (0).swap(b_MS);
  vector<bam1_t> (0).swap(b_SM);
  
  int old_minOverlap = msc::minOverlap;
  msc::minOverlap += msc::minOverlapPlus;  
  for(int i=0; i<(int)mcbp.size(); ++i) {
    pairinfo_st ibp=mcbp[i];
    
    stat_region(ibp, FASTA, msc::bam_l_qseq);
    assess_rd_rp_sr_infomation(ibp);
    mcbp[i]=ibp;
    
    // good signal pass
    if ( mcbp[i].rdscore>=2 && mcbp[i].rpscore>=2 ) continue;
    // if ( mcbp[i].rpscore>=2 && mcbp[i].sr_count>=9 ) continue;
    // too short for reads matching
    if ( abs(mcbp[i].F2-mcbp[i].R1)<msc::bam_l_qseq/2 ) continue;
    
    string cnv=cnv_format1(mcbp[i]);
    for(size_t t=0; t<cnv.size(); ++t) if (cnv[t]=='\t') cnv[t]=' ';
    cerr << "checking: " << cnv << endl;
    match_reads_for_pairs(ibp, FASTA, -1, true);
    assess_rd_rp_sr_infomation(ibp);
    mcbp[i]=ibp;
    
    if ( ibp.F2_sr>1 && ibp.R1_sr>1 ) {
      // if MR found, checked new matched
      ibp.F2=ibp.MS_F2;
      ibp.R1=ibp.MS_R1;
      stat_region(ibp, FASTA, msc::bam_l_qseq);
      assess_rd_rp_sr_infomation(ibp);
    }
    
    cnv=mr_format1(ibp);
    for(size_t t=0; t<cnv.size(); ++t) if (cnv[t]=='\t') cnv[t]=' ';
    cerr << "matched:  " << cnv << endl;
    cnv=cnv_format1(ibp);
    for(size_t t=0; t<cnv.size(); ++t) if (cnv[t]=='\t') cnv[t]=' ';
    cerr << "matched:  " << cnv << endl;
    
    if ( ibp.srscore>0 && 
	 ibp.rdscore>=mcbp[i].rdscore && 
	 ibp.rpscore>=mcbp[i].rpscore  &&
	 bool(ibp.MS_F2>ibp.MS_R1 ) == bool(ibp.F2>ibp.R1 ) ) mcbp[i]=ibp;
    if ( ibp.srscore>=2 && 
	 (ibp.rdscore>=mcbp[i].rdscore || ibp.rpscore>=mcbp[i].rpscore)  &&
	 ibp.un<msc::minOverlap &&
	 bool(ibp.MS_F2>ibp.MS_R1 ) == bool(ibp.F2>ibp.R1 ) ) mcbp[i]=ibp; 
    
    cnv=cnv_format1(mcbp[i]);
    for(size_t t=0; t<cnv.size(); ++t) if (cnv[t]=='\t') cnv[t]=' ';
    cerr << "updated:  " << cnv << "\n" << endl;
  }
  msc::minOverlap = old_minOverlap ;
  
  return;
}

void exhaustive_search(int ref, int beg, int end, 
		       int min_pair_length, string& FASTA,
		       vector<pairinfo_st>& mcbp) 
{
  mcbp.clear();
  vector<bam1_t> b_MS(0), b_SM(0);
  get_softclip_reads(ref, beg, end, FASTA, b_MS, b_SM) ;
  exhaustive_search(b_MS, b_SM, min_pair_length, FASTA,  mcbp) ;
  return;
}
