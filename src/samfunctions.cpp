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
#include <unistd.h>
#include <assert.h>
//#include "functions.h"
//#include "matchreads.h"
using namespace std;

#include "samfunctions.h"

template <class T>
static inline std::string to_string (const T& t)
{
  std::stringstream ss;
  ss << t;
  return ss.str();
}



static inline int is_overlap(uint32_t beg, uint32_t end, const bam1_t *b)
{
  uint32_t rbeg = b->core.pos;
  uint32_t rend = b->core.n_cigar ? 
    bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + 1;
  return (rend > beg && rbeg < end);
}

/*!
  @abstract     Duplicate an alignment; but only keep cigar and qseq
  @param  src   source alignment struct
  @return       pointer to the destination alignment struct
  @comment      adapted from bam_dup1()
*/
bam1_t *bam_shorten2(const bam1_t *src)
{
  bam1_t *b;
  b = bam_init1();
  *b = *src;
  b->core.l_qname=0;
  b->l_aux=0;
  b->data_len=(b)->core.n_cigar*4 + ((b)->core.l_qseq + 1)/2;
  b->m_data = b->data_len;
  b->data = (uint8_t*)calloc(b->data_len, 1);
  memcpy(b->data, src->data+src->core.l_qname, b->data_len);
  return b;
}
bam1_t *bam_shorten2(const bam1_t *src, uint8_t* bdata)
{
  bam1_t *b;
  b = bam_init1();
  *b = *src;
  b->core.l_qname=0;
  b->l_aux=0;
  b->data_len=(b)->core.n_cigar*4 + ((b)->core.l_qseq + 1)/2;
  b->m_data = b->data_len;
  //  b->data = (uint8_t*)calloc(b->data_len, 1);
  b->data = bdata;
  memcpy(b->data, src->data+src->core.l_qname, b->data_len);
  return b;
}
/*!
  @abstract     Duplicate an alignment; but only keep cigar qseq and qual
  @param  src   source alignment struct
  @return       pointer to the destination alignment struct
  @comment      adapted from bam_dup1()
*/
bam1_t *bam_shorten(const bam1_t *src)
{
  bam1_t *b;
  b = bam_init1();
  *b = *src;
  b->core.l_qname=0;
  b->l_aux=0;
  b->data_len=(b)->core.n_cigar*4 + (b)->core.l_qseq + ((b)->core.l_qseq + 1)/2;
  b->m_data = b->data_len;
  b->data = (uint8_t*)calloc(b->data_len, 1);
  memcpy(b->data, src->data+src->core.l_qname, b->data_len);
  
  return b;
}


/*! callback function 
 *  check map quality and pairs that support the break points
 */
static int fetch_func(const bam1_t *b, void *data)  
{  
  bam_pileup_api_wrapper_t *tmp = (bam_pileup_api_wrapper_t*)data;  
  bam_plbuf_push(b, tmp->buf);  
  if ( is_overlap(tmp->bp1-1, tmp->bp2-1, b) ) {
    if ( (int)b->core.qual <= 1  ) ++(tmp->c0);
    else ++(tmp->c1);
  }
  
  //*! BAM_FPROPER_PAIR may flag long/negative inserts as not properly paired
  if ( (b->core.flag & BAM_FPAIRED) && (b->core.mtid == b->core.tid) ) {
    int r1,r2,ISIZE;
    if ( (b->core.flag & BAM_FREVERSE) ==0 ) {   // read is Forward 
      r1= b->core.pos+1;
      r2= b->core.mpos+b->core.l_qseq+1;       // only approximate
      ISIZE=b->core.isize;
    }
    else {  // read is Reverse
      r1= b->core.mpos+1; 
      r2=b->core.pos+b->core.l_qseq+1;         // only approximate
      ISIZE=-b->core.isize;
    }
    bool is_support=false;
    if ( tmp->T=='D' && r1<tmp->bp1 && r2>tmp->bp2 ) is_support=true;
    if ( tmp->T=='A' && ISIZE<0  ) is_support=true;
    tmp->pair += is_support;
  }
  
  return 0;  
}  

/*! callback function */
static int pileup_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data)  
{  
  bam_pileup_api_wrapper_t *tmp = (bam_pileup_api_wrapper_t*)data;  
  if ((int)pos >= tmp->beg && (int)pos < tmp->end) { 
    //printf("%s\t%d\t%d\n", tmp->in->header->target_name[tid], pos + 1, n);  
    tmp->n[pos-tmp->beg]=n;
  }
  return 0;  
}

void get_qseq(const bam1_t *b,  string& seq)  
{
  seq.clear();
  uint8_t *s = bam1_seq(b);
  // cerr << get_cigar(b) << endl;
  // cerr << bam_format1(msc::fp_in->header ,b) << endl;
  for (int i=0; i<b->core.l_qseq; ++i) 
    seq.push_back(bam_nt16_rev_table[bam1_seqi(s, i)]);
  // cerr << seq << "\n" << endl;
  return;
}
string get_qseq(const bam1_t *b)
{  
  string seq;
  get_qseq(b, seq);
  return seq;
}

void get_rname(const bam_header_t *header, const bam1_t *b,  string& rname)
{
  if (header) rname=string(header->target_name[b->core.tid]);
  else rname=to_string(b->core.tid);
  return;
}  

void get_cigar(const bam1_t *b,  string& cigar)  
{
  cigar.clear();
  if (b->core.n_cigar == 0) cigar="*";
  else {
    for (unsigned int i = 0; i < b->core.n_cigar; ++i) {
      cigar+=to_string( bam1_cigar(b)[i]>>BAM_CIGAR_SHIFT );
      cigar.push_back( BAM_CIGAR_STR[ bam1_cigar(b)[i]&BAM_CIGAR_MASK] );
    }
  }
  return;
}
string get_cigar(const bam1_t *b)  
{
  string cigar;
  get_cigar(b,  cigar) ;
  return cigar;
}
string get_cigar(const POSCIGAR_st& m)  
{
  string cigar="*";
  if ( m.pos<m.base || m.op.size()==0 ) return cigar;
  cigar="";
  for(int i=0; i<(int)m.op.size(); ++i) {
    cigar+=to_string( m.nop[i] )+string(1, BAM_CIGAR_STR[m.op[i]] );
  }
  return cigar;
}


/*! 
  @abstract resolve CIGAR
  note: pos, cpos is 1-based, anchor, iclip, qpos is 0-based
  
  @param  b     pointer to a bam
  @return m.anchor cigar[m.anchor] is first match(including del)   
  @return m.iclip  cigar[m.iclip] has the longest soft clipped bases   
  @return m.pos   match position of anchor on reference
  @return m.base base;
  @return m.op   operators, 
  @return m.nop  length of operatoration
  @return m.cop  reference position of operators, 1-based
  @return m.qop  position of operators on seq, 0-based
*/
void resolve_cigar_pos(const bam1_t *b,  POSCIGAR_st& m, int base)  
{  
  m.base=base;
  if ( base!=1 && base!=0 ) {
    m.pos=0-base; // 1-based
    cerr << "resolve_cigar_pos(): base must be 0 or 1\n";
    exit(0);
    return;
  }
  
  if ( b->core.n_cigar<=0 || b->core.pos<0 ) {
    m.pos=base-1;
    return;
  }
  
  int k;
  int ncigar=b->core.n_cigar;
  uint32_t *cigar = bam1_cigar(b);
  uint32_t end;
  
  m.l_qseq=b->core.l_qseq;
  m.tid=b->core.tid;
  m.pos=b->core.pos+m.base;
  m.qual=b->core.qual;
  m.anchor=-1;
  m.iclip=-1;
  m.op.resize(ncigar,0);  
  m.nop.resize(ncigar,0);  
  m.cop.resize(ncigar,0);  
  m.qop.resize(ncigar,0);
  
  //get position on qseq 0-based
  int ns=0;
  end=0;  
  for (k = 0; k < ncigar; ++k) {
    int op = bam_cigar_op(cigar[k]);
    int l = bam_cigar_oplen(cigar[k]);
    m.op[k]=op;
    m.nop[k]=l;
    m.qop[k]=end;
    if ( op == BAM_CMATCH || 
	 op == BAM_CINS || 
	 op == BAM_CSOFT_CLIP || 
	 op == BAM_CEQUAL || 
	 op == BAM_CDIFF ) end+=l;
    
    if (  op == BAM_CSOFT_CLIP && l > ns) { ns=l; m.iclip=k; }
    
    if ( m.op[k] == BAM_CMATCH || 
	 m.op[k] == BAM_CDEL || 
	 m.op[k] == BAM_CEQUAL || 
	 m.op[k] == BAM_CDIFF) { if ( m.anchor<0 ) m.anchor=k; }
    
  }
  if ( m.anchor<0 ) { m.pos=0-base; return; }
  
  // position on reference, base is given by parameter
  // BAM_CSOFT_CLIP is added because we need positions of S 
  end = m.pos;
  for (k = m.anchor; k < ncigar; ++k) {
    m.cop[k]=end;
    if (  m.op[k] == BAM_CMATCH || 
	  m.op[k] == BAM_CDEL || 
	  m.op[k] == BAM_CREF_SKIP || 
	  m.op[k] == BAM_CSOFT_CLIP ) end += m.nop[k] ;
  }
  end = m.pos;
  for (k = m.anchor-1; k >= 0; --k) {
    if (  m.op[k] == BAM_CMATCH || 
	  m.op[k] == BAM_CDEL || 
	  m.op[k] == BAM_CREF_SKIP || 
	  m.op[k] == BAM_CSOFT_CLIP )  end -= m.nop[k] ;
    m.cop[k]=end;
  }
  
  if ( m.qop.back() + m.nop.back()*(m.op.back()!=BAM_CHARD_CLIP) 
       != m.l_qseq ) {
    cerr << "[resolve_cigar_pos] Something wrong reading CIGAR. " 
	 << "length in bam_t b: " << m.l_qseq << " , "
	 << "length calculated: " << m.qop.back()+m.nop.back() << "\t"
	 << endl;
    cerr << get_qseq(b) << "\n"
	 << get_cigar(b) << endl;
    //cerr << bam_format1(msc::fp_in->header, b)  << endl;
	
  }
  
  return;  
}  
void resolve_cigar_pos(const bam1_t *b,  POSCIGAR_st& m)  {
  resolve_cigar_pos(b,  m, 1);
}

void resolve_cigar_pos(int POS, string& CIGAR, POSCIGAR_st& m)  
{  
  if ( CIGAR=="*" ) { m.pos=0; return; }
  
  char numchar[]="99999999999999";
  const string ops=BAM_CIGAR_STR;
  
  m.base=1;
  m.pos=POS;
  m.anchor=-1;
  m.iclip=-1;
  m.op.resize(0);
  m.nop.resize(0);
  m.l_qseq=0;
  
  for(int i=0, k=0; i<(int)CIGAR.size(); ++i) {
    if ( CIGAR[i]>='0' && CIGAR[i]<='9'  ) {
      numchar[k]=CIGAR[i];
      ++k;
      if ( k>9 ) { cerr << "CIGAR error" << CIGAR << endl; exit(0); }
    }
    else {
      if ( CIGAR[i]=='*' ) { m.pos=0; return; }
      size_t cigar_num=ops.find(CIGAR[i]);
      if ( cigar_num==string::npos) { 
	cerr << "CIGAR error" << CIGAR << endl; 
	exit(0); 
      }
      numchar[k]=0;
      m.op.push_back(cigar_num);
      m.nop.push_back(atoi(numchar));
      k=0;
    }
  }
  
  m.cop.resize(m.op.size());
  m.qop.resize(m.op.size());
  
  int ns=0;
  uint32_t end=0; //position on seq 0-based
  for (int k = 0; k < (int)m.op.size(); ++k) {
    
    int op=m.op[k];
    int l=m.nop[k];
    m.qop[k]=end;
    
    if ( op == BAM_CMATCH || 
	 op == BAM_CINS || 
	 op == BAM_CSOFT_CLIP || 
	 op == BAM_CEQUAL || 
	 op == BAM_CDIFF ) end+=l;
    
    if (  op == BAM_CSOFT_CLIP && l > ns) { ns=l; m.iclip=k; }
    
    if ( m.op[k] == BAM_CMATCH || 
	 m.op[k] == BAM_CDEL || 
	 m.op[k] == BAM_CEQUAL || 
	 m.op[k] == BAM_CDIFF) { if ( m.anchor<0 ) m.anchor=k; }
    
  }
  if ( m.anchor<0 ) { m.pos=0; return; }
  
  //position on reference 1-based
  // BAM_CSOFT_CLIP is added because we need positions of S 
  end = m.pos;  //position on reference
  for (int k = m.anchor; k < (int)m.op.size(); ++k) {
    m.cop[k]=end;
    if (  m.op[k] == BAM_CMATCH || 
	  m.op[k] == BAM_CDEL || 
	  m.op[k] == BAM_CREF_SKIP || 
	  m.op[k] == BAM_CSOFT_CLIP ) end += m.nop[k] ;
  }
  end = m.pos;
  for (int k = m.anchor-1; k >= 0; --k) {
    if (  m.op[k] == BAM_CMATCH || 
	  m.op[k] == BAM_CDEL || 
	  m.op[k] == BAM_CREF_SKIP || 
	  m.op[k] == BAM_CSOFT_CLIP )  end -= m.nop[k] ;
    m.cop[k]=end;
  }
  
  m.l_qseq=m.qop.back()+m.nop.back()*(m.op.back()!=BAM_CHARD_CLIP);
  return;  
}  

//! p 0 based
//! return position on reference for p base in a read
int get_pos_for_base(const POSCIGAR_st& m, int p)
{
  if ( m.pos<m.base ) return -1;
  
  bool braced=false;
  int i=0;
  for(i=0; i<(int)m.qop.size(); ++i) {
    if ( m.op[i] == BAM_CMATCH || 
	 m.op[i] == BAM_CINS || 
	 m.op[i] == BAM_CSOFT_CLIP || 
	 m.op[i] == BAM_CEQUAL || 
	 m.op[i] == BAM_CDIFF ) {
      if ( (int)m.qop[i] <= p && (int)m.qop[i]+(int)m.nop[i] > p ) {
	braced=true;
	break;
      }
    }
  }
  if ( ! braced ) {
    cerr << "get_pos_for_base(): Error finding the position " 
	 << p << " in " << get_cigar(m) << endl;
    return -1;
  }
  
  return m.cop[i]+p-m.qop[i]-m.base;
}
int get_pos_for_base(const bam1_t *b , int p)
{
  POSCIGAR_st m;
  resolve_cigar_pos(b, m, 0);  
  int pos=get_pos_for_base(m, p);
  if ( pos<0 ) {
    cerr << "get_pos_for_base(): Error finding the position "
	 << p << " in " << get_cigar(b) << endl;
  }
  return pos;
}

//! expand cigar string
void expand_cigar(const POSCIGAR_st& m, vector<int>& e_cigar )
{
  e_cigar.clear();
  if ( m.pos<m.base ) return;
  
  for(int i=0; i<(int)m.op.size(); ++i) {
    int op=m.op[i];
    if ( op == BAM_CMATCH || 
	 op == BAM_CINS || 
	 op == BAM_CSOFT_CLIP || 
	 op == BAM_CEQUAL || 
	 op == BAM_CDIFF ) {
      e_cigar.insert(e_cigar.end(), m.nop[i], op);
    }
  }

  if ( (int)e_cigar.size() != m.l_qseq ) {
    cerr << "expaned cigar size does not match l_qseq\t"
	 << e_cigar.size() << "\t" << m.l_qseq
	 << endl;
  } 
  return;
}
void expand_cigar(const bam1_t *b, vector<int>& e_cigar)
{
  e_cigar.clear();
  if ( b->core.tid<0 || b->core.n_cigar<=0 || b->core.pos<0 ) return;
  
  uint32_t *cigar = bam1_cigar(b);
  for (int k = 0; k < (int)b->core.n_cigar ; ++k) {
    int op = bam_cigar_op(cigar[k]);
    int l = bam_cigar_oplen(cigar[k]);
    if ( op == BAM_CMATCH || 
	 op == BAM_CINS || 
	 op == BAM_CSOFT_CLIP || 
	 op == BAM_CEQUAL || 
	 op == BAM_CDIFF ) {
      e_cigar.insert(e_cigar.end(), l, op);
    }
  }
  
  if ( (int)e_cigar.size() != b->core.l_qseq ) {
    cerr << "expaned cigar size does not match l_qseq\t"
	 << e_cigar.size() << "\t" << b->core.l_qseq
	 << endl;
  } 
  return;
}

//! get reference position for each base 
//! CINS  are given POS_BAM_CINS (-1)
//! CPADS are given POS_BAM_CPAD (-2)
//! CDEL CREF_SKIP CHARD_CLIP are ignroed
void expand_pos(const POSCIGAR_st& m, vector<int>& e_pos )
{
  e_pos.clear();
  if ( m.pos < m.base ) return;
  
  for(int i=0; i<(int)m.op.size(); ++i) {
    if ( m.op[i]==BAM_CDEL || 
	 m.op[i]==BAM_CREF_SKIP || 
	 m.op[i]==BAM_CHARD_CLIP ) continue;
    for( int k=0; k<(int)m.nop[i]; ++k) {
      int pos= m.cop[i]+k;
      if ( m.op[i]==BAM_CINS ) pos=POS_BAM_CINS;
      if ( m.op[i]==BAM_CPAD ) pos=POS_BAM_CPAD;
      e_pos.push_back(pos);
    }
  }
  
  if ( (int)e_pos.size() !=  m.l_qseq ) {
    cerr << "expand_pos(): error\t" 
	 <<  get_cigar(m) << "\t" << m.l_qseq << "\t" << e_pos.size()
	 << endl;
  }
  
  return;
}
void expand_pos(const bam1_t *b, vector<int>& e_pos)
{
  e_pos.clear();
  if ( b->core.tid<0 || b->core.n_cigar<=0 || b->core.pos<0 ) return;
  POSCIGAR_st m;
  resolve_cigar_pos(b, m, 0);
  if ( m.pos < m.base ) return;
  expand_pos(m, e_pos);
  if ( (int)e_pos.size() !=  b->core.l_qseq ) {
    cerr << "expand_pos(): error\t" 
	 <<  get_cigar(b) << "\t" << b->core.l_qseq << "\t" << e_pos.size()
	 << endl;
  }
  return;
}

//! reference projected onto read, taking bases from expanded_pos
//! CINS CPADS are taken from read
//! CDEL CREF_SKIP CHARD_CLIPS are ignored
string ref_projected_onto_qseq(const bam1_t *b, const string& FASTA)
{
  string e_ref="";
  
  if ( b->core.tid<0 || b->core.n_cigar<=0 || b->core.pos<0 ) return e_ref;
  vector<int> e_pos(0);
  expand_pos(b, e_pos);
  
  uint8_t *s = bam1_seq(b);
  for(int i=0; i<(int)e_pos.size(); ++i) {
    char base_ref;
    if ( (size_t)e_pos[i]>=FASTA.size() ) base_ref='N';
    else if ( e_pos[i]==POS_BAM_CPAD ) base_ref=bam_nt16_rev_table[bam1_seqi(s, i)];
    else if ( e_pos[i]==POS_BAM_CINS ) base_ref='I';
    else base_ref = FASTA[ e_pos[i] ] ;
    e_ref.push_back(base_ref);
  }
  
  if ( (int)e_ref.size() !=  b->core.l_qseq ) {
    cerr << "expand_pos(): error\t" 
	 <<  get_cigar(b) << "\t" << b->core.l_qseq << "\t" << e_pos.size()
	 << endl;
  }
  
  return e_ref;
}

string qseq_projected_onto_ref(const bam1_t *b, const string& FASTA)
{
  string e_ref="";
  
  if ( b->core.tid<0 || b->core.n_cigar<=0 || b->core.pos<0 ) return e_ref;
  vector<int> e_pos(0);
  expand_pos(b, e_pos);
  
  uint8_t *s = bam1_seq(b);
  for(int i=0; i<(int)e_pos.size(); ++i) {
    char base_ref;
    if ( (size_t)e_pos[i]>=FASTA.size() ) base_ref='N';
    else if ( e_pos[i]==POS_BAM_CPAD ) base_ref=bam_nt16_rev_table[bam1_seqi(s, i)];
    else if ( e_pos[i]==POS_BAM_CINS ) base_ref='I';
    else base_ref = FASTA[ e_pos[i] ] ;
    e_ref.push_back(base_ref);
  }
  
  if ( (int)e_ref.size() !=  b->core.l_qseq ) {
    cerr << "expand_pos(): error\t" 
	 <<  get_cigar(b) << "\t" << b->core.l_qseq << "\t" << e_pos.size()
	 << endl;
  }
  
  return e_ref;
}



/*! 
  @abstract calibrate the S parts of CIGAR by comparing FASTA and SEQ
  note: pos, cpos is 1-based, anchor, iclip, qpos is 0-based
  
  @param  FASTA   reference string
  @param  SEQ     read string
  @param  m       resolved mapping by resolve_cigar_pos()

  @return m.anchor cigar[m.anchor] is first match(including del)   
  @return m.iclip  cigar[m.iclip] has the longest soft clipped bases   
  @return m.pos   match position of anchor on reference, 1-based
  @return m.op   operators, 
  @return m.nop  reference position of operators, 1-based
  @return m.cop  reference position of operators, 1-based
  @return m.qop  position of operators on seq, 0-based
*/
int calibrate_resolved_cigar_pos(string& FASTA, string& SEQ, POSCIGAR_st& m)
{
  int nChanged=0;
  if ( m.op.size()<=1 ) return nChanged;
  
  if ( m.cop[0]+m.nop[0] < m.cop[0] ) return nChanged;
  if ( m.cop[0]+m.base < 1 ) return nChanged;
  if ( m.cop.back()+m.nop.back() >= (int)FASTA.size() ) return nChanged;

  //  cerr << m.cop[0] << "\t" << m.cop[1] << "\t" << m.cop[2] << endl;
  
  // remove hard clip at ends
  if ( m.op[0]== BAM_CHARD_CLIP ) {
    m.op.erase(m.op.begin());
    m.nop.erase(m.nop.begin());
    m.cop.erase(m.cop.begin());
    m.qop.erase(m.qop.begin());
    --m.anchor;
    --m.iclip;
  }
  if ( m.op.back()== BAM_CHARD_CLIP ) {
    m.op.pop_back();
    m.nop.pop_back();
    m.cop.pop_back();
    m.qop.pop_back();
  }
  if ( m.op.size()<=1 ) return nChanged;
  
  // adjust S at the 3' end, m.pos don't change
  if ( m.op.back() == BAM_CSOFT_CLIP && 
       ( m.op[m.op.size()-2] == BAM_CMATCH || m.op[m.op.size()-2] == BAM_CEQUAL ) ) {
    string CLIPPEDSEQ=SEQ.substr(m.qop.back(), m.nop.back());
    string REFCLIP=FASTA.substr(m.cop.back()-m.base, m.nop.back());
    
    int iclip=m.op.size()-1;
    int p1misMatch=REFCLIP.size();
    for(int i=0;i<(int)REFCLIP.size();++i) {
      if ( REFCLIP[i]==CLIPPEDSEQ[i] ) continue;
      p1misMatch=i;
      break;
    }
    int dx=p1misMatch;
    nChanged+=dx;
    //if ( dx>0 ) cerr << CLIPPEDSEQ << "\n" << REFCLIP << "\t" << dx << "\t" << iclip << endl;
    
    m.nop[iclip-1]+=dx;
    m.nop[iclip]-=dx;
    m.cop[iclip]+=dx;
    m.qop[iclip]+=dx;
  }
  
  // adjust S at the 5' end, m.pos change too 
  if ( m.op[0] == BAM_CSOFT_CLIP && 
       ( m.op[1] == BAM_CMATCH || m.op[1] == BAM_CEQUAL ) ) {
    string CLIPPEDSEQ=SEQ.substr(m.qop[0], m.nop[0]);
    string REFCLIP=FASTA.substr(m.cop[0]-m.base, m.nop[0]);
    
    int p2misMatch=0;
    for(int i=REFCLIP.size()-1;i>=0;--i) {
      if ( REFCLIP[i]==CLIPPEDSEQ[i] ) continue;
      p2misMatch=i+1;
      break;
    }
    int dx=REFCLIP.size()-p2misMatch;
    nChanged+=dx;
    //if ( dx>0 ) cerr << CLIPPEDSEQ << "\n" << REFCLIP << "\t" << dx << "\t0" << endl;
    
    m.nop[0]-=dx;
    m.nop[1]+=dx;
    m.cop[1]-=dx;
    m.qop[1]-=dx;
    m.pos=m.cop[1];
  }
  
  m.anchor=-1;
  m.iclip=-1;
  int ns=0;
  for (int k = 0; k < (int)m.op.size(); ++k) {
    if ( m.op[k] == BAM_CSOFT_CLIP && (int)m.nop[k] >= ns) { ns=m.nop[k]; m.iclip=k; }
    if ( m.op[k] == BAM_CMATCH || 
	 m.op[k] == BAM_CDEL || 
	 m.op[k] == BAM_CEQUAL || 
	 m.op[k] == BAM_CDIFF) { if ( m.anchor<0 ) m.anchor=k; }
  }
  if ( m.anchor<0 ) m.pos=0;
  
  if ( m.qop.back()+m.nop.back() != m.l_qseq ) {
    cerr << "[calibrate_resolved_cigar_pos] Something wrong calibrating CIGAR. " 
	 << "length in bam: " << m.l_qseq << " , "
	 << "length in calculated: " << m.qop.back()+m.nop.back()
	 << endl;
  }
  
  return nChanged;
}

int calibrate_cigar_pos(const string& FASTA, bam1_t *b)
{
  int nChanged=0;
  if ( b->core.n_cigar<=0 || b->core.pos<0 ) return nChanged;
  
  POSCIGAR_st m;
  resolve_cigar_pos(b, m, 0);
  
  if ( m.op.size()<=1 ) return nChanged;
  if ( m.cop[0]+m.nop[0] < m.cop[0] ) return nChanged;
  if ( m.cop[0]+m.base < 1 ) return nChanged;
  if ( m.cop.back()+m.nop.back() >= (int)FASTA.size() ) return nChanged;
  
  uint32_t *cigar = bam1_cigar(b);
  
  int num_op=0;
  for(int i=0; i<(int)m.op.size(); ++i) num_op += m.op[i]!=BAM_CHARD_CLIP;
  //cerr << "CIGAR OPs\t" << get_cigar(b) << "\t" << num_op << endl; 
  if ( num_op<2 ) return nChanged;
  
  string SEQ=get_qseq(b);
  string precigar=get_cigar(b);
  // int pre_pos=b->core.pos;
  int S_idx;
  
  vector<int> pre_e_pos(0);
  expand_pos(b, pre_e_pos);
  
  // adjust S at the 3' end, m.pos don't change
  S_idx = m.op.size() - 1 - ( m.op.back()==BAM_CHARD_CLIP );
  if ( m.op[S_idx] == BAM_CSOFT_CLIP && 
       ( m.op[S_idx-1] == BAM_CMATCH || m.op[S_idx-1] == BAM_CEQUAL ) ) {
    string CLIPPEDSEQ=SEQ.substr(m.qop[S_idx], m.nop[S_idx]);
    string REFCLIP=FASTA.substr(m.cop[S_idx]-m.base, m.nop[S_idx]);
    
    int p1misMatch=REFCLIP.size();
    for(int i=0;i<(int)REFCLIP.size();++i) {
      if ( REFCLIP[i]==CLIPPEDSEQ[i] ) continue;
      p1misMatch=i;
      break;
    }
    int dx=p1misMatch;
    nChanged+=dx;
    //if ( dx>0 ) cerr << CLIPPEDSEQ << "\n" << REFCLIP << "\t" << dx << "\t" << iclip << endl;
    
    m.nop[S_idx-1]+=dx;
    m.nop[S_idx]-=dx;
    m.cop[S_idx]+=dx;
    m.qop[S_idx]+=dx;
    
    // adjust cigar at 3' end MS
    cigar[S_idx]=bam_cigar_gen(m.nop[S_idx], m.op[S_idx] );
    cigar[S_idx-1]=bam_cigar_gen(m.nop[S_idx-1], m.op[S_idx-1]);
  }
  
  // adjust S at the 5' end, m.pos change too 
  S_idx = 0 + (m.op[0]==BAM_CHARD_CLIP);
  if ( m.op[S_idx] == BAM_CSOFT_CLIP && 
       ( m.op[S_idx+1] == BAM_CMATCH || m.op[S_idx+1] == BAM_CEQUAL ) ) {
    string CLIPPEDSEQ=SEQ.substr(m.qop[S_idx], m.nop[S_idx]);
    string REFCLIP=FASTA.substr(m.cop[S_idx]-m.base, m.nop[S_idx]);
    
    int p2misMatch=0;
    for(int i=REFCLIP.size()-1;i>=0;--i) {
      if ( REFCLIP[i]==CLIPPEDSEQ[i] ) continue;
      p2misMatch=i+1;
      break;
    }
    int dx=REFCLIP.size()-p2misMatch;
    nChanged+=dx;
    //if ( dx>0 ) cerr << CLIPPEDSEQ << "\n" << REFCLIP << "\t" << dx << "\t0" << endl;
    
    m.nop[S_idx]-=dx;
    m.nop[S_idx+1]+=dx;
    m.cop[S_idx+1]-=dx;
    m.qop[S_idx+1]-=dx;
    m.pos=m.cop[S_idx+1];
    
    // adjust cigar at 5' end SM
    cigar[S_idx]=bam_cigar_gen(m.nop[S_idx], m.op[S_idx]);
    cigar[S_idx+1]=bam_cigar_gen(m.nop[S_idx+1], m.op[S_idx+1] );
    b->core.pos=m.cop[S_idx+1];
  }
  
  vector<int> after_e_pos(0);
  expand_pos(b, after_e_pos);
  
  if ( (int)after_e_pos.size() != b->core.l_qseq ||
       (int)pre_e_pos.size() != b->core.l_qseq ) {
    cerr << "size error after expanding pos\t" 
	 << pre_e_pos.size() << "\t" 
	 << after_e_pos.size() << "\t" 
	 << get_cigar(b) << "\t" 
	 << precigar
	 << endl;
    exit(0);
  }
  for(int i=0; i<(int)after_e_pos.size(); ++i) {
    if ( after_e_pos[i] != pre_e_pos[i] ) {
      cerr << "pos error after calibrating\t" 
	   << after_e_pos[i] << "\t"
	   << pre_e_pos[i] << "\t" << get_cigar(b) << "\t" << precigar
	   << endl;
    }
  }
  
  //if ( nChanged>0 ) {
  //  cerr << "CIGAR update\t" 
  //	 << pre_pos << " " << precigar << "\t" 
  //	 << b->core.pos << " " << get_cigar(b) 
  //	 << endl;
  //}
  
  return nChanged;
}



/*! 
  @abstract calculate read depth on [p1,p2] and save them in tmp.n[];
  note: all reads are counted, equivalent to samtools mpileup -A -QB0; 
  
  @param  tmp    BAM file container; tmp.in and tmp.bamidx musted be loaded
  @param  RNAME  chromosome name
  @param  p1     1-based 5' position for pileup
  @param  p2     1-based 3' position for pileup
  @param  tmp.bp1  1-based 5' position for break point inside pileup
  @param  tmp.bp2  1-based 3' position for break point inside pileup
  @return tmp.c0   number of reads with mapq<=1
  @return tmp.c1   number of reads with mapq>1
  @return tmp.n  vector that hodlds readdepth for each position
*/
void bampileup(bam_pileup_api_wrapper_t& tmp, string RNAME, int p1, int p2)
{
  if ( tmp.bamidx==NULL ) {
    cerr << "Bam index not loaded" << endl;
    exit(0);
  }
  tmp.buf=NULL;  
  tmp.c0=0;
  tmp.c1=0;
  tmp.pair=0;
  string bamRegion=RNAME+":"+to_string(p1)+"-"+to_string(p2);
  bam_parse_region(tmp.in->header, bamRegion.c_str(), 
		   &tmp.ref, &tmp.beg, &tmp.end);
  if ( tmp.ref<0 ) {
    cerr << "Region not understood\t" << bamRegion << endl;
    exit(0);
  }
  
  tmp.n.resize(tmp.end-tmp.beg,0);
  tmp.buf = bam_plbuf_init(pileup_func, &tmp); // initialize pileup  
  //  bam_fetch(tmp.in->x.bam, tmp.bamidx, tmp.ref, tmp.beg, tmp.end, buf, fetch_func);  
  bam_fetch(tmp.in->x.bam, tmp.bamidx, tmp.ref, tmp.beg, tmp.end, &tmp, fetch_func);  
  bam_plbuf_push(0, tmp.buf); // finalize pileup  
  
  bam_plbuf_destroy(tmp.buf);  
}
