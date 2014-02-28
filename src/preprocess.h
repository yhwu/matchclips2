#ifndef _PRE_PROCESS_H
#define _PRE_PROCESS_H

// save read b's qseq and cigar in buffer v
// note since v is a vector, while v is growing, the address of v may change
// hence, ibam.data become invalid
// so, initially, ibam.data saves the relative distance to the head of v
// after v is done growing, ibam.data will be displaced by &v[0]
#define _save_read_in_vector(b, ibam, v) {				\
    ibam = *b;								\
    ibam.core.l_qname=ibam.l_aux=0;					\
    ibam.data_len=(b)->core.n_cigar*4 + ((b)->core.l_qseq + 1)/2;	\
    ibam.m_data = ibam.data_len;					\
    ibam.data = (uint8_t*) ( v.size() ) ;				\
    v.insert(v.end(),							\
	     b->data+b->core.l_qname,					\
	     b->data+b->core.l_qname+ibam.data_len);			\
  }

// save read b's complete data in buffer v
#define _save_read_in_vector_all(b, ibam, v) {				\
    ibam = *b;								\
    ibam.m_data = ibam.data_len;					\
    ibam.data = (uint8_t*) ( v.size() ) ;				\
    v.insert(v.end(), b->data, b->data+b->data_len);			\
  }

bool is_read_count_for_depth(const bam1_t *b);
bool is_read_count_for_depth(const bam1_t *b, int qual);
bool is_read_count_for_pair(const bam1_t *b);

void check_map_quality(int ref, int beg, int end, double& q0, double& q1);
int mean_readdepth(int ref, int beg, int end) ;
int median_readdepth(int ref, int beg, int end) ;
void check_cnv_readdepth(int ref, int beg, int end, int dx, 
			 int& d1, int& d2, int& din);

void check_cnv_readdepth_2(int ref, int beg, int end, int dx, 
			   int& d1, int& d2, int& din); 

void find_displacement(string& FASTA, int F2, int R1, 
		       int& dx_F2, int& dx_R1);

bool string_overlap(const string& readMS, const string& readSM, 
		    const int minOver, const int maxErr,
		    int& p1, vector<int>& p_err);

void get_break_points(const string& FASTA, bam1_t *bF2, bam1_t *bR1, int p1, vector<int>& p_err, int& F2, int& R1, int& e_dis);

bool is_keep_read(const bam1_t *b, string& FASTA, RSAI_st& iread );

void prepare_pairend_matchclip_data(int ref, int beg, int end, 
				    int min_pair_length,
				    string& FASTA,
				    vector<intpair_st>& pairs,
				    vector<bam1_t>& b_MS, vector<bam1_t>& b_SM) ;

void stat_region(pairinfo_st& ibp, string& FASTA, int dx) ;

//void stat_region(pairinfo_st& bp);

//void stat_regions(vector<pairinfo_st>& bp, string& FASTA);

void assess_rd_rp_sr_infomation(pairinfo_st& ibp) ;

void assess_rd_rp_sr_infomation(vector<pairinfo_st>& bp);

//void finalize_output(vector<pairinfo_st>& bp, 
//		     vector<pairinfo_st>& strong, 
//		     vector<pairinfo_st>& weak) ;
//void preprocess_bam_file(string bamFile, string bamRegion, int minS, int minq, int minQ, string fastaFile);

#endif
