#ifndef _EXHAUSTIVE_H
#define _EXHAUSTIVE_H

struct exhaustive_search_thread_data_t {
  int thread_id;
  int NUM_THREADS;
  int check_length;
  string* FASTA;
  vector<bam1_t>* b_MS;
  vector<bam1_t>* b_SM;
  vector<uint8_t>* bdata;
  vector<ED_st>* bp;
};

void reduce_matched_break_points(vector<ED_st>& bp, vector<ED_st>& reduced);

void match_reads_for_exhaustive_search(int thread_id,
				       int NUM_THREADS,
				       vector<bam1_t>& b_MS,
				       vector<bam1_t>& b_SM,
				       string& FASTA, 
				       int check_length,
				       vector<ED_st>& bp);

void multithreads_read_matching(vector<bam1_t>& b_MS,
				vector<bam1_t>& b_SM,
				string& FASTA, 
				int check_length,
				vector<ED_st>& bp);

void exhaustive_search(vector<bam1_t>& b_MS, vector<bam1_t>& b_SM,
		       int min_pair_length, string& FASTA,
		       vector<pairinfo_st>& mcbp) ;

void exhaustive_search(int ref, int beg, int end, int min_pair_length, string& FASTA,
		       vector<pairinfo_st>& pairbp) ;


#endif
