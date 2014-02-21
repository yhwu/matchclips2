#ifndef _PAIRGUIDE_H
#define _PAIRGUIDE_H

void check_inner_pair_ends(const bam1_t *b, 
			   int& r1, bool& ia1, int& r2, bool& ia2);

void check_outer_pair_ends(const bam1_t *b, 
			   int& r1, bool& ia1, int& r2, bool& ia2);

int check_normalpairs_cross_pos(int ref, int end);

int check_abnormalpairs_cross_region(int ref, int F2, int R1);

void check_normal_and_abnormalpairs_cross_region(int ref, int F2, int R1,
						 int& p_F2, int& p_R1, int& p_F2R1);

void check_pair_group(vector<intpair_st>& pairs, vector<pairinfo_st>& bpinfo);

void get_break_points(const string& FASTA, bam1_t *bF2, bam1_t *bR1, int p1, vector<int>& p_err, int& F2, int& R1, int& e_dis);

void match_reads_for_pairs(pairinfo_st& ipairbp, string& FASTA, int dx, bool pointmode);
void match_reads_for_pairs(vector<pairinfo_st>& pairbp, string& FASTA, int dx, bool pointmode);

void stat_pair_group(vector<pairinfo_st>& bp);

void stat_pair_group(vector<intpair_st>& bp, vector<pairinfo_st>& bpinfo) ;

void pair_guided_search(vector<intpair_st>& pairs, string& FASTA, 
			vector<pairinfo_st>& pairbp) ;
void pair_guided_search(int ref, int beg, int end, int min_pair_length, string& FASTA,
			vector<pairinfo_st>& pairbp) ;

#endif
