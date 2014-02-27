#ifndef _READ_CIGAR_H
#define _READ_CIGAR_H

void read_POS_CIGAR(int POS, string& CIGAR,
		    size_t& anchor, size_t& iclip, 
		    vector<char>& op, vector<int>& opnum,
		    vector<int>& opposseq, vector<int>& oppos );

int calibrate_cigar(string& SEQ, int& POS, string& FASTA,
		    size_t& anchor, size_t& iclip, 
		    vector<char>& op, vector<int>& opnum,
		    vector<int>& opposseq, vector<int>& oppos );

#endif
