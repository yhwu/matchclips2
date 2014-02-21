#ifndef _READ_REF_H
#define _READ_REF_H

bool read_fasta(string fastaFile, string chr, string& ref);
void get_N_regions(string& fasta, vector<int>& beg, vector<int>& end);
void load_reference(string fastaFile, string chr, string& ref);

#endif
