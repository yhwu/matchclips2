#ifndef _MATCH_READS_H
#define _MATCH_READS_H

#define MAX_THREADS 64

#define CNVTYPE "DAU"
#define TYPE_DEL 0
#define TYPE_DUP 1
#define TYPE_UNKNOWN 9

class msc {
public: 
  static int verbose;
  static long int BUFFERSIZE;
  static int pid;
  static long seed;
  static bool debug;
  static int start;
  static int end;
  static ofstream fout;
  static string mycommand;
  static string execinfo;
  static string bamFile;
  static vector<string> bamRegion;
  static string refFile;
  static string cnvFile;
  static string outFile;
  static string logFile;
  static string function;
  static int numThreads;
  static int maxMR;
  static int errMatch;
  static int minSNum;
  static int minOverlap;
  static int minOverlapPlus;
  static int minMAPQ;
  static int minBASEQ;
  static bool search_length_set_by_user;
  static int maxDistance;
  static int minClusterSize;
  static bool noSecondary;
  static bool dumpBam;
  static bool header;
  static int bam_ref;
  static int bam_l_qseq;
  static bool bam_is_paired;
  static bool bam_pe_disabled;
  static bool bam_pe_set_by_user;
  static int bam_pe_insert;
  static int bam_pe_insert_sd;
  static int bam_rd;
  static int bam_rd_sd;
  static int bam_tid;
  static vector<uint8_t> bdata;
  static vector<uint8_t> bpdata;
  static vector<uint32_t> rd;
  static vector<string> bam_target_name;
  static string chr;
  static samfile_t *fp_in;
  static bam_index_t *bamidx;
  static samfile_t *fp_out;
  ~msc(){};
};

struct RSAI_st {        // read suffic array index
  int tid;              // target id as in TARGET
  int pos;              // 1-based position on reference
  size_t p1;          // 0-based position in memory buffer
  int q1;
  int len;              // length of read
  int len_cigar;        // length of CIGAR
  int M;                // number of M bases
  int Mrpos;            // position of first matched base in the read
  int S;                // number of S bases
  int sbeg;             // 1-based pos begin of S, sbeg=Spos
  int send;             // 1-based pos end of S
  int pos_beg;          // 1-based pos of first base of read on reference
  int pos_end;          // 1-based pos of last base of read on reference
  int mms;              // mismatched wrt REF in S part of read
  int mmm;              // mismatched wrt REF in M part of read
  RSAI_st():tid(-1),
	    pos(0),
	    p1(0),
	    q1(0),
	    len(0), 
	    len_cigar(0), 
	    M(0),
	    Mrpos(0),
	    S(0),
	    pos_beg(0),
	    pos_end(0),
	    mms(0),
	    mmm(0){};
};

struct BREAKSEEK_st {        // describes break point discovery
  static string BUFFER;      // buffer to save information
  int tid;              // target id as in TARGET
  int P1;               // first break point position in chromosome
  int P2;               // second break point position in chromosome
  int L1;               // length of readMS
  int S1;               // clipped bases of readMS
  int n1;               // mismatched bases in clipped port of readMS
  int L2;               // length of readSM
  int S2;               // clipped bases of readSM
  int n2;               // mismatched bases in clipped port of readMS
  int q1;               // map quality of read1
  int q2;               // map quality of read2
  int CL;               // common string length of reads
  int len;              // length of variation
  int ED;               // edit_distance between matched and reference
  int ED5p;             // edit_distance between merged and reference at 5'
  int ED3p;             // edit_distance between merged and reference at 3'
  int UN;               // uncertainty of break point
  char T;               // D for DEL, A for ADD
  char SVSEQ[51];       // 30 chars befroe and after the SV, includes BP 
  int count;            // number of supportive pairs
  bool isit;            // number of supportive pairs
  string INSEQ;         // inserted sequence for Insertion only
  string MERGE;         // sequence of merged reads
  string TAG;           // comment on cnv
  BREAKSEEK_st():tid(-9),
		 P1(0),
		 P2(0),
		 L1(0),
		 n1(0),
		 L2(0),
		 n2(0),
		 CL(0),
		 len(0),
		 ED(0),
		 UN(0),
		 count(0),
		 isit(true),
		 INSEQ("."),
		 MERGE(200,'C'),
		 TAG(30,'d'){};
};
//  static int ctid;      // current tid
//  string RNAME( ) { return ( msc::bam_target_name[tid] ); }
//map<string, int> BREAKSEEK_st::RNAMEINT;
//int BREAKSEEK_st::ctid;

struct intpair_st {     
  int F2;
  bool F2_acurate;
  int R1;
  bool R1_acurate;
  int pair_count;
  intpair_st(): F2(0), F2_acurate(0), R1(0), R1_acurate(0), pair_count(0) {};
};
bool sort_pair(const intpair_st& p1, const intpair_st& p2);
bool sort_pair_len(const intpair_st& p1, const intpair_st& p2);

// edit distance structure information
struct ED_st {     
  static int square;
  static double linc;
  static double rinc;
  int iL;  // index of b in left array
  int iR;  // index of b in right array
  int F2;  // calculated break point at 5' end
  int R1;  // calculated break point at 3' end
  int ED;  // edit distance       
  int count;
  ED_st(): iL(-1), iR(-1), F2(-1), R1(-1), ED(100000), count(1) {};
};


struct pairinfo_st {     
  int tid;   // ref id
  int F2;    // right position of Forward read
  int F2_rp; // num of normal pairs cross F2
  int F2_rd; // read depth within length of a long pair left of F2
  bool F2_acurate; // 1: by calend; 0: by mpos+l_qseq
  int R1;    // left position of Reverse Complement read
  int R1_rp; // num of normal pairs cross R1
  int R1_rd; // read depth within length of a long pair right of R1
  bool R1_acurate; // 1, always given by pos or mpos
  int un;          // break point uncertainty
  int MS_F2;       // F2 calculated by matching
  int MS_F2_rd;       // F2 calculated by matching
  int MS_R1;       // R1 ...
  int MS_R1_rd;       // R1 ...
  int MS_ED;       // edit distance
  int MS_ED_count; // count of pairs with min ED
  int MS_S_count;   // pairs discovered by pure soft clip reads
  int F2_sr;       // reads across F2 match with those across R1
  int R1_sr;       // reads across R1 ..                      F2
  int sr_ed;       // total such
  int sr_count;    // total such
  int rpscore; // score of 2 is considered convincing
  int rdscore; // score of 2 is considered convincing
  int srscore; // score of 2 is considered convincing
  int pair_count; // at least >=6 (3 templates) to be considered 
  int rd;         // read depth between F2 and R1 
  pairinfo_st(): 
    tid(-1), F2(-1), F2_rp(-1), F2_rd(-1), F2_acurate(0), 
    R1(-1), R1_rp(-1), R1_rd(-1), R1_acurate(0), un(-1), 
    MS_F2(-1), MS_F2_rd(-1), MS_R1(-1), MS_R1_rd(-1), MS_ED(-1), MS_ED_count(-1), MS_S_count(-1), 
    F2_sr(-1), R1_sr(-1), sr_ed(-1), sr_count(-1),
    rpscore(-1), rdscore(-1), srscore(-1), pair_count(-1), rd(-1) {};
};


struct thread_data_t {
  int thread_id;
  int NUM_THREADS;
  string* SEQBUFFER;
  vector<RSAI_st>* MSREAD;
  vector<RSAI_st>* SMREAD;
  string* FASTA;
  string* RNAME;
  vector<BREAKSEEK_st>* BPSEEK;
};

string tmpfile(int thread_id);

void write_ED_st_to_file(vector<ED_st>& ed, string fn);

string cnv_format1(pairinfo_st &bp);
string mr_format1(pairinfo_st &bp);

void match_MS_SM_reads(int argc, char* argv[]);

#endif
