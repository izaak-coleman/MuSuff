#ifndef GENOMEMAPPER_H
#define GENOMEMAPPER_H

#include <string>
#include <vector>
#include <fstream>

#include "BranchPointGroups.h"
#include "Reads.h"

struct snv_aln_info {
 std::vector<int> SNV_pos;
 int flag;
 int chr;
 int position;
 std::string reference_mismatch;
 std::string non_mutated_cns;
 std::string mutated_cns;
};

struct single_snv {
 int flag;
 int chr;
 int position;
 char mutation_base;
 char healthy_base;
};

class GenomeMapper {

private:
  struct mutation_classes{
    std::vector<int> SNV_pos;
  };

  struct consensus_pair {
    std::string mutated;
    std::string non_mutated;
    std::string read_freq_m;
    std::string read_freq_nm;
    int mut_offset;
    int nmut_offset;
    mutation_classes mutations;

    // used to compare with cancer cns without trimming
    int healthy_ohang_right;    
    int healthy_ohang_left;
  };


  BranchPointGroups *BPG; // access to breakpoint groups
  ReadsManipulator *reads;
  std::vector<consensus_pair> consensus_pairs;


  void buildConsensusPairs();
  // Function fills the consensus_pairs vector with 
  // consensus pairs generated from the breakpoint blocks
  // of BPG

  void countSNVs();

  std::string generateParseString(mutation_classes &m);
  // Function generates a string of the following string (mutation string)
  // [SNV:a;b;c][SSV:e;f;g][LSV:x;y;z]
  // this string is stored as the header for each
  // aligned healthy read. This allows identification
  // of the mutation indexes directly from the SAM file

  void constructSNVFastqData();
  // generates the fastq file of non_mutated reads to align
  // to the refrence genome. 
  // a fastq element has format 

  void callBWA();
  // Function first calls bwa aln:
  // cmd ./bwa/bwa aln -t 16 hg19.fa cns_pairs.fastq > cns_pairs.sai
  // Then cals to bwa samse
  // ./bwa/bwa samse hg19.fa cna_pairs.sai cns_pair.fastq > cns_pairs.sam


  void call_SNV_variants(std::vector<snv_aln_info> &alignments, std::string filename);
  void printAllAlignments(std::vector<snv_aln_info> &alignments);
  void printSingleAlignment(snv_aln_info &snv);
  std::string reverseComplementString(std::string s);
  void correctReverseCompSNV(std::vector<snv_aln_info> &alignments);
  static bool compareSNVLocations(const single_snv &a, const single_snv &b);
  void outputSNVToUser(std::vector<snv_aln_info> &alignments, std::string report_filename);


public:
    GenomeMapper(BranchPointGroups &bpgroups, ReadsManipulator &reads,
        std::string outfile);


    void printConsensusPairs();
    // print out each mutated and non mutated string
    void printMutation(char healthy, char cancer, std::ofstream &mut_file);
};
#endif
