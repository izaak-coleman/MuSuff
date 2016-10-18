#ifndef GENOMEMAPPER_H
#define GENOMEMAPPER_H

#include <string>
#include <vector>

#include "BranchPointGroups.h"
#include "Reads.h"

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






public:
    GenomeMapper(BranchPointGroups &bpgroups, ReadsManipulator &reads);


    void printConsensusPairs();
    // print out each mutated and non mutated string
};
#endif
