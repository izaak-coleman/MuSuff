#include <iostream>
#include <vector>
#include <string>

#include <sys/time.h>
#include <sys/resource.h>

#include "SuffixArray.h"
#include "BranchPointGroups.h"
#include "util_funcs.h"
#include "Reads.h"
#include "GenomeMapper.h"
#include "benchmark.h"


using namespace std;

int main(int argc, char **argv) {

  // Construct the reads buffer (stores FASTA strings)
  ReadsManipulator reads(argc, argv);

  //SA construction called in constructor
  SuffixArray SA(reads, reads.getMinSuffixSize());

  // BG construction called in constructor
  BranchPointGroups BG(SA, reads, reads.getEcont());

  GenomeMapper mapper(BG, reads, reads.outFile());
  return 0;

}

