#include <iostream>
#include <vector>
#include <string>

#include <sys/time.h>
#include <sys/resource.h>

#include "SuffixArray.h"
#include "BranchPointGroups.h"
#include "helper_functions.h"
#include "Reads.h"
#include "GenomeMapper.h"
#include "benchmark.h"


using namespace std;

int main(int argc, char **argv) {

  struct rusage main_res;
  struct rusage reads_res;
  struct rusage suf_res;
  struct rusage bpg_res;
  struct rusage gm_res;

  getrusage(RUSAGE_SELF, &main_res);
  getrusage(RUSAGE_SELF, &reads_res);
  getrusage(RUSAGE_SELF, &suf_res);
  getrusage(RUSAGE_SELF, &bpg_res);
  getrusage(RUSAGE_SELF, &gm_res);

  START(main);

  command_line_params params;
  parseCommandLine(params, argc, argv);

  // Construct the reads buffer (stores FASTA strings)

  START(reads_manip);
  ReadsManipulator reads(params);
  END(reads_manip);
  TIME(reads_manip);
  PRINT(reads_manip);
  getrusage(RUSAGE_SELF, &reads_res);
  cerr << "rss, reads, " << reads_res.ru_maxrss << endl;


  //SA construction called in constructor
  START(suf_array);
  SuffixArray SA(reads, params.minimum_suffix_size);
  END(suf_array);
  TIME(suf_array);
  PRINT(suf_array);
  getrusage(RUSAGE_SELF, &suf_res);
  cerr << "rss, suf, " << suf_res.ru_maxrss << endl;

  //reads.printReads();
  //SA.printSuffixes();
  // BG construction called in constructor
  START(bpg);
  BranchPointGroups BG(SA, reads, params.econt);
  END(bpg);
  TIME(bpg);
  PRINT(bpg);
  getrusage(RUSAGE_SELF, &bpg_res);
  cerr << "rss, bpg, " << bpg_res.ru_maxrss << endl;

  START(gm);
  GenomeMapper mapper(BG, reads);
  END(gm);
  TIME(gm);
  PRINT(gm);
  getrusage(RUSAGE_SELF, &gm_res);
  cerr <<"rss, gm, " << gm_res.ru_maxrss << endl;
 // mapper.printConsensusPairs();



 END(main);
 TIME(main);
 PRINT(main);
 getrusage(RUSAGE_SELF, &main_res);
 cerr << "rss, main, " << main_res.ru_maxrss << endl;
 


  return 0;

}

