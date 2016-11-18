#include <iostream>
#include <string>
#include <vector>

#include "VariantCall.h"
#include "xa_as_lines.h"

using namespace std;

int main(int argc, char **argv) {
  if( argc != 3) {
    cout << "usage: <exe> <samfile> <outfilename>" << endl;
  }
  string samfile = argv[1];
  string outfile = argv[2];


  vector<snv_aln_info> alignments;
  string sam_alt_alns;
  xaToSeparateLine(samfile, sam_alt_alns);
  parseSamFile(alignments, sam_alt_alns);
  identifySNVs(alignments);
  outputSNVToUser(alignments, outfile);

  return 0;
}
