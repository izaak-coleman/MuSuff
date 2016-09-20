#include <zlib.h>
#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <cstring>

#include "kseq.h"

KSEQ_INIT(gzFile, gzread);
using namespace std;


struct fastq_t {
  string name, seq, qual;
};


void copyFromKSeq(fastq_t &to, kseq_t *from) {
  // copy id
  to.name = from->name.s;
  // copy sequence
  to.seq = from->seq.s;
  // copy quality string

  if(from->qual.l) {    // check that qual is present
    to.qual = from->seq.s;
  }
}



int main(int argc, char **argv){

  gzFile data_file;
  vector<fastq_t> fastq_reads;

  data_file = gzopen(argv[1], "r");      // link stream to GZ file
  kseq_t *seq = kseq_init(data_file);// initialize memory alloc

  int eof_check;
  fastq_t next_read;
  while ((eof_check = kseq_read(seq)) >= 0) {

    copyFromKSeq(next_read, seq); 

    fastq_reads.push_back(next_read);
  }
  cout << "FROM VECTOR" << endl << endl;
  cout << fastq_reads.size() << endl;

  for (fastq_t seq : fastq_reads) {
    cout << "Sequence Header: " << seq.name << endl;

    cout << "Sequence: " << seq.seq << endl;

    if(seq.qual.size()) {
      cout << "Sequence Quality: " << seq.qual << endl;
    }

    cout << "Size of kseq_t: " << sizeof(seq) << endl;

  }


  kseq_destroy(seq);
  gzclose(data_file);

  return 0;

}
