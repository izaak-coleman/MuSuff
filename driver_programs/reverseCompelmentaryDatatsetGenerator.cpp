#include <iostream>

#include <fstream>
#include <string>

#include <zlib.h>
#include "kseq.h"


KSEQ_INIT(gzFile, gzread);
using namespace std;

string reverseComplementString(string s);


int main(int argc, char **argv) {

  gzFile data_file;
  data_file = gzopen(argv[1], "r");   // link to .fastq.gz

  kseq_t *seq = kseq_init(data_file);   


  int eof_check;

  string of_name = argv[1];
  of_name += ".rev";
  ofstream ofile;
  ofile.open(of_name);

  while ((eof_check = kseq_read(seq)) >=0) {

    // As we only want high quality reads, 
    // only reads with a quality score have potential to be added

    ofile << seq->name.s << endl; // print name

    string dna = seq->seq.s;      // print string in reverse complement
    ofile << reverseComplementString(dna) << endl;



    if (seq->comment.l){          // if comment print
      ofile << seq->comment.s << endl;
    }
    else {                        // else print "+"
      ofile << "+" << endl;
    }


    if(seq->qual.l) {             // print quality score
      ofile << seq->qual.s << endl;
    }
  }

  kseq_destroy(seq);
  gzclose(data_file);

  return 0;
}

string reverseComplementString(string s){
  string revcomp = "";

  for(int i = s.size()-1; i >= 0; i--) {
  // travel in reverse and switch for complementary
    switch(s[i]) {

      case 'A':{
        revcomp += "T";
        break;
       }

      case 'T':{
        revcomp += "A";
        break;
      }

      case 'C':{
        revcomp += "G";
        break;
      }

      case 'G':{
        revcomp += "C";
        break;
      }
    }
  }

  return revcomp;

}
