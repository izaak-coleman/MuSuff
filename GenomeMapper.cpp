#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>

#include "GenomeMapper.h"
#include "BranchPointGroups.h"
#include "Reads.h"
#include "helper_functions.h"

using namespace std;

GenomeMapper::GenomeMapper(BranchPointGroups &bpgroups, ReadsManipulator &reads){

  this->reads = &reads;
  this->BPG = &bpgroups;


  buildConsensusPairs();

  countSNVs();
  constructSNVFastqData();
}

void callBWA() {
  string command_aln = 
    "./bwa/bwa aln -t 16 hg19.fa cns_pairs.fastq > cns_pairs.sai";

  string command_sampe = 
    "./bwa/bwa sampe hg19.fa cns_pairs.sai cns_pairs.fastq > cns_pairs.sam";

  string command_sam_parser = 
    "./sam_parser cns_pairs.sam mutations.snv.txt";

  system(command_aln.c_str());
  system(command_sampe.c_str());
  system(command_sam_parser.c_str());
}


void GenomeMapper::buildConsensusPairs() {

  consensus_pairs.reserve(BPG->getSize()); // make room

  // generate consensus pair for each breakpoint block
  // and then add starting gaps to align sequence pair
  for (int i=0; i < BPG->getSize(); ++i) {

    // generate consensus sequences
    consensus_pair next_pair;
    next_pair.mutated = BPG->generateConsensusSequence(i,
      next_pair.mut_offset, TUMOUR, next_pair.read_freq_m);

    next_pair.non_mutated = BPG->generateConsensusSequence(i,
        next_pair.nmut_offset, HEALTHY, next_pair.read_freq_nm);

    // discard sequences that do not contain both a non-mutated
    // and mutated cns pair
    if(next_pair.mutated == "\0" || next_pair.non_mutated == "\0") {
      continue;
    }


    // cleave consensus sequences, so they
    // align (by cleaving start), and remove hanging tails (cleaving end)

    if(next_pair.mut_offset < next_pair.nmut_offset) {
      //next_pair.mutated.insert(0, BPG->addGaps(next_pair.nmut_offset -
      //                    next_pair.mut_offset));

      next_pair.non_mutated.erase(0, next_pair.nmut_offset -
          next_pair.mut_offset);
      next_pair.read_freq_nm.erase(0, next_pair.nmut_offset -
          next_pair.mut_offset);
    }
    else {
      //next_pair.non_mutated.insert(0, BPG->addGaps(next_pair.mut_offset -
      //                    next_pair.nmut_offset));
      next_pair.mutated.erase(0, next_pair.mut_offset - next_pair.nmut_offset);
      next_pair.read_freq_m.erase(0, next_pair.mut_offset - next_pair.nmut_offset);
    }

    // cleave tails
    if(next_pair.mutated.size() > next_pair.non_mutated.size()) {
      int pos_to_cleave = next_pair.mutated.size() -
        next_pair.non_mutated.size();

      next_pair.mutated.erase(next_pair.mutated.size()-1-pos_to_cleave,
          pos_to_cleave);
    }
    else {
      int pos_to_cleave = next_pair.non_mutated.size() -
        next_pair.mutated.size();

      next_pair.non_mutated.erase(next_pair.non_mutated.size()-1-pos_to_cleave,
          pos_to_cleave);
    }


    consensus_pairs.push_back(next_pair);
  }

}

void GenomeMapper::countSNVs() {

  for(consensus_pair &p : consensus_pairs) {
//    for( int i=0; i < p.mutated.size();i++) {
//      if(p.mutated[i] != p.non_mutated[i]) {
//        p.mutations.SNV_pos.push_back(i);
//      }  
//    }
    bool discard = false; 
    // remove any consensus pairs that
    // do not contain SNV style sequence differences
    // these could result from indels

    for(int i=0; i < p.mutated.size() - 1; i++){
      if (p.mutated[i] != p.non_mutated[i] && 
	  p.mutated[i+1] != p.non_mutated[i+1]) {
          discard = true;
      }
    }
    
    if(discard) {
      continue;
    }


    // identify SNV occuring at the start of the string
    if (p.mutated[0] != p.non_mutated[0] && 
        p.mutated[1] == p.non_mutated[1]) {

      p.mutations.SNV_pos.push_back(0);
    }

    // identify SNVs occuring within the string
    for(int i=1; i < p.mutated.size() - 1; i++) {
      if(p.mutated[i] != p.non_mutated[i] &&        // char diffrent but
          p.mutated[i+1] == p.non_mutated[i+1] &&   // ...next char same
          p.mutated[i-1] == p.non_mutated[i-1]) {   // ...prev char same

        // then count as SNV
        p.mutations.SNV_pos.push_back(i);   // store index of variants
      }
    }

    // identify SNV occuring at the very end of the string
    if (p.mutated[p.mutated.size()-1] !=
        p.non_mutated[p.non_mutated.size()-1] &&
        p.mutated[p.mutated.size()-2] == 
        p.non_mutated[p.non_mutated.size()-2]){
      
      p.mutations.SNV_pos.push_back(p.mutated.size()-1);
    }
  }
}


void GenomeMapper::printConsensusPairs() {
  for(consensus_pair cns_pair : consensus_pairs) {
    cout << "Mutation Sequence:" << endl;
    cout << cns_pair.mutated << endl;
    cout << "Healthy Sequence:" << endl;
    cout << cns_pair.non_mutated << endl;
    cout << "FreqMut String:" << endl;
    cout << cns_pair.read_freq_m << endl;
    cout << "FreqHeal String:" << endl;
    cout << cns_pair.read_freq_nm << endl;

    cout << "SNV locations" << endl;
    for(int pos : cns_pair.mutations.SNV_pos) {
      cout << pos << ", ";
    }
    cout << endl << endl;
  }
}


void GenomeMapper::constructSNVFastqData() {
  ofstream snv_fq;
  snv_fq.open("cns_pairs.fastq");

  for (consensus_pair &cns_pair : consensus_pairs) {
    if(cns_pair.mutations.SNV_pos.size() == 0) {
      continue;
    }
    
    // otherwise, write healthy consensus as a fastq
    // with the SNVs stored in the header parse string
    
    snv_fq << "@" + generateParseString(cns_pair.mutations) + "[" + 
    cns_pair.mutated + "]"
    << endl;
    snv_fq << cns_pair.non_mutated << endl;
    snv_fq << "+" << endl;
    string qual(cns_pair.non_mutated.size(), '~'); // set quality to highest, as dummy
    snv_fq << qual << endl;
  }
}

string GenomeMapper::generateParseString(mutation_classes &m) {
  // this codes the mutations found in a consensus sequence
  // the format is [SNV:a;b;c][SSV:e;f;g][LSV:x;y;z]
  // lower case characters represent the index values of
  // various mutations. 

  string mutation_string;
  mutation_string = "[SNV:";
  // code all the SNV positions into the string
  // at the same time, converting the positions 
  // from 0-based to 1-based

  for (int snv_index : m.SNV_pos) {
    mutation_string += to_string(snv_index+1) + ";";
  }
  mutation_string.pop_back();	// remove trailing ";"
  mutation_string += "][SSV:][LSV:]";
  return mutation_string;
}
