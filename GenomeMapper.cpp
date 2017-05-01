// GenomeMapper.cpp

#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <boost/regex.hpp>

#include "GenomeMapper.h"
#include "BranchPointGroups.h"
#include "Reads.h"
#include "util_funcs.h"
#include "string.h"
#include "SamEntry.h"
#include "benchmark.h"

using namespace std;

static const int MUT_CODE = 0;
static const int FLAG = 1;
static const int CHR = 2;
static const int AL_INDEX = 3;
static const int AL_CNS = 9;
static const int MISMATCHES = 18;
static const int SNV = 1;
static const int SSV = 2;
static const int LSV = 3;
static const int MUT_CNS = 1;

static const int MAX_LOW_CONF_POSITIONS = 10;


static const int REVERSE_FLAG = 16;
static const int FORWARD_FLAG = 0;


GenomeMapper::GenomeMapper(BranchPointGroups &bpgroups, ReadsManipulator &reads,
    string outfile){

  this->reads = &reads;
  this->BPG = &bpgroups;


  cout << "Building consensus pairs" << endl;
  START(buildcns);
  buildConsensusPairs();
  END(buildcns);
  TIME(buildcns);
  PRINT(buildcns);
  cout << "Writing fastq" << endl;
  constructSNVFastqData();
  callBowtie2();

  vector<SamEntry> alignments;
  cout << "Parsing sam" << endl;
  parseSamFile(alignments, "/data/ic711/result/cns_pairs.sam");
  cout << "Identifying SNV" << endl;
  identifySNVs(alignments);

//  printAlignmentStructs(alignments);

  outputSNVToUser(alignments, outfile);
}

void GenomeMapper::callBowtie2() {
  cout << "Calling Bowtie2" << endl;
  string command_aln = 
    "./bowtie2-2.3.1/bowtie2 -x /data/ic711/insilico_data/bowtie_index/hg19 -U /data/ic711/result/cns_pairs.fastq -S /data/ic711/result/cns_pairs.sam";
  system(command_aln.c_str());
}

void GenomeMapper::callBWA() {
  cout << "Calling bwa..." << endl;

  string command_aln = 
    "./bwa/bwa aln -t 16 /data/ic711/insilico_data/smufin_provided_data/ref_genome/hg19.fa /data/ic711/result/cns_pairs.fastq > /data/ic711/result/cns_pairs.sai";

  string command_sampe = 
    "./bwa/bwa samse /data/ic711/insilico_data/smufin_provided_data/ref_genome/hg19.fa /data/ic711/result/cns_pairs.sai /data/ic711/result/cns_pairs.fastq > /data/ic711/result/cns_pairs.sam";

  system(command_aln.c_str());  
  system(command_sampe.c_str());

  cout << "Finished bwa call." << endl;
}

void GenomeMapper::buildConsensusPairs() {
  // generate consensus pair for each breakpoint block
  // and then add starting gaps to align sequence pair

  consensus_pairs.reserve(BPG->getSize()); // make room
  int continued{0};
  for (int i=0; i < BPG->getSize(); ++i) {
    consensus_pair pair;
    pair.left_ohang = pair.right_ohang = 0;   // default to no overhang

    bool skip_mutated{false}, skip_non_mutated{false};
    skip_mutated = BPG->generateConsensusSequence(i, pair.mut_offset, TUMOUR, pair.pair_id, pair.mutated, pair.mqual);
    skip_non_mutated = BPG->generateConsensusSequence(i, pair.nmut_offset, HEALTHY, pair.pair_id, pair.non_mutated, pair.nqual);

    // discard sequences that do not contain both a non-mutated
    // and mutated cns pair
    if(skip_mutated == true  || skip_non_mutated == true) {
      continued++;
      continue;
    }

    //cout << "Cancer: " << endl;
    //cout << pair.mutated << endl
    //     << pair.mqual << endl;
    //cout << "Healthy: " << endl;
    //cout << pair.non_mutated << endl
    //     << pair.nqual << endl;
    trimCancerConsensus(pair);                // trim extra cancer sequence

    // TURN LQF OFF MASK
    bool low_quality_block = false;
    maskLowQualityPositions(pair, low_quality_block);
    if (low_quality_block) {
      continue;
    }
    consensus_pairs.push_back(pair);
  }
  cout << "Skipped " << continued << endl;
  consensus_pairs.shrink_to_fit();
}

void GenomeMapper::maskLowQualityPositions(consensus_pair & pair, bool &low_quality) {
  int num_low_quality_positions{0};
  for (int pos=0; pos < pair.mutated.size(); pos++) {
    if (pair.mqual[pos] != '-' || pair.nqual[pos + pair.left_ohang] != '-') {
        pair.mutated[pos] = pair.non_mutated[pos + pair.left_ohang];
        num_low_quality_positions++;
    }
  }

  if (num_low_quality_positions > MAX_LOW_CONF_POSITIONS) low_quality = true;
}

//void GenomeMapper::buildConsensusPairs() {
//  // generate consensus pair for each breakpoint block
//  // and then add starting gaps to align sequence pair
//
//  consensus_pairs.reserve(BPG->getSize()); // make room
//  int continued{0};
//  for (int i=0; i < BPG->getSize(); ++i) {
//
//    vector< vector<int> > tumour_base_frequency, healthy_base_frequency;
//    consensus_pair pair;
//    pair.left_ohang = pair.right_ohang = 0;   // default to no overhang
//
//    pair.mutated = BPG->generateConsensusSequence(i,
//      pair.mut_offset, TUMOUR, pair.pair_id, tumour_base_frequency);
//
//    pair.non_mutated = BPG->generateConsensusSequence(i,
//        pair.nmut_offset, HEALTHY, pair.pair_id, healthy_base_frequency);
//
//    // discard sequences that do not contain both a non-mutated
//    // and mutated cns pair
//    if(pair.mutated == "\0" || pair.non_mutated == "\0") {
//      continued++;
//      continue;
//    }
//
//    trimCancerConsensus(pair);                // trim extra cancer sequence
//
//    //// TURN OFF MASK

//    bool low_quality_block = false;

//    maskLowConfidencePositions(pair, healthy_base_frequency, 
//        tumour_base_frequency, low_quality_block);
//
//    if (low_quality_block) {
//      continue;
//    }
//    consensus_pairs.push_back(pair);
//  }
//  cout << "Skipped " << continued << endl;
//  consensus_pairs.shrink_to_fit();
//}

//void GenomeMapper::maskLowConfidencePositions(consensus_pair &pair,
//                                vector< vector<int> > &healthy_base_freq,
//                                vector< vector<int> > &tumour_base_freq,
//                                bool &discard) {
//  discard = false;
//  unsigned int start_h= 0, start_t= 0;
//
//  // start points do not take into account the updated
//  // matrix - which should have been cleaved
//  for(int i=0; i < healthy_base_freq[0].size(); i++) {
//    if(healthy_base_freq[0][i] != -1) {
//      start_h = i; // + left_ohang
//      break;    
//    }
//  }
//  for(int i=0; i < tumour_base_freq[0].size(); i++) {
//    if(tumour_base_freq[0][i] != -1) {
//      start_t = i; // should be corrected based on -1 updated matrix
//      break;
//    }
//  }
//
//  // mask based on tumour cns
//  for(int pos = start_t; pos < pair.mutated.size() + start_t; pos++) {
//    int n_tumour_bases_above_err_freq = 0;
//
//    if(tumour_base_freq[0][pos] == -1) {    // is this squared to the final base
//      break;
//    }
//
//    // get total reads
//    double total_reads = 0;
//    for(int base=0; base < 4; base++) {
//      total_reads += tumour_base_freq[base][pos];
//    }
//
//    // calc number of bases over the error frequency
//    for(int base=0; base < 4; base++) {
//      if(tumour_base_freq[base][pos] / total_reads > ALLELIC_FREQ_OF_ERROR) {
//        n_tumour_bases_above_err_freq++;
//      }
//    }
//
//    // if the number of bases with a high allelic frequency is above
//    // one, then the position is of low condifence, so mask
//    if (n_tumour_bases_above_err_freq > 1) {
//      pair.mutated[pos - start_t] = pair.non_mutated[pos - start_t + pair.left_ohang];
//    }
//  }
//
//  int number_of_low_conf_positions = 0; 
//  // mask based on healthy cns
//  for(int pos = start_h + pair.left_ohang; pos < pair.mutated.size() +
//      start_h + pair.left_ohang; pos++) {
//
//    int n_healthy_bases_above_err_freq = 0;
//    if(healthy_base_freq[0][pos] == -1) {
//      break;
//    }
//
//    // get total reads
//    double total_reads = 0;
//    for(int base = 0; base < 4; base++) {
//      // am I not counting over the non-updated start position (should be 
//      // pos + left_ohang, so i'm counting the cancer positions)
//      total_reads += healthy_base_freq[base][pos];
//    }                                                 
//
//    // cals number of bases over the error frequency
//    for(int base=0; base < 4; base++) {
//      if(healthy_base_freq[base][pos] / total_reads > ALLELIC_FREQ_OF_ERROR) {
//        n_healthy_bases_above_err_freq++;
//      }
//    }
//
//    // if n bases with high allelic freq. is above one, then 
//    // position is low confidence so mask
//    if(n_healthy_bases_above_err_freq > 1) {
//      pair.mutated[pos - pair.left_ohang - start_h] = pair.non_mutated[pos -
//        start_h];
//      number_of_low_conf_positions++;
//    }
//  }
//
//
//  // if there were too many low condidence positions, 
//  // the block is low quality, so discard
//  if(number_of_low_conf_positions > MAX_LOW_CONF_POSITIONS) {
//    discard = true;
//  }
//  cout << "number of low conf: " << number_of_low_conf_positions << endl;
//}


void GenomeMapper::trimCancerConsensus(consensus_pair & pair) {
  // Trim portions of the cancer consensus sequence that are
  // longer than heathy. Leave healthy if longer.
  // left cleave
  if (pair.mut_offset > pair.nmut_offset) {
    pair.mutated.erase(0, pair.mut_offset - pair.nmut_offset);
    pair.mqual.erase(0, pair.mut_offset - pair.nmut_offset);
  }
  else if (pair.mut_offset < pair.nmut_offset) {
    pair.left_ohang = pair.nmut_offset - pair.mut_offset;
  }

  // right cleave
  if (pair.mutated.size() > (pair.non_mutated.size() - pair.left_ohang)) {
    int dist = pair.mutated.size() - (pair.non_mutated.size() - pair.left_ohang);
    pair.mutated.erase(pair.mutated.size() - dist, dist);
    pair.mqual.erase(pair.mqual.size() - dist, dist);
  }
  else if (pair.mutated.size() < (pair.non_mutated.size() - pair.left_ohang)) {
    pair.right_ohang = (pair.non_mutated.size() - pair.left_ohang) - pair.mutated.size();
  }
  // else, equal length. Do nothing
}


void GenomeMapper::printMutation(char healthy, char cancer, ofstream &mut_file) {
  mut_file << healthy <<  ", " << cancer << endl;
}

void GenomeMapper::printAlignmentStructs(vector<SamEntry> &alignments) {
  for (SamEntry &entry: alignments) {
    cout << "Cancer seq:" << endl;
    if(entry.get<int>(SamEntry::FLAG) == FORWARD_FLAG) {
      printGaps(entry.get<int>(SamEntry::LEFT_OHANG));
    }
    else {
      printGaps(entry.get<int>(SamEntry::RIGHT_OHANG));
    }

    cout << entry.get<string>(SamEntry::HDR) << endl;
    cout << "Healthy seq:" << endl;
    cout << entry.get<string>(SamEntry::SEQ) << endl;
    cout << "Left ohang:  " << entry.get<int>(SamEntry::LEFT_OHANG) << endl;
    cout << "Right ohang: " << entry.get<int>(SamEntry::RIGHT_OHANG) << endl;
    cout << "Pair id: " << entry.get<int>(SamEntry::BLOCK_ID) << endl;
    for(int i = 0; i < entry.snvLocSize(); i++) {
      cout << entry.snvLocation(i) << ", ";
    }
    cout << endl << endl;
  }
}

void GenomeMapper::printGaps(int gaps) {
  for(; gaps > 0; gaps--) {
    cout << "-";
  }
}


void GenomeMapper::printConsensusPairs() {
  for(consensus_pair cns_pair : consensus_pairs) {
    cout << "Mutation Sequence:" << endl;
    cout << cns_pair.mutated << endl;
    cout << "Healthy Sequence:" << endl;
    cout << cns_pair.non_mutated << endl;

    cout << "SNV locations" << endl;
    for(int pos : cns_pair.mutations.SNV_pos) {
      cout << pos << ", ";
    }
    cout << endl << endl;
  }
}


void GenomeMapper::constructSNVFastqData() {
  ofstream snv_fq;
  snv_fq.open("/data/ic711/result/cns_pairs.fastq");

  for (consensus_pair &cns_pair : consensus_pairs) {
//    if(cns_pair.mutations.SNV_pos.size() == 0) {    // REMOVE IF on 14/11/16
//      continue;
//    }
    
    // otherwise, write healthy consensus as a fastq
    // with the SNVs stored in the header parse string
    
    snv_fq << "@" + cns_pair.mutated + 
    "[" + to_string(cns_pair.left_ohang) + ";" + 
    to_string(cns_pair.right_ohang) + ";" + 
    to_string(cns_pair.pair_id) + "]" << endl;

    snv_fq << cns_pair.non_mutated << endl;
    snv_fq << "+" << endl;
    string qual(cns_pair.non_mutated.size(), '!'); // set quality to highest, as dummy
    snv_fq << qual << endl;
  }
}



void GenomeMapper::parseSamFile(vector<SamEntry> &alignments, string filename) {
  ifstream snv_sam(filename);	// open alignment file
  
  boost::regex rgx_header("(@).*");
  string line;
  while(getline(snv_sam, line)) {

    if (boost::regex_match(line, rgx_header)) {	// skip past headers
      continue;
    }
   
    SamEntry entry(line);     // parse the entry into a SamEntry
    if (entry.get<string>(SamEntry::RNAME) != "22") {
      continue;
    }
    alignments.push_back(entry);
  }
}

void GenomeMapper::printAllAlignments(vector<SamEntry> &alignments){
  for(SamEntry &entry: alignments) {
    cout << "FLAG  :" << entry.get<int>(SamEntry::FLAG) << endl;
    cout << "CHR  :" << entry.get<string>(SamEntry::RNAME) << endl;
    cout << "POS :" << entry.get<int>(SamEntry::POS) << endl;
    cout << "HEALTHY: " << entry.get<string>(SamEntry::SEQ) << endl;
    cout << "TUMOUR : " << entry.get<string>(SamEntry::HDR) << endl;
    for(int i = 0; i < entry.snvLocSize(); i++) {
      cout << entry.snvLocation(i)  << ", ";
    }
    cout << endl << endl;
  }
}

void GenomeMapper::printSingleAlignment(SamEntry &entry) {
  cout << "FLAG  :" << entry.get<int>(SamEntry::FLAG) << endl;
  cout << "CHR  :" << entry.get<string>(SamEntry::RNAME) << endl;
  cout << "POS :" << entry.get<int>(SamEntry::POS) << endl;
  cout << "HELATHY :" << entry.get<string>(SamEntry::SEQ) << endl;
  cout << "TUMOUR :" << entry.get<string>(SamEntry::HDR) << endl;
  for(int i = 0; i < entry.snvLocSize(); i++) {
    cout << entry.snvLocation(i)  << ", ";
  }
  cout << endl << endl;
}

string GenomeMapper::reverseComplementString(string s){
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

void GenomeMapper::correctReverseCompSNV(vector<SamEntry> &alignments) {

  for(SamEntry &entry : alignments) {

    if(entry.get<int>(SamEntry::FLAG) == FORWARD_FLAG) {
      for(int i = 0; i < entry.snvLocSize(); i++) {
        entry.setSNVLocation(i, entry.snvLocation(i) - 1);
      }
    }
    else if(entry.get<int>(SamEntry::FLAG) == REVERSE_FLAG) { // convert indecies to rev comp and rev comp cns
      entry.get<string>(SamEntry::HDR) = reverseComplementString(entry.get<string>(SamEntry::HDR));
      for(int i = 0; i < entry.snvLocSize(); i++) {
        entry.setSNVLocation(i, entry.get<string>(SamEntry::HDR).size() - entry.snvLocation(i));
      }
    }
  }
}

void GenomeMapper::identifySNVs(vector<SamEntry> &alignments) {
  for (SamEntry & entry : alignments) {
      if(entry.get<int>(SamEntry::FLAG) == FORWARD_FLAG) {
        countSNVs(entry, entry.get<int>(SamEntry::LEFT_OHANG));
      }
      else if (entry.get<int>(SamEntry::FLAG) == REVERSE_FLAG) {
        entry.set(SamEntry::HDR, reverseComplementString(entry.get<string>(SamEntry::HDR))); 
        countSNVs(entry, entry.get<int>(SamEntry::RIGHT_OHANG)); // invert overhangs due to rev comp
      }
  }
}
void GenomeMapper::countSNVs(SamEntry &alignment, int ohang) {
  string mutated = alignment.get<string>(SamEntry::HDR);
  string non_mutated = alignment.get<string>(SamEntry::SEQ);
  // indel signature
  for(int i=0; i < mutated.size() - 1; i++) {
    if (mutated[i] != non_mutated[i + ohang] &&
        mutated[i+1] != non_mutated[i+1 + ohang]) {
      return;
    }
  }

  // SNV at start
  if (mutated[0] != non_mutated[0 + ohang] &&
      mutated[1] == non_mutated[1 + ohang]) {
      alignment.snv_push_back(0); 
  }

  // SNV at end 
  int cns_len = mutated.size();
  if (
      mutated[cns_len-1] != 
      non_mutated[cns_len-1 + ohang] &&
      mutated[cns_len-2] == 
      non_mutated[cns_len-2 + ohang]
      ) {
      alignment.snv_push_back(cns_len-1);
  }

  // SNV in body

  for (int i=1; i < mutated.size() - 1; i++) {
    if (mutated[i-1] == non_mutated[i-1 + ohang] &&
        mutated[i] != non_mutated[i + ohang] &&
        mutated[i+1] == non_mutated[i+1 + ohang] ) {
      alignment.snv_push_back(i);
    }
  }
}



bool GenomeMapper::compareSNVLocations(const single_snv &a, const single_snv &b) {
  return a.position < b.position;
}


void GenomeMapper::outputSNVToUser(vector<SamEntry> &alignments, string report_filename) {


  // load each snv into a separate struct, so each can be easily sorted

  vector<single_snv> separate_snvs;
  for(SamEntry & entry : alignments) {
    for(int i=0; i < entry.snvLocSize(); i++) {
      int snv_index = entry.snvLocation(i);
      int overhang = 0;
      if (entry.get<int>(SamEntry::FLAG) == FORWARD_FLAG) {
        overhang = entry.get<int>(SamEntry::LEFT_OHANG);
      }
      else if (entry.get<int>(SamEntry::FLAG) == REVERSE_FLAG) {
        overhang = entry.get<int>(SamEntry::RIGHT_OHANG);
      }

      single_snv snv;
      snv.chr = entry.get<string>(SamEntry::RNAME);
      snv.position = (entry.get<int>(SamEntry::POS) + snv_index + overhang); // location of snv
      snv.healthy_base = entry.get<string>(SamEntry::SEQ)[snv_index + overhang];
      snv.mutation_base = entry.get<string>(SamEntry::HDR)[snv_index];
      snv.pair_id = entry.get<int>(SamEntry::BLOCK_ID);

      separate_snvs.push_back(snv);
    }
  }

  // sort the snvs 
  std::sort(separate_snvs.begin(), separate_snvs.end(), compareSNVLocations);
  
  ofstream report(report_filename);
  report << "Mut_ID\tType\tChr\tPos\tNormal_NT\tTumor_NT" << endl;

  unsigned int i=0;
  for(single_snv &snv : separate_snvs) {
    report << i << "\t" << "SNV\t" << snv.chr << "\t"
           << snv.position << "\t"
           << snv.healthy_base << "\t" 
           << snv.mutation_base << "\t"
           << snv.pair_id
           << endl;
    i++;
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

//void GenomeMapper::countSNVs() {
//
//
//
//  for(consensus_pair &p : consensus_pairs) {
//    bool discard = false; 
//    // remove any consensus pairs that
//    // do not contain SNV style sequence differences
//    // these could result from indels
//
//    for(int i=0; i < p.mutated.size() - 1; i++){
//      if (p.mutated[i] != p.non_mutated[i] && 
//	  p.mutated[i+1] != p.non_mutated[i+1]) {
//          discard = true;
//      }
//    }
//    
//    if(discard) {
//      continue;
//    }
//
//
//    // identify SNV occuring at the start of the string
//    if (p.mutated[0] != p.non_mutated[0] && 
//        p.mutated[1] == p.non_mutated[1]) {
//
//      p.mutations.SNV_pos.push_back(0);
//    }
//
//    // identify SNVs occuring within the string
//    for(int i=1; i < p.mutated.size() - 1; i++) {
//      if(p.mutated[i] != p.non_mutated[i] &&        // char diffrent but
//          p.mutated[i+1] == p.non_mutated[i+1] &&   // ...next char same
//          p.mutated[i-1] == p.non_mutated[i-1]) {   // ...prev char same
//
//        // then count as SNV
//        p.mutations.SNV_pos.push_back(i);   // store index of variants
//      }
//    }
//
//    // identify SNV occuring at the very end of the string
//    if (p.mutated[p.mutated.size()-1] !=
//        p.non_mutated[p.non_mutated.size()-1] &&
//        p.mutated[p.mutated.size()-2] == 
//        p.non_mutated[p.non_mutated.size()-2]){
//      
//      p.mutations.SNV_pos.push_back(p.mutated.size()-1);
//    }
//  }
//}
