#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <boost/regex.hpp>

#include "GenomeMapper.h"
#include "BranchPointGroups.h"
#include "Reads.h"
#include "helper_functions.h"

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
static const int MUT_CNS = 4;


static const int REVERSE_FLAG = 16;
static const int FORWARD_FLAG = 0;

GenomeMapper::GenomeMapper(BranchPointGroups &bpgroups, ReadsManipulator &reads){

  this->reads = &reads;
  this->BPG = &bpgroups;


  buildConsensusPairs();

  countSNVs();
  constructSNVFastqData();
  callBWA();
  callSamParser();

  vector<snv_aln_info> alignments;
  call_SNV_variants(alignments, "cns_pairs.sam");
  correctReverseCompSNV(alignments);
  outputSNVToUser(alignments, "reported_SNV.txt");
}

void GenomeMapper::callBWA() {
  cout << "Calling bwa..." << endl;

  string command_aln = 
    "./bwa/bwa aln -t 16 ./dataset/ref_genome/hg19.fa cns_pairs.fastq > cns_pairs.sai";

  string command_sampe = 
    "./bwa/bwa samse ./dataset/ref_genome/hg19.fa cns_pairs.sai cns_pairs.fastq > cns_pairs.sam";

  system(command_aln.c_str());  
  system(command_sampe.c_str());

  cout << "Finished bwa call." << endl;
}

void GenomeMapper::callSamParser() {

  cout << "Begining sam_parser call..." << endl;

  string command_sam_parser = 
    "./sam_parser cns_pairs.sam mutations.snv.txt";

  system(command_sam_parser.c_str());
  cout << "Finishing sam_parser call." << endl;
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

// samparser content

void GenomeMapper::call_SNV_variants(vector<snv_aln_info> &alignments, string filename) {
  ifstream snv_sam(filename);	// open alignment file
  
  boost::regex header("(@).*");		// matches header

  boost::regex mut_code_parser("^\\[(.*)\\]\\[(.*)\\]\\[(.*)\\]\\[(.*)\\]$");
  boost::regex position_parser("[A-Z][A-Z][A-Z]:(.*)$");
  boost::regex md_field_parser("MD:Z:([^\t]+)");  // with bug

  string line;

  while(getline(snv_sam, line)) {
    if (boost::regex_match(line,header)) {	// skip past headers
      continue;
    }
   
  
    vector<string> fields;
    string token("");
    while(token != line) {
      token = line.substr(0, line.find_first_of("\t"));
      line = line.substr(line.find_first_of("\t") + 1);
      fields.push_back(token);
    }

    // check that algorithm aligned read to chromosome 22
    if (fields[CHR] != "22") {
      continue;
    }

    // load relevant fields into each align info struct that dont require parsing
    snv_aln_info al_info;
    al_info.flag = stoi(fields[FLAG]);
    al_info.chr = stoi(fields[CHR]);
    al_info.position = stoi(fields[AL_INDEX]);
    al_info.non_mutated_cns = fields[AL_CNS];



    // parse the mutation code string into the four fields
    boost::smatch mut_code_fields;
    boost::regex_match(fields[MUT_CODE], mut_code_fields, mut_code_parser);

    // load the mutated cns sequence
    al_info.mutated_cns = mut_code_fields[MUT_CNS];

    // parse the SNV string to positions
    boost::smatch snv_string; 
    string snv_field = mut_code_fields[SNV];
    boost::regex_match(snv_field, snv_string, position_parser);
  
    string snv_positions = snv_string[1];
    string pos("");
    while(pos !=  snv_positions) {
      pos = snv_positions.substr(0, snv_positions.find_first_of(";"));
      snv_positions = snv_positions.substr(snv_positions.find_first_of(";") + 1);
      al_info.SNV_pos.push_back(stoi(pos));
    }
    alignments.push_back(al_info);
  }

}

void GenomeMapper::printAllAlignments(vector<snv_aln_info> &alignments){
  for(snv_aln_info snv: alignments) {
    cout << "FLAG  :" << snv.flag << endl;
    cout << "CHR  :" << snv.chr << endl;
    cout << "POS :" << snv.position << endl;
    cout << "HELATHY: " << snv.non_mutated_cns << endl;
    cout << "TUMOUR : " << snv.mutated_cns << endl;
    for(int m : snv.SNV_pos) {
      cout << m << ", ";
    }
    cout << endl << endl;
  } 
}

void GenomeMapper::printSingleAlignment(snv_aln_info &snv) {
  cout << "FLAG  :" << snv.flag << endl;
  cout << "CHR  :" << snv.chr << endl;
  cout << "POS :" << snv.position << endl;
  cout << "HELATHY :" << snv.non_mutated_cns << endl;
  cout << "TUMOUR :" << snv.mutated_cns << endl;
  for(int m : snv.SNV_pos) {
    cout << m << ", ";
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

void GenomeMapper::correctReverseCompSNV(vector<snv_aln_info> &alignments) {

  for(snv_aln_info &al : alignments) {

    if(al.flag == FORWARD_FLAG) {
      for(int &snv : al.SNV_pos) {
        snv--;			// return to 0 index for addition to aln val
      }
    }
    else if(al.flag == REVERSE_FLAG) { // convert indecies to rev comp and rev comp cns
      al.mutated_cns = reverseComplementString(al.mutated_cns);
      for(int &snv : al.SNV_pos) {
        snv = (al.mutated_cns.size() - snv);
      } 
    }
  }
}

bool GenomeMapper::compareSNVLocations(const single_snv &a, const single_snv &b) {
  return a.position < b.position;
}


void GenomeMapper::outputSNVToUser(vector<snv_aln_info> &alignments, string report_filename) {


  // load each snv into a separate struct, so each can be easily sorted

  vector<single_snv> separate_snvs;
  for(snv_aln_info &al : alignments) {
    for(int snv_index : al.SNV_pos) {

      single_snv snv;
      snv.chr = al.chr;
      snv.position = (al.position + snv_index); // location of snv
      snv.healthy_base = al.non_mutated_cns[snv_index];
      snv.mutation_base = al.mutated_cns[snv_index];

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
           << endl;
    i++;
  }
}
