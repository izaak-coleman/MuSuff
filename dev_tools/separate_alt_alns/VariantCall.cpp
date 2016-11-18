// GenomeMapper.cpp
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <boost/regex.hpp>

#include "VariantCall.h"
#include "xa_as_lines.h"
#include "string.h"

using namespace std;

static const int MUT_CODE = 0;
static const int AL_INDEX = 3;
static const int AL_CNS = 9;
static const int MISMATCHES = 18;
static const int SNV = 1;
static const int SSV = 2;
static const int LSV = 3;
static const int MUT_CNS = 4;


static const int REVERSE_FLAG = 16;
static const int FORWARD_FLAG = 0;



std::vector<consensus_pair> consensus_pairs;

void callBWA() {
  cout << "Calling bwa..." << endl;

  string command_aln = 
    "./bwa/bwa aln -t 16 /data/insilico_data/smufin_provided_data/ref_genome/hg19.fa cns_pairs.fastq > cns_pairs.sai";

  string command_sampe = 
    "./bwa/bwa samse /data/insilico_data/smufin_provided_data/ref_genome/hg19.fa cns_pairs.sai cns_pairs.fastq > cns_pairs.sam";

  system(command_aln.c_str());  
  system(command_sampe.c_str());

  cout << "Finished bwa call." << endl;
}


void printMutation(char healthy, char cancer, ofstream &mut_file) {
  mut_file << healthy <<  ", " << cancer << endl;
}


void printConsensusPairs() {
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


void constructSNVFastqData() {
  ofstream snv_fq;
  snv_fq.open("cns_pairs.fastq");

  for (consensus_pair &cns_pair : consensus_pairs) {
    if(cns_pair.mutations.SNV_pos.size() == 0) {    // REMOVE IF on 14/11/16
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

string generateParseString(mutation_classes &m) {
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

void parseSamFile(vector<snv_aln_info> &alignments, string filename) {
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

    // DEV RULE: clearing al_info.SNV_pos to build from scrach using
    // post snv identification design
    al_info.SNV_pos.clear(); // should be empty anyways as no call to countSNVs()


    alignments.push_back(al_info);
  }

}

void printAllAlignments(vector<snv_aln_info> &alignments){
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

void printSingleAlignment(snv_aln_info &snv) {
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


void correctReverseCompSNV(vector<snv_aln_info> &alignments) {

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

void identifySNVs(vector<snv_aln_info> &alignments) {
  for (snv_aln_info &a : alignments) {
      if(a.flag == FORWARD_FLAG) {
        countSNVs(a);
      }
      else if (a.flag == REVERSE_FLAG) {
        a.mutated_cns = reverseComplementString(a.mutated_cns); 
        countSNVs(a);
      }
  }
}
void countSNVs(snv_aln_info &alignment) {
  
  // indel signature
  for(int i=0; i < alignment.mutated_cns.size() - 1; i++) {
    if (alignment.mutated_cns[i] != alignment.non_mutated_cns[i] &&
        alignment.mutated_cns[i+1] != alignment.non_mutated_cns[i+1]) {
      return; // no not count 
    }
  }

  // SNV at start
  if (alignment.mutated_cns[0] != alignment.non_mutated_cns[0] &&
      alignment.mutated_cns[1] == alignment.non_mutated_cns[1]) {
      alignment.SNV_pos.push_back(0);   
  }

  // SNV at end 
  int cns_len = alignment.mutated_cns.size();
  if (alignment.mutated_cns[cns_len-1] != alignment.non_mutated_cns[cns_len-1]
      && alignment.mutated_cns[cns_len-2] == alignment.non_mutated_cns[cns_len-2]) {
      alignment.SNV_pos.push_back(cns_len-1);
  }

  // SNV in body

  for (int i=1; i < alignment.mutated_cns.size() - 1; i++) {
    if (alignment.mutated_cns[i-1] == alignment.non_mutated_cns[i-1] &&
        alignment.mutated_cns[i] != alignment.non_mutated_cns[i] &&
        alignment.mutated_cns[i+1] == alignment.non_mutated_cns[i+1] ) {
      alignment.SNV_pos.push_back(i);
    }
  }
}



bool compareSNVLocations(const single_snv &a, const single_snv &b) {
  return a.position < b.position;
}


void outputSNVToUser(vector<snv_aln_info> &alignments, string report_filename) {


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
