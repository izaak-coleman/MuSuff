#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <algorithm>

#include <boost/regex.hpp>
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
static const int MUT_CNS = 4;


static const int REVERSE_FLAG = 16;
static const int FORWARD_FLAG = 0;



struct snv_aln_info {
 vector<int> SNV_pos;
 int flag;
 int chr;
 int position;
 string reference_mismatch;
 string non_mutated_cns;
 string mutated_cns;
};

struct single_snv {
 int flag;
 int chr;
 int position;
 char mutation_base;
 char healthy_base;
};

void call_SNV_variants(vector<snv_aln_info> &alignments, string filename) {
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

int main(int argc, char **argv) {
  START(sp);

  if(argc != 3) {
    cout << "usage: <file>.sam <report_filename> " << endl;
    exit(1);
  }
  string sam_filename = argv[1];
  string report_filename = argv[2];



  vector<snv_aln_info> alignments;
  call_SNV_variants(alignments, sam_filename);
  correctReverseCompSNV(alignments);
  outputSNVToUser(alignments, report_filename);
  printAllAlignments(alignments);
  END(sp);
  TIME(sp);
  PRINT(sp);
  return 0;
}
