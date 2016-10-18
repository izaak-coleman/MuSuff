// gsa_file_gen.cpp
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

#include "../Suffix_t.h"
#include "gsa_file_gen.h"
#include "string.h"   // split_string()

static const std::string EXT = ".gsa";
static const int READ_ID_IDX = 0;
static const int OFFSET_IDX  = 1;
static const int TYPE_IDX = 2;
static const int NUM_FIELDS = 3;


using namespace std;

void buildGSAFile(vector<Suffix_t> &GSA, string filename) {
  if (filename.substr(filename.size() - EXT.size()) != EXT){
    filename += EXT;
  }
  ofstream gsa_file;
  gsa_file.open(filename);

  for(int i = 0; i < GSA.size(); i++) {
    string next_suf = "";
    next_suf += (to_string(GSA[i].read_id) + ",");
    next_suf += (to_string(GSA[i].offset)  + ",");
    next_suf += to_string(GSA[i].type);
    gsa_file << next_suf << endl;
  }
}

void constructGSAFromFile(string filename, vector<Suffix_t> &GSA) {
  ifstream gsa_file;
  gsa_file.open(filename);
  string next_line;

  while (getline(gsa_file, next_line)) {
    vector<string> elements;
    split_string(next_line, ",", elements); // split.gsa line into suffix fields
    Suffix_t suf;

    if (elements.size() == NUM_FIELDS) {
      suf.read_id = stoi(elements[READ_ID_IDX]);    // load data
      suf.offset = stoi(elements[OFFSET_IDX]);
      suf.type = stoi(elements[TYPE_IDX]);
      GSA.push_back(suf);
    }
  }
}
