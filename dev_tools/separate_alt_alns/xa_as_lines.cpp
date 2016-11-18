#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <boost/regex.hpp>

#include "xa_as_lines.h"
#include "string.h"

using namespace std;

void xaToSeparateLine(string filename, string &out_filename) {
  ifstream samXA; ofstream samSepLines;
  samXA.open(filename.c_str());
  out_filename = filename.substr(0,filename.size()-SAM_EXT) +
    "_alt_aln.sam";
  samSepLines.open(out_filename.c_str());
  skipSamHeader(samXA);

  // begin moving XA to seperate line
  string sam_entry;
  while (getline(samXA, sam_entry)) {
    vector<string> entry_fields;
    split_string(sam_entry, "\t", entry_fields);

    if (xaPresent(entry_fields)) {
      vector<string> xa_colon_fields;
      vector<string> alt_alns;
      split_string(getXA(entry_fields), ":", xa_colon_fields);
      split_string(xa_colon_fields[2], ";", alt_alns);
      for(int i = 0; i < alt_alns.size(); i++) {
        samSepLines << constructAlternateSamEntry(entry_fields, alt_alns[i]) 
                    << endl;
      }
    }
    samSepLines << sam_entry << endl; // output entry as before
  }
}

string constructAlternateSamEntry(vector<string> const &entry_fields, string
    alignment) {

  vector<string> aln_fields;
  split_string(alignment, ",", aln_fields);   // get alignment fields

  string orientation;
  int alt_pos = stoi(aln_fields[1]);
  string query = entry_fields[SEQ];

  if (entry_fields[FLAG] == FWD && alt_pos >= 0) {
    orientation = FWD;
  } 

  else if (entry_fields[FLAG] == FWD && alt_pos < 0) {
    orientation = REV;   // alt aligned in reverse
    alt_pos *= -1;
    query = reverseComplementString(query);
  }

  else if (entry_fields[FLAG] == REV && alt_pos >= 0) {
    orientation = FWD;
    query = reverseComplementString(query);
  }
  else if (entry_fields[FLAG] == REV && alt_pos < 0) {
    orientation = REV;
    alt_pos *= -1;
  }

  string sam_entry = "";
  for(int i=0; i < entry_fields.size(); i++) {

    switch(i) {
      case FLAG:
        sam_entry += orientation + "\t";
        break;
      case POS:
        sam_entry += to_string(alt_pos) + "\t";
        break;
      case SEQ:
        sam_entry += query + "\t";
        break;
      case CIGAR:
        sam_entry += aln_fields[2] + "\t";
        break;
      case CHR:
        sam_entry += aln_fields[0] + "\t";
        break;
      default: {
        if (i == entry_fields.size() - 1) {
          sam_entry += entry_fields[i]; // skip tab concat
        }
        else {
          sam_entry += entry_fields[i] + "\t";
        }
      }
    }

  }
  return sam_entry;
}

bool xaPresent(vector<string> &fields) {
  for (int i=0; i < fields.size(); i++) {
    if (fields[i].substr(0,2) == XA_PREFIX) {
      return true;
    }
  }
  return false;
}

string getXA(vector<string> &fields) {
  for (int i=0; i < fields.size(); i++) {
    if (fields[i].substr(0,2) == XA_PREFIX) {
      return fields[i];
    }
  }
  return "\0";
}


void skipSamHeader(ifstream &samfile) {
  boost::regex header("(@).*");   // regex differentiates header / non-header
  string line;
  int end_of_header_line_num = 1;
  while (getline(samfile, line)) {
    if (!boost::regex_match(line, header)) {
      break;
    }
    end_of_header_line_num++;
  }

  samfile.clear();
  samfile.seekg(0, ios::beg); // return to start
  for(; end_of_header_line_num > 1; end_of_header_line_num--) {
    getline(samfile, line);
  }
}
