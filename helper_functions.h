#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H

#include <vector>
#include <string>
#include <utility>



#include "Suffix_t.h"

class ReadsManipulator;       // forward decl.

enum {TUMOUR, HEALTHY, SWITCHED};
enum {LEFT, RIGHT};         // distinguish, paired end type

// the reason for file_and_type weird field names is this struct used to 
// be a pair, with just first and second. So, to save myself changing
// all the field names in code, I just extended namesystem to struct names

struct file_and_type {
  std::string first; // file name
  bool second;  // tissue type
};


struct command_line_params {
  std::vector<file_and_type> datafiles;
  uint8_t minimum_suffix_size;
  double econt;
};




void parseCommandLine(command_line_params &params, int argc, char **argv);
// reads in char **argv[] command line options and sets them as
// struct datafeilds of the correct type


void loadWordsFromFile(std::vector<std::string> &reads, std::string filname);
// load file words into text


std::string::iterator returnStartIterator(Suffix_t &suf, std::vector<std::string> &reads);
// Function locates the read corresponding to suf.read_id and 
// Sets a pointer in that read starting at suf.offset

std::string::iterator returnEndIterator(Suffix_t &suf, std::vector<std::string> &reads);
// Function returns a pointer to the end of the read corresponding to 
// suf.read_id


std::string returnSuffix(std::vector<std::string> &reads, Suffix_t suf);
// Function returns the suffix that the suffix_t represents


int computeLongestCommonPrefix(Suffix_t &isuf, Suffix_t &jsuf,
                               ReadsManipulator &reads);
// Function returns the length of the lcp between isuf and jsuf



#endif
