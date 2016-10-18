// gsa_file_gen.h
#ifndef GSA_FILE_GEN 
#define GSA_FILE_GEN 

#include <vector>
#include <string>

#include "../Suffix_t.h"

void buildGSAFile(std::vector<Suffix_t> &GSA, std::string filename);
// Function outputs a GSA file for persistent use;
// A GSA can be constructed directly from this output file

void constructGSAFromFile(std::string filename, std::vector<Suffix_t> &GSA);
// Function constructs a GSA from  a .gsa file

#endif
