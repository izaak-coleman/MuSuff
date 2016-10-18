#ifndef GENERATE_SUFFIX_ARRAY
#define GENERATE_SUFFIX_ARRAY

#include <vector>
#include <string>

#include "../Suffix_t.h"

void generateGenSuffixArray(std::string filename, std::vector<std::string> &text,
                            std::vector<Suffix_t> &GSA);
// USAGE:
// string fin; vector<Suffix_t> GSA;
// generateGenSuffixArray(fin, GSA);
// DES:
// fin is datafile containing strings. This data is first extracted, and
// loaded into the text vector. Suffix_t elements are build of the 
// data in text, and then lexicographically sorted to build GSA.


void loadTextFromFile(std::string filename, std::vector<std::string> &text);
// load toy words into Text                                   

void loadUnsortedSuffixes(std::vector<std::string> &text,
                          std::vector<Suffix_t> &GSA);
//  DES:
//  Generates suffix_t objects for each suffix, upto MIN_SUF_SIZE suffixes
//  from each string in the text

void lexMergeSort(std::vector<std::string> &text, 
                  std::vector<Suffix_t> &GSA);
// Lexicographical merge sorts the suffix_t elements in GSA

void sort(std::vector<std::string> &text,
          std::vector<Suffix_t> &GSA, int from, int to);
// recursive divide function or lexMergeSort

void merge(std::vector<std::string> &text,
           std::vector<Suffix_t> &GSA, int from, int mid, int to);
// conquer function of lexMergeSort

bool lexCompare(Suffix_t &lhs, Suffix_t &rhs, std::vector<std::string> &text);
// Function lexiographically compares two suffixes, returning true if 
// lhs is before rhs, false otherwise. 

std::string::iterator returnStartIterator(Suffix_t &suf, 
                                          std::vector<std::string> &text);
// Function locates the read corresponding to suf.read_id and 
// Sets a pointer in that read starting at suf.offset

std::string::iterator returnEndIterator(Suffix_t &suf, 
                                        std::vector<std::string> &text);
// Function returns a pointer to the end of the read corresponding to 
// suf.read_id

std::vector<Suffix_t>* copyOf(std::vector<Suffix_t> &GSA, int from, int to);


#endif

