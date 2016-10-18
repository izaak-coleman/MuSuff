#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "../Suffix_t.h" 
#include "generate_suffix_array.h"

static const int MIN_SUF_SIZE = 1;

using namespace std;


void generateGenSuffixArray(string filename, 
                       vector<string> &text, 
                       vector<Suffix_t> &GSA) {

  loadTextFromFile(filename, text);
  loadUnsortedSuffixes(text, GSA);
  lexMergeSort(text, GSA);
}

void loadTextFromFile(string filename, vector<string> &text) {
  ifstream fin;
  fin.open(filename);
  string next_line;
  while (getline(fin, next_line)) {
    text.push_back(next_line + "$");
  }
}

void loadUnsortedSuffixes(vector<string> &text, vector<Suffix_t> &GSA){
  
  for (int read_id = 0; read_id < text.size(); read_id++) { 
    for(int offset = 0; 
        offset <= text[read_id].size() - MIN_SUF_SIZE;
        offset++) {                                 // for each suf in each read

      Suffix_t suf;                                 // generate its suffix_t
      suf.read_id = read_id;
      suf.offset = offset;
      suf.type= true;
      GSA.push_back(suf);                           // load data
    }
  }
}

void lexMergeSort(vector<string> &text, vector<Suffix_t> &GSA) {
   int l = GSA.size();  
   sort(text, GSA, 0, l);      // start recursive mergesort
}

void sort(vector<string> &text, vector<Suffix_t> &GSA, int from, int to) {
  if ((to - from) <= 1) {return;} // base case

  int mid = (from + to) / 2;      // keep dividing
  sort(text, GSA, from, mid);
  sort(text, GSA, mid, to);
  merge(text, GSA, from, mid, to);
}

void merge(vector<string> &text, vector<Suffix_t> &GSA, 
           int from, int mid, int to) {

  // make out of place copies of SA section 
  vector<Suffix_t> *left = copyOf(GSA, from, mid);
  vector<Suffix_t> *right = copyOf(GSA, mid, to);


  int left_ptr=0, right_ptr=0, GSA_ptr=from;
  bool end_of_left = false;      // stop range bound errors 
  bool end_of_right = false;

  // add suffix_t's back to GSA in lexicographical order
  while (GSA_ptr < to) {

    // right finished... keep adding lefts elems
    if (!end_of_left && end_of_right) {
      GSA[GSA_ptr++] = (*left)[left_ptr++];
    }

    // left finished... keep adding rights elems
    else if (!end_of_right && end_of_left) {
      GSA[GSA_ptr++] = (*right)[right_ptr++];
    }

    // left lexiocographically before right element, so add left next
    else if (!end_of_left &&
             lexCompare((*left)[left_ptr], (*right)[right_ptr], text)
             ) {

      GSA[GSA_ptr++] = (*left)[left_ptr++];

    }
    else if(!end_of_right){   // right lexicographcially before left
      GSA[GSA_ptr++] = (*right)[right_ptr++];
    }
    else{
      cout << "MERGE SORT ERROR" << endl;
    }


    // check bounds
    if (left_ptr == left->size()) {
      end_of_left = true;
    }
    if (right_ptr == right->size()) {
      end_of_right = true;
    }
  }

  delete left;
  delete right;
}

bool lexCompare(Suffix_t &lhs, Suffix_t &rhs, vector<string> &text) {
  // Generate pointers to lhs and rhs suffixes in text
  string::iterator lhs_iter = returnStartIterator(lhs, text);
  string::iterator lhs_end   = returnEndIterator(lhs, text);
  string::iterator rhs_iter = returnStartIterator(rhs, text);
  string::iterator rhs_end   = returnEndIterator(rhs, text);

  for( ; (lhs_iter != lhs_end && rhs_iter != rhs_end); lhs_iter++, rhs_iter++){
    // lex compare character
    if (*lhs_iter < *rhs_iter) { return true; }
    if (*rhs_iter < *lhs_iter) { return false; }
    // equiv char so move to next...
    }
    // One suf is prefix of other. Return true if left is prefix of right
    return (lhs_iter == lhs_end) && (rhs_iter != rhs_end); 
}

string::iterator returnStartIterator(Suffix_t &suf, vector<string> &text) {
  // Use suf.read_id to locate the read, and then set an iterator
  // pointing to the suf's offset
  string::iterator iter = text[suf.read_id].begin() + suf.offset;
  return iter;
}

string::iterator returnEndIterator(Suffix_t &suf, vector<string> &text) {
  // use suf.read_id to locate the read, then return an iterator to the 
  // end of that read
  string::iterator iter = text[suf.read_id].end();
  return iter;
}

vector<Suffix_t> *copyOf(vector<Suffix_t> &GSA, int from, int to){
  // make new section of SA on heap
  vector<Suffix_t> *GSA_section_ptr;
  GSA_section_ptr = new vector<Suffix_t>;

  // load with section
  for (int i = from; i < to; i++) {
    GSA_section_ptr->push_back(GSA[i]);
  }

  return GSA_section_ptr;
}

// end of file
