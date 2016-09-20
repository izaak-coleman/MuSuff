// toyMops.cpp
#include <iostream>
#include <vector>
#include <string>
#include <fstream>


using namespace std;

/*----------Decls----------*/


struct suffix_t{
  int read_id;
  uint8_t offset;
  bool is_common;
  bool in_readgroup;
};

// deals with command line input
struct command_line_params {
  string filename;  // datafile name NOTE: will have a tumour and healthy
  uint8_t minimum_suffix_size;  // define the lower bound length of prefixes
};

typedef vector<suffix_t> sarray;
// sarray is SA for each prefix

typedef vector<int> readgroup;
// Defines groups of reads that have an idential prefix, 
// but differ at some branchpoint. Stores the index in the SA of each
// suffix_t in the readgroup. This is more memory efficient

typedef vector<readgroup>  readgroups;
// Vector contains all the readgroups

typedef vector<int> lcp_array;
// lpc_array contains an array of the longest common prefix for adjacent
// suffixes in the SA
// An lcp_array only needs to be computed once


// ARRAY INIT AND GENERIC HELPERS


void parseCommandLine(command_line_params &params, 
                                     int argc, char **argv);
// reads in char **argv[] command line options and sets them as
// struct datafeilds of the correct type

void loadWordsFromFile(vector<string> &reads, string filname);
// load toy words into Text

void printSuffixData(sarray &SA);
// Function prints out data for each element of the suffix array

string returnSuffix(vector<string> &reads, suffix_t suf);
// Function returns the suffix that the suffix_t represents

void loadUnsortedSuffixes(sarray &SA, 
                          const vector<string> &reads,
                          uint8_t min_suffix);
// Function reads through Text (DNA reads) generating a suffix_t for 
// each suffix of each read, and loads into SA




// MERGE SORT ---------------------------------

void lexMergeSort(sarray &SA, vector<string> &reads);
// Functin lexicographically sorts the suffixes in SA, by the end of this
// function SA has been transformed into an actual suffix Array
// NOTE: probably best to return a sorted SA, deleting a heap input

void sort(sarray &SA, vector<string> &reads, int from, int to);
// recursive divide function or lexMergeSort

void merge(sarray &SA, vector<string> &reads, int from, int mid, int to);
// conquer function of lexMergeSort

sarray* copyOf(sarray &SA, int from, int to);
// SA copy function used by merge. Returns a pointer to a subsection
// of SA spanning from indicies [from, to). On the heap

bool leftHigherThanRight(suffix_t lhs, suffix_t rhs, vector<string> &reads);
// Function compares two suffix_t to see if the lhs suffix is lexiographically
// higher than the rhs suffix, returning true if this is the case
// return the iterator coresponding to the suffix_t

string::iterator returnStartIterator(suffix_t &suf, vector<string> &reads);
// Function locates the read corresponding to suf.read_id and 
// Sets a pointer in that read starting at suf.offset

string::iterator returnEndIterator(suffix_t &suf, vector<string> &reads);
// Function returns a pointer to the end of the read corresponding to 
// suf.read_id

bool lexCompare(suffix_t &lhs, suffix_t &rhs, vector<string> &reads);
// Function lexiographically compares two suffixes, returning true if 
// lhs is before rhs, false otherwise. 
// Also, sets lhs/rhs.is_common to true
// NOTE: THIS FUNCTION WILL BE CHANGED TO ALLOW COUNTING OF COMMON READS
// NOTE: IS_COMMON IS FOR TOYMAPS EXAMPLE



// LCP GENERATION -------------------------------


int computeLongestCommonPrefix(suffix_t &isuf, suffix_t &jsuf, vector<string>
    &reads);
// Function computes the lcp shared between isuf and jsuf and returns the


// BRANCHPOINT GROUP GENERATION ---------------------



void generateReadgroups(readgroups &groups, sarray &SA,
                        vector<string> &reads);
// Function makes a pass through SA, grouping suffixes that share the same
// LCP. Such suffixes share the same branchpoint, and thus represent
// a toy example of a healthy > cancer mutation

void makeReadGroup(int index_of_group_seed, readgroup &group, 
                   sarray &SA, vector<string> &reads);
// Function is directly called by generateReadGroups. It takes the 
// index of a unique suffix_t in the SA, and uses this index to locate
// all the adjacent suffixes that share the same branch point. 
// This is done by searching searching through the LCP
// All the sequences with the same LCP max LCP to the suffixes adjacent
// to the seed are added to group


void getSuffixesFromLeft(int start_point, int lcp,
                         sarray &SA, readgroup &group, vector<string> &reads);
// Function caled directly by makeReadGroup if max LCP is between seed and 
// seed-1 in LCP. Adds all the suffixes that have lcp == lcp(seed-1, seed)


void getSuffixesFromRight(int start_point, int lcp,
                          sarray &SA, readgroup &group, vector<string> &reads);
// Function called directly by makeReadGroup if max LCP is between seed
// and seed + 1 in LCP. Adds all the suffixes that have lcp == lcp(seed+1, seed)


/*----------Definitions----------*/


int main(int argc, char **argv) {


  // process c-line args
  command_line_params params;
  parseCommandLine(params, argc, argv);
 
  /* Load Texts into the text vector*/
  vector<string> reads;    // Text
  loadWordsFromFile(reads, params.filename);        // Add toy words

  /* Generate the unsorted SA*/
  sarray SA;
  loadUnsortedSuffixes(SA, reads, params.minimum_suffix_size);
  printSuffixData(SA);    // check loaded

  // print suffixes to screen
  for(suffix_t suffix : SA) {
    cout << returnSuffix(reads, suffix) << endl;  // confirm suffixes
  }
 
  // generate SA!!
  lexMergeSort(SA, reads);

  // print suffix data and string to screen
  printSuffixData(SA);
  for(suffix_t suffix : SA) {
    cout << returnSuffix(reads, suffix) << endl;  // confirm suffixes
  }

  // pull off unique sequences and associate readgroups

  readgroups branchpoint_groups;
  generateReadgroups(branchpoint_groups, SA, reads);

  cout << "All good with read group generation... " << endl;
  
  int i=1;
  for(readgroup g : branchpoint_groups) {
    cout << i << "st group: " << endl;
    for(int suf_index : g) {
      cout << "Suffix_t index in array " << suf_index << endl;
      cout << returnSuffix(reads, SA[suf_index]) << endl;
    }
    cout << "End of group" << endl << "---------------------------------" <<
      endl << endl;
    i++;
  }



  return 0;
}


void parseCommandLine(command_line_params &params, 
                                     int argc, char **argv){
  if (argc != 3) {
    cout << "Usage: <exec> <datafile> <minimum_prefix_size>" << endl;
    exit(1);
  }
  params.filename = static_cast<string>(argv[1]);
  params.minimum_suffix_size = static_cast<uint8_t>(std::stoi(argv[2]) + 1);
}




void loadWordsFromFile(vector<string> &reads, string filename){
  ifstream datafile;
  datafile.open(filename);
  string next_read;
  while (datafile.good()) {
    getline(datafile, next_read);
    reads.push_back(next_read);
  }
}


void printSuffixData(sarray &SA) {
  for(suffix_t s : SA) {
    cout << "read_id: " << s.read_id <<  " --- " 
         << "offset: "  << unsigned(s.offset)  <<  " --- "
         << "is_common: " << ((s.is_common) ? "duplicate" : "unique") 
         << endl;
  }
}

string returnSuffix(vector<string> &reads, suffix_t suf){
  return reads[suf.read_id].substr(suf.offset); 
}


void loadUnsortedSuffixes(sarray &SA, 
                          const vector<string> &reads, 
                          uint8_t min_suffix) {

  // Loop through each read, and all its suffixes
  for(int read_id = 0; read_id < reads.size(); read_id++){

    for(int offset = 0; 
        offset < reads[read_id].size() && 
        offset <= (reads[read_id].size() - min_suffix); // stop at min suf length
        offset++) {

      // load the next suffix with relevant data
      suffix_t suf;
      suf.read_id = read_id;
      suf.offset = offset;
      suf.is_common = false;
      suf.in_readgroup = false;


      // add to SA
      SA.push_back(suf);
    }
  }
}

void lexMergeSort(sarray &SA, vector<string> &reads) {
   int l = SA.size();  
   sort(SA, reads, 0, l);      // start recursive mergesort
}

void sort(sarray &SA, vector<string> &reads, int from, int to) {
  if ((to - from) <= 1) {return;} // dividing reaches single elem, hit base case

  // otherwise...keep dividing
  int mid = (from + to) / 2;
  sort(SA, reads, from, mid);
  sort(SA, reads, mid, to);
  merge(SA, reads, from, mid, to);
}


void merge(sarray &SA, vector<string> &reads, int from, int mid, int to) {
  // make out of place copies of SA section 
  sarray *left = copyOf(SA, from, mid);
  sarray *right = copyOf(SA, mid, to);


  int left_ptr=0, right_ptr=0, SA_ptr=from;
  bool end_of_left = false;      // stop range bound errors 
  bool end_of_right = false;

  // add suffix_t's back to SA in lexicographical order
  while (SA_ptr < to) {

    // right finished... keep adding lefts elems
    if (!end_of_left && end_of_right) {
      SA[SA_ptr++] = (*left)[left_ptr++];
    }

    // left finished... keep adding rights elems
    else if (!end_of_right && end_of_left) {
      SA[SA_ptr++] = (*right)[right_ptr++];
    }

    // left lexiocographically before right element, so add left next
    else if (!end_of_left &&
             lexCompare((*left)[left_ptr], (*right)[right_ptr], reads)
             ) {

      SA[SA_ptr++] = (*left)[left_ptr++];

    }
    else if(!end_of_right){   // right lexicographcially before left
      SA[SA_ptr++] = (*right)[right_ptr++];
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


sarray* copyOf(sarray &SA, int from, int to){
  // make new section of SA on heap
  sarray *SA_section_ptr;
  SA_section_ptr= new sarray;

  // load with section
  for ( int i = from; i < to; i++) {
    SA_section_ptr->push_back(SA[i]);
  }

  return SA_section_ptr;
}

bool leftHigherThanRight(suffix_t lhs, suffix_t rhs, vector<string> &reads) {
  // get suffixes from text, lhs and rhs represent
  string lhs_string = returnSuffix(reads, lhs);
  string rhs_string = returnSuffix(reads, rhs);

  // define string iterators for lexicographical_compare
  string::iterator lhs_iter = lhs_string.begin();
  string::iterator lhs_end  = lhs_string.end();
  string::iterator rhs_iter = rhs_string.begin();
  string::iterator rhs_end  = rhs_string.end();

  return lexicographical_compare(lhs_iter, lhs_end, rhs_iter, rhs_end);
  // returns TRUE if lhs, higher then rhs
}

bool lexCompare(suffix_t &lhs, suffix_t &rhs, vector<string> &reads) {
  // Generate pointers to lhs and rhs suffixes in reads
  string::iterator lhs_iter = returnStartIterator(lhs, reads);
  string::iterator lhs_end   = returnEndIterator(lhs, reads);
  string::iterator rhs_iter = returnStartIterator(rhs, reads);
  string::iterator rhs_end   = returnEndIterator(rhs, reads);

  for( ; (lhs_iter != lhs_end && rhs_iter != rhs_end); lhs_iter++, rhs_iter++){
    // lex compare character
    if (*lhs_iter < *rhs_iter) { return true; }
    if (*rhs_iter < *lhs_iter) { return false; }
    // equiv char so move to next...
    }

    // Identical suffixes, set is_common to true
    if (lhs_iter == lhs_end && rhs_iter == rhs_end) {
      lhs.is_common = true; rhs.is_common = true;
    } 
    // One is prefix of other, return the prefix as higher suffix
    return (lhs_iter == lhs_end) && (rhs_iter != rhs_end);
}

string::iterator returnStartIterator(suffix_t &suf, vector<string> &reads) {
  // Use suf.read_id to locate the read, and then set an iterator
  // pointing to the suf's offset
  string::iterator iter = reads[suf.read_id].begin() + suf.offset;
  return iter;
}

string::iterator returnEndIterator(suffix_t &suf, vector<string> &reads) {
  // use suf.read_id to locate the read, then return an iterator to the 
  // end of that read
  string::iterator iter = reads[suf.read_id].end();
  return iter;
}


void computeLongestCommonPrefixArray(lcp_array &LCP, sarray &SA, 
                                     vector<string> &reads) {

  // Reserver space in lcp_array for efficieny
  LCP.reserve(SA.size() - 1);

  // iterate through adjacent pairs in SA, computing the Lcp
  for(int i=0, j=1; j < SA.size(); i++, j++) {
    LCP.push_back( computeLongestCommonPrefix(SA[i], SA[j], reads) );
  }

}

int computeLongestCommonPrefix(suffix_t &isuf, suffix_t &jsuf, 
                               vector<string> &reads) {

  // Get suffix pointers in reads
  string::iterator isuf_iter = returnStartIterator(isuf, reads);
  string::iterator isuf_end   = returnEndIterator(isuf, reads);
  string::iterator jsuf_iter = returnStartIterator(jsuf, reads);
  string::iterator jsuf_end   = returnEndIterator(jsuf, reads);

  // computes lcp
  int lcp = 0;

  while (*isuf_iter == *jsuf_iter &&
         isuf_iter != isuf_end    &&
         jsuf_iter != jsuf_end      ) {

    lcp++;      // matched next char, so extend lcp

    isuf_iter++; jsuf_iter++;
  }

  return (lcp);   
}

void generateReadgroups(readgroups &groups, sarray &SA, 
                        vector<string> &reads) {

  // Pass through SA locating unique sequences
  for(int i=0; i < SA.size(); i++) {

    if (!SA[i].is_common) { // not common, so unique

      readgroup newgroup;     
      newgroup.push_back(i);

      // extend left and/or right of the unique suffix
      // gathering the branchpoint group
      makeReadGroup(i, newgroup, SA, reads);

      groups.push_back(newgroup); // store the new group
    }
  }
}

void makeReadGroup(int seed_index, readgroup &group, 
                   sarray &SA, vector<string> &reads) {

  // seed_index is the index of the unique suffix_t this function
  // was called with. 


  // assume no lcp
  int lcp_right_adjacent = 0, lcp_left_adjacent = 0;

  // left bound SEG FAULT PROTECT
  if (seed_index == 0) { 
    // Calculat the lcp between the seed and the right adj. suffix in SA
    lcp_right_adjacent = computeLongestCommonPrefix(SA[seed_index],
        SA[seed_index+1], reads);
  }

  // right bound SEG FAULT PROTECT
  else if (seed_index == SA.size()-1) {
    // Calculate the lcp between the seed and left adj. suffix in SA
    lcp_left_adjacent = computeLongestCommonPrefix(SA[seed_index-1],
        SA[seed_index], reads);
  }

  // Safe to acess left and right
  else {
    lcp_right_adjacent = computeLongestCommonPrefix(SA[seed_index],
        SA[seed_index+1], reads);
    lcp_left_adjacent = computeLongestCommonPrefix(SA[seed_index-1],
        SA[seed_index], reads);
  }
  


  // Compute the position of the branchpoint by identifying the max lcp
  // of adjacent sequence. The branchpoint will always be the longest
  // of the LCPs.

  if (lcp_right_adjacent == 0 && lcp_left_adjacent == 0) {
    return;     // the suffix is completely unique, so a suffix of a 
                // shorter unique branch, do not use.
  }

  // So, branchpoint exists to the left, gather sequences
  else if (lcp_left_adjacent > lcp_right_adjacent){
    getSuffixesFromLeft(seed_index, lcp_left_adjacent, SA, group, reads);
  }
  
  // So, branchpoint exists to right, gather sequences
  else if (lcp_left_adjacent < lcp_right_adjacent) {
    getSuffixesFromRight(seed_index, lcp_right_adjacent, SA, group, reads);
  }

  else{ // already in the center of the branchpoint, so gather left and right
    getSuffixesFromLeft(seed_index, lcp_left_adjacent, SA, group, reads);
    getSuffixesFromRight(seed_index, lcp_right_adjacent, SA, group, reads);
  }
}

void getSuffixesFromLeft(int seed_index, int lcp, 
                         sarray &SA, readgroup &group, vector<string> &reads) {

  int left_arrow = seed_index-1;
  // While lexicographally adjacent suffixes share the same lcp value
  // they have the same branchpoint, thus they are in the same group,
  // so add them

  while( // lcps are same AND not out of bounds AND not already in group...
      SA[left_arrow].in_readgroup == false &&
      left_arrow > 0           
      && computeLongestCommonPrefix(SA[left_arrow], SA[seed_index], reads) == lcp
      ) {

    // ...add to seed group
    SA[left_arrow].in_readgroup = true;
    group.push_back(left_arrow);
    left_arrow--;
  }
}

void getSuffixesFromRight(int seed_index, int lcp, 
                          sarray &SA, readgroup &group, vector<string> &reads) {

  int right_arrow = seed_index+1;
  // While lexicographically adjacent suffixes share the same lcp val
  // they have the same branchpoint, thus they are in the same group
  // so add them!!

  while ( // lcps are the same AND not out of bounds AND not already in group...
      SA[right_arrow].in_readgroup == false &&
      right_arrow < (SA.size()-1)         // max LCP size is one less than SA
      && computeLongestCommonPrefix(SA[right_arrow], SA[seed_index], reads) == lcp
      ) {

    // ...add to seed group
    SA[right_arrow].in_readgroup = true;
    group.push_back(right_arrow);
    right_arrow++;
  }
}


// end of file
