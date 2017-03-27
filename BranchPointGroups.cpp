// BranchPointGroups.cpp
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <set>
#include <utility>
#include <iterator>
#include <map>
#include <algorithm>
#include <limits>
#include <thread>
#include <mutex>
#include <cstdlib> // exit

#include "radix.h"
#include "util_funcs.h"
#include "BranchPointGroups.h"
#include "SuffixArray.h"
#include "Suffix_t.h"
#include "GenomeMapper.h"

#include "benchmark.h"


using namespace std;
static const int N_THREADS = 64;
static const int TRIM_VALUE = 4;
static const int COVERAGE_UPPER_THRESHOLD = 80;
static const int READ_LENGTH = 100;
static const int CTR_GSA1 = 4;
static const int CTR_GSA2 = 4;


BranchPointGroups::BranchPointGroups(SuffixArray &_SA, 
                                     ReadsManipulator &_reads,
                                     double _econt) {
  econt = _econt;
  reads = &_reads;  
  SA = &_SA;    

  // Generate branchpoint groups
  cout << "Extracting cancer-specific reads..." << endl;
   extractCancerSpecificReads(); 
  //outputExtractedCancerReads("/data/ic711/point3.txt");
  cout << "No of extracted reads: " << CancerExtraction.size() << endl;

  // Group blocks covering same mutations in both orientations
  cout << "Generating breakpoint blocks..." << endl;
//  makeBreakPointBlocks();
  seedBreakPointBlocks();
//  outputFromBPB("/data/ic711/point4.txt");

  cout << "made " << BreakPointBlocks.size() << " blocks. " << endl;
  cout << "Adding non-mutated alleles to blocks." << endl;
  extractNonMutatedAlleles();
  //outputFromBPB("/data/ic711/point5.txt");
}


void BranchPointGroups::outputExtractedCancerReads(std::string const& filename) {
  ofstream ofHandle(filename.c_str());
  for (unsigned int read_idx : CancerExtraction) {
    ofHandle << "(" << read_idx << ",T)" << std::endl;
  }
  ofHandle.close();
}

void BranchPointGroups::outputFromBPB(std::string const& filename) {
  ofstream ofHandle(filename.c_str());
  for (bp_block const& b : BreakPointBlocks) {
    for (read_tag const& tag : b.block) {
      ofHandle << "(" << tag.read_id << ","
               << ((tag.tissue_type) ? "H" : "T")
               << ")" << std::endl;
    }
  }
  ofHandle.close();
}


BranchPointGroups::~BranchPointGroups() {
}

void BranchPointGroups::extractCancerSpecificReads() {

  unsigned int elements_per_thread = (SA->getSize()  / N_THREADS);
  cout << "Elements per thread" << elements_per_thread << endl;
  cout << "GSA size: " << SA->getSize() << endl;
  vector<thread> workers;

  unsigned int from=0, to= elements_per_thread;
  for(int i=0; i < N_THREADS; i++) {
    workers.push_back(
      std::thread(&BranchPointGroups::extractionWorker, 
        this, from, to)
    );

    from = to;
    if(i == N_THREADS - 2) {
      to = SA->getSize();
    }
    else {
      to += elements_per_thread;
    }
  }

  for(auto &thread : workers) {
      thread.join();
  }
}

void BranchPointGroups::extractionWorker(unsigned int seed_index, unsigned int to) {
  set<unsigned int> threadExtr;
  unsigned int extension {seed_index + 1};


  while (seed_index < to && seed_index != SA->getSize() - 1) {   // CONFIRM EFFECT OF THIS
    double c_reads{0}, h_reads{0};    // reset counts

    // Assuming that a > 2 group will form, start counting from seed_index
    if (SA->getElem(seed_index).type == HEALTHY) h_reads++;
    else c_reads++;

    while (::computeLCP(SA->getElem(seed_index), SA->getElem(extension), *reads)
        >= reads->getMinSuffixSize()) {
      // tally tissue types of group
      if(SA->getElem(extension).type == HEALTHY) h_reads++;
      else c_reads++;
      extension++;
      if (extension == SA->getSize()) break;    // bound check GSA
    }

    // Group size == 1 and group sizes of 1 permitted and groups is cancer read
    if (extension - seed_index == 1 && CTR_GSA1  == 1 && SA->getElem(seed_index).type == TUMOUR) {
      threadExtr.insert(SA->getElem(seed_index).read_id);           // extract read
    }
    else if (c_reads >= CTR_GSA1 && (h_reads / c_reads) <= econt)  {
      for (unsigned int i = seed_index; i < extension; i++) {
        if (SA->getElem(i).type == TUMOUR) {
          threadExtr.insert(SA->getElem(i).read_id);
        }
      }
    }
    seed_index = extension++; 
  }

  // On very rare cases, the very last value in the GSA is unique. 
  // In such a case, this value will not have been analysed 
  // by the loop. In these cases, seed_index == SA.size()-1, rather
  // than SA.size() - which is the case if the extension proceeds to the
  // end.  Therefore, if the case condition is true, we need to check
  // it separately.
  if (seed_index == SA->getSize() -1) {
    if (CTR_GSA1 == 1 && SA->getElem(seed_index).type == TUMOUR) {
      threadExtr.insert(SA->getElem(seed_index).read_id);
    }
  }

  // Load to CancerExtraction, avoiding thread interference
  std::lock_guard<std::mutex> lock(cancer_extraction_lock);
  for(unsigned int extracted_cancer_read : threadExtr) {
      CancerExtraction.insert(extracted_cancer_read);
  }
}


////Old version, does not correctly take into account CTR, and 
//// also suffers from rare case bug
//void BranchPointGroups::extractionWorker(unsigned int seed_index, 
//                                                        unsigned int to) {
//
//  set<unsigned int> localThreadExtraction;
//  double cancer_sequences = 0, healthy_sequences = 0;
//  unsigned int extension;
//
//
//  while (seed_index < to-1) {
//   
//    // Calc econt of suffixes sharing same genomic location (assumption: 
//    // LCP >= 30). If econt below user specified value, extract group
//    if(::computeLCP(SA->getElem(seed_index),
//                                  SA->getElem(seed_index+1), *reads) >= 30) {
//
//      // tally up the tissue types of seed suffix
//      if (SA->getElem(seed_index).type == HEALTHY) {
//        healthy_sequences++;
//      } else {  // TUMOUR
//        cancer_sequences++;
//      }
//
//      extension = seed_index+1;
//      while ((::computeLCP(SA->getElem(seed_index), SA->getElem(extension), *reads) >= 30)
//            ) {
//
//        if (SA->getElem(extension).type == HEALTHY) {
//          healthy_sequences++;
//        } else { // TUMOUR 
//          cancer_sequences++;
//        }
//        extension++;
//        if(extension == SA->getSize()) {break;}
//      }
//
//      // check below econt
//      if (cancer_sequences >= CTR_GSA1  && 
//         ((healthy_sequences / cancer_sequences) <= econt )) { // permit group
//        for (unsigned int i = seed_index; i < extension; i++) {
//          if(SA->getElem(i).type == TUMOUR) {     // only extr. cancer reads
//            localThreadExtraction.insert(SA->getElem(i).read_id);
//            //{
//            //  std::lock_guard<std::mutex> cout_guard(cout_lock);
//            //  cout << reads->returnSuffix(SA->getElem(i)) << endl;
//            //}
//          }
//        }
//      }
//      seed_index = extension; // jump scan to end of group
//    }
//
//    else { // if CTR == 1 {
//      // In this case, there is no > 2 group of reads containing 
//      // this 30bp stretch of DNA. However, the read is still unique
//      // it could be unique due to a mutation, and a group could not
//      // form due to MTT, mutation distribution or low coverage 
//      // therefore, if it is a cancer read, extract it
//      //if (SA->getElem(seed_index).type == TUMOUR) {
//      //  localThreadExtraction.insert(seed_index);
//      //}
//      seed_index++;
//    }
//
//    healthy_sequences = 0; 
//    cancer_sequences = 0; 
//
//  }//end while
//
//  // After extraction, threads need to load their data into global CancerExtractionSet
//
//
//  std::lock_guard<std::mutex> lock(cancer_extraction_lock); // avoid thread interference
//  for(unsigned int read_with_mutation : localThreadExtraction) {
//      CancerExtraction.insert(read_with_mutation);
//    }
//}

void BranchPointGroups::seedBreakPointBlocks() {

  string concat("");
  vector<pair<unsigned int, unsigned int>> binary_search_array;

  // load the first element
  set<unsigned int>::iterator it = CancerExtraction.begin();
  pair<unsigned int, unsigned int> first_read(*it,0);
  binary_search_array.push_back(first_read);    // first read starts at zero
  concat += reads->getReadByIndex(*it, TUMOUR);
  concat += reverseComplementString( 
      reads->getReadByIndex(*it, TUMOUR)
  ) + "$";

  cout << "Cancer Extraction size: " <<  CancerExtraction.size() << endl;
  unsigned int concat_idx = 0;
  it++;
  for (; it != CancerExtraction.end(); it++) {
    // build concat
    concat += reads->getReadByIndex(*it, TUMOUR);
    concat += reverseComplementString(reads->getReadByIndex(*it, TUMOUR)) + "$";
    // add bsa values
    concat_idx += reads->getReadByIndex(*std::prev(it), TUMOUR).size() * 2;
    pair<unsigned int, unsigned int> read_concat_pair (*it, concat_idx);
    binary_search_array.push_back(read_concat_pair);
  }

  // Build SA
  cout << "Building cancer specific sa" << endl;
  unsigned long long *radixSA = 
    Radix<unsigned long long>((uchar*) concat.c_str(), concat.size()).build();
  cout << "Size of cancer specific sa: " << concat.size() << endl;


  cout << "Transforming to cancer specfic gsa" << endl;

  // transform to GSA
  vector<read_tag> gsa;
  for (unsigned long long i=0; i < concat.size(); i++) {
    pair <unsigned int, unsigned int> read_concat_pair = 
      SA->binarySearch(binary_search_array, radixSA[i]);

    unsigned int offset = radixSA[i] - read_concat_pair.second;
    int read_size = reads->getReadByIndex(read_concat_pair.first, TUMOUR).size();

    bool orientation = RIGHT;
    if (offset >= read_size) {
      orientation = LEFT;
      offset -= read_size;
    }
    // remove suffixes that are too short
    if (read_size - offset <= reads->getMinSuffixSize())  continue;

    // read is stored in forward orientation, convert if LEFT
    if (orientation == LEFT) {
      offset = read_size - offset - reads->getMinSuffixSize() - 1; // dollar -1
    }

    read_tag tag;
    tag.read_id = read_concat_pair.first;
    tag.orientation = orientation;
    tag.offset = offset;
    tag.tissue_type = TUMOUR;
    gsa.push_back(tag);       // gsa should be built
  } 

  cout << "Extracting groups from cancer specific gsa" << endl;
  delete [] radixSA;

  //// PRINT gsa
  //for (read_tag const& tag : gsa) {
  //  std::cout << readTagToString(tag) 
  //            << ((tag.orientation) ? " -- R" : " -- L") << endl;
  //}

  extractGroups(gsa);
}

string BranchPointGroups::readTagToString(read_tag const& tag) {
  string read = reads->getReadByIndex(tag.read_id, tag.tissue_type);
  int offset = tag.offset;
  string dollar = "";
  if (tag.orientation == LEFT) {
    offset = read.size() - (tag.offset + reads->getMinSuffixSize() + 1);
    read = reverseComplementString(read);
    dollar = "$";
  }
  return read.substr(offset) + dollar;
}

int BranchPointGroups::computeLCP(read_tag const& a, read_tag const& b) {
  string a_str = readTagToString(a);
  string b_str = readTagToString(b);
  int lcp = 0;
  while (lcp < a_str.length() &&
         lcp < b_str.length() &&
         a_str[lcp] == b_str[lcp]) lcp++;
  return lcp;
}

//int BranchPointGroups::computeLCP(read_tag a, read_tag b) {
//  string::const_iterator a_it, b_it, a_end, b_end;
//  const string &astr = reads->getReadByIndex(a.read_id, TUMOUR);
//  const string &bstr = reads->getReadByIndex(b.read_id, TUMOUR);
//
//  int a_inc, b_inc;
//  bool a_rc = false, b_rc = false;
//
//  if (a.orientation == RIGHT) {
//    a_it = astr.begin() + a.offset;
//    a_end = astr.end();
//    a_inc = 1;
//  }
//  else { 
//    a_it = astr.begin() + a.offset + reads->getMinSuffixSize() - 1;
//    a_end = astr.begin();
//    a_inc = -1;
//    a_rc = true;
//
//  }
//
//  if (b.orientation == RIGHT) {
//    b_it = bstr.begin() + b.offset;
//    b_end = bstr.end();
//    b_inc = 1;
//  }
//  else {
//    b_it = bstr.begin() + b.offset + reads->getMinSuffixSize() - 1;
//    b_end = bstr.begin();
//    b_inc = -1;
//    b_rc = true;
//  }
//
//  for(; revCompCharacter(*a_it, a_rc) == revCompCharacter(*b_it, b_rc) &&
//        b_it != b_end &&
//        a_it != a_end; 
//        a_it += a_inc, b_it += b_inc, lcp++);
//
//}


//int BranchPointGroups::computeLCP(read_tag a, read_tag b) {
//  string::const_iterator a_it, b_it, a_end, b_end;
//  const string &astr = reads->getReadByIndex(a.read_id, TUMOUR);
//  const string &bstr = reads->getReadByIndex(b.read_id, TUMOUR);
//  int a_inc, b_inc;
//  bool a_rc = false, b_rc = false;
//
//  if (a.orientation == RIGHT) {
//    a_it = astr.begin() + a.offset;
//    a_end = astr.end();
//    a_inc = 1;
//  }
//  else {
//    a_it = astr.end() - 2 - a.offset; // 2 avoids dollar symbol
//    a_end = astr.begin();
//    a_inc = -1;
//    a_rc = true;
//  }
//
//  if (b.orientation == RIGHT) {
//    b_it = bstr.begin() + b.offset;
//    b_end = bstr.end();
//    b_inc = 1;
//  }
//  else {
//    b_it = bstr.end() - 2 - b.offset; // 2 avoids dollar symbol
//    b_end = bstr.begin();
//    b_inc = -1;
//    b_rc = true;
//  }
//
//  int lcp = 0;
//  for(; revCompCharacter(*a_it, a_rc) == revCompCharacter(*b_it, b_rc) &&
//        b_it != b_end &&
//        a_it != a_end; 
//        a_it += a_inc, b_it += b_inc, lcp++);
//
//  return lcp;
//}

char BranchPointGroups::revCompCharacter(char ch, bool rc) {
  if (!rc) return ch;
  switch (ch) {
    case 'A': return 'T';
    case 'T': return 'A';
    case 'G': return 'C';
    case 'C': return 'G';

    default: {cout << "Hit Dollar" << endl; exit(1);}
  }
}

//// New version: Takes into account CTR and does not have rare case bug
void BranchPointGroups::extractGroups(vector<read_tag> const& gsa) {
  // Adding from, to parameters, in order to make logic compatible
  // with later multithreading modifications.
  unsigned int to{gsa.size()}, seed_index{0};
  unsigned int extension{seed_index + 1};
  unsigned int block_id;

  // < to is thread safe, != gsa.size() - 1 is rare case bound safe
  while (seed_index < to && seed_index != gsa.size() - 1) {
    bp_block block;
    // compute group size - avoid inserting here to avoid unnecessary mallocs
    while (computeLCP(gsa[seed_index], gsa[extension]) >= reads->getMinSuffixSize()) {
      extension++;
      if (extension == gsa.size()) break;
    }
    // Make alloc if group above CTR
    if (extension - seed_index >= CTR_GSA2) {
      for (int i=seed_index; i < extension; i++) block.block.insert(gsa[i]); 
    }
    else {    // continue, discarding group
      seed_index = extension++;
      continue;
    }
  
    // Load group into break point blocks
    block.id = block_id;
    BreakPointBlocks.push_back(block);
    block_id++;
    seed_index = extension++;
  }

  if (seed_index == gsa.size() - 1 && CTR_GSA2 == 1) {
    bp_block block;
    block.block.insert(gsa[seed_index]);
    block.id = block_id;
    BreakPointBlocks.push_back(block);
  }
}

//// Old version: Does not take into account CTR, and rare case bug
//void BranchPointGroups::extractGroups(vector<read_tag> &gsa) {
//
//  unsigned int seed_index{0};
//  unsigned int extension{0};
//  unsigned int block_id{0};
//  while (seed_index < gsa.size()-1) {
//   
//    // seed
//    if (computeLCP(gsa[seed_index], gsa[seed_index+1]) >= reads->getMinSuffixSize()) {
//      bp_block block;
//      block.block.insert(gsa[seed_index]); // add seed to block
//      extension = seed_index+1;
//      // extend
//      while ((computeLCP(gsa[seed_index], gsa[extension])
//           >= reads->getMinSuffixSize())) {
//        block.block.insert(gsa[extension]);
//        extension++;
//        if (extension == gsa.size()) break;
//      }
//      if (block.block.size() < 4)  {
//        seed_index = extension;
//        continue;
//      }
//
//
//      block.id = block_id;
//      BreakPointBlocks.push_back(block);
//      block_id++;
//
//      seed_index = extension;
//    }
//    else {
//      seed_index++;
//    }
//  }
//  
//  // PRINT groups
//  for (bp_block const& b : BreakPointBlocks) {
//    cout << "Block id: " << b.id << endl;
//    for (read_tag const& tag : b.block) {
//      cout << readTagToString(tag) << endl;
//    }
//  }
//}

void BranchPointGroups::makeBreakPointBlocks(){
  if(CancerExtraction.size() == 0) {
    cout << "No mutations were identified " << endl;
    exit(1);
  }
 

  std::hash<std::string> hash;
  multimap<string, read_tag> mutation_grouper;

  for(unsigned int read_index : CancerExtraction) {
    string cancer_read = reads->getReadByIndex(read_index, TUMOUR);
    for(int i=0; i < cancer_read.size()-30; i++) {    // make 30bp fwd/rev

      pair<string, read_tag> fwd_key_val, rev_key_val;

      // forward and reverse pair of each 30pb window
      string sub_str_fwd = cancer_read.substr(i, 30);
      string sub_str_rev = reverseComplementString(sub_str_fwd);
      unsigned long long fwd_hash, rev_hash;
      fwd_hash = hash(sub_str_fwd);
      rev_hash = hash(sub_str_rev);
    

      // init forward tag
      read_tag fwd_tag;
      fwd_tag.read_id = read_index;
      fwd_tag.offset = i;
      fwd_tag.orientation = RIGHT;
      fwd_tag.tissue_type = TUMOUR;

      // init reverse tag
      read_tag rev_tag;
      rev_tag.read_id = read_index;
      rev_tag.offset = i;
      rev_tag.orientation = LEFT;
      rev_tag.tissue_type = TUMOUR;

      // make pairs and load into map
      fwd_key_val.first = sub_str_fwd;
      fwd_key_val.second = fwd_tag;

      rev_key_val.first = sub_str_rev;
      rev_key_val.second = rev_tag;

      mutation_grouper.insert(fwd_key_val);
      mutation_grouper.insert(rev_key_val);

    }
 }

  // now loaded into map, extract blocks
  bp_block block;

  multimap<string, read_tag>::iterator it = mutation_grouper.begin();  

  unsigned int block_id = 0;
  while(it != std::prev(mutation_grouper.end())){
    bool left_mate = false, right_mate = false;

    if(it->first == std::next(it)->first) { // then same genomic location

      string seed = it->first;
      block.insert(it->second);           

      if(it->second.orientation) {
        right_mate = true;
      }
      else {
        left_mate = true;  
      }

      // extend group, all equal to seed
      it++;  
      while(seed == it->first) {
        block.insert(it->second);

        if(it->second.orientation) {
          right_mate = true;
        }
        else {
          left_mate = true;
        }

        it++;
        if(it == mutation_grouper.end()) {break;}
      }    
    
      // only extract double orientation groups with CTR 4
      if(/*left_mate & right_mate && */ block.size() >= 4) {
        block.id = block_id;
        BreakPointBlocks.push_back(block);
        block_id++;
      }

      block.clear();
    }

    else {
      it++;
    }

    if(it == mutation_grouper.end()) {break;}
  }

}



string BranchPointGroups::reverseComplementString(string s){
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


//void BranchPointGroups::extractNonMutatedAlleles() {
//// NB: When we perform integer grouping, this will have be be performed
//// in forward and reverse orientation
//
//  // only need to seach overlap in each orientation, rather than whole set
//  // because the sequence is identical
//  bool searched_forward = false, searched_reverse = false;
//
//  for(bp_block &block : BreakPointBlocks) {
//
//    if(searched_forward && searched_reverse) {  // finished search
//      break;
//    }
//
//    for(read_tag tag : block.block) {
//      if(searched_forward && (tag.orientation == RIGHT)) { // done fwd search
//        continue;
//      }
//      if(searched_reverse && (tag.orientation == LEFT)) {
//        continue;
//      }
//
//      // if passed here, we will begin extracting the non mutated reads
//      // that also share the same 30bp alignment with the group in 
//      // both orientations, 
//
//      // get the read the tag represents
//      string read = reads->getReadByIndex(tag.read_id, tag.tissue_type);
//
//      // get the alligned section
//      string aligned_section = read.substr(tag.offset, 30);
//
//      long long int index = binarySearch(aligned_section); // MAYBE SWITCH FOr REGULAR BS
//      if(index == -1) {
//        cout << "Binary search failed for block: " << block.id << endl;
//        cout << "Searching on: " << ((tag.orientation == RIGHT) ?  "RIGHT" :
//          "LEFT") << endl;
//      }
//      else {
//        extendBlock(index, block.block, tag.orientation); // h reads assigned to same or
//      }
//
//      // update search switches
//      if(tag.orientation == RIGHT) {
//        searched_forward = true;
//      }
//      else {
//        searched_reverse = true;
//      }
//    }
//
//    searched_forward = false;
//    searched_reverse = false;
//
//  }
//}

void BranchPointGroups::extractNonMutatedAlleles() {

  for (bp_block &block : BreakPointBlocks) {
    read_tag tag = *block.block.begin();  // copy first element
    string read = reads->getReadByIndex(tag.read_id, tag.tissue_type);
    read = read.substr(tag.offset, reads->getMinSuffixSize());
    string rev_read = reverseComplementString(read);

    long long int fwd_index = binarySearch(read);
    long long int rev_index = binarySearch(rev_read);

    if (fwd_index != -1) {
      extendBlock(fwd_index, block.block, tag.orientation);
    }
    else {
      //std::cout << "Search failed" << std::endl;
    }

    if (rev_index != -1) {
      extendBlock(rev_index, block.block, !tag.orientation);
    }
    else {
      //std::cout << "Search failed" << std::endl;
    }
  }
}

bool BranchPointGroups::generateConsensusSequence(unsigned int block_idx, 
    int &cns_offset, bool tissue_type, unsigned int &pair_id, 
    string & cns, string & qual) {

  // DEBUG BOOL

  // all seqs get converted to RIGHT orientation, before consensus

  if (BreakPointBlocks[block_idx].size() > COVERAGE_UPPER_THRESHOLD) {
    return true;
  } 


  // select only one tissue type
  vector<read_tag> type_subset;
  for(read_tag tag : BreakPointBlocks[block_idx].block) {
    if (tissue_type == HEALTHY  && 
       (tag.tissue_type == HEALTHY || tag.tissue_type == SWITCHED)) {
      type_subset.push_back(tag);
    }
    else if(tissue_type == TUMOUR && tag.tissue_type == TUMOUR) {
      type_subset.push_back(tag);
    }
  }

  pair_id = BreakPointBlocks[block_idx].id;

  if(type_subset.size() == 0) {   // no seq. of tissue type, cannot be mapped
    cout << "SKIPPING FROM HERE" << endl;
    return true;
   // DEBUG_BOOL = true;
  }

  // perform offset conversion, converting LEFT offsets to RIGHT
  int max_offset = 0;
  int min_offset = numeric_limits<int>::max();

  // convert the offset indexes from LEFT to equivalent index in RIGHT
  for(read_tag &tag : type_subset) {
    if (tag.orientation == LEFT) {
      int read_size = reads->getReadByIndex(tag.read_id,
          tag.tissue_type).size();
      tag.offset = (read_size - (tag.offset + reads->getMinSuffixSize() + 1));
    }
  //for (read_tag &tag : type_subset) {
    if(tag.offset > max_offset) {   // also record max and min offset
      max_offset = tag.offset;
    }
    if(tag.offset < min_offset) {
      min_offset = tag.offset;
    }
  }

  // value did not change and will cause std::bad_alloc, so set to 0
  if(min_offset == numeric_limits<int>::max()) {
    min_offset = 0;
  }

  // decl alignment block
  vector<string> aligned_block;
  // initialize align_counter.
  vector< vector<int> > align_counter;
  for (int n_vectors=0; n_vectors < 4; n_vectors++) {
    vector<int> v(max_offset + READ_LENGTH - min_offset, 0);
    align_counter.push_back(v);
  }


  for(read_tag tag : type_subset) {


    // if left, convert to right

    string read;
    if(tag.orientation == LEFT) {
     read = reverseComplementString( 
         reads->getReadByIndex(tag.read_id, tag.tissue_type)
         );
    }
    else {
      read = reads->getReadByIndex(tag.read_id, tag.tissue_type);
      read.pop_back();    // remove dollar symbol
    }

    for(int i=0; i < read.size(); i++) {
      switch(read[i]) {
        case 'A':
          align_counter[0][(max_offset - tag.offset) + i]++;
          break;

        case 'T':
          align_counter[1][(max_offset - tag.offset) + i]++;
          break;

        case 'C':
          align_counter[2][(max_offset - tag.offset) + i]++;
          break;

        case 'G':
          align_counter[3][(max_offset - tag.offset) + i]++;
          break;
      }
    }
    aligned_block.push_back(addGaps(max_offset - tag.offset) + read);
  }


  cns = ""; // seed empty consensus sequence

  int n_skipped_start_pos=0;
  bool hit_consensus_min = false;
  for (int pos=0; pos < align_counter[0].size(); pos++) {

    // only start when the number of reads aligned at that position
    /// is >= TRIM_VALUE

    int nreads=0;
    for (int base=0; base < 4; base++) {
      nreads += align_counter[base][pos];
    }
    if (nreads < TRIM_VALUE) {
      if (!hit_consensus_min) {
        n_skipped_start_pos++;		// need to nibble off skipped positions from cns offset
      }
      invalidatePosition(align_counter, pos);   // invalidate rather than trim
      continue;                   // continue performs the trimming
    }
    else {
      hit_consensus_min = true;
    }


    int maxVal = 0, maxInd = 0;

    // find highest freq. base
    for(int base=0; base < 4; base++) {
      if (align_counter[base][pos] > maxVal) {
        maxVal = align_counter[base][pos];
        maxInd = base;
      }
    }

    switch(maxInd) {
      case 0:
        cns += "A";
        break;
      case 1:
        cns += "T";
        break;
      case 2:
        cns += "C";
        break;
      case 3:
        cns += "G";
        break;
    }
  }// end for
  
  // if empty string, return that we need to skip block
  if (cns == "") {
    return true;
  }
  //if (cns == "" && type_subset.size() < 4) DEBUG_BOOL = true;


  //// DEBUG
  //std::cout << ((tissue_type) ? "Healthy" : "Cancer") << " sub-block below" <<
  //std::endl;
  //cout << "Block id: " << BreakPointBlocks[block_idx].id  << endl;
  //for (int i=0; i < aligned_block.size(); i++) { // SHOW ALIGNED BLOCK
  //  cout << aligned_block[i]  << ((type_subset[i].orientation == RIGHT) ? ", R" : ", L") 
  //       << ((type_subset[i].tissue_type % 2) ? ", (H," : 
  //           ((type_subset[i].tissue_type == SWITCHED) ? ", (S," : ", (T, ridx: ")) 
  //       << type_subset[i].read_id << ")" << endl;
  //}
  //cout << "CONSENSUS AND CNS LEN" <<  cns.size() << endl;
  //cout << cns << endl << endl;
  //cout << "QSTRING" << endl;
  //cout << buildQualityString(align_counter, cns,  tissue_type) << endl;


  //for(vector<int> v : align_counter) {
  //  for(int i : v) {
  //    cout << i << ",";
  //  }
  //  cout << endl;
  //}
  //cout << endl;
  

  // set the offset value
  cns_offset = (max_offset - n_skipped_start_pos);

  // set the quality string
  qual = buildQualityString(align_counter, cns, tissue_type);
  return false;
  //return DEBUG_BOOL;
}

string BranchPointGroups::buildQualityString(vector< vector<int> > const& freq_matrix,
    string const& cns, bool tissue) {
  // Function steps through each position of the string and
  // determines whether a position should be masked if:
  //  -- Cancer: if the number of bases contributing to the 
  //             consensus base is < CTR then mask. Or, if the 
  //             number of bases with frequency above the error threshold
  //             (ALLELIC_ERROR_THRESH) is > 1 then mask.
  //  -- Health: if the number of bases with frequency above the error threshold
  //             (ALLELIC_ERROR_THRESH) is > 1 then mask.
  
  string q_str("");
  int m_start{0};
  for (int i=0; i < freq_matrix[0].size(); i++) {
    if (freq_matrix[0][i] != -1) {
      m_start = i;
      break;
    }
  }

  //if (cns.size() == 0) {
  //  cout << "No consensus sequence! " << endl;
  //  exit(1);
  //}
  for (int pos=0; pos < cns.size(); pos++) {

    // Determine the number of bases above the error frequency
    // by the calulation:
    //      bases / total bases = freq.
    // bases is equivalent to the number of reads with a base of type given
    // base. 
    if (tissue == HEALTHY) {

      double total_bases = 0;                     // get total_bases
      for (int base=0; base < 4; base++) {
        total_bases += freq_matrix[base][pos + m_start];
      }
      int n_bases_above_err_freq{0};
      for(int base=0; base < 4; base++) {
        if((freq_matrix[base][pos + m_start] / total_bases) > ALLELIC_FREQ_OF_ERROR) {
          n_bases_above_err_freq++;
        }
      }
      if (n_bases_above_err_freq > 1) { // then mask
        q_str += "L";
        continue;
      }
    }

    // Unique masking logic to cancer reads
    // The supporting evidence for the chosen consensus
    // position must be above 4
    if (tissue == TUMOUR) {
      int base{0};
      switch (cns[pos]) {
        case 'A': base = 0; break;
        case 'T': base = 1; break;
        case 'C': base = 2; break;
        case 'G': base = 3; break;
      }
      if (freq_matrix[base][pos + m_start] < CTR_GSA2) {
        q_str += "M";
        continue;
      }
    }

    q_str += "-";
  }
  return q_str;
}


string BranchPointGroups::addGaps(int n_gaps) {
  string gaps = "";

  for(int i=0; i < n_gaps; i++) {
    gaps += "-";
  }

  return gaps;
}

void BranchPointGroups::invalidatePosition(vector< vector<int> > &alignment_counter, int pos) {
  for (int base=0; base < 4; base++) {
    alignment_counter[base][pos] = -1;
  }
}


void BranchPointGroups::extendBlock(int seed_index, 
    set<read_tag, read_tag_compare> &block, bool orientation) {

  // seed_index is the index of the unique suffix_t this function
  // was called with. 

  int left_of_seed = 0, right_of_seed = 0;

  if (seed_index != 0){
    left_of_seed = ::computeLCP(SA->getElem(seed_index),
                     SA->getElem(seed_index-1), *reads);
  }

  if(seed_index != SA->getSize()-1) {
    right_of_seed = ::computeLCP(SA->getElem(seed_index), 
                     SA->getElem(seed_index+1), *reads);
  }


  if (left_of_seed >= 30 && seed_index != 0){
    getSuffixesFromLeft(seed_index, block, orientation);
  }
  
  if (right_of_seed >= 30 && seed_index != (SA->getSize() -1)) {
    getSuffixesFromRight(seed_index, block, orientation);
  }
}

void BranchPointGroups::getSuffixesFromLeft(int seed_index,
  set<read_tag, read_tag_compare> &block, bool orientation) {

  int left_arrow = seed_index-1;
  // While lexicographally adjacent suffixes share the same lcp value
  // they have the same branchpoint, thus they are in the same group,
  // so add them

  while( // lcps are same AND not out of bounds AND not already in group...
      left_arrow > 0           
      && ::computeLCP(SA->getElem(left_arrow),
                                    SA->getElem(seed_index), *reads) >= 30) {

    // ...add read pointed to by suffix to the block
    // however now add all reads as healthy. 
    // because of the way the set performs comparison, identical reads
    // now labled healthy will be rejected
    read_tag next_read;
    next_read.read_id = SA->getElem(left_arrow).read_id;
    next_read.offset = SA->getElem(left_arrow).offset;
    next_read.orientation = orientation;

    if (SA->getElem(left_arrow).type == HEALTHY) {
      next_read.tissue_type = HEALTHY;
    }
    else {    // we need to switch the type of TUMOUR to SWITCHED
      next_read.tissue_type = SWITCHED;
    }

    // insert tag into block
    pair<std::set<read_tag, read_tag_compare>::iterator, bool> correct_ins = block.insert(next_read);
    left_arrow--;
  }
}

void BranchPointGroups::getSuffixesFromRight(int seed_index,
    set<read_tag, read_tag_compare> &block, bool orientation) {

  int right_arrow = seed_index+1;
  // While lexicographically adjacent suffixes share the same lcp val
  // they have the same branchpoint, thus they are in the same group
  // so add them

  while ( // lcps are the same AND not out of bounds AND not already in group...
      right_arrow < (SA->getSize()-1)         // max LCP size is one less than SA
      && ::computeLCP(SA->getElem(right_arrow), 
                                    SA->getElem(seed_index), *reads) >= 30 ) {

    // ...add read pointed to by suffix to the block
    read_tag next_read;
    next_read.read_id = SA->getElem(right_arrow).read_id;
    next_read.offset = SA->getElem(right_arrow).offset;
    next_read.orientation = orientation;

    if (SA->getElem(right_arrow).type == HEALTHY) {
      next_read.tissue_type = HEALTHY;
    }
    else {    // we need to switch the type of TUMOUR to SWITCHED
      next_read.tissue_type = SWITCHED;
    }

    block.insert(next_read);
    right_arrow++;
  }
}


int BranchPointGroups::minVal(int a, int b) {
  return (a > b) ? b : a;
}

bool BranchPointGroups::lexCompare(string l, string r, unsigned int min_lr) {
  // return true if l < r

  // Generate pointers to lhs and rhs suffixes in reads

  // * min_lr avoids redundant searches droping bound to, in practice
  // O(n + log m)

  string::iterator l_iter  = l.begin() + min_lr;    
  string::iterator l_end   = l.end();
  string::iterator r_iter  = r.begin() + min_lr;
  string::iterator r_end   = r.end();

  for( ; (l_iter != l_end && r_iter != r_end); l_iter++, r_iter++){
    // lex compare character
    if (*l_iter < *r_iter) { return true; }
    if (*r_iter < *l_iter) { return false; }
    // equiv char so move to next...
    }

    // One is prefix of other, return the prefix as higher suffix
    return (l_iter == l_end) && (r_iter != r_end);
}

long long int BranchPointGroups::binarySearch(string query) {
  // Search bounds
  unsigned int right{SA->getSize() - 1};   // start at non-out of bounds
  unsigned int left{0};
  unsigned int mid;

  // prefix lengths with query
  unsigned int min_left_right;
  unsigned int lcp_left_query;
  unsigned int lcp_right_query;

  // find minimum prefix length of left and right bounds with query
  lcp_left_query = lcp(reads->returnSuffix(SA->getElem(left)), query, 0);
  lcp_right_query = lcp(reads->returnSuffix(SA->getElem(right)), query, 0);
  min_left_right = minVal(lcp_left_query, lcp_right_query);

  while (left <= right) {
    bool left_shift{false};
    mid = (left + right) / 2;
    if(lcp(reads->returnSuffix(SA->getElem(mid)), query,  min_left_right) == query.size()) {
      // 30bp stretch covered. Arrived at genomic location. Return.
      return mid; // backUpToFirstMatch(mid, query);
    }
    if(lexCompare(reads->returnSuffix(SA->getElem(mid)), query, min_left_right)) {
      // then query lexicographically lower (indexed > mid) (higher ranked
      // characters) so move left bound towards right
      left = mid+1;
      left_shift = true;
    }
    else {
      // then query is  lexicographically higher (indexed < mid due to lower
      // ranking characters) than mid. Therefore, need to move right bound
      // towards left
       right = mid-1;
    }

    // only recompute the moved bound
    if (left_shift) {
      lcp_left_query  = lcp(reads->returnSuffix(SA->getElem(left)), query, min_left_right);
    }
    else {  // must be right_shift
      lcp_right_query = lcp(reads->returnSuffix(SA->getElem(right)), query, min_left_right);
    }
    min_left_right = minVal(lcp_left_query, lcp_right_query);
  }

  return -1; // no match
}

int BranchPointGroups::lcp(string l, string r, unsigned int mlr) {
  while (l[mlr] == r[mlr]) {
    mlr++;

  }
  return mlr;
}

long long int BranchPointGroups::backUpToFirstMatch(long long int bs_hit, string query) {
  while (bs_hit >= 0) {
    if (lcp(reads->returnSuffix(SA->getElem(bs_hit)), query, 0) != query.size()){
      return bs_hit+1;
    }
    else {
      bs_hit--;
    }
  }
  return bs_hit + 1; // 0
}

unsigned int BranchPointGroups::getSize() {
  return BreakPointBlocks.size();
}



void BranchPointGroups::printBreakPointBlocks() {
  int i=1;
  for(bp_block block : BreakPointBlocks) {
    cout << i << "th block: " << endl;
    cout << "block size" << block.size() << endl;
    for (read_tag tag: block.block) {
      cout << ((tag.tissue_type) ? "Healthy" : "Cancer") 
        << endl << "read_id: " << tag.read_id << endl
        << "offset: " << tag.offset << endl
        << ((tag.orientation) ? "RIGHT" : "LEFT" ) << endl;
      cout << reads->getReadByIndex(tag.read_id, tag.tissue_type) << endl;
    }

    cout << "End of block: -----------------------------------------" << endl;
    i++;
  }
}


//void BranchPointGroups::extractMutationSites() {
//
//  unsigned int elements_per_thread = (SA->getSize()  / N_THREADS);
//  vector<thread> workers;
//
//  unsigned int from=0, to= elements_per_thread;
//  for(int i=0; i < N_THREADS; i++) {
//    workers.push_back(
//      std::thread(&BranchPointGroups::generateBranchPointGroupsWorker, 
//        this, from, to)
//    );
//
//    from = to;
//    if(i == N_THREADS -2) {
//      to = SA->getSize();
//    }
//    else {
//      to += elements_per_thread;
//    }
//  }
//
//  for(auto &thread : workers) {
//      thread.join();
//  }
//}
//
//void BranchPointGroups::generateBranchPointGroupsWorker(unsigned int from, 
//                                                        unsigned int to) {
//
//  set<set<unsigned int>> localThreadExtraction;
//  double cancer_sequences = 0, healthy_sequences = 0;
//  unsigned int extension;
//
//
//  // Make a pass through SA, extracting reads with cancer specific mutations
//  unsigned int seed_index = from;
//  //for (int seed_index = 0; seed_index < SA->getSize()-1; seed_index++) {
//  while (seed_index < to-1) {
////      cout << reads->returnSuffix(SA->getElem(seed_index)) << endl;
//   
//    // If lcp == 30, the reads are from the same genomic location compute
//    // the boundaries for the genomic location, and calc econt to see
//    // if location is cancer specific
//    if(computeLCP(SA->getElem(seed_index),
//                                  SA->getElem(seed_index+1), *reads) >= 30) {
//
//      // tally up the tissue types of seed suffix
//      if (SA->getElem(seed_index).type) {
//        healthy_sequences++;
//      }
//      else {  // cancer
//        cancer_sequences++;
//      }
//
//      // extend branchpoint group
//      extension = seed_index+1;
//      while((computeLCP(SA->getElem(seed_index),
//                                       SA->getElem(extension), *reads) >= 30)
//            ) {
//
//        // tally tissue type of the suffixes in the group
//        if (SA->getElem(extension).type) {
//          healthy_sequences++;
//        }
//        else { // cancer
//          cancer_sequences++;
//        }
//
//        extension++;
//        if(extension == SA->getSize()) {break;}
//      }
//
//
//      // determine if the group is cancer specific
//      // To do this, calculate contamination ration and determine
//      // if it is higher than the user defined level
//      // also, make sure the number of tumour specific reads is 
//      // > 4 (CTR>4)
//
//      //cout << "HealthySequences: " << healthy_sequences << endl;
//      //cout << "cancer_sequences: " << cancer_sequences << endl;
//      //cout << "Ratio: " << healthy_sequences / cancer_sequences << endl;
//      //cout << boolalpha <<  ((healthy_sequences / cancer_sequences)  < econt) <<
//      //  endl;
//      //cout << "ECONT: " << econt << endl;
//
//      if (cancer_sequences >= 4 && 
//         ((healthy_sequences / cancer_sequences) < econt )) { // then make group
//
//
//          // set up break point group
//          set<unsigned int>  block;
//
//          // store all the indicies in the group (pointers to suffixes)
//          unsigned int next_read;
//          for (unsigned int i = seed_index; i < extension; i++) {
//
//            if(SA->getElem(i).type == TUMOUR) {
//              next_read = SA->getElem(i).read_id;
//              block.insert(next_read);   // store reads in block
//            }
//          }
//          localThreadExtraction.insert(block); // store block
//      }
//
//      seed_index = extension; // jump scan to end of group
//    }//end if
//
//    else {  //   there was no match found, so we need to increment seed_index + 1
//      seed_index++;
//    }
//
//    healthy_sequences = 0; // reset
//    cancer_sequences = 0; // reset
//
//  }//end while
//
//  // After extraction, threads need to load their data into global CancerExtractionSet
//
//
//  std::lock_guard<std::mutex> lock(cancer_extraction_lock); // avoid thread interference
//  for(set<unsigned int> mutation_site  : localThreadExtraction) {
//    for(unsigned int read_with_mutation : mutation_site) {
//      CancerExtraction.insert(read_with_mutation);
//    }
//  }
//}

//void BranchPointGroups::printCancerExtractionGroups() {
//  int i=1;
//  for(set<unsigned int> pre_block : *CancerExtraction) {
//    cout << i << "th pre_block: " << endl;
//
//    set<unsigned int>::iterator start = pre_block.begin();
//    set<unsigned int>::iterator end   = pre_block.end();
//
//    while (start != end) {
//      cout << "read: " << *start << endl;
//      start++;
//    }
//
//    cout << "End of pre_block -  --------------------------------------" <<
//      endl;
//    i++;
//  }
//
//}


//void BranchPointGroups::unifyBreakPointBlocks(){ O(N^2)
//  pair<bool, unsigned int> test(false, 100);
//  pair<bool, unsigned int> test1(false, 54);
//  pair<bool, unsigned int> test2(false, 88);
//
//  (*BreakpointBlocks)[10].insert(test);
//  (*BreakpointBlocks)[14].insert(test1);
//  (*BreakpointBlocks)[25].insert(test2);
//  for(set<pair<bool, unsigned int>> mutation_set : *BreakpointBlocks){
//    for(pair<bool, unsigned int> read : mutation_set) {
//      cout << read.first << " --- " << read.second  << endl;
//    }
//  }
//
//
//  for(int i=0; i < BreakpointBlocks->size(); i++) { // compare each block
//    for(int j=i; j < BreakpointBlocks->size(); j++) { // agains all other blocks
//
//      if(i == j){ continue;}    // skip self comparisons
//
//      std::set<pair<bool, unsigned int>>::iterator i_iter, j_iter, i_end, j_end;

//      i_iter = (*BreakpointBlocks)[i].begin();    // point it to start of set
//      i_end  = (*BreakpointBlocks)[i].end();
//
//      j_iter = (*BreakpointBlocks)[j].begin();
//      j_end  = (*BreakpointBlocks)[j].end();
//
//      // quick compare.. if largest element is smaller than smallest in other
//      // no overlap
//      i_end--; j_end--;
//      if( (*(i_end) < *j_iter)   // no elements can be same
//          || 
//          (*(j_end) < *i_iter)   // again no elems can be same
//        ) { continue;}
//
//      i_end++; j_end++;
//
//      // make pass through set, comparing elems
//      while ( (i_iter != i_end) && (j_iter != j_end) ) {
//
//        if (*i_iter > *j_iter) {
//          j_iter++;
//        }
//        else if (*j_iter > *i_iter) {
//          i_iter++;
//        }
//        else { // *i_iter == *j_iter so unify groups
//
//          // unify to i
//          std::set<pair<bool, unsigned int> >::iterator unify;
//          unify = (*BreakpointBlocks)[j].begin();
//          for (;unify != (*BreakpointBlocks)[j].end(); unify++) {
//            (*BreakpointBlocks)[i].insert(*unify);  // add unique elems j to i
//          }
//
//          // unify to j
//          unify = (*BreakpointBlocks)[i].begin();
//          for (;unify != (*BreakpointBlocks)[i].end(); unify++) {
//            (*BreakpointBlocks)[j].insert(*unify);  // add unique elems i to j
//          }
//        }
//
//        i_iter++; j_iter++;
//
//      } // end while
//
//
//    } // end j for 
//  }   // end i for
//
//}
//void BranchPointGroups::extractNonMutatedAllelesDeep() {
//  for(set<read_tag, read_tag_compare> &block : BreakPointBlocks) {
//
//    for(read_tag tag : block) {
//      string read = reads->getReadByIndex(tag.read_id, tag.tissue_type);
//      for(int i=0; i < reads->getReadByIndex(tag.read_id,
//            tag.tissue_type).size()-30; i++) {
//
//        string aligned_section = read.substr(i, 30);
//
//        long long int index = binarySearch(aligned_section);
//        if(index == -1) {
//          cout << "Error, a previously identified sequence no longer exists" <<
//            " in suffix array... program terminating. "  << endl;
//          exit(1);
//        }
//        else {
//          extendBlock(index, block, tag.orientation);
//        }
//      }
//    }
//  }
//}






//bool hitDollar(Suffix_t &lhs, Suffix_t &rhs) {
//
//  // Already know that the 30 positions match, so start computation from 
//  // 29th position in string
//  string::iterator lhs_string_iter = returnStartIterator(&lhs, *reads) + 29;
//  string::iterator rhs_string_iter = returnStartIterator(&rhs, *reads) + 29;
//
//  string::iterator lhs_end =         returnEndIterator(&lhs, *reads); 
//  string::iterator rhs_end =         returnEndIterator(&rhs, *reads);
//
//  // extend further and compare
//  while (*lhs_string_iter  == *rhs_string_iter) {
//    lhs_string_iter++; rhs_string_iter++;
//  }
//
//  // hit a mismatch, check if the mismatch is a dollar symbol
//  if(*lhs_string_iter == '$' || *rhs_string_iter == '$') {
//    return true;
//  }
//
//  // else
//  return false;
//}
//
//
//


//void BranchPointGroups::generateBranchpointGroupa() {
//
//  // Pass through SA locating unique sequences
//  for(int i=0; i < SA->getSize(); i++) {
//
//    if (!SA->getUnique(i)) { // not common, so unique
//
//      cout << "CALLED TOO OFTEN " << i << endl;
//      vector<unsigned int> newgroup;
//      newgroup.push_back(i);
//
//      // extend left and/or right of the unique suffix
//      // gathering the branchpoint group
//      makeReadGroup(i, newgroup);
//      BPG->push_back(newgroup); // store the new group
//    }
//  }
//}



// end of file






//void BranchPointGroups::mergeHash() {
//
//  struct hashtag{
//    unsigned long long int hash_id;
//    unsigned int group;
//    unsigned int read_id;
//    bool orienation;
//  };
//
//  unsigned long long int gmin, gmax;
//  gmin = std::numerical_limits<unsigned long long int>::max();
//  gmax = 0;
//
//  std::hash<std::string> h;
//  vector<vector<hashtag>> hash_list;
//  hashtag next_forward_hash;
//  hashtag next_reverse_hash;
//  for(int i=0; i < ComplementaryUnified->size(); i++){         // M
//    set<unsigned int> group = (*ComplementaryUnified)[i];    
//    vector<hashtag> group_hashes;
//
//    for(unsigned int read : group) {                          // M*C
//
//      string read_string = reads->getReadByIndex(read, TUMOUR);  
//      for(int i=0; i < read_string.size()-30; i++) {            // M*C*n
//        next_forward_hash.hash_id = h(read_string.substr(i, 30));
//        next_forward_hash.group = i;
//        next_forward_hash.orientation = true;
//        next_forward_hash.read_id = read;
//      }
//
//      string reverse_read = reverseComplementString(read_string);
//      for(int i=0; i < reverse_read.size()-30; i++) {
//        next_reverse_hash.hash_id = h(reverse_read.substr(i, 30));
//        next_reverse_hash.group = i;
//        next_reverse_hash.orientation = false;
//        next_reverse_hash.read_id = read;
//      }
//
//      group_hashes.push_back(next_forward_hash);
//      group_hashes.push_back(next_reverse_hash);
//    }
//
//    hash_list.push_back(group_hashes);
//  }
//
//
//  // find global min max
//  for(hashtag hashed_substring : hash_list) {      // M.C.n
//    if(hashed_substring.hash_id < gmin) {
//      gmin = hashed_substring.hash_id;
//    }
//    if(hashed_substring.hash_id > gmax) {
//      gmax = hashed_substring.hash_id;
//    }
//  }
//
//  struct group_unif_info {
//    unsigned int group;
//    unsigned int prev_unification_group;
//    bool unified;
//  };
//
//
//  // make unification table
//  vector<group_unif_info> unif_table;
//  unif_table.reserve(CompelmentaryUnified->size());
//  for(unsigned int i=0; i < ComplemenratyUnified->size(); i++) {
//    group_unif_unfo info;
//    info.group = i;
//    info.unified = false;
//    unif_table.push_back(info);
//  }
//
//
//  // Make tagArray
//
//  typedef pair<vector<group_unif_info>::iterator, vector<hashtag>::iterator> 
//  lut_hash_pointer;
//
//  lut_hash_pointer next_ptr;
//  vector<vector<lut_hash_pointer>> tagArray;    // tagArrau
//
//  for(unsigned int i=0; i < ComplemenratyUnified->size(); i++) {
//    for(int j=0; j < hash_list[i].size(); j++) {
//      lut_hash_pointer.first = unif_table.begin() + i;
//      lut_hash_pointer.second = hash_list[i].begin() + j;
//
//      tagArray[(*lut_hash_pointer.second).hash_id] = lut_hash_pointer;
//    }
//  }
//
//
//
//  // begin unif
//
//  for(int i; i < (gmax - gmin); i++) {
//
//    if(tagArray[i].size() ==  1) {continue;}
//
//    bool forward = false, reverse = false;
//    for(lut_hash_pointer tag: tagArray[i]) {
//      if(tag.second->orientation = true) {
//        forward = true;
//      }
//      if(tag.second->orientation = false) {
//        reverse = true;
//      }
//    }
//
//    // group not single orientation
//    if (forward == false ||  reverse == false) {continue;}
//
//    // all reads in should share 30bp in common, and if no collisions
//    // counter should thus equal size()-1
//    vector<unsigned int> collisions;
//
//    // for every tag in hashpoint, try to merge to base
//    for(int j=0; j< tagArray[i].size(); j++) {
//
//        if(tagArray[i][0].second->orientation == 
//           tagArray[i][j].second->orientation) {counter++; continue;}
//
//        // collision check, i.e validate with actuall slow test
//        // should be an exact match!! so change for seq matcher
//        if(thirtyBasePairOverlap( tagArray[i][0].second->read_id, 
//                                  tagArray[i][j].second->read_id) ) {
//          // merge
//
//          // if the base group has already been unified with group x, then
//          // unify all in list, with group x
//          int unify_upon_this_group;
//
//
//          // check unif_table to see if the base group has been unified...
//          if(tagArray[i][0].first->unified) {  // if so, unify to that group
//            unify_upon_this_group = tagArray[i][0].first->prev_unification_group;  
//          }
//          else {  // if not, unify to the base
//            unify_upon_this_group = tagArray[i][0].first->group; 
//          }
//
//          // begin merging
//        
//          // point to group to merge
//          set<unsigned int>::iterator j_it = 
//          (*ComplementaryUnfied)[tagArray[i][j].first->group].begin();
//
//          set<unsigned int::iterator j_end = 
//          (*ComplementaryUnfied)[tagArray[i][j].first->group].end();
//
//          // merge with base group
//          while (j_it != j_end) {
//            (*CompelementaryUnfied)[unify_upon_this_group]->insert(*j_it);  
//          }
//
//        }
//        else {
//          collisions.push_back(j);
//        }
//      }
//      for(int collision : collisions) {
//        cout << "Collision with: " << tagArray[i][collision].second->read_id <<
//        endl;
//      }
//    
//  }
//}

//void BranchPointGroups::discardSingleOrientationBlocks() {
//
//  bool right_pair = false, left_pair = false;
//
//  vector<set<pair<bool, unsigned int>>>::iterator block = BreakpointBlocks->begin();
//
//
//  while(block != BreakpointBlocks->end()) {
//
//
//    for (pair<bool, unsigned int> read : *block) {
//
//      // check if the block contains two different read orientations
//      if ((*reads->TumourMateOrder)[read.second]) {
//        right_pair = true;
//      }
//      else {
//        left_pair = true;
//      }
//    }
//
//    // delete if group has one orientation
//    if (right_pair == false || left_pair == false) {
//      BreakpointBlocks->erase(block);
//    }
//    else {
//      block++;
//    }
//
//    right_pair = false; left_pair = false;    // reset
//
//  }
//}

//void BranchPointGroups::makeHashGroups() {
//  enum {LEFT, RIGHT};      // mate pair enum
//
//  struct hashtag {
//    unsigned int hash_id;
//    unsigned int read_id;
//    unsigned int offset;
//    bool orientation;
//    bool inblock;
//  };
//
//
//  // Make all the hashes...
//
//  // this will contain all 100 hashes from each extracted read and its rev comp
//  vector<hashtag> hash_list;
//  hash_list.reserve(CancerExtraction->size() * 30 * 100);
//
//  std::hash<std::string> h; // the glorious h() function
//
//  // the forward and reverse hashes for each string
//  hashtag fwd_hash;
//  hashtag rev_hash;
//  
//  // for each read from each extracted group..
//  for(set<unsigned int> mutation_site: *CancerExtraction) {
//    for (unsigned int read_index : mutation_site) {  
//    // .. hash all reads and rev comps by sliding a 30pb window
//
//      string read = reads->getReadByIndex(read_index, TUMOUR);
//      string rev_read = reverseComplementString(read);
//
//      for(int i=0; i < read.size()-30; i++) {
//
//        // iterate a hash over a 30bp window
//        fwd_hash.hash_id = h(read.substr(i, 30)); 
//        rev_hash.hash_id = h(rev_read.substr(i, 30));
//
//        fwd_hash.read_id = read_index;
//        rev_hash.read_id = read_index;
//
//        fwd_hash.offset = i;
//        rev_hash.offset = i;
//
//        fwd_hash.inblock= false;
//        rev_hash.inblock= false;
//
//        fwd_hash.orientation = RIGHT;
//        rev_hash.orientation = LEFT;    // call LEFT reverse
//
//        hash_list.push_back(fwd_hash);
//        hash_list.push_back(rev_hash);
//
//      }
//    }
//  }
//
//
//  hash_list.shrink_to_fit();      // neaten up vector
//
//
//  // find global min and max hash value...
//  unsigned int gmin, gmax;
//  gmin = std::numeric_limits<unsigned int>::max();
//  gmax = 0;
//
//  for(hashtag hashed_substring : hash_list) {      // M.C.n
//    if(hashed_substring.hash_id < gmin) {
//      gmin = hashed_substring.hash_id;
//    }
//    if(hashed_substring.hash_id > gmax) {
//      gmax = hashed_substring.hash_id;
//    }
//  }
//
//  // using its hashvalue as an index, load every hashtag into an array
//  // common 30bp seqs from different reads will have hashed to the same value
//
//  cout << "make past here" << endl;
//  cout << "gmax - gmin" << gmax - gmin << endl;
//  vector<vector<hashtag>> commonReadsArray(gmax+1-gmin);
//  cout << commonReadsArray.size();
//  for (hashtag tag : hash_list) {
//    cout << tag.hash_id << endl;
//    commonReadsArray[tag.hash_id - gmin].push_back(tag);  // place tag in array
//  }
//
//  for(vector<hashtag> commonReads: commonReadsArray) {
//    if (commonReads.size() < 2) { continue; } 
//
//
//    set<set<unsigned int>> blocks;
//    set<unsigned int> block;
//
//    bool block_made = true;
//    bool left = false;
//    bool right = false;
//
//    for(int i=0; i < commonReads.size(); i++) {
//      if(commonReads[i].inblock) {continue;}
//
//      block.insert(commonReads[i].read_id);
//      commonReads[i].inblock = true;
//
//      for(int j=0; j < commonReads.size(); j++) {
//        string block_seed = 
//        reads->getReadByIndex(commonReads[0].read_id,
//        TUMOUR).substr(commonReads[0].offset, 30);
//
//        string potential_member = 
//        reads->getReadByIndex(commonReads[i].read_id,
//        TUMOUR).substr(commonReads[i].offset, 30);
//
//
//        // sequences should match, if not, its a collision, this is 
//        // not a member, so continue
//        if (!sequenceMatch(block_seed, potential_member)) { continue; }
//
//        // group!!
//        block_made = true;
//        if(commonReads[i].orientation == RIGHT) {
//          right = true;
//        }
//        else {
//          left = true;
//        }
//
//        if(commonReads[j].orientation== RIGHT) {
//          right = true;
//        }
//        else {
//          left = true;
//        }
//        block.insert(commonReads[j].read_id);
//        commonReads[j].inblock = true;
//      }
//  
//      if(left && right && block_made) {
//        blocks.insert(block);
//      }
//      
//      // reset
//      block.clear();
//      block_made = left = right = false;
//    }
//
//    for (set<unsigned int> block : blocks) {
//      if(block.size() >= 4) {
//        HashedPreBlocks->push_back(block);
//      }
//    }
//  }
//
//}

//bool BranchPointGroups::sequenceMatch(string right, string left) {
//
//  for(int i=0; i < right.size(); i++) {
//    if(right[i] != left[i]) {return false;}
//  }
//
//  return true;
//}



//void BranchPointGroups::unifyBreakpointBlocks() {
//  unsigned int gmin, gmax;
//  gmin = reads->n_tumour_reads;
//  gmax = 0;
//
//  // scan all blocks and identify global min and max
//  for(set<unsigned int> block : *HashedPreBlocks) {
//    set<unsigned int>::iterator it, end;
//    it  = block.begin();
//    end = block.end();
//
//    for (; it != end; it++) {
//      if (*it < gmin) {gmin = *it;}   // update global minimum
//      if (*it > gmax) {gmax = *it;}   // update global max
//    }
//  }
//
//  // unif_table inst.
//  // unif_table keeps track of which groups have been unified to which others
//
//  struct block_unif_info {    // the elements of the table
//    unsigned int id;          // the id of the group
//    unsigned int prev_unification_group;  // what group it was unified to
//    bool unified;             // whether the group has been unified
//  };
//
//  vector<block_unif_info> unif_table;
//  unif_table.reserve(HashedPreBlocks->size());
//  //cout << "unif table size: " << CancerExtraction->size() << endl;
//
//  for(unsigned int i=0; i < HashedPreBlocks->size(); i++) {
//    block_unif_info b;
//    b.id = i;
//    b.unified = false;
//
//    unif_table.push_back(b);          // intialize table
//  }
//
//
//
//  // make tagArray:
//  // tagArray is a 2d array of the range gmin-gmax. 
//  // The idea is that, each block id "tag" is added to this 2d array, 
//  // such that, the read_id is used as an index. 
//  // For example, given block_1 = {1, 2, 4, 8}. A tag to block_1
//  // will be placed at positions 1, 2, 4, 8 in tag array. 
//  // the tag is accompanied by a pointer to other blocks, and whether the block
//  // has been merged. This allows effective merging
//
//  // decl. tagArray as an array of vectors of pointers that will be used
//  // to look up the unification status of a group in unif_table
//
//  bool unification_occured = false;
//  do {
//    vector<vector<vector<block_unif_info>::iterator>> tagArray(gmax+1 -gmin);
//
//    //cout << "unif_table status: " << endl;
//    //for (block_unif_info block: unif_table) {
//    //  cout << block.id << endl;
//    //  cout << block.unified << endl;
//    //  cout << block.prev_unification_group << endl;
//    //}
//    unification_occured = false;    // reset after loop
//
//    // load read_id's
//    for (unsigned int tag=0; tag < HashedPreBlocks->size(); tag++) {
//      if (unif_table[tag].unified) {continue;}  // dont add info if merged
//
//      // generate pointer to block
//      vector<block_unif_info>::iterator ptr_to_block_info = 
//        unif_table.begin() + tag;
//
//
//      // set up iterators to next block, to read the read_id's
//      set<unsigned int>::iterator read_id, end;
//      read_id = (*HashedPreBlocks)[tag].begin();
//      end = (*HashedPreBlocks)[tag].end();
//
//      for (; read_id != end; read_id++) {   // for each read in block
//
//        tagArray[*read_id- gmin].push_back(ptr_to_block_info);     
//        // load the pointer to the group in
//        // unif_table into the vector att he read_id'th index in tagArray
//      }
//    }
//    //for(vector<vector<block_unif_info>::iterator> v : tagArray) {
//    //  cout << "Tag elem size: " << v.size() << endl;
//    //}
//
//    // begin unification
//    for (unsigned int i=0; i < (gmax+1-gmin); i++) {
//
//
//      if (tagArray[i].size() < 2) {
//          continue;
//      } // nothing in list to unify
//
//      else {
//
//        // if the base group has already been unified with group x, then
//        // unify all in list, with group x
//        int unify_upon_this_group;
//
//        // make syntax easier on the eye ;-)
//
//        // check unif_table to see if the base group has been unified...
//        if(tagArray[i][0]->unified) {  // if so, unify to that group
//          unify_upon_this_group = tagArray[i][0]->prev_unification_group;  
//        }
//        else {  // if not, unify to the base
//          unify_upon_this_group = tagArray[i][0]->id; 
//        }
//
//        // working backwards, unify the groups in block_tags, to 
//        // unify_upon_this_group
//        for (int k=tagArray[i].size()-1; k > 0; k--) {
//
//          if (tagArray[i][k]->unified == true) {continue;} // already merged!
//
//          unification_occured = true;   // going to unify groups
//
//
//          set<unsigned int>::iterator it = 
//            (*HashedPreBlocks)[tagArray[i][k]->id].begin(); // start of group
//
//          set<unsigned int>::iterator end = 
//            (*HashedPreBlocks)[tagArray[i][k]->id].end(); // end of group
//
//          while(it != end) {
//            (*HashedPreBlocks)[unify_upon_this_group].insert(*it); it++;
//          }
//          tagArray[i][k]->unified = true;
//          tagArray[i][k]->prev_unification_group = unify_upon_this_group;
//        }
//      }
//    }
//  } while(unification_occured);
//
//
//  // now, just extract merged groups and load to ComplementryUnified
//  for (unsigned int i=0; i < unif_table.size(); i++) {
//    if (unif_table[i].unified == false) {
//
//      // pass the group into ComplementaryUnified
//      FinalUnify->push_back((*HashedPreBlocks)[i]);
//    }
//  }
//
//  // done with this data as processed into ComplementyUnfied
//  //cout << "PRINTING NUMERICAL GROUPS" << endl;
//  //printCancerExtractionGroups();
//  delete CancerExtraction;
//  CancerExtraction = nullptr;
//
//  // finally, merging has completed. So, load the merged reads into 
//  // ComplementaryGroups
//  //for (int i=0; i < unif_table.size(); i++) {
//
//  //  if(unif_table[i].unified == false) { // only add groups that were not unified
//
//  //    // need to convert each read from each group into 
//  //    // a tissue type, read_id pair, of a "block"
//
//  //    set<unsigned int>::iterator it = (*CancerExtraction)[i].begin();
//  //    set<unsigned int>::iterator end = (*CancerExtraction)[i].end();
//
//  //    vector<pair<bool, unsigned int> > block;
//  //    pair<bool, unsigned int> next_read;   
//  //    while (it != end) {
//  //      next_read.first = TUMOUR;
//  //      next_read.second = *it;
//  //      block.push_back(next_read);
//
//  //      it++;
//  //    }
//
//  //    ComplementaryUnfied->push_back(block);
//  //  }
//  //}
//
//}

//void BranchPointGroups::unifyComplementaryGroups() {
//  unsigned int gmin, gmax;
//  gmin = reads->n_tumour_reads;
//  gmax = 0;
//
//  // scan all blocks and identify global min and max
//  for(set<unsigned int> block : *CancerExtraction) {
//    set<unsigned int>::iterator it, end;
//    it  = block.begin();
//    end = block.end();
//
//    for (; it != end; it++) {
//      if (*it < gmin) {gmin = *it;}   // update global minimum
//      if (*it > gmax) {gmax = *it;}   // update global max
//    }
//  }
//
//  // unif_table inst.
//  // unif_table keeps track of which groups have been unified to which others
//
//  struct block_unif_info {    // the elements of the table
//    unsigned int id;          // the id of the group
//    unsigned int prev_unification_group;  // what group it was unified to
//    bool unified;             // whether the group has been unified
//  };
//
//  vector<block_unif_info> unif_table;
//  unif_table.reserve(CancerExtraction->size());
//  //cout << "unif table size: " << CancerExtraction->size() << endl;
//
//  for(unsigned int i=0; i < CancerExtraction->size(); i++) {
//    block_unif_info b;
//    b.id = i;
//    b.unified = false;
//
//    unif_table.push_back(b);          // intialize table
//  }
//
//
//
//  // make tagArray:
//  // tagArray is a 2d array of the range gmin-gmax. 
//  // The idea is that, each block id "tag" is added to this 2d array, 
//  // such that, the read_id is used as an index. 
//  // For example, given block_1 = {1, 2, 4, 8}. A tag to block_1
//  // will be placed at positions 1, 2, 4, 8 in tag array. 
//  // the tag is accompanied by a pointer to other blocks, and whether the block
//  // has been merged. This allows effective merging
//
//  // decl. tagArray as an array of vectors of pointers that will be used
//  // to look up the unification status of a group in unif_table
//
//  bool unification_occured = false;
//  do {
//    vector<vector<vector<block_unif_info>::iterator>> tagArray(gmax+1 -gmin);
//
//    //cout << "unif_table status: " << endl;
//    //for (block_unif_info block: unif_table) {
//    //  cout << block.id << endl;
//    //  cout << block.unified << endl;
//    //  cout << block.prev_unification_group << endl;
//    //}
//    unification_occured = false;    // reset after loop
//
//    // load read_id's
//    for (unsigned int tag=0; tag < CancerExtraction->size(); tag++) {
//      if (unif_table[tag].unified) {continue;}  // dont add info if merged
//
//      // generate pointer to block
//      vector<block_unif_info>::iterator ptr_to_block_info = 
//        unif_table.begin() + tag;
//
//
//      // set up iterators to next block, to read the read_id's
//      set<unsigned int>::iterator read_id, end;
//      read_id = (*CancerExtraction)[tag].begin();
//      end = (*CancerExtraction)[tag].end();
//
//      for (; read_id != end; read_id++) {   // for each read in block
//
//        tagArray[*read_id- gmin].push_back(ptr_to_block_info);     
//        // load the pointer to the group in
//        // unif_table into the vector att he read_id'th index in tagArray
//      }
//    }
//    //for(vector<vector<block_unif_info>::iterator> v : tagArray) {
//    //  cout << "Tag elem size: " << v.size() << endl;
//    //}
//
//    // begin unification
//    for (unsigned int i=0; i < (gmax+1-gmin); i++) {
//
//
//      if (tagArray[i].size() < 2) {
//          continue;
//      } // nothing in list to unify
//
//      else {
//
//        // if the base group has already been unified with group x, then
//        // unify all in list, with group x
//        int unify_upon_this_group;
//
//        // make syntax easier on the eye ;-)
//
//        // check unif_table to see if the base group has been unified...
//        if(tagArray[i][0]->unified) {  // if so, unify to that group
//          unify_upon_this_group = tagArray[i][0]->prev_unification_group;  
//        }
//        else {  // if not, unify to the base
//          unify_upon_this_group = tagArray[i][0]->id; 
//        }
//
//        // working backwards, unify the groups in block_tags, to 
//        // unify_upon_this_group
//        for (int k=tagArray[i].size()-1; k > 0; k--) {
//
//          if (tagArray[i][k]->unified == true) {continue;} // already merged!
//
//          unification_occured = true;   // going to unify groups
//
//
//          set<unsigned int>::iterator it = 
//            (*CancerExtraction)[tagArray[i][k]->id].begin(); // start of group
//
//          set<unsigned int>::iterator end = 
//            (*CancerExtraction)[tagArray[i][k]->id].end(); // end of group
//
//          while(it != end) {
//            (*CancerExtraction)[unify_upon_this_group].insert(*it); it++;
//          }
//          tagArray[i][k]->unified = true;
//          tagArray[i][k]->prev_unification_group = unify_upon_this_group;
//        }
//      }
//    }
//  } while(unification_occured);
//
//
//  // now, just extract merged groups and load to ComplementryUnified
//  for (unsigned int i=0; i < unif_table.size(); i++) {
//    if (unif_table[i].unified == false) {
//
//      // pass the group into ComplementaryUnified
//      ComplementaryUnfied->push_back((*CancerExtraction)[i]);
//    }
//  }
//
//  // done with this data as processed into ComplementyUnfied
//  //cout << "PRINTING NUMERICAL GROUPS" << endl;
//  //printCancerExtractionGroups();
//  delete CancerExtraction;
//  CancerExtraction = nullptr;
//
//  // finally, merging has completed. So, load the merged reads into 
//  // ComplementaryGroups
//  //for (int i=0; i < unif_table.size(); i++) {
//
//  //  if(unif_table[i].unified == false) { // only add groups that were not unified
//
//  //    // need to convert each read from each group into 
//  //    // a tissue type, read_id pair, of a "block"
//
//  //    set<unsigned int>::iterator it = (*CancerExtraction)[i].begin();
//  //    set<unsigned int>::iterator end = (*CancerExtraction)[i].end();
//
//  //    vector<pair<bool, unsigned int> > block;
//  //    pair<bool, unsigned int> next_read;   
//  //    while (it != end) {
//  //      next_read.first = TUMOUR;
//  //      next_read.second = *it;
//  //      block.push_back(next_read);
//
//  //      it++;
//  //    }
//
//  //    ComplementaryUnfied->push_back(block);
//  //  }
//  //}
//
//}
//
//
//
//void BranchPointGroups::unifyReverseComplementaryReads() {
//
//  for(unsigned int i = 0; i < ComplementaryUnfied->size(); i++ ) {
//
//    // if the first read of the ith block of ComplementryUnified 
//    // is oposite pair of the first read from the jth block of ComplementryUnified
//
//    for(unsigned int j = i; j < ComplementaryUnfied->size(); j++) {
//      if (i == j) {continue;} // dont self compare
//
//      // if pair type of j and i is the same, we dont want to compare them
//      std::set<unsigned int>::iterator i_it, j_it;
//      i_it = (*ComplementaryUnfied)[i].begin();
//      j_it = (*ComplementaryUnfied)[j].begin();
//
//      if( (*reads->TumourMateOrder)[*i_it] ==   // skip if common mate 
//          (*reads->TumourMateOrder)[*j_it]) {
//        continue;
//      }
//
//
//      // groups are opposite pair type, so need rev comp comparison
//
//      if(thirtyBasePairOverlap(*i_it, *j_it)) {
//
//        // unify groups
//        // unify to i
//        for (;j_it != (*ComplementaryUnfied)[j].end(); j_it++) {
//          (*ComplementaryUnfied)[i].insert(*j_it);  // add unique elems j to i
//        }
//
//        // unify to j
//        for (;i_it != (*ComplementaryUnfied)[i].end(); i_it++) {
//          (*ComplementaryUnfied)[j].insert(*i_it);  // add unique elems i to j
//        }
//      }
//    }
//  }
//
//
//  // one compared reverse complementary sequences, make unique groups
//  // by passing the groups through a set
//
//  set<set<unsigned int>> filter_uniques;
//  for(unsigned int i=0; i < ComplementaryUnfied->size(); i++) {
//    filter_uniques.insert((*ComplementaryUnfied)[i]);
//  }
//
//  cout << "FILTER UNIQUES" << endl;
//  for(set<unsigned int> s : filter_uniques) {
//    for(unsigned int n : s) {
//      cout << n << ", ";
//    }
//    cout << endl;
////  }
//
//
//
//
//  // after passed through set groups unique, so place into vector
//
//  // pass through each group
//  for(set<unsigned int> group : filter_uniques) {
//
//    set<pair<bool, unsigned int>> block;
//    pair<bool, unsigned int> tumour_read;
//
//    // load each read into next block
//    for(unsigned int read_index : group) {
//      tumour_read.first = TUMOUR;
//      tumour_read.second = read_index;
//
//      block.insert(tumour_read);
//    }
//    // load block into BreakPointBlocks 
//    BreakpointBlocks->push_back(block);
//  }
//}


//bool BranchPointGroups::thirtyBasePairOverlap(unsigned int lhs, 
//                           unsigned int rhs) {
//
//  // string as rev comp of each other, by getting the rev comp of 
//  // one, the two string should represent the same strand, and
//  // will have a match of 30bp if they share the same genomic location
//
//  string rhs_read_rev_comp = reverseComplementString(      // rev comp of right
//                     reads->getReadByIndex(rhs, TUMOUR)
//                     );
//  string lhs_read = reads->getReadByIndex(lhs, TUMOUR);
//  
//
//
//  // window along 30bp along sequence looking for match
//  for(unsigned int i=0; i < lhs_read.size() - 30; i++) {
//    string lhs_sub = lhs_read.substr(i, 30);
//
//    // if we find a 30bp match...
//    if(rhs_read_rev_comp.find(lhs_sub) != std::string::npos) {
//      return true;    // report match found
//    }
//  }
//
//
//  return false;   // didnt find a match
//}
