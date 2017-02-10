#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <string>
#include <mutex>
#include <thread>
#include <algorithm>
#include <cstdlib>


// reading tar.gz
#include <zlib.h>
#include "kseq.h"

// sa construction
#include "radix.h"
#include <sys/types.h>
#include <cstring>
#include <cstdlib>

#include "string.h"
#include "build_suffix_array.h" 
std::mutex quality_processing_lock;

using namespace std;

KSEQ_INIT(gzFile, gzread);
static const int NUM_ARGS = 4;
static const int HEADER_FILE_IDX = 1;
static const int ECONT_IDX = 2;
static const int OFILE_IDX = 3;
static const int TERM_CHAR_CORRECTION = 1;
static const int N_THREADS = 1;
static const int MIN_SUFFIX_SIZE = 30;  // remove user definable variable

static const double QUALITY_THRESH = 0.1; // 10% 
static const char PHRED_20 = '5';   // lowest high quality phred score

static const string REMOVED_TOKENS = "N"; // remove N from fastq
static const string TERM_CHAR = "$";      // suffix termination character

enum {MUTANT, NON_MUTANT};

// main

//int main(int argc, char **argv) {
//  if (argc != 2) {
//    cout << "Usage: <exe> <datafilelist>" << endl;
//    return -1;
//  }
//
//  vector<string> cancerFileNames, healthyFileNames;
//  splitFileNamesOnDataType(healthyFileNames,cancerFileNames, argv[1]);
//
//  vector<string> cancerReads, healthyReads;
//  vector<gsaTuple> cancerGSA, healthyGSA;
//  if (!healthyFileNames.empty()) {
//    healthyGSA = buildGSA(healthyFileNames, healthyReads);
//  }
//  if (!cancerFileNames.empty()) {
//    cancerGSA  = buildGSA(cancerFileNames, cancerReads);
//  }
//
//  cout << "Size of cancer gsa: " << cancerGSA.size() << endl;
//  cout << "Size of cancer reads: " << cancerReads.size() << endl;
//
//  for (gsaTuple const& tup : cancerGSA) {
//    cout << cancerReads[tup.read_idx].substr(tup.offset) << endl;
//  }
//  
//  string read = "HelloMrPostmanLookAtMeWoahYeah";
//  set<unsigned int> found_reads = findReadsCoveringLocation(cancerReads,
//      cancerGSA, read);
//  cout << endl << endl << "Identfied reads: " << endl;
//  for (unsigned int read_idx : found_reads) {
//    cout << cancerReads[read_idx] << endl;
//  }
//}

string rc(string const& s) {
  string rc_s = "";
  for(int i = s.size()-1; i >=0; i--) {
    switch (s[i]) {
      case 'A': rc_s.push_back('T'); break;
      case 'T': rc_s.push_back('A'); break;
      case 'C': rc_s.push_back('G'); break;
      case 'G': rc_s.push_back('C'); break;
    }
  }
  return rc_s;
}
vector<snippetData>
extractReadsCoveringSnippets(vector<coordinateData> const& coords,
                             vector<string> const& reads,
                             vector<gsaTuple> const& gsa,
                             TissueType tissue) {
  vector<snippetData> results;
  for (coordinateData const& coord : coords) {
    // set up snippets
    set<gsaTuple, gsaTupleCompare> fwdNonMutated, revNonMutated, fwdMutated, revMutated;
    string fwd_non_mut(coord.sequence), fwd_mut(coord.sequence);
    fwd_mut[coord.flanking_dist] = coord.cBase;
    string rev_non_mut(rc(fwd_non_mut)), rev_mut(rc(fwd_mut));


    // search for reads 
    fwdNonMutated =  
      findReadsCoveringLocation(reads, gsa, fwd_non_mut, Orientation::fwd, Covers::nonMut);
    fwdMutated    =  
      findReadsCoveringLocation(reads, gsa, fwd_mut,     Orientation::fwd, Covers::mut);
    revNonMutated =  
      findReadsCoveringLocation(reads, gsa, rev_non_mut, Orientation::rev, Covers::nonMut);
    revMutated    =  
      findReadsCoveringLocation(reads, gsa, rev_mut,     Orientation::rev, Covers::mut);

    // load data int snippetCoord
    snippetData entry;
    entry.header           = coord.header; 
    entry.snippet          = coord.sequence;
    entry.healthy          = coord.hBase; 
    entry.cancer           = coord.cBase;
    entry.mutationLocation = coord.coordinate; 
    entry.tissue           = tissue;
    entry.fd               = coord.flanking_dist;

    // first insert the reads covering non mutated sequence in both orientations
    entry.coveringReads.insert(fwdNonMutated.begin(), fwdNonMutated.end())
    entry.coveringReads.insert(revNonMutated.begin(), revNonMutated.end())

    // Then insert reads covering mutated sequence. Note, most of these
    // will be rejected, as only a reads uniquely cover the mutated sequence
    entry.coveringReads.insert(fwdMutated.begin(), fwdMutated.end())
    entry.coveringReads.insert(revMutated.begin(), revMutated.end())

    // done filling entry data
    results.push_back(entry);
  }
  return results;
}

void printSnippetData(ostream & out, vector<snippetData> const& data) {
  for (snippetData const& entry : data) {
    out << "Header: "  << entry.header 
              << " - " << entry.healthy << ":" << entry.cancer
              << " - " << ((entry.tissue == TissueType::healthy) ? "H":"T")
              << endl;

    out << "Snippt: "  << entry.snippet << endl;
    out << "Coord : "  << entry.mutationLocation << endl;
    out << "Reads : " << endl;
    for (int i = 0; i < entry.reads.size(); i++) {
      out << entry.reads[i] << " :: "
                << entry.read_idx[i]  << " :: "
                << ((entry.found_on[i]) ? "N" : "M") << endl;
    }
  }
}

void printAlignedRead(string const& snippet, string const& read, 
    int read_idx

set<gsaTuple, gsaTupleCompare> 
findReadsCoveringLocation(vector<string> const& reads, 
  vector<gsaTuple> const& gsa, string const& query, Orientation ori,
  Covers cov) {

  set<gsaTuple, compareGSATuple> readsCoveringLocation;
  for(int i=0; i < query.size - 30; i++) {
    string querySubstr = query.substr(i, 30);
    vector<gsaTuple>::const_iterator result = binarySearch(reads, gsa,
        querySubstr);

    if (result != gsa.end()) {

      readsCoveringLocation.insert(*result); // input the hit element

      // search back
      vector<gsaTuple>::const_iterator leftArrow = result-1;
      if (leftArrow >= gsa.begin()) {
        while (
            lcp(reads[leftArrow->read_idx]
              .substr(leftArrow->offset), querySubstr)
              >= MIN_SUFFIX_SIZE) {

          // try and insert into set
          gsaTuple entry(leftArrow->read_idx, leftArrow->offset, ori, cov, i);
          readsCoveringLocation.insert(entry);
          leftArrow--;
          if (leftArrow < gsa.begin()) break;
        }
      }
      // search forward
      vector<gsaTuple>::const_iterator rightArrow = result + 1;
      if (rightArrow < gsa.end()) {
        while (lcp(reads[rightArrow->read_idx]
              .substr(rightArrow->offset), querySubstr)
              >= MIN_SUFFIX_SIZE) {
          
          // try and insert into set
          gsaTuple entry(rightArrow->read_idx, rightArrow->offset, ori, cov, i);
          readsCoveringLocation.insert(entry);
          rightArrow++;
          if (rightArrow >= gsa.end()) break;
        }
      }
    }
  }
  return readsCoveringLocation;
}

//set<unsigned int> findReadsCoveringLocation(vector<string> const& reads,
//    vector<gsaTuple> const& gsa, string const& query) {

//  set<unsigned int> readsCoveringLocation;
//
//  for (int i=0; i <= query.size() - 30; i++) {
//    string querySubstr = query.substr(i, 30);
//    vector<gsaTuple>::const_iterator result = binarySearch(reads, gsa, querySubstr);
//    if (result != gsa.end()) {
//
//      // input result
//      readsCoveringLocation.insert(result->read_idx);
//      // search back
//      vector<gsaTuple>::const_iterator leftArrow = result-1;
//      if (leftArrow >= gsa.begin()) {
//        while (
//            lcp(reads[leftArrow->read_idx]
//              .substr(leftArrow->offset), querySubstr)
//              >= MIN_SUFFIX_SIZE) {
//          readsCoveringLocation.insert(leftArrow->read_idx);
//          leftArrow--;
//          if (leftArrow < gsa.begin()) break;
//        }
//      }
//      // search forward
//      vector<gsaTuple>::const_iterator rightArrow = result + 1;
//      if (rightArrow < gsa.end()) {
//        while (lcp(reads[rightArrow->read_idx]
//              .substr(rightArrow->offset), querySubstr)
//              >= MIN_SUFFIX_SIZE) {
//            readsCoveringLocation.insert(rightArrow->read_idx);
//            rightArrow++;
//            if (rightArrow >= gsa.end()) break;
//        }
//      }
//    }
//  }
//  return readsCoveringLocation;
//}

vector<gsaTuple>::const_iterator  binarySearch(vector<string> const& reads, 
    vector<gsaTuple> const& gsa, string const& query)
{
  long unsigned int left{0};
  long unsigned int right{gsa.size()};
  long unsigned int mid;
  while (left < right) {
    mid = left + ((right - left) / 2);
    string suffix = reads[gsa[mid].read_idx].substr(gsa[mid].offset);
    if (lcp(query, reads[gsa[mid].read_idx].substr(gsa[mid].offset)) >=
        MIN_SUFFIX_SIZE) {
      return gsa.begin() + mid;
    }
    else if (suffix < query) left = mid+1;
    else right = mid;
  }
  return gsa.end();
}

int lcp(string const& a, string const& b) {
  int lcp = 0;
  for(; a[lcp] == b[lcp]; lcp++);
  return lcp;
}

void printReadsAndId(int from, int to, int step, vector<string> const&
    reads, string const& nameOfContainer) {
  ofstream ofile("/data/ic711/idEquivQueryFastq.txt", ios_base::app);

  ofile << nameOfContainer << endl;
  if(from < 0 || to > reads.size()) {
    cout << "Out of range" << endl;
    exit(1);
  }
  for(; from < to; from += step) {
    ofile << reads[from] << " : " << from << endl;
  }
  ofile.close();
}


vector<gsaTuple> buildGSA(vector<string> const& fileNames,
    vector<string> & reads, string const& nameOfContainer) {

  // extract reads from tar.gz file
  for (string const& filename  : fileNames) {
    loadFastq(filename, reads);
  }

  cout << "Loaded fastq data." << endl;

  // construct bsa
  vector<pair<unsigned int, unsigned int > > bsa = constructBSA(reads);
  cout << "Constructed bsa " << endl;

  // construct sa
  string concat;
  for(auto r : reads) concat += r;
  unsigned long long * sa;
  sa = Radix<unsigned long long>((uchar*) concat.c_str(), concat.size()).build();
  cout << "Constructed sa" << endl;

  // construct gsa
  vector<gsaTuple> gsa;
  for(unsigned long long int i=0; i < concat.size(); i++) {
    pair <unsigned int, unsigned int> read_concat_pair;
    read_concat_pair = binarySearch(bsa, sa[i]);

    unsigned int offset = sa[i] - read_concat_pair.second;
    int read_size = reads[read_concat_pair.first].size();

    // remove suffixes that are too short
    if (read_size - offset <= MIN_SUFFIX_SIZE)  continue;

    gsaTuple gsaElem(read_concat_pair.first, offset);
    gsa.push_back(gsaElem); 
  }
  delete [] sa;
  cout << "Constructed gsa" << endl;
  return gsa;
}

pair<unsigned int, unsigned int> binarySearch(
                  vector<pair<unsigned int, unsigned int> > const& BSA, 
                  unsigned int suffix_index) {

  unsigned int right = BSA.size();
  unsigned int left = 0;
  unsigned int mid;

  // binary search to home in on read
  while(left < right) {
    mid = left + ((right - left) / 2);

    if(suffix_index == BSA[mid].second) {
        return BSA[mid];
    }
    else if (suffix_index > BSA[mid].second) {
      left = mid+1;
    }
    else {
      right = mid;
    }
  }
  left--;
  return BSA[left]; // left should be on the seq
}

vector<pair<unsigned int, unsigned int> > constructBSA(
    vector<string> const & reads) {

  vector<pair<unsigned int, unsigned int> > bsa;
  unsigned int index_in_concat = 0;
  pair<unsigned int, unsigned int> firstRead(0, 0);
  bsa.push_back(firstRead);

  for (unsigned int i=1; i < reads.size(); i++) {
    index_in_concat += reads[i-1].size();
    pair<unsigned int, unsigned int> bsa_element(i, index_in_concat);
    bsa.push_back(bsa_element);
  }

  return bsa;
}

void splitFileNamesOnDataType(vector<string> & hFiles, 
    vector<string> & cFiles, string const& dataFile) {
  ifstream dataFileHandle(dataFile);
  string line;
  while (getline(dataFileHandle, line)) {
    vector<string> fields;
    string token = ",";
    split_string(line, token, fields);
    if (fields[1] == "H") {         // store filenames
      hFiles.push_back(fields[0]);
    }
    else if (fields[1] == "T"){
      cFiles.push_back(fields[0]);
    }
  }
}


void loadFastq(string filename, vector<string> &p_data) {

  cout << "Loading reads from " << filename << endl;
  gzFile data_file;
  data_file = gzopen(filename.c_str(), "r");    // open stream to next fastq.gz 
  kseq_t *seq = kseq_init(data_file);           // init parser

  vector<fastq_t> fastq_elements;

  // Load data from file into array
  int eof_check;
  fastq_t next_read;
  while ((eof_check = kseq_read(seq)) >=0) {

    // As we only want high quality reads, 
    // only reads with a quality score have potential to be added

    if (seq->qual.l) {
      next_read.id   = seq->name.s;
      // copy sequence
      next_read.seq  = seq->seq.s;
      // copy quality string
      next_read.qual = seq->qual.s;

      fastq_elements.push_back(next_read);     // load into global vector
    }
  }
  kseq_destroy(seq);
  gzclose(data_file);
 

  // MULTITHREADED SECTION
  vector<thread> thread_group;                      // set up thread store
  thread_group.reserve(N_THREADS);

  vector<fastq_t> *fastq_elements_p = &fastq_elements;

  int elements_per_thread = (fastq_elements.size() / N_THREADS);
  int from = 0, to = elements_per_thread;

  // run threads...

  for(int i=0; i < N_THREADS; ++i) {
    // run thread, processing a chunk of the raw data
    thread_group.push_back(
         thread(&qualityProcessRawData, 
                     fastq_elements_p, &p_data, from, to, i));


    // set up next chunk
    from = to;
    if (i == N_THREADS-1) {// last thread
      cout << "All reads input" << endl;
      to = fastq_elements.size();
    }
    else {
      to += elements_per_thread;
    }
  }

  // wait for threads to finish task, terminating parallel section
  for (auto & thread : thread_group) {
    thread.join();
  }

}

void qualityProcessRawData(vector<fastq_t> *r_data, 
                           vector<string> *p_data,
                           int from,
                           int to, int tid){

  // from, to define the range this thread will process
  vector<string> localThreadsStore;
  localThreadsStore.reserve(to - from);

  cout << "Thread " << tid << " processing section: " << from  
       << " - " << to << endl;


  // Search through string, first determining quality, then 
  // removing substrings
  double n_low_qual_bases = 0.0;
  for(int i = from; i < to; i++) {

    // Begin quality processing using quality sequence. 
    // Using the substring coordinate pairs, search the corresponding
    // positions in the quality sequence. Only keep coordinate pairs
    // where the substring they represent has <10% of positions with score
    // under 20 pthread, which is ascii char '5'


    // link iterators to quality read of the next fastq elem
    string::iterator iter  = (*r_data)[i].qual.begin();
    string::iterator end   = (*r_data)[i].qual.end();
    n_low_qual_bases = 0.0;                           // count low qual
    while (iter != end) {
      if (*iter < PHRED_20) {
        n_low_qual_bases++;
      }
      iter++;
    }

    // if number of low quality bases above QUALITY_THRESH reject it 
    // (skip over the read in the for loop without adding to localThreadStore)
    if( (n_low_qual_bases / (*r_data)[i].qual.size()) >  QUALITY_THRESH) {
      continue;   // skip remaining for iteration
    }

    // else, deemed high quality


    vector<string> tokenless_set;
    split_string((*r_data)[i].seq, REMOVED_TOKENS, tokenless_set); // remove tok
    for (int i=0; i < tokenless_set.size(); i++) {
      if(tokenless_set[i].length() >= MIN_SUFFIX_SIZE - TERM_CHAR_CORRECTION) {
        localThreadsStore.push_back(tokenless_set[i] + TERM_CHAR);
      }
    }
  }

  lock_guard<mutex> lock(quality_processing_lock); // coordinate threads
  for (string accepted_read : localThreadsStore) {
    p_data->push_back(accepted_read);
  }

    // Link iterators to string
    //string::iterator left = (*r_data)[i].seq.begin();
    //string::iterator right = (*r_data)[i].seq.begin();

    //while (right != (*r_data)[i].seq.end()) {    // while not at end...

    //  // if the character 'N' is hit, and the distance from left to 
    //  // right is more than 30 store substring
    //  if(*right == 'N' && (distance(left, right) >= 30)) { 

    //    localThreadsStore.push_back(
    //        (*r_data)[i].seq.substr( (left - (*r_data)[i].seq.begin()), // from left
    //          distance(left, right) ) + "$"
    //    );

    //    // right now pointing at 'N', so move one char forward
    //    right++;
    //    // we just started a new substring sequence, so set left to right
    //    left = right;
    //  }

    //  else if (*right == 'N'){    // hit an 'N', but subseq too small
    //    right++;
    //    left = right;
    //  }

    //  else {                      // normal character
    //    right++;
    //  }
    //}

    //// Whole string was 'N'-less
    //if (distance(left, right) >= 30) {
    //    localThreadsStore.push_back(
    //        (*r_data)[i].seq.substr( (left - (*r_data)[i].seq.begin()), // from left
    //          distance(left, right) ) + "$"
    //    );
    //}
}
