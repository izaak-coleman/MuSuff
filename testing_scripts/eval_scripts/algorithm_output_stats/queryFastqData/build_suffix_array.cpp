#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <mutex>
#include <thread>
#include <algorithm>


// reading tar.gz
#include <zlib.h>
#include "kseq.h"

// sa construction
#include "radix.h"
#include <sys/types.h>
#include <cstring>
#include <cstdlib>

#include "string.h"

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

static const std::string REMOVED_TOKENS = "N"; // remove N from fastq
static const std::string TERM_CHAR = "$";      // suffix termination character
std::mutex quality_processing_lock;

struct fastq_t {  
  std::string id, seq, qual;
};

struct gsaTuple {
  unsigned int read_idx;
  unsigned int offset;
  gsaTuple(unsigned int r, unsigned int o):read_idx(r), offset(o) {}

  friend ostream & operator<<(ostream & cout, gsaTuple const& tup) {
    std::cout << "(" << tup.read_idx << "," << tup.offset << ")" << std::endl;
  }
};

void splitFileNamesOnDataType(std::vector<std::string> & hFiles, 
    std::vector<std::string> & cFiles, std::string const& dataFile);

void loadFastq(std::string filename, std::vector<std::string> &p_data);

void qualityProcessRawData(std::vector<fastq_t> *r_data, 
                           std::vector<std::string> *p_data,
                           int from,int to, int tid);
std::vector<gsaTuple> buildGSA(std::vector<std::string> fileNames);

std::pair<unsigned int, unsigned int> binarySearch(
                  std::vector<std::pair<unsigned int, unsigned int> > const& BSA, 
                  unsigned int suffix_index);

std::vector<std::pair<unsigned int, unsigned int> > constructBSA(
    std::vector<std::string> const & reads);

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cout << "Usage: <exe> <datafilelist>" << std::endl;
    return -1;
  }

  std::vector<std::string> cancerFileNames, healthyFileNames;
  splitFileNamesOnDataType(healthyFileNames,cancerFileNames, argv[1]);

  std::vector<std::string> cancerReads, healthyReads;
  if (!healthyFileNames.empty()) {
    std::vector<gsaTuple> healthyGSA = buildGSA(healthyFileNames);
  }
  if (!cancerFileNames.empty()) {
    std::vector<gsaTuple> cancerGSA  = buildGSA(cancerFileNames);
  }
}

gsaTuple binarySearch(std::vector<std::string> const& reads, 
    std::vector<gsaTuples> const& gsa, std::string query)
{
  unsigned int left{0}, right{reads.size()}, mid;

  while (left < right) {
    mid = left + ((right - left) / 2);
    std::string suffix = reads[gsa[mid].read_idx].substr(gsa[mid].offset);

    if (suffix < query) left = mid+1;
    else if (suffix > query) right = mid;
    else {
      return gsa[mid]
    }
  }
}


std::vector<gsaTuple> buildGSA(std::vector<std::string> fileNames) {

  // extract reads from tar.gz file
  std::vector<std::string> reads;
  for (std::string const& filename  : fileNames) {
    loadFastq(filename, reads);
  }
  std::cout << "Loaded fastq data." << std::endl;

  // construct bsa
  std::vector<std::pair<unsigned int, unsigned int > > bsa = constructBSA(reads);
  std::cout << "Constructed bsa " << std::endl;

  // construct sa
  std::string concat;
  for(auto r : reads) concat += r;
  unsigned long long * sa;
  sa = Radix<unsigned long long>((uchar*) concat.c_str(), concat.size()).build();
  std::cout << "Constructed sa" << std::endl;

  // construct gsa
  std::vector<gsaTuple> gsa;
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
  std::cout << "Constructed gsa" << std::endl;
  return gsa;
}

std::pair<unsigned int, unsigned int> binarySearch(
                  std::vector<std::pair<unsigned int, unsigned int> > const& BSA, 
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

std::vector<std::pair<unsigned int, unsigned int> > constructBSA(
    std::vector<std::string> const & reads) {

  std::vector<std::pair<unsigned int, unsigned int> > bsa;
  unsigned int index_in_concat = 0;
  std::pair<unsigned int, unsigned int> firstRead(0, 0);
  bsa.push_back(firstRead);

  for (unsigned int i=1; i < reads.size(); i++) {
    index_in_concat += reads[i-1].size();
    std::pair<unsigned int, unsigned int> bsa_element(i, index_in_concat);
    bsa.push_back(bsa_element);
  }

  return bsa;
}

void splitFileNamesOnDataType(std::vector<std::string> & hFiles, 
    std::vector<std::string> & cFiles, std::string const& dataFile) {
  std::ifstream dataFileHandle(dataFile);
  std::string line;
  while (std::getline(dataFileHandle, line)) {
    std::vector<std::string> fields;
    std::string token = ",";
    split_string(line, token, fields);
    if (fields[1] == "H") {         // store filenames
      hFiles.push_back(fields[0]);
    }
    else if (fields[1] == "T"){
      cFiles.push_back(fields[0]);
    }
  }
}


void loadFastq(std::string filename, std::vector<std::string> &p_data) {

  gzFile data_file;
  data_file = gzopen(filename.c_str(), "r");    // open stream to next fastq.gz 
  kseq_t *seq = kseq_init(data_file);           // init parser

  std::vector<fastq_t> fastq_elements;

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
  std::vector<std::thread> thread_group;                      // set up thread store
  thread_group.reserve(N_THREADS);

  std::vector<fastq_t> *fastq_elements_p = &fastq_elements;

  int elements_per_thread = (fastq_elements.size() / N_THREADS);
  int from = 0, to = elements_per_thread;

  // run threads...

  for(int i=0; i < N_THREADS; ++i) {
    // run thread, processing a chunk of the raw data
    thread_group.push_back(
         std::thread(&qualityProcessRawData, 
                     fastq_elements_p, &p_data, from, to, i));


    // set up next chunk
    from = to;
    if (i == N_THREADS-1) {// last thread
      std::cout << "All reads input" << std::endl;
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

void qualityProcessRawData(std::vector<fastq_t> *r_data, 
                           std::vector<std::string> *p_data,
                           int from,
                           int to, int tid){

  // from, to define the range this thread will process
  std::vector<std::string> localThreadsStore;
  localThreadsStore.reserve(to - from);

  std::cout << "Thread " << tid << " processing section: " << from  
       << " - " << to << std::endl;


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
    std::string::iterator iter  = (*r_data)[i].qual.begin();
    std::string::iterator end   = (*r_data)[i].qual.end();
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


    std::vector<std::string> tokenless_set;
    split_string((*r_data)[i].seq, REMOVED_TOKENS, tokenless_set); // remove tok
    for (int i=0; i < tokenless_set.size(); i++) {
      if(tokenless_set[i].length() >= MIN_SUFFIX_SIZE - TERM_CHAR_CORRECTION) {
        localThreadsStore.push_back(tokenless_set[i] + TERM_CHAR);
      }
    }
  }

  std::lock_guard<std::mutex> lock(quality_processing_lock); // coordinate threads
  for (std::string accepted_read : localThreadsStore) {
    p_data->push_back(accepted_read);
  }

    // Link iterators to string
    //string::iterator left = (*r_data)[i].seq.begin();
    //string::iterator right = (*r_data)[i].seq.begin();

    //while (right != (*r_data)[i].seq.end()) {    // while not at end...

    //  // if the character 'N' is hit, and the distance from left to 
    //  // right is more than 30 store substring
    //  if(*right == 'N' && (std::distance(left, right) >= 30)) { 

    //    localThreadsStore.push_back(
    //        (*r_data)[i].seq.substr( (left - (*r_data)[i].seq.begin()), // from left
    //          std::distance(left, right) ) + "$"
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
    //if (std::distance(left, right) >= 30) {
    //    localThreadsStore.push_back(
    //        (*r_data)[i].seq.substr( (left - (*r_data)[i].seq.begin()), // from left
    //          std::distance(left, right) ) + "$"
    //    );
    //}
}
