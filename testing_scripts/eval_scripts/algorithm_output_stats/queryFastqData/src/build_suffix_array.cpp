#include <iostream>
#include <fstream>
#include <vector>
#include <set>
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
#include "build_suffix_array.h" 
std::mutex quality_processing_lock;

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

// main

//int main(int argc, char **argv) {
//  if (argc != 2) {
//    std::cout << "Usage: <exe> <datafilelist>" << std::endl;
//    return -1;
//  }
//
//  std::vector<std::string> cancerFileNames, healthyFileNames;
//  splitFileNamesOnDataType(healthyFileNames,cancerFileNames, argv[1]);
//
//  std::vector<std::string> cancerReads, healthyReads;
//  std::vector<gsaTuple> cancerGSA, healthyGSA;
//  if (!healthyFileNames.empty()) {
//    healthyGSA = buildGSA(healthyFileNames, healthyReads);
//  }
//  if (!cancerFileNames.empty()) {
//    cancerGSA  = buildGSA(cancerFileNames, cancerReads);
//  }
//
//  std::cout << "Size of cancer gsa: " << cancerGSA.size() << std::endl;
//  std::cout << "Size of cancer reads: " << cancerReads.size() << std::endl;
//
//  for (gsaTuple const& tup : cancerGSA) {
//    std::cout << cancerReads[tup.read_idx].substr(tup.offset) << std::endl;
//  }
//  
//  std::string read = "HelloMrPostmanLookAtMeWoahYeah";
//  std::set<unsigned int> found_reads = findReadsCoveringLocation(cancerReads,
//      cancerGSA, read);
//  std::cout << std::endl << std::endl << "Identfied reads: " << std::endl;
//  for (unsigned int read_idx : found_reads) {
//    std::cout << cancerReads[read_idx] << std::endl;
//  }
//}


std::vector<snippetData>
extractReadsCoveringSnippets(std::vector<coordinateData> const& coords,
                             std::vector<std::string> const& reads,
                             std::vector<gsaTuple> const& gsa,
                             TissueType tissue) {

  std::vector<snippetData> results;
  for (coordinateData const& coord : coords) {
    std::set<unsigned int> readsCoveringCoord;
    readsCoveringCoord = findReadsCoveringLocation(reads, gsa, coord.sequence);
    snippetData entry; 
    entry.header = coord.header;
    entry.snippet = coord.sequence;
    entry.healthy = coord.hBase;
    entry.cancer = coord.cBase;
    entry.mutationLocation  = coord.coordinate;
    entry.tissue = tissue;
    for (unsigned int read_idx : readsCoveringCoord) {
      entry.reads.push_back(reads[read_idx]);
      entry.read_idx.push_back(read_idx);
    }

    // done filling entry data
    results.push_back(entry);
  }
  return results;
}

void printSnippetData(std::ostream & out, std::vector<snippetData> const& data) {
  for (snippetData const& entry : data) {
    out << "Header: "  << entry.header 
              << " - " << entry.healthy << ":" << entry.cancer
              << " - " << ((entry.tissue == TissueType::healthy) ? "H":"T")
              << std::endl;

    out << "Snippt: "  << entry.snippet << std::endl;
    out << "Coord : "  << entry.mutationLocation << std::endl;
    out << "Reads : " << std::endl;
    for (int i = 0; i < entry.reads.size(); i++) {
      out << entry.reads[i] << " :: "
                << entry.read_idx[i] 
                << std::endl;
    }
  }
}


std::set<unsigned int> findReadsCoveringLocation(std::vector<std::string> const& reads,
    std::vector<gsaTuple> const& gsa, std::string const& query) {
  std::set<unsigned int> readsCoveringLocation;

  for (int i=0; i <= query.size() - 30; i++) {
    std::string querySubstr = query.substr(i, 30);
    std::vector<gsaTuple>::const_iterator result = binarySearch(reads, gsa, querySubstr);
    if (result != gsa.end()) {
      // input result
      readsCoveringLocation.insert(result->read_idx);
      // search back
      std::vector<gsaTuple>::const_iterator leftArrow = result-1;
      if (leftArrow >= gsa.begin()) {
        while (
            lcp(reads[leftArrow->read_idx]
              .substr(leftArrow->offset), querySubstr)
              >= MIN_SUFFIX_SIZE) {
          readsCoveringLocation.insert(leftArrow->read_idx);
          leftArrow--;
          if (leftArrow < gsa.begin()) break;
        }
      }
      // search forward
      std::vector<gsaTuple>::const_iterator rightArrow = result + 1;
      if (rightArrow < gsa.end()) {
        while (lcp(reads[rightArrow->read_idx]
              .substr(rightArrow->offset), querySubstr)
              >= MIN_SUFFIX_SIZE) {
            readsCoveringLocation.insert(rightArrow->read_idx);
            rightArrow++;
            if (rightArrow >= gsa.end()) break;
        }
      }
    }
  }
  return readsCoveringLocation;
}

std::vector<gsaTuple>::const_iterator  binarySearch(std::vector<std::string> const& reads, 
    std::vector<gsaTuple> const& gsa, std::string const& query)
{
  long unsigned int left{0};
  long unsigned int right{gsa.size()};
  long unsigned int mid;
  while (left < right) {
    mid = left + ((right - left) / 2);
    std::string suffix = reads[gsa[mid].read_idx].substr(gsa[mid].offset);
    if (lcp(query, reads[gsa[mid].read_idx].substr(gsa[mid].offset)) >=
        MIN_SUFFIX_SIZE) {
      return gsa.begin() + mid;
    }
    else if (suffix < query) left = mid+1;
    else right = mid;
  }
  return gsa.end();
}

int lcp(std::string const& a, std::string const& b) {
  int lcp = 0;
  for(; a[lcp] == b[lcp]; lcp++);
  return lcp;
}


std::vector<gsaTuple> buildGSA(std::vector<std::string> const& fileNames,
    std::vector<std::string> & reads) {

  // extract reads from tar.gz file
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

  std::cout << "Loading reads from " << filename << std::endl;
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
