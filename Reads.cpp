// Reads.cpp
#include <string>
#include <vector>
#include <iostream>
#include <utility>
#include <fstream>

// For parallel data processing 
#include <thread>
#include <mutex>

#include <zlib.h>   // gunzip parser
#include "kseq.h"   // fastq parser


#include "helper_functions.h"
#include "Reads.h"

KSEQ_INIT(gzFile, gzread);    // initialize .gz parser


using namespace std;



ReadsManipulator::ReadsManipulator(command_line_params &params) {

  // Inititalize Read stores
  HealthyReads = new vector<string>;
  TumourReads  = new vector<string>;


  // start reading files
  cout << "About to parse " << params.datafiles.size() << " files..." << endl;

  for(int i=0; i < params.datafiles.size(); i++) {

    vector<fastq_t> raw_data; // hold next file raw data

    cout << "Loading data from " << params.datafiles[i].first << " ..." << endl;

    if (params.datafiles[i].second) { // == HEALTHY
      cout << params.datafiles[i].first << "Healllthy" << endl;
      loadFastqRawDataFromFile(params.datafiles[i].first, raw_data, HealthyReads);
    }

    else{ // filetype == TUMOUR
      cout << params.datafiles[i].first << "tumour" << endl;
      loadFastqRawDataFromFile(params.datafiles[i].first, raw_data, TumourReads);
    }
  }

  n_healthy_reads = HealthyReads->size();
  n_tumour_reads  = TumourReads->size();

}

void ReadsManipulator::loadFastqRawDataFromFile(string filename, 
                              vector<fastq_t> &r_data, 
                              vector<string> *p_data) {

  gzFile data_file;
  data_file = gzopen(filename.c_str(), "r");    // link to next fastq.gz 
  kseq_t *seq = kseq_init(data_file);           // init parser

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

      r_data.push_back(next_read);     // load into global vector
    }
  }
  

  cout << "Size of raw_data: " << r_data.size() << endl;
  //cout << "printing raw_data: " << endl;
  //for(fastq_t f: r_data) {
  //  cout << f.id << endl << f.seq << endl << endl;
  //}


  kseq_destroy(seq);
  gzclose(data_file);
 

  // MULTITHREADED SECTION
  int nthreads = 1;
  vector<thread> thread_group;                      // set up thread store
  thread_group.reserve(nthreads);

  vector<fastq_t> *r_data_p = &r_data;

  int elements_per_thread = (r_data.size() / nthreads);
  int from = 0, to = elements_per_thread;

  // run threads...

  for(int i=0; i < nthreads; ++i) {
    // run thread, processing a chunk of the raw data
    thread_group.push_back(
         std::thread(&ReadsManipulator::qualityProcessRawData, this, 
                     r_data_p, p_data, from, to, i));


    // set up next chunk
    from = to;
    if (i == nthreads-1) {// last thread
      cout << "All reads input" << endl;
      to = r_data.size();
    }
    else {
      to += elements_per_thread;
    }
  }

  // wait for threads to finish task, terminating parallel section
  for (auto &thread : thread_group) {
    thread.join();
  }

}

void ReadsManipulator::qualityProcessRawData(vector<fastq_t> *r_data, 
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


    n_low_qual_bases = 0.0; 
      // count bases with quality < 20
    while (iter != end) {
      if (*iter < '5') {
        n_low_qual_bases++;
      }
      iter++;
    }

    // if number of low quality bases above 10%, reject it 
    // (skip over the read in the for loop without adding to localThreadStore)
    if( (n_low_qual_bases / (*r_data)[i].qual.size()) >  0.1) {
      continue;   // skip remaining for iteration
    }


    
    // If quality is high enough, now spliting strings with N
    // characters


    // Link iterators to string
    string::iterator left = (*r_data)[i].seq.begin();
    string::iterator right = (*r_data)[i].seq.begin();

    while (right != (*r_data)[i].seq.end()) {    // while not at end...

      // if the character 'N' is hit, and the distance from left to 
      // right is more than 30 store substring
      if(*right == 'N' && (std::distance(left, right) >= 30)) { 

        localThreadsStore.push_back(
            (*r_data)[i].seq.substr( (left - (*r_data)[i].seq.begin()), // from left
              std::distance(left, right) ) + "$"
        );

        // right now pointing at 'N', so move one char forward
        right++;
        // we just started a new substring sequence, so set left to right
        left = right;
      }

      else if (*right == 'N'){    // hit an 'N', but subseq too small
        right++;
        left = right;
      }

      else {                      // normal character
        right++;
      }
    }

    // Whole string was 'N'-less
    if (std::distance(left, right) >= 30) {
        localThreadsStore.push_back(
            (*r_data)[i].seq.substr( (left - (*r_data)[i].seq.begin()), // from left
              std::distance(left, right) ) + "$"
        );
    }
  }


  // Now the reads have been processed, each reads section of processed
  // data is stored in the thread local variable localThreadStore. 
  // This data needs to be passed to the global store, and to stop interference
  // (Multiple threads writing to the same location - data race)
  // Mutex is needed. 

  // Copy local thread store to global store


  std::lock_guard<std::mutex> lock(quality_processing_lock); // coordinate threads
  for (string accepted_read : localThreadsStore) {
    p_data->push_back(accepted_read);
  }

  cout << "reach end of function" << endl;

  // when lock goes out of scope... quality_processing_lock is released
  // And p_data can be accessed by another thread
}



ReadsManipulator::~ReadsManipulator() {
  delete HealthyReads;
  delete TumourReads;
}


void ReadsManipulator::printReads(){
  cout << "Healthy file reads: " << endl;
  cout << "Size of HealthyReads: " << getSize(HEALTHY) << endl;
  for(string s : *HealthyReads) {
    cout << s << endl;
  }
  cout << endl << endl;

  cout << "Tumour file reads: " << endl;
  cout << "Size of TumourReads: " << getSize(TUMOUR) << endl;
  for(string s : *TumourReads) {
    cout << s << endl;
  }
}

string::iterator ReadsManipulator::returnStartIterator(Suffix_t &suf) {
  // Use suf.type and suf.read_id to locate the read, and then set an iterator
  // pointing at suf.offset dist from begining

  string::iterator iter;
  if(suf.type == HEALTHY) { 
    iter = (*HealthyReads)[suf.read_id].begin() + suf.offset;
  }
  else {  // suf.type == TUMOUR
    iter = (*TumourReads)[suf.read_id].begin() + suf.offset;
  }

  return iter;
}

string::iterator ReadsManipulator::returnEndIterator(Suffix_t &suf) {
  // Use suf.type and suf.read_id to locate the read, then return an iterator to the 
  // end of that read

  string::iterator iter;
  if (suf.type == HEALTHY) {
    iter = (*HealthyReads)[suf.read_id].end();
  }
  else {  // suf.type == TUMOUR
    iter = (*TumourReads)[suf.read_id].end();
  }

  return iter;
}

string ReadsManipulator::returnSuffix(Suffix_t &suf){
  // return the string assoc. with suf
  if (suf.type == HEALTHY) {
    return (*HealthyReads)[suf.read_id].substr(suf.offset);
  }
  else { // suf.type == TUMOUR
    return (*TumourReads)[suf.read_id].substr(suf.offset);
  }
}

unsigned int ReadsManipulator::getSize(bool tissueType) {
  if (tissueType == HEALTHY) {
    return HealthyReads->size();
  }
  else {    // == TUMOUR
    return TumourReads->size();
  }
}

string & ReadsManipulator::getReadByIndex(int index, int tissue) {
  if(tissue == HEALTHY) {
    return (*HealthyReads)[index];
  }
  else {  // tissue == TUMOUR || tissue == SWITCHED
    return (*TumourReads)[index];
  }
}


// end of file

