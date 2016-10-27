// Reads.cpp
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <utility>
#include <cstring>

// For parallel data processing 
#include <thread>
#include <mutex>

#include <zlib.h>   // gunzip parser
#include "kseq.h"   // fastq parser


#include "util_funcs.h"
#include "Reads.h"

KSEQ_INIT(gzFile, gzread);    // initialize .gz parser


using namespace std;
static const int NUM_ARGS = 4;
static const int HEADER_FILE_IDX = 1;
static const int MIN_SUF_IDX = 2;
static const int ECONT_IDX = 3;
static const int TERM_CHAR_CORRECTION = 1;
static const int N_THREADS = 1;

ReadsManipulator::ReadsManipulator(int argc, char **argv) {

  vector<file_and_type> datafiles;
  parseCommandLine(argc, argv, datafiles);

  // Inititalize Read stores
  HealthyReads = new vector<string>;
  TumourReads  = new vector<string>;


  // start reading files
  cout << "Extracting fastq data from " << datafiles.size() << " files..." << endl;

  for(int i=0; i < datafiles.size(); i++) {
    cout << "Loading data from " << datafiles[i].first << "..." << endl;
    if (datafiles[i].second) { // == HEALTHY
      cout << datafiles[i].first << endl;
      loadFastqRawDataFromFile(datafiles[i].first, HealthyReads);
    }
    else{ // filetype == TUMOUR
      cout << datafiles[i].first << endl;;
      loadFastqRawDataFromFile(datafiles[i].first, TumourReads);
    }
  }
}

void ReadsManipulator::parseCommandLine(int argc, char** argv,
    vector<file_and_type> &datafiles) {

  // Only NUM_ARGS input. If fail count, explain to user options and
  // file format
  if (argc != NUM_ARGS) {
    cout << "Usage: <exec> <datafile_names.txt>"
         << " <minimum_suffix_size> <contamination_ratio>" << endl;

    cout << endl << endl;
    cout << "Please note the the format of headerfile.txt:"
         << endl 
         << "<file1>\t<H/T> <1/2>" << endl
         << "<file2>\t<H/T> <1/2>" << endl
         << "..." << endl;
    exit(1);
  }

  // save params
  minimum_suffix_size = std::stoi(argv[MIN_SUF_IDX]) + TERM_CHAR_CORRECTION;
  econt = std::stod(argv[ECONT_IDX]);

  // store data file names
  ifstream headerfile;
  headerfile.open(argv[HEADER_FILE_IDX]);
  cout << "Parsing inputs..." << endl;
  cout << "Reading filenames from " << argv[HEADER_FILE_IDX] << " ... " << endl;

  string file_string;
  file_and_type file_info;

  while (getline(headerfile, file_string)) {

    // momentarily work with c strings, to make use of the glorious strtok()...

    // load the filename into file_and_type...
    char * c_file = const_cast<char*>(file_string.c_str()); 
    file_info.first = strtok(c_file, ",\t ");

    // load the type bool into file_and_type...
    char * c_datatype = strtok(NULL, ",\t ");      // get next feild
    string datatype = c_datatype;                    // convert to string

    if (datatype == "H") {              // then bool
      file_info.second = HEALTHY;
    }
    else if (datatype == "T") {
      file_info.second = TUMOUR;
    }
    else {
      cout << datatype << " is not a valid datatype, either H or T. " << endl
           << "Program terminating..." << endl;
      exit(1);
    }


    cout << "Storing " << file_info.first << " as "
         << ((file_info.second) ? "healthy" : "tumour") << " data "
         << endl;

    datafiles.push_back(file_info);      // store in params
  }

}

void ReadsManipulator::loadFastqRawDataFromFile(string filename, 
                              vector<string> *p_data) {

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
         std::thread(&ReadsManipulator::qualityProcessRawData, this, 
                     fastq_elements_p, p_data, from, to, i));


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

int ReadsManipulator::getMinSuffixSize() {
  return minimum_suffix_size;
}
double ReadsManipulator::getEcont() {
  return econt;
}


// end of file

