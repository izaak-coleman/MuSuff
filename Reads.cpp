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
#include "string.h" // split_string()
#include "Reads.h"

KSEQ_INIT(gzFile, gzread);    // initialize .gz parser


using namespace std;
static const int HEADER_FILE_IDX = 1;
static const int ECONT_IDX = 2;
static const int OFILE_IDX = 3;
static const int TERM_CHAR_CORRECTION = 1;
static const int MIN_SUFFIX_SIZE = 30;
static const int DISTAL_TRIM = 0;

static const double QUALITY_THRESH = 0.1; // 10% 
static const char PHRED_20 = '5';   // lowest high quality phred score
static const string REMOVED_TOKENS = "N"; // remove N from fastq
static const string TERM_CHAR = "$";      // suffix termination character



ReadsManipulator::ReadsManipulator(int n_threads, string const& inputFile):
N_THREADS(1){
  minimum_suffix_size = MIN_SUFFIX_SIZE;
  distal_trim_len = DISTAL_TRIM;

  vector<file_and_type> datafiles;
  parseInputFile(inputFile, datafiles);
  cout << "DISTAL_TRIM: " << DISTAL_TRIM << endl;

  cout << "Loaded " << datafiles.size() << " data files." << endl;

  for(int i=0; i < datafiles.size(); i++) {
    cout << "Extracting data from " << datafiles[i].first << "..." << endl;
    if (datafiles[i].second) { // == HEALTHY
      loadFastqRawDataFromFile(datafiles[i].first, HealthyReads, HealthyPhreds);
    }
    else{ // filetype == TUMOUR
      loadFastqRawDataFromFile(datafiles[i].first, TumourReads, TumourPhreds);
    }
  }

//  printRemainingReads("/data/ic711/point1.txt");
//  printAllReads();
  cout << "End of ReadsManipulator constructor " << endl;
//  writeContainer(HealthyReads, "/data/ic711/HealthyReads.txt");
//  writeContainer(TumourReads, "/data/ic711/TumourReads.txt");
//  writeContainer(HealthyPhreds, "/data/ic711/HealthyPhreds.txt");
//  writeContainer(TumourPhreds, "/data/ic711/TumourPhreds.txt");
}

void ReadsManipulator::writeContainer(vector<string> const& c, string const& fname) {
  ofstream of(fname.c_str());
  for (string const& s : c) {
    of << s << endl;
  }
  of.close();
}

void ReadsManipulator::printAllReads() {
  ofstream osock("/data/ic711/checking_read_trim_len.txt");
  for (string s : HealthyReads) {
    s.pop_back();
    osock << s << endl;
  }
  for (string s : TumourReads) {
    s.pop_back();
    osock << s << endl;
  }
  osock.close();
}

void ReadsManipulator::parseInputFile(string const& inputFile, 
                                        vector<file_and_type> &datafiles) {
  ifstream sock;
  sock.open(inputFile.c_str());
  cout << "Gathering datafiles from " << inputFile << "." << endl;
  string file_string;
  while (getline(sock, file_string)) {
    file_and_type file_info;
    vector<string> fields;
    split_string(file_string, ",\t ", fields);

    // load data
    file_info.first = fields[0];
    if (fields[1] == "H") file_info.second = HEALTHY;
    else if (fields[1] == "T") file_info.second = TUMOUR;
    else {
      cout << fields[1] 
           << " is not a valid datatype, either H or T." 
           << endl << "Program terminating." << endl;
      exit(1);
    }
    cout << "Input " << file_info.first << " as "
         << ((file_info.second) ? "healthy" : "tumour") << " data "
         << endl;
    datafiles.push_back(file_info);
  }
  sock.close();
}

string ReadsManipulator::performDistalTrim(string & s) {
  s.erase(0, distal_trim_len);
  s.erase(s.length() - distal_trim_len);
  return s;
}

void ReadsManipulator::loadFastqRawDataFromFile(string filename, 
                              vector<string> &processed_reads, 
                              vector<string> & processed_phreds) {

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
                     fastq_elements_p, &processed_reads, &processed_phreds, from, to, i));


    // set up next chunk
    from = to;
    if (i == N_THREADS-1) {// last thread
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
                           vector<string> *processed_reads,
                           vector<string> *processed_phreds,
                           int from,
                           int to, int tid){

  // from, to define the range this thread will process
  vector<string> readThreadStore;
  vector<string> phredThreadStore;
  readThreadStore.reserve(to - from);
  phredThreadStore.reserve(to - from);

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


// REMEMBER TO REMOVE PERFORM DISTAL TRIM FOR PRODUCTION
    vector<string> read_substrs;  // remove N chars
    vector<string> phred_substrs;
    int left_arrow{0}, right_arrow{0};
    while((right_arrow = (*r_data)[i].seq.find(REMOVED_TOKENS, left_arrow)) != string::npos) {
      if (left_arrow == right_arrow) {
        left_arrow++;
        continue;
      }
      read_substrs.push_back((*r_data)[i].seq.substr(left_arrow, right_arrow - left_arrow));
      phred_substrs.push_back((*r_data)[i].qual.substr(left_arrow, right_arrow - left_arrow));
      left_arrow = right_arrow + 1;
    }
    read_substrs.push_back((*r_data)[i].seq.substr(left_arrow));
    phred_substrs.push_back((*r_data)[i].qual.substr(left_arrow));

    for (int i=0; i < read_substrs.size(); i++) {
      if(read_substrs[i].length() >= MIN_SUFFIX_SIZE /*- TERM_CHAR_CORRECTION */) {
        readThreadStore.push_back(performDistalTrim(read_substrs[i]) + TERM_CHAR);
        phredThreadStore.push_back(performDistalTrim(phred_substrs[i]));
      }
    }
  }

  std::lock_guard<std::mutex> lock(quality_processing_lock); // coordinate threads
  for (string accepted_read : readThreadStore) {
    processed_reads->push_back(accepted_read); // load trimmed read
  }
  for (string phred : phredThreadStore) {
    processed_phreds->push_back(phred);
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



ReadsManipulator::~ReadsManipulator() {
}


void ReadsManipulator::printReads(){
  std::cout << "Healthy file reads: " << endl;
  std::cout << "Size of HealthyReads: " << getSize(HEALTHY) << endl;
  for(string s : HealthyReads) {
    std::cout << s << endl;
  }
  std::cout << endl << endl;

  std::cout << "Tumour file reads: " << endl;
  std::cout << "Size of TumourReads: " << getSize(TUMOUR) << endl;
  for(string s : TumourReads) {
    std::cout << s << endl;
  }
}


void ReadsManipulator::printReadsAndId(int from, int to, int step) const {
  ofstream ofile("/data/ic711/readIdEquivICSmuFin.txt");
  if (from < 0 || to > HealthyReads.size() || to > TumourReads.size()) {
    cout << "Out of range" << endl;
    exit(1);   // should throw
  }
  ofile << "Healthy Reads" << endl;
  for (int it = from; it < to; it += step) {
    ofile << HealthyReads[it] << " : " << it << endl;
  }

  ofile << "Cancer Reads" << endl;
  for (int it = from; it < to; it += step) {
    ofile << TumourReads[it] << " : " << it << endl;
  }

  ofile.close();
}

string::iterator ReadsManipulator::returnStartIterator(Suffix_t &suf) {
  // Use suf.type and suf.read_id to locate the read, and then set an iterator
  // pointing at suf.offset dist from begining

  string::iterator iter;
  if(suf.type == HEALTHY) { 

    if (suf.read_id >= HealthyReads.size() || suf.read_id < 0) {
      cout << "returnStartIterator() out of bounds " << endl;
      exit(1);
    }
    iter = HealthyReads[suf.read_id].begin() + suf.offset;
  }
  else {  // suf.type == TUMOUR
    if (suf.read_id >= TumourReads.size() || suf.read_id < 0) {
      cout << "returnStartIterator() out of bounds " << endl;
      exit(1);
    }
    iter = TumourReads[suf.read_id].begin() + suf.offset;
  }

  return iter;
}

string::iterator ReadsManipulator::returnEndIterator(Suffix_t &suf) {
  // Use suf.type and suf.read_id to locate the read, then return an iterator to the 
  // end of that read

  string::iterator iter;
  if (suf.type == HEALTHY) {
    if (suf.read_id >= HealthyReads.size() || suf.read_id < 0) {
      cout << "returnEndIterator() out of bounds " << endl;
      exit(1);
    }
    iter = HealthyReads[suf.read_id].end();
  }
  else {  // suf.type == TUMOUR
    if (suf.read_id >= TumourReads.size() || suf.read_id < 0) {
      cout << "returnEndIterator() out of bounds " << endl;
      exit(1);
    }
    iter = TumourReads[suf.read_id].end();
  }

  return iter;
}

string ReadsManipulator::returnSuffix(Suffix_t &suf){
  // return the string assoc. with suf
  if (suf.type == HEALTHY) {
    if (suf.read_id >= HealthyReads.size() || suf.read_id < 0) {
      cout << "returnSuffix() out of bounds " << endl;
      exit(1);
    }
    return HealthyReads[suf.read_id].substr(suf.offset);
  }
  else { // suf.type == TUMOUR
    if (suf.read_id >= TumourReads.size() || suf.read_id < 0) {
      cout << "returnSuffix() out of bounds " << endl;
      exit(1);
    }
    return TumourReads[suf.read_id].substr(suf.offset);
  }
}

unsigned int ReadsManipulator::getSize(bool tissueType) {
  if (tissueType == HEALTHY) {
    return HealthyReads.size();
  }
  else {    // == TUMOUR
    return TumourReads.size();
  }
}

string & ReadsManipulator::getReadByIndex(int index, int tissue) {
  if(tissue == HEALTHY) {
    if (index >= HealthyReads.size() || index < 0) {
      cout << "getReadByIndex() out of bounds" << endl;
      exit(1);
    }
    return HealthyReads[index];
  }
  else {  // tissue == TUMOUR || tissue == SWITCHED
    if (index >= TumourReads.size() || index < 0) {
      cout << "getReadsByIndex() out of bounds" << endl;
      exit(1);
    }
    return TumourReads[index];
  }
}
string & ReadsManipulator::getPhredString(int index, int tissue) {
  if(tissue == HEALTHY) {
    if (index >= HealthyPhreds.size() || index < 0) {
      cout << "getReadByIndex() out of bounds" << endl;
      exit(1);
    }
    return HealthyPhreds[index];
  }
  else {  // tissue == TUMOUR || tissue == SWITCHED
    if (index >= TumourPhreds.size() || index < 0) {
      cout << "getReadsByIndex() out of bounds" << endl;
      exit(1);
    }
    return TumourPhreds[index];
  }
}

char ReadsManipulator::baseQuality(int index, int tissue, int pos) {
  if (tissue == HEALTHY) {
    return HealthyPhreds[index][pos];
  }
  else { // tissue == TUMOUR || tissue == SWITCHED
    return HealthyPhreds[index][pos];
  }
}

int ReadsManipulator::getMinSuffixSize() {
  return minimum_suffix_size;
}


void ReadsManipulator::printRemainingReads(std::string const& filename) {
  ofstream fileHandle(filename.c_str());

  for (unsigned int i=0; i < HealthyReads.size(); i++) {
    fileHandle << "(" << i << ",H)"  << std::endl;
  }
  for (unsigned int i=0; i < TumourReads.size(); i++) {
    fileHandle << "(" << i << ",T)"  << std::endl;
  }
  fileHandle.close();
} 


// end of file


