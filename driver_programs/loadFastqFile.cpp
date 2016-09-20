#include <iostream>
#include <vector>
#include <string>
#include <iterator>

// Thread libs
#include <thread>
#include <mutex>

#include <zlib.h> // need zlib.h for parsing of .gz files
#include "kseq.h" // nice fastq parser


KSEQ_INIT(gzFile, gzread);
using namespace std;




struct fastq_t {
  string id, seq, qual;
};


void loadFastqRawDataFromFile(string filename, 
                              vector<fastq_t> &r_data,
                              vector<string> *p_data);
// Function first reads in data from the filename.fasta.gz file
// storing it in the global raw data store
// r_data global buffer. 
// loadFastqRawDataFromFile then threads off, calling the main processing
// thread routine, threads store processed data in p_data;


void qualityProcessRawData(vector <fastq_t> *r_data, 
                           vector <string> *p_data,
                           int from,
                           int to,
                           int tid);
// This function is called by each thread. Use of threads
// is for true parallelism, aiming to speed up quality checking.
// Q.check. consists of 1) discaring if >10% of bases have phread score
// of < q20. 2) spliting a given string, if 'N' exists.


void deleteMemory(vector<string> *data) {
  delete data;
}



int main(int argv, char **argc) {
  vector<string> *TumourReads;      // Simulates ReadManipulator.reads storage

  vector<fastq_t> raw_data;         // holds the rawdata as it comes in

  // init stuff...
  TumourReads = new vector<string>;

  loadFastqRawDataFromFile(argc[1], raw_data, TumourReads);


  cout << "RAW DATA --------------------------------------" << endl << endl;
  cout << "The size of raw_data is: ";
  cout << raw_data.size() << endl;
  cout << sizeof(raw_data);


  for(int i=0; i < raw_data.size(); i++) {
      cout << "READ " << i << endl;
      cout << "id: " << raw_data[i].id << endl;
      cout << "seq: " << raw_data[i].seq << endl;
      cout << "qual: " << raw_data[i].qual << endl;

      cout << "Size of fastq_t is: " << sizeof(raw_data[i]) << endl;
  }

  cout << "ACCEPTED DATA ---------------------------------" << endl << endl;
  cout << "The size of the processed data is: " << TumourReads->size() << endl;
  for(string read : *TumourReads) {
    cout << read << endl;
  }

  delete TumourReads;
  
  return 0;
}

void loadFastqRawDataFromFile(string filename, 
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
         std::thread(qualityProcessRawData, r_data_p, p_data, from, to, i));


    // set up next chunk
    from = to;
    if (i == nthreads-1) {// last thread
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


    n_low_qual_bases = 0.0; 
      // count bases with quality < 20
    while (iter != end) {
      if (*iter < '5') {
        n_low_qual_bases++;
      }
      iter++;
    }
    cout << "NUMBER OF LOW QUAL BASES: "<< n_low_qual_bases << endl;
    cout << "DENOM: " << (*r_data)[i].qual.size() << endl;
    cout << "Ratio: " <<  n_low_qual_bases / (*r_data)[i].qual.size() << endl;

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
              std::distance(left, right) )
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
              std::distance(left, right) )
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

  // when lock goes out of scope... quality_processing_lock is released
  // And p_data can be accessed by another thread
}

// End of file
