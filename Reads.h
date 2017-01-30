#ifndef READS_H
#define READS_H

#include <string>
#include <iostream>
#include <vector>
#include <utility>

#include <mutex>  // lock

#include "util_funcs.h"

struct fastq_t {      // Struct only read needs to know about
  std::string id, seq, qual;
};

struct file_and_type {
  std::string first;    // filename
  bool second;          // data set
};


class ReadsManipulator{
  // This class stores all the reads for the dataset. 
  // Functions are allowed direct access to reads. 

private:

  int minimum_suffix_size;
  double econt;
  std::string ofile;
  std::vector<std::string> HealthyReads;  // Container for healthy dataset
  std::vector<std::string> TumourReads;   // Container for cancer dataset 
  std::mutex quality_processing_lock;  // lock for thread copy to H/T.Reads


  void parseCommandLine(int argc, char **argv, 
       std::vector<file_and_type> &datafiles);
  // extracts data file information from header file
  // extracts econt and minimum suffix length

  void loadFastqRawDataFromFile(std::string filename, 
                              std::vector<std::string> &p_data);
  // Function loads DNA reads from filename.fasta.gz
  // either HealthyReads/TumourReads dep. on passes pointer

  void qualityProcessRawData(std::vector<fastq_t> *r_data, 
                            std::vector<std::string> *p_data,
                            int from, int to, int tid);
  // Function deployed on threads. Function acts to:
  // 1) Discard reads where number positions in a read with a value
  // less than '5' is over 10% (QUALITY_THRESH)
  // 2) Removes N characters spliting the read. Fragments with size < 30 are
  // min_suffix_size

public:
 // these arrays have a 1:1 mapping with the HealthyReads, TumourReads arrays
 // where their values express whether the read at the same index in 
 // the Healthy/TumourReads arrays if of LEFT or RIGHT type

 ReadsManipulator(int argc, char **argv);
 // Constructor for loading and processing reads


 ~ReadsManipulator();
 // Deletes the reads


 void printReads();
 // Prints the loaded reads REMOVE THIS FUNCTION

 unsigned int getSize(bool tissueType);
 // returns the size HealthyReads, or TumourReads dep. on tissueType


 std::string::iterator returnStartIterator(Suffix_t &suf);
// Function locates the read corresponding to suf.read_id and 
// Sets a pointer in that read starting at suf.offset
 
 std::string::iterator returnEndIterator(Suffix_t &suf);
// Function returns a pointer to the end of the read corresponding to 
// suf.read_id

 std::string returnSuffix(Suffix_t &suf);
// Function returns the suffix that the suffix_t represents

 std::string & getReadByIndex(int index, int tissue);
 // This function returns the read by adress from the given 
 // Reads vector

 int getMinSuffixSize();
 // returns the min suffix size of suffixes in the suffix array
 // as specified by the user

 double getEcont();
 // returns the econt value as specified by the user at cmd-line

 std::string outFile();
 // returns the name of the final output file

 void printRemainingReads(std::string const& filename);

};
#endif
