#ifndef BUILD_SUFFIX_ARRAY
#define BUILD_SUFFIX_ARRAY

#include <vector>
#include <set>
#include <string>
#include <mutex>
#include <iostream>
#include <thread>
#include <algorithm>


#include "string.h"
#include "parse_coordinate_data.h"

enum class TissueType {tumour, healthy};
enum class Orientation {fwd, rev};
enum class Covers{mut, nonMut};


struct fastq_t {  
  std::string id, seq, qual;
};



struct gsaTuple {
  unsigned int read_idx;
  unsigned int offset;
  Orientation orientation;
  Covers covers;
  int relative_to;

  gsaTuple(unsigned int idx, unsigned int off):
  read_idx(idx), offset(off) {
    relative_to = -1;
  }
  gsaTuple(unsigned int idx, unsigned int off, Orientation ori, Covers cov, 
      int rel)
    :read_idx(idx), offset(off), orientation(ori), covers(cov), relative_to(rel) {
  }

  friend std::ostream & operator<<(std::ostream & cout, gsaTuple const& tup) {
    cout << "(" << tup.read_idx << "," << tup.offset << ")" << std::endl;
    return cout;
  }
};

struct compareGSATuple {
  bool operator() (gsaTuple const& a, gsaTuple const& b) {
    return a.read_idx < b.read_idx;
  }
};

struct snippetData {
  std::string snippet;
  std::string header;
  unsigned int mutationLocation;
  char healthy;
  char cancer;
  unsigned int fd;
  TissueType tissue; 

  std::set<gsaTuple, compareGSATuple> coveringReads;
  //std::vector<std::string> reads;
  //std::vector<unsigned int> read_idx;
  //std::vector<bool> found_on;
};

void splitFileNamesOnDataType(std::vector<std::string> & hFiles, 
    std::vector<std::string> & cFiles, std::string const& dataFile);

void loadFastq(std::string filename, std::vector<std::string> &p_data);

void qualityProcessRawData(std::vector<fastq_t> *r_data, 
                           std::vector<std::string> *p_data,
                           int from,int to, int tid);
std::vector<gsaTuple> buildGSA(std::vector<std::string> const& fileNames,
    std::vector<std::string> & reads, std::string const& nameOfContainer);

std::pair<unsigned int, unsigned int> binarySearch(
                  std::vector<std::pair<unsigned int, unsigned int> > const& BSA, 
                  unsigned int suffix_index);

std::vector<std::pair<unsigned int, unsigned int> > constructBSA(
    std::vector<std::string> const & reads);


std::vector<gsaTuple>::const_iterator binarySearch(std::vector<std::string> const& reads,
    std::vector<gsaTuple> const& gsa, std::string const& query);

int lcp(std::string const& a, std::string const& b);

std::set<gsaTuple, compareGSATuple> findReadsCoveringLocation(std::vector<std::string> const&
    reads, std::vector<gsaTuple> const& gsa, std::string const& query,
    Orientation ori, Covers cov);

std::vector<snippetData>
extractReadsCoveringSnippets(std::vector<coordinateData>
    const& coords, std::vector<std::string> const& reads, 
    std::vector<gsaTuple> const& gsa,
    TissueType tissue);

void printSnippetData(std::ostream & out, std::vector<snippetData> const& data,
     std::vector<std::string> const& healthyReads, std::vector<std::string>
     const& cancerReads);

void printReadsAndId(int from, int to, int step, std::vector<std::string> const&
    reads, std::string const& nameOfContainer);

std::string rc(std::string const& s);
// returns the reverse complement of a string
#endif
