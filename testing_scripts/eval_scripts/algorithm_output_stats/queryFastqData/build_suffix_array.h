#ifndef BUILD_SUFFIX_ARRAY
#define BUILD_SUFFIX_ARRAY
#endif

#include <vector>
#include <set>
#include <string>
#include <mutex>
#include <thread>
#include <algorithm>


#include "string.h"
#include "parse_coordinate_data.h"

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

struct snippetData {
  std::string snippet;
  std::string header;
  unsigned int mutationLocation;
  char healthy;
  char cancer;
  std::vector<std::string> reads;
  std::vector<unsigned int> read_idx;
};

void splitFileNamesOnDataType(std::vector<std::string> & hFiles, 
    std::vector<std::string> & cFiles, std::string const& dataFile);

void loadFastq(std::string filename, std::vector<std::string> &p_data);

void qualityProcessRawData(std::vector<fastq_t> *r_data, 
                           std::vector<std::string> *p_data,
                           int from,int to, int tid);
std::vector<gsaTuple> buildGSA(std::vector<std::string> const& fileNames,
    std::vector<std::string> & reads);

std::pair<unsigned int, unsigned int> binarySearch(
                  std::vector<std::pair<unsigned int, unsigned int> > const& BSA, 
                  unsigned int suffix_index);

std::vector<std::pair<unsigned int, unsigned int> > constructBSA(
    std::vector<std::string> const & reads);

std::vector<gsaTuple>::const_iterator binarySearch(std::vector<std::string> const& reads,
    std::vector<gsaTuple> const& gsa, std::string const& query);

int lcp(std::string const& a, std::string const& b);

std::set<unsigned int> findReadsCoveringLocation(std::vector<std::string> const&
    reads, std::vector<gsaTuple> const& gsa, std::string const& query);

std::vector<snippetData> extractreadsCoveringSnippets(std::vector<coordinateData
    const& coords, std::vector<std::string> const& reads, 
    std::vector<gsaTuple> const& gsa);

