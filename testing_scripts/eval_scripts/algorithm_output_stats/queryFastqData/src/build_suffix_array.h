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


struct fastq_t {  
  std::string id, seq, qual;
};

struct gsaTuple {
  unsigned int read_idx;
  unsigned int offset;
  gsaTuple(unsigned int r, unsigned int o):read_idx(r), offset(o) {}

  friend std::ostream & operator<<(std::ostream & cout, gsaTuple const& tup) {
    cout << "(" << tup.read_idx << "," << tup.offset << ")" << std::endl;
    return cout;
  }
};

struct snippetData {
  std::string snippet;
  std::string header;
  unsigned int mutationLocation;
  char healthy;
  char cancer;
  TissueType tissue; 
  std::vector<std::string> reads;
  std::vector<unsigned int> read_idx;
  std::vector<bool> mut_type;
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

std::vector<snippetData>
extractReadsCoveringSnippets(std::vector<coordinateData>
    const& coords, std::vector<std::string> const& reads, 
    std::vector<gsaTuple> const& gsa,
    TissueType tissue);

void printSnippetData(std::ostream & out, std::vector<snippetData> const& data);
#endif
