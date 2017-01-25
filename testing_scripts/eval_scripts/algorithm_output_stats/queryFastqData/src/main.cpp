// main.cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include <algorithm>

#include "build_suffix_array.h"
#include "parse_coordinate_data.h"

using namespace std;

int main(int argc, char **argv) {

  if (argc != 4) {
    std::cout << "Usage: <exe> <coord_data> <readpaths> <ofile>" << std::endl;
    return -1;
  }

  std::vector<coordinateData> coordData = parseCoordinateData(argv[1]);
  std::vector<std::string> cancerFileNames, healthyFileNames;
  splitFileNamesOnDataType(healthyFileNames, cancerFileNames, argv[2]);


  // load reads and build gsa's
  std::vector<std::string> cancerReads, healthyReads;
  std::vector<gsaTuple> healthyGSA, cancerGSA;
  if (!healthyFileNames.empty()) {
    healthyGSA = buildGSA(healthyFileNames, healthyReads);
  }
  if (!cancerFileNames.empty()) {
    cancerGSA = buildGSA(cancerFileNames, cancerReads);
  }

  // extract the cancer and healthy reads covering each coordinate
  // load into a snippetData struct
  vector<snippetData> cancerExtraction =
    extractReadsCoveringSnippets(coordData, 
                                 cancerReads,
                                 cancerGSA,
                                 TissueType::tumour);
  vector<snippetData> healthyExtraction = 
    extractReadsCoveringSnippets(coordData, 
                                 healthyReads, 
                                 healthyGSA,
                                 TissueType::healthy);

  // concat healthy and cancer entries into cancer
  cancerExtraction.insert(cancerExtraction.end(), healthyExtraction.begin(),
      healthyExtraction.end());   // concat vector
  // sort all the entries on mutationLocation - same indexes will be 
  // next to each other
  std::sort(cancerExtraction.begin(), cancerExtraction.end(), 
      [] (snippetData const& a, snippetData const& b) {
        return a.mutationLocation < b.mutationLocation;
      });

  ofstream oFileHandle(argv[3]);
  printSnippetData(oFileHandle, cancerExtraction);

  return 0;
}
