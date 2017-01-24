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
  std::vector<std::string> cancerReads, healthyReads;

  std::vector<gsaTuple> healthyGSA, cancerGSA;
  if (!healthyFileNames.empty()) {
    healthyGSA = buildGSA(healthyFileNames, healthyReads);
  }
  if (!cancerFileNames.empty()) {
    cancerGSA = buildGSA(cancerFileNames, cancerReads);
  }

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

  cancerExtraction.insert(cancerExtraction.end(), healthyExtraction.begin(),
      healthyExtraction.end());   // concat vector

  std::sort(cancerExtraction.begin(), cancerExtraction.end(), 
      [] (snippetData const& a, snippetData const& b) {
        return a.mutationLocation < b.mutationLocation;
      });

  ofstream oFileHandle(argv[3]);
  printSnippetData(oFileHandle, cancerExtraction);

  return 0;
}
