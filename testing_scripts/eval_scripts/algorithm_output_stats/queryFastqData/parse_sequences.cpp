// parse_coordinate_data.cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "parse_coordinate_data.h"

using namespace std;

void printCoordinateData(std::vector<coordinateData> const& coordData) {
  for (coordinateData const& c : coordData) {
    std::cout << "Header: " << c.header << std::endl;
    std::cout << "Coord : " << c.coordinate << std::endl;
    std::cout << "Seq   : " << c.sequence << std::endl;
    std::cout << std::endl;
  }
}


vector<coordinateData> parseCoordinateData(std::string filename) {
  ifstream coordDataHandle;
  coordDataHandle.open(filename);
  std::vector<coordinateData> coordData;

  std::string line;
  std::string token = " --- ";
  while (std::getline(coordDataHandle, line)) {
    std::size_t endOfHeader = line.find(token);
    std::string header = line.substr(0, endOfHeader);
    unsigned int coordinate = 
      std::stoi(line.substr(endOfHeader+token.size()));
    std::string sequence;
    std::getline(coordDataHandle, sequence);    // advance to sequence line
    coordinateData coordEntry(header, coordinate, sequence);
    coordData.push_back(coordEntry);

    std::getline(coordDataHandle, line);        // advance to line break
  }

  coordDataHandle.close();
  return coordData;                             // move vector
}
