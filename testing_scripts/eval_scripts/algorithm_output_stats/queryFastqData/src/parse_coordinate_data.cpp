// parse_coordinate_data.cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "parse_coordinate_data.h"
#include "string.h"

using namespace std;
static int const COORDINATE = 0;
static int const HBASE = 1;
static int const CBASE = 2;
static int const FLANKING_DIST = 3;

void printCoordinateData(std::vector<coordinateData> const& coordData) {
  for (coordinateData const& c : coordData) {
    std::cout << "Header: " << c.header << std::endl;
    std::cout << "Coord : " << c.coordinate << std::endl;
    std::cout << "Seq   : " << c.sequence << std::endl;
    std::cout << "Health: " << c.hBase << std::endl;
    std::cout << "Cancer: " << c.cBase << std::endl;
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

    // parse header info
    std::size_t endOfHeader = line.find(token);
    std::string header = line.substr(0, endOfHeader);
    std::string rest = line.substr(endOfHeader+token.size());
    std::vector<std::string> otherFields;
    split_string(rest, ":", otherFields);
    unsigned int coordinate = std::stoi(otherFields[COORDINATE]);
    char hBase = otherFields[HBASE][0];
    char cBase = otherFields[CBASE][0];
    unsigned int flanking_dist = std::stoi(otherFields[FLANKING_DIST])


    std::string sequence;
    std::getline(coordDataHandle, sequence);    // advance to sequence line
    coordinateData coordEntry(header, 
                              coordinate, 
                              sequence, 
                              hBase, 
                              cBase
                              flanking_dist);
    coordData.push_back(coordEntry);

    std::getline(coordDataHandle, line);        // advance to line break
  }

  coordDataHandle.close();
  return coordData;                             // move vector
}
