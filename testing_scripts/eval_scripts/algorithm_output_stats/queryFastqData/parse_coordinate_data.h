// parse_coordinate_data.h
#ifndef PARSE_COORDINATE_DATA
#define PARSE_COORDINATE_DATA

#include <string>
#include <vector>

struct coordinateData {
  std::string header;
  unsigned int coordinate;
  std::string sequence;
  coordinateData(std::string h, unsigned int c, std::string s): header{h},
    coordinate{c}, sequence{s} {}
};

vector<coordinateData> parseCoordinateData(std::string filename);
void printCoordinateData(std::vector<coordinateData> const& coordData);

#endif
