// parse_coordinate_data.h
#ifndef PARSE_COORDINATE_DATA
#define PARSE_COORDINATE_DATA

#include <string>
#include <vector>

struct coordinateData {
  std::string header;
  unsigned int coordinate;
  std::string sequence;
  char hBase;
  char cBase;
  coordinateData(std::string h, unsigned int c, std::string s,
                 char hb, char cb): 
    header{h}, coordinate{c}, sequence{s}, hBase(hb), cBase(cb) {}
};

std::vector<coordinateData> parseCoordinateData(std::string filename);
void printCoordinateData(std::vector<coordinateData> const& coordData);

#endif
