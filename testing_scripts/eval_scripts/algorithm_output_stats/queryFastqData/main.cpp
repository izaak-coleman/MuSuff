// main.cpp
#include <iostream>
#include <vector>
#include "parse_coordinate_data.h"

using namespace std;
int main(int argc, char **argv) {

  if (argc != 2) {
    std::cout << "Usage: <exe> <coord_data>" << std::endl;
    return -1;
  }

  std::vector<coordinateData> coordData = parseCoordinateData(argv[1]);
  printCoordinateData(coordData);
  return 0;
}
