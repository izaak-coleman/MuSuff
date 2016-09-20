#include <vector>
#include <string>
#include <iostream>

#include "helper_functions.h"
#include "Reads.h"

using namespace std;


bool tissueCheck(string type) {
  if (type == "H") {
    return true;
  }
  if (type == "T") {
    return false;
  }
}
int main(int argc, char **argv){

  file_and_type file1(argv[1], tissueCheck(argv[2]));
  file_and_type file2(argv[3], tissueCheck(argv[4]));
  file_and_type file3(argv[5], tissueCheck(argv[6]));

  command_line_params params;
  params.datafiles.push_back(file1);
  params.datafiles.push_back(file2);
  params.datafiles.push_back(file3);

  cout << std::boolalpha << "File1: " << params.datafiles[0].first << 
    " - " << params.datafiles[0].second << endl;
  cout << std::boolalpha << "File2: " << params.datafiles[1].first << 
    " - " << params.datafiles[1].second << endl;
  cout << std::boolalpha << "File3: " << params.datafiles[2].first << 
    " - " << params.datafiles[2].second << endl;
  cout << "Params size: " <<  params.datafiles.size() << endl;
  params.minimum_suffix_size = static_cast<uint8_t>(std::stoi(argv[7]) + 1);

  ReadsManipulator reads(params);
  reads.printReads();


  return 0;
}
