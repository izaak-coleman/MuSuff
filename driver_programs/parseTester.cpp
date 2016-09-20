// parseTester.cpp
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <utility>
#include <vector>
#include <cstring>


using namespace std;

enum {TUMOUR, HEALTHY};

typedef std::pair<std::string, bool> file_and_type;

struct command_line_params {
  std::vector<file_and_type> datafiles;
  uint8_t minimum_suffix_size;
  float econt;
};

void parseCommandLine(command_line_params &params, int argc, char **argv) {

  // Only four param input. If fail count, explain to user options and
  // file format
  if (argc != 4) {
    cout << "Usage: <exec> <datafile_names.txt>"
         << " <minimum_suffix_size> <contamination_ratio>" << endl;

    cout << endl << endl;
    cout << "Please note the the format of headerfile.txt:"
         << endl 
         << "<file1>\t<H/T>" << endl
         << "<file2>\t<H/T>" << endl
         << "..." << endl;
    exit(1);
  }

    // Open the header file that contains the filenames and types
    // of the syntax <filename.fasta> \t <H/T>
    // Parse the file, and load into file_and_type pair

    ifstream headerfile;
    headerfile.open(argv[1]);   // link stream to file

    // file_string holds each line (file metadata)
    string file_string;
    file_and_type file_info;

    while (getline(headerfile, file_string)) {
      cout << file_string << endl;
      
      // momentarily work with c strings, so I can use the glorious strtok()...

      // load the filename into file_and_type...
      char * c_file = const_cast<char*>(file_string.c_str()); // un const
      file_info.first = strtok(c_file, ",\t ");
      cout << file_info.first << endl;

      // load the type bool into file_and_type...
      char * c_datatype = strtok(NULL, ",\t ");   // remainder of cstring
      string datatype = c_datatype;                    // convert to string
      cout << datatype << endl;

      if (datatype == "H") {              // then bool
        file_info.second = HEALTHY;
        cout << file_info.second << endl;
      }
      else if (datatype == "T") {
        file_info.second = TUMOUR;
        cout << file_info.second << endl;
      }
      else {
        cout << datatype << " is not a valid datatype, either H or T. " << endl
             << "Program terminating..." << endl;
        exit(1);
      }

      params.datafiles.push_back(file_info);      // store in params
  }

  // Parse the minimum suffix size parameter, setting it as an int, 
  // also the econt, setting the ratio as a float
  params.minimum_suffix_size = static_cast<uint8_t>(std::stoi(argv[2]) + 1);
  params.econt = std::stof(argv[3]);
}

int main(int argc, char **argv) {

  command_line_params params;

  parseCommandLine(params, argc, argv);

  cout << "minimum_suffix_size: " << (int) params.minimum_suffix_size << endl;
  cout << "ECONT: " << params.econt << endl;


  for(int i=0; i < params.datafiles.size(); i++) {
    cout << "Filename:  " << params.datafiles[i].first << std::boolalpha
         << " ---Type: " << params.datafiles[i].second << endl;
  }
  return 0;
}
