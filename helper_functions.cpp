#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cstring>

#include "Suffix_t.h"
#include "helper_functions.h"
#include "Reads.h"

using namespace std;


void parseCommandLine(command_line_params &params, int argc, char **argv) {

  // Only four param input. If fail count, explain to user options and
  // file format
  if (argc != 4) {
    cout << "Usage: <exec> <datafile_names.txt>"
         << " <minimum_suffix_size> <contamination_ratio>" << endl;

    cout << endl << endl;
    cout << "Please note the the format of headerfile.txt:"
         << endl 
         << "<file1>\t<H/T> <1/2>" << endl
         << "<file2>\t<H/T> <1/2>" << endl
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

  cout << "Parsing inputs..." << endl;
  cout << "Reading filenames from " << argv[1] << " ... " << endl;

  while (getline(headerfile, file_string)) {

    // momentarily work with c strings, to make use of the glorious strtok()...

    // load the filename into file_and_type...
    char * c_file = const_cast<char*>(file_string.c_str()); // un const
    file_info.first = strtok(c_file, ",\t ");

    // load the type bool into file_and_type...
    char * c_datatype = strtok(NULL, ",\t ");      // get next feild
    string datatype = c_datatype;                    // convert to string

    if (datatype == "H") {              // then bool
      file_info.second = HEALTHY;
    }
    else if (datatype == "T") {
      file_info.second = TUMOUR;
    }
    else {
      cout << datatype << " is not a valid datatype, either H or T. " << endl
           << "Program terminating..." << endl;
      exit(1);
    }


    cout << "Storing " << file_info.first << " as "
         << ((file_info.second) ? "healthy" : "tumour") << " data "
         << endl;

    params.datafiles.push_back(file_info);      // store in params
  }

  // Parse the minimum suffix size parameter, setting it as an int, 
  // also the econt, setting the ratio as a float
  params.minimum_suffix_size = static_cast<uint8_t>(std::stoi(argv[2]) + 1);
  params.econt = std::stod(argv[3]);
}

int computeLongestCommonPrefix(Suffix_t &isuf, Suffix_t &jsuf, 
                               ReadsManipulator &reads) {

  // Get suffix pointers in reads
  string::iterator isuf_iter  = reads.returnStartIterator(isuf);
  string::iterator isuf_end   = reads.returnEndIterator(isuf);
  string::iterator jsuf_iter  = reads.returnStartIterator(jsuf);
  string::iterator jsuf_end   = reads.returnEndIterator(jsuf);

  // computes lcp
  int lcp = 0;

  while (*isuf_iter == *jsuf_iter &&
         isuf_iter != isuf_end    &&
         jsuf_iter != jsuf_end      ) {

    lcp++;      // matched next char, so extend lcp

    isuf_iter++; jsuf_iter++;
  }

  return (lcp);   
}
