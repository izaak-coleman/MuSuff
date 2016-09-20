#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <fstream>



using namespace std;


// deals with command line input
struct command_line_params {
  string filename;  // datafile name NOTE: will have a tumour and healthy
  uint8_t minimum_suffix_size;  // define the lower bound length of prefixes
};



void parseCommandLine(command_line_params &params, 
                                     int argc, char **argv);
// reads in char **argv[] command line options and sets them as
// struct datafeilds of the correct type

void loadWordsFromFile(vector<string> &reads, string filname);
// load toy words into Text

struct suffix_t{
  int read_id;
  uint8_t offset;
  bool is_common = false;
};
typedef vector<suffix_t> sarray;




bool lexCompare(suffix_t &lhs, suffix_t &rhs, vector<string> &reads);



int main(int argc, char **argv){
  std::cout << std::boolalpha;


  // process c-line args
  command_line_params params;
  parseCommandLine(params, argc, argv);
 
  /* load texts into the text vector*/
  vector<string> reads;    // text
  loadWordsFromFile(reads, params.filename);        // add toy words

  suffix_t test_suffix_one, test_suffix_two;

  test_suffix_one.read_id = 1;
  test_suffix_one.offset = 1;

  test_suffix_two.read_id = 1;
  test_suffix_two.offset = 1;

  cout << "test_suffix_one.is_common: " << test_suffix_one.is_common << endl;
  cout << "test_suffix_two.is_common: " << test_suffix_two.is_common << endl;
  cout << "test_suffix_two.is_common: " << test_suffix_one.count<< endl;
  cout << "test_suffix_two.is_common: " << test_suffix_two.count<< endl;

  cout << lexCompare(test_suffix_one, test_suffix_two, reads) << endl;

  cout << "test_suffix_one.is_common: " << test_suffix_one.is_common << endl;
  cout << "test_suffix_two.is_common: " << test_suffix_two.is_common << endl;
  cout << "test_suffix_two.is_common: " << test_suffix_one.count<< endl;
  cout << "test_suffix_two.is_common: " << test_suffix_two.count<< endl;

  string::iterator str1 = returnStartIterator(test_suffix_one, reads);
  string::iterator str2 = returnEndIterator(test_suffix_one, reads);
  while(str1 != str2) {
    cout << *str1 << endl;
    str1++;
  }



  /*----------------------- CPP REF lexocomp -------------------- */
  cout << endl << endl << "CPP REF lexocomp" << endl << endl;

  string a = reads[test_suffix_one.read_id];
  string b = reads[test_suffix_two.read_id];



	std::string::iterator a_it = a.begin();
	std::string::iterator aend = a.end();
	std::string::iterator b_it = b.begin()+1;
	std::string::iterator bend = b.end();

	cout << lexicographical_compare(a_it, aend, b_it, bend) << endl;

	return 0;
}



bool lexCompare(suffix_t &lhs, suffix_t &rhs, vector<string> &reads) {
  // Generate pointers to lhs and rhs suffixes in reads
  string::iterator lhs_iter = returnStartIterator(lhs, reads);
  string::iterator lhs_end   = returnEndIterator(lhs, reads);
  string::iterator rhs_iter = returnStartIterator(rhs, reads);
  string::iterator rhs_end   = returnEndIterator(rhs, reads);

  for( ; (lhs_iter != lhs_end && rhs_iter != rhs_end); lhs_iter++, rhs_iter++){
    // lex compare character
    if (*lhs_iter < *rhs_iter) { return true; }
    if (*rhs_iter < *lhs_iter) { return false; }
    // equiv char so move to next...
    }

    // Identical suffixes, set is_common to true
    if (lhs_iter == lhs_end && rhs_iter == rhs_end) {
      lhs.is_common = true; rhs.is_common = true;
      lhs.count++; rhs.count++;
    } 
    // One is prefix of other, return the prefix as higher suffix
    return (lhs_iter == lhs_end) && (rhs_iter != rhs_end);
}

string::iterator returnStartIterator(suffix_t &suf, vector<string> &reads) {
  // Use suf.read_id to locate the read, and then set an iterator
  // pointing to the suf's offset
  string::iterator iter = reads[suf.read_id].begin() + suf.offset;
  return iter;
}
string::iterator returnEndIterator(suffix_t &suf, vector<string> &reads) {
  // use suf.read_id to locate the read, then return an iterator to the 
  // end of that read
  string::iterator iter = reads[suf.read_id].end();
  return iter;
}

void parseCommandLine(command_line_params &params, 
                                     int argc, char **argv){
  if (argc != 3) {
    cout << "Usage: <exec> <datafile> <minimum_prefix_size>" << endl;
    exit(1);
  }
  params.filename = static_cast<string>(argv[1]);
  params.minimum_suffix_size = static_cast<uint8_t>(std::stoi(argv[2]));
}

void loadWordsFromFile(vector<string> &reads, string filename){
  ifstream datafile;
  datafile.open(filename);
  string next_read;
  while (datafile.good()) {
    getline(datafile, next_read);
    reads.push_back(next_read);
  }
}
