#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include "string.h"

using namespace std;


void split_string(string s, string tokens, vector<string> &split_strings) {
  char *c_s = const_cast<char*>(s.c_str());
  char *c_tokens = const_cast<char*>(tokens.c_str());
  char *c_split = strtok(c_s, c_tokens); // split string into delimited cstrings
  while (c_split != NULL) {
    split_strings.push_back(c_split);
    c_split = strtok(NULL, c_tokens);
  }
}
