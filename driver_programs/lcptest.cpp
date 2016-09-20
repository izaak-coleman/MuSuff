#include <string>
#include <vector>
#include <iostream>


using namespace std;
void computeLongestCommonPrefixArray(vector<int> &lcp, vector<string> &strings);

int main() {
  
  string a = "apple";
  string b = "appied";
  string c = "dangerous";
  string d = "danger";

  vector<string> strings;
  strings.reserve(4);
  strings.push_back(a);
  strings.push_back(b);
  strings.push_back(c);
  strings.push_back(d);

  vector<int> lcp;
  lcp.reserve(strings.size() - 1);   // lcp is one less than SA

  computeLongestCommonPrefixArray(lcp, strings);

  for(int v : lcp){
    cout << "Next lcp: " << v << endl;
  }




  return 0;
}

void computeLongestCommonPrefixArray(vector<int> &lcp, vector<string> &strings) {

  for(int i=0, j=1; j < strings.size(); i++, j++){
    int next_lcp = 0;
    // setup iterators
    string::iterator i_iter = strings[i].begin();
    string::iterator i_end  = strings[i].end();
    string::iterator j_iter = strings[j].begin();
    string::iterator j_end  = strings[j].end();

    while(*i_iter == *j_iter &&
          i_iter != i_end &&
          j_iter != j_end ) {
      next_lcp++;
      i_iter++;  j_iter++;
    }
    lcp.push_back(next_lcp);
  }
}
