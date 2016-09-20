#include <iostream>
#include <string>
#include <vector>

using namespace std;


long long int binarySearch(vector<string> &SA, string query);
// Function uses the Myers and Manber 
// O(n * log m) binary search 
// that is in practice, basically O(n + log m)

int lcp(string l, string r, unsigned int min_lr);

bool lexCompare(string l, string r, unsigned int min_lr);

int minVal(int a, int b);

int main() {

  vector<string> suffixes;
  suffixes.push_back("$"                );
  suffixes.push_back("a$"               );
  suffixes.push_back("delta$"           );
  suffixes.push_back("elta$"            );
  suffixes.push_back("idelta$"          );
  suffixes.push_back("ipidelta$"        );
  suffixes.push_back("issipidelta$"     );
  suffixes.push_back("ississipidelta$"  );
  suffixes.push_back("lta$"             );
  suffixes.push_back("mississipidelta$" );
  suffixes.push_back("pidelta$"      );
  suffixes.push_back("sipidelta$"    );
  suffixes.push_back("sissipidelta$" );
  suffixes.push_back("ssipidelta$"   );
  suffixes.push_back("ssissipidelta$");
  suffixes.push_back("ta$");

  string query;
  cout << "enter next search: ";
  cin >> query;

  cout << query << endl;
  while (query != "end") {
    long long int test =  binarySearch(suffixes, query);

    if (test == -1) {   // does not exist in tree
      cout << "substring not in text" << endl;
    }
    else {
      unsigned int index = (unsigned int) test;
      cout << "Index of element is at: " << index << endl;
      cout << suffixes[index] << endl;
    }

    cout << "next query: ";
    cin >> query;
  }


  return 0;
}

long long int binarySearch(vector<string> &SA, string query) {

  // binary search bound pointers
  unsigned int right = SA.size();     // left right BS pointers
  unsigned int left = 0;
  unsigned int  mid;

  // lcp bound pointers
  unsigned int min_left_right;
  unsigned int lcp_left_query;
  unsigned int lcp_right_query;


  // find left and right lcps before loop
  lcp_left_query = lcp(SA[left], query, 0);
  lcp_right_query = lcp(SA[right], query, 0);

  // get the min
  min_left_right = min(lcp_left_query, lcp_right_query);

  // only when both right and left pointers have shifted, do 
  // we need to recalculate the right and left lcp with query
  bool right_shift = false, left_shift = false;
  
  while(left <  right) {

    mid = left + ((right - left) / 2);
    cout << "left: " << left << endl;
    cout << "right: " << right << endl;
    cout << "mid : " << mid << endl;
    cout << endl << endl;

    if(lcp(SA[mid], query, min_left_right) == query.size()) {
      // then we have covered query, by 30bp, and so the suffix is in the
      // correct genomic location
      return mid;
    }

    if(lexCompare(SA[mid], query, min_left_right)) {

      // then query is lower, so move left bound down
      left = mid+1;
      left_shift = true;
    }

    else {
      // then query is higher, so move right bound up

      right = mid;
      right_shift = true;
    }



    if (right_shift && left_shift) {
      // all chars from [0] to min_left_right must already be identical
      // so dont check
      lcp_left_query = lcp(SA[left], query, min_left_right);
      lcp_right_query = lcp(SA[right], query, min_left_right);

      min_left_right = minVal(lcp_left_query, lcp_right_query);
    }
  }


  return -1;                        // no match
}

bool lexCompare(string l, string r, unsigned int min_lr) {
  // Generate pointers to lhs and rhs suffixes in reads

  // * min_lr avoids redundant searches droping bound to, in practice
  // O(n + log m)

  string::iterator l_iter  = l.begin() + min_lr;    
  string::iterator l_end   = l.end();
  string::iterator r_iter  = r.begin() + min_lr;
  string::iterator r_end   = r.end();

  for( ; (l_iter != l_end && r_iter != r_end); l_iter++, r_iter++){
    // lex compare character
    if (*l_iter < *r_iter) { return true; }
    if (*r_iter < *l_iter) { return false; }
    // equiv char so move to next...
    }

    // One is prefix of other, return the prefix as higher suffix
    return (l_iter == l_end) && (r_iter != r_end);
}

int lcp(string l, string r, unsigned int mlr) {
  while (l[mlr] == r[mlr]) {
    mlr++;
  }
  return mlr;
}

int minVal(int a, int b) {
  if(a > b){
    return b;
  }
  return a;
}

