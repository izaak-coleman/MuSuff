#include <iostream>
#include <string>

using namespace std;

string reverseComplementString(string s){
  string revcomp = "";

  for(int i = s.size()-1; i >= 0; i--) {
  // travel in reverse and switch for complementary
    switch(s[i]) {

      case 'A':{
        revcomp += "T";
        break;
       }

      case 'T':{
        revcomp += "A";
        break;
      }

      case 'C':{
        revcomp += "G";
        break;
      }

      case 'G':{
        revcomp += "C";
        break;
      }
    }
  }

  return revcomp;

}

int main() {
  string s;
  cout << "Enter string: " ;
  cin >> s;
  cout << reverseComplementString(s);
  return 0;
}
