#include <string>
#include <iostream>
#include <sys/types.h>
#include "radix.h"


using namespace std;



int main() {
  string str = "mapsmopsmopsbanana";

  unsigned long long *SA;


  //*data = const_cast<char*>(str.c_str());

  SA = Radix<unsigned long long>(str.c_str(), str.size()).build();

  for(int i=0; i < str.size(); i++) {
    cout << SA[i] << endl;
    cout << str.substr(SA[i]) << endl;
  }

  delete [] SA;

}
