#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main() {
  ifstream fin;

  fin.open("lexfourTupel.txt");

  string next_tuple;
  int i=0;
  while(fin.good()) {
    getline(fin, next_tuple);

    cout << "if(!strcmp(\"" << next_tuple << "\"," << " seq)) {" << endl;
    cout << "  return " << i << ";" << endl;
    cout << "}" << endl;
    i++;
  }
  return 0;
}



