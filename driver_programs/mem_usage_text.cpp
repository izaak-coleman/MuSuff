#include <sys/time.h>
#include <sys/resource.h>

#include <iostream>

#include <vector>

using namespace std;

int main() {

 struct rusage sys_usage;

 getrusage(RUSAGE_SELF, &sys_usage);
 vector<int> stuff;

 for(int i=0; i < 1000000000; i++) {
  stuff.push_back(i);
 }
 for(int &s : stuff) {
   s = s*199;
 }
 cout << endl;
 getrusage(RUSAGE_SELF, &sys_usage);

 cout << "Max res set: " << sys_usage.ru_maxrss << endl;

 

 return 0;
}
