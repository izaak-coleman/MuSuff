#include <string>
#include <vector>
#include <fstream>
#include <cstring>
#include <sys/types.h>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>

#include <iostream>
#include <algorithm>

#include <thread>

#include "radix.h"

using namespace std;

static const int LUT_DISTANCE = 250000;


struct suffix{
  unsigned int read_id;
  int offset;
  bool type;
};

void generateTissueSA(vector<suffix> *TissueSA, vector<string> *reads, 
    bool type);

void generateSA(vector<suffix> &SA, vector<string> &t_reads, 
                vector<string> &h_reads);

string concatenateReads(vector<string> &reads);

void generateBSA(vector<pair<unsigned int, unsigned int> > &BSA, vector<string>
    &reads);

pair<unsigned int, unsigned int> binarySearch(vector<pair<unsigned int, 
    unsigned int> > &BSA, unsigned int suffix_index);




int main(int argc, char** argv) {
  ifstream hfin, tfin;
  tfin.open(argv[1]);

  hfin.open(argv[2]);

  vector<string> tumour_reads, healthy_reads;
  string tread;
  while(tfin.good()) {
    getline(tfin, tread);
    tread += "$";
    tumour_reads.push_back(tread);
  }


  string hread;
  while(hfin.good()) {
    getline(hfin, hread);
    hread += "$";
    healthy_reads.push_back(hread);
  }

  healthy_reads.pop_back();
  tumour_reads.pop_back();




  vector<suffix> genSA;

  generateSA(genSA, tumour_reads, healthy_reads);


  for(suffix s : genSA) {
    if(s.type) {
      cout << healthy_reads[s.read_id].substr(s.offset) << endl;
    cout << "read id: " << s.read_id << endl;
    cout << "offset: " << s.offset << endl;
    cout << "type:" << ((s.type) ? " H" : " T" ) << endl << endl;
    }
    else {
      cout << tumour_reads[s.read_id].substr(s.offset) << endl;
      cout << "read id: " << s.read_id << endl;
    cout << "offset: " << s.offset << endl;
    cout << "type:" << ((s.type) ? " H" : " T" ) << endl << endl;
    }

  }


  return 0;
}

void generateSA(vector<suffix> &SA, vector<string> &t_reads, vector<string>
    &h_reads) {

  vector<suffix> healthy_SA;
  vector<suffix> tumour_SA;

  // build local suffix arrays on two threads in parallel
  vector<thread> tissue_SA_threads;
  tissue_SA_threads.push_back(
      std::thread(generateTissueSA, &healthy_SA, &h_reads, true));

  tissue_SA_threads.push_back(
      std::thread(generateTissueSA, &tumour_SA, &t_reads, false));

  for(auto &thread : tissue_SA_threads) {
    thread.join();
  }

  

  cout << "size of healhty SA: " << healthy_SA.size() << endl;
  cout << "size of tumour SA: " << tumour_SA.size() << endl;


  SA.reserve(50 * (t_reads.size() + h_reads.size())); // make room

  bool end_of_healthy = false;
  bool end_of_tumour = false; // not for long...lets hope ;)

  unsigned int tind=0, hind=0;

  cout << "Merging tissue general suffix arrays to build final suffix " << 
    "array " << endl;


  while ((tind < tumour_SA.size()) || (hind < healthy_SA.size())) {

    // bound checks...
    if(hind == healthy_SA.size()) {
      end_of_healthy = true;
    }

    if(tind == tumour_SA.size()) {
      end_of_tumour = true;
    }


    // one has reached end, so add all of other
    if(end_of_healthy) {
      SA.push_back(tumour_SA[tind]);
      tind++;
    }

    else if(end_of_tumour) {
      SA.push_back(healthy_SA[hind]);
      hind++;
    }

    else {  // noone reached end so add based on lex order
      string::iterator h_start, h_end, t_start, t_end;

      t_start = t_reads[tumour_SA[tind].read_id].begin() +
        tumour_SA[tind].offset;
      t_end= t_reads[tumour_SA[tind].read_id].end();

      h_start = h_reads[healthy_SA[hind].read_id].begin() + 
        healthy_SA[hind].offset;

      h_end = h_reads[healthy_SA[hind].read_id].end();

      if(lexicographical_compare(t_start, t_end, h_start, h_end)) {
        SA.push_back(tumour_SA[tind]);
        tind++;
      }
      else {
        SA.push_back(healthy_SA[hind]);
        hind++;
      }

    }

  }

  cout << "Done with suffix array build" << endl;

}

void generateTissueSA(vector<suffix> *TissueSA, vector<string> *reads, bool type) {

  // concatenate all reads to one giant read
  string concat = concatenateReads(*reads);


  // generate the Binary Search Array used to transform output from 
  // suffix array (one string) to generalized suffix array
  // (multi string suffix array)
  vector<pair<unsigned int, unsigned int>> readToConcatMap;
  generateBSA(readToConcatMap, *reads);


  // Builds SA
  cout << "Generating " << ((type) ? "healthy" : "tumour" ) << " suffix" << 
    " array. " << endl;

  unsigned long long *SA;
  SA = Radix<unsigned long long>((uchar*) concat.c_str(), concat.size()).build();



  cout << "Generating " << ((type) ? "healthy" : "tumour" ) << 
    " general suffix array. " << endl;

  TissueSA->reserve(reads->size() * 50);


  // generate generalizedSA
  for(unsigned int i=0; i < concat.size(); i++) {

    pair<unsigned int, unsigned int> read_concat_tup = 
                       binarySearch(readToConcatMap, SA[i]);



    // suffixes with length >= 30 only
    if((SA[i] - read_concat_tup.second) < ((*reads)[read_concat_tup.first].size() -
        30)) {
      suffix s;
      s.read_id = read_concat_tup.first;
      s.offset = SA[i] - read_concat_tup.second; 
      s.type = type;

      TissueSA->push_back(s);
    }
    else {    // add suffix oo generalizedSA
      continue;
    }
  }

  TissueSA->shrink_to_fit();

  delete [] SA;
}


pair<unsigned int, unsigned int> 
            binarySearch(vector<pair<unsigned int, unsigned int> > &BSA, 
                          unsigned int suffix_index) {

  unsigned int right = BSA.size();
  unsigned int left = 0;
  unsigned int mid;





  while(left < right) {
    mid = left + ((right - left) / 2);

    if(suffix_index == BSA[mid].second) {
        return BSA[mid];
    }
    else if (suffix_index > BSA[mid].second) {
      left = mid+1;
    }
    else {
      right = mid;
    }
  }
  left--;
  return BSA[left]; // left should be on the seq
}




void generateBSA(vector<pair<unsigned int, unsigned int>> &BSA, vector<string> &reads) {

  
  unsigned int index_in_concat = 0;

  BSA.reserve(reads.size());
  pair<unsigned int, unsigned int> zeroPos(0,0);
  BSA.push_back(zeroPos);

  for(unsigned int i=1; i < reads.size(); i++) { 
    pair<unsigned int, unsigned int> read_total;
    index_in_concat += reads[i-1].size();

    read_total.first = i;
    read_total.second = index_in_concat;

    BSA.push_back(read_total);
  }

}


string concatenateReads(vector<string> &reads) {
  string concat = "";
  for(string &read : reads) {
    concat += read; 
  }

  return concat;
}

