#include <iostream>
#include <string>
#include <vector>

using namespace std;

void unifyReverseComplementaryReads();

bool thirtyBasePairOverlap(string lhs, string rhs);

string reverseComplementString(string s);

int main(){

  string dna =
    "AGTATTCAATATTCATAAAAAATAGTAGAAAATGTAGATTACCGTTCGGCAAAATTTGCCGGCTTGTCGA";
  string dna_rev =
    "TTCTACTATTTTTTATGAATATTGAATACT";

  cout << reverseComplementString(dna_rev) << endl;

  if(thirtyBasePairOverlap(dna, dna_rev)){
    cout << "overlap" << endl;
  }
  return 0;
}

//void unifyReverseComplementaryReads() {
//
//  for(unsigned int i = 0; i < ComplementaryUnfied->size(); i++ ) {
//
//    // if the first read of the ith block of ComplementryUnified 
//    // is oposite pair of the first read from the jth block of ComplementryUnified
//
//    for(unsigned int j = i; j < ComplementaryUnfied->size(); j++) {
//      if (i == j) {continue;} // dont self compare
//
//      // if pair type of j and i is the same, we dont want to compare them
//      std::set<pair<unsigned int>>::iterator i_it, j_it;
//      i_it = (*ComplementaryUnfied)[i].begin();
//      j_it = (*ComplementaryUnfied)[j].begin();
//
//      if( (*reads->TumourMateOrder)[*i_it] ==
//          (*reads->TumourMateOrder)[*j_it]) {
//        continue;
//      }
//
//
//      // groups are opposite pair type, so need rev comp comparison
//
//      if(thirtyBasePairOverlap(*i_it, *j_it) {
//
//        // unify groups
//        // unify to i
//        for (;j_it != (*ComplementaryUnfied)[j].end(); j_it++) {
//          (*ComplementaryUnfied)[i].insert(*j_it);  // add unique elems j to i
//        }
//
//        // unify to j
//        for (;i_it != (*ComplementaryUnfied)[i].end(); i_it++) {
//          (*ComplementaryUnfied)[j].insert(*i_it);  // add unique elems i to j
//        }
//      }
//    }
//  }
//
//
//  // one compared reverse complementary sequences, make unique groups
//  // by passing the groups through a set
//
//  set<set<unsigned int>> filter_uniques;
//  for(int i=0; i < ComplementaryUnfied->size(); i++) {
//    filter_uniques.insert(ComplementaryUnfied->push_back(i));
//  }
//
//
//
//
//  // after passed through set groups unique, so place into vector
//
//  std::set<set<unsigned int>>::iterator group_it;
//  group_it = ComplementaryUnfied->begin();
//
//  vector<pair<bool, unsigned int>> next_read;
//  set<unsigned int>::iterator read_it;
//
//  for(; group_it != ComplementaryUnfied->end(); group_it++) {
//    read_it = group_it->begin();
//
//    for(; read_it != group_it->end(); read_it++) {
//      next_read.first = TUMOUR;     // pass tissue type
//      next_read.second = *read_it;  // pass index
//
//      BreakpointBlocks->push_back;
//    }
//  }
//
//}

bool thirtyBasePairOverlap(string lhs, string rhs) {

  string rhs_read_rev_comp =  reverseComplementString(rhs);

  for(int i=0; i < lhs.size() - 30; i++) {
    string lhs_sub = lhs.substr(i, 30);

    // if we find a 30bp match...
    if(rhs_read_rev_comp.find(lhs_sub) != std::string::npos) {
      return true;    // report match found
    }
  }

  return false;   // didnt find a match
}

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



