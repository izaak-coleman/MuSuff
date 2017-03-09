#ifndef SAMENTRY_H
#define SAMENTRY_H

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <boost/any.hpp>

class SamEntry {
  // Class takes a sam entry string from the file,
  // and extracts the fields. 
  // It also stores the information pertaining to the number
  // mutations idenified in the particular sam entry

  // Fields
public:

  // COMPULSORY SAM FIELDS
  static const int QNAME;  // <int, string> K, V
  static const int FLAG;  // <int, int> K, V
  static const int RNAME;  // <int, string> K, V
  static const int POS;  // <int, int> K, V
  static const int MAPQ;  // <int, int> K, V
  static const int CIGAR;  // <int, string> K, V
  static const int RNEXT;  // <int, string> K, V
  static const int PNEXT;  // <int, string> K, V
  static const int TLEN;  // <int, string> K, V
  static const int SEQ;  // <int, string> K, V
  static const int QUAL; // <int, string> K, V

  // BWA OPT FIELDS
  static const int NM;    // <int, string> K, V 
  static const int MD;    // <int, string> K, V 
  static const int AS;    // <int, string> K, V 
  static const int BC;    // <int, string> K, V 
  static const int X0;    // <int, string> K, V 
  static const int X1;    // <int, string> K, V 
  static const int XN;    // <int, string> K, V 
  static const int XM;    // <int, string> K, V 
  static const int XO;    // <int, string> K, V 
  static const int XG;    // <int, string> K, V 
  static const int XT;    // <int, string> K, V 
  static const int XA;    // <int, string> K, V 



  static const int XS;    // <int, string> K, V 
  static const int XF;    // <int, string> K, V 
  static const int XE;    // <int, string> K, V 

  // ICSMuFin FIELDS
  static const int LEFT_OHANG;   // <int, int> K,V
  static const int RIGHT_OHANG;  // <int, int> K, V
  static const int BLOCK_ID;     // <int, int> K, V
  static const int CANCER_SEQ;   // <int, string> K, V

  SamEntry(std::string const& entry); // information parsed from this




  // substript access wrapper for fields (map<int, boost::any>)
  template <typename RT>
  RT get(int key) {

    try {
      return boost::any_cast<RT> (fields[key]);
    }
    catch(...) {
      std::cout << "Exception occured" << std::endl
                << "Likely a cast to incorrect type, or" << std::endl
                << "you tried to access an absent key." << std::endl;

    }
  }


private:
  // Fields
  std::vector<int> SNVLocations; // SNV index relative to QUAL string
  unsigned int pair_id;          // not sure if this is nec.

  std::map<int, boost::any> fields; // contains all the sam data fields
                                    // accessed via the static consts

  std::string startsWith(std::string const& tok, std::vector<std::string> const&
      fields);


};
#endif
