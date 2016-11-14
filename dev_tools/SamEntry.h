#ifndef SAMENTRY_H
#define SAMENTRY_H

#include <vector>
#include <string>

static const int NUM_BWA_FIELDS = 15;
static const int REVERSE_FLAG = 16;
static const int FORWARD_FLAG = 0;
static const int SNV = 1;
static const int SSV = 2;
static const int LSV = 3;
static const int MUT_CNS = 4;

// sam compulsary indexes
static const int MUT_CODE = 0;
static const int FLAG = 1;
static const int CHR = 2;
static const int ALIGNMENT = 3;
static const int MAPQ = 4;
static const int CIGAR = 5;
static const int INS_SIZE = 8;
static const int AL_CNS = 9;
static const int QSTRING = 10;
static const int MISMATCHES = 18;


// sam augmented (bwa uniques) indexes
static const int NM = 0;
static const int MD = 1;
static const int AS = 2;
static const int BC = 3;
static const int X0 = 4;
static const int X1 = 5;
static const int XN = 6;
static const int XM = 7;
static const int XO = 8;
static const int XG = 9;
static const int XT = 10;
static const int XA = 11;
static const int XS = 12;
static const int XF = 13;
static const int XE = 14;

// bwa field prefixes
static const std::string PREFIX_NM = "NM";
static const std::string PREFIX_MD = "MD";
static const std::string PREFIX_AS = "AS";
static const std::string PREFIX_BC = "BC";
static const std::string PREFIX_X0 = "X0";
static const std::string PREFIX_X1 = "X1";
static const std::string PREFIX_XN = "XN";
static const std::string PREFIX_XM = "XM";
static const std::string PREFIX_XO = "XO";
static const std::string PREFIX_XG = "XG";
static const std::string PREFIX_XT = "XT";
static const std::string PREFIX_XA = "XA";
static const std::string PREFIX_XS = "XS";
static const std::string PREFIX_XF = "XF";
static const std::string PREFIX_XE = "XE";

class SamEntry {
private:

  // whole record
  std::string record;
  // sam mandatory fields
  int flag;
  int aln_pos;
  int mapq;
  int ins_size;
  std::string chr_num;
  std::string mutation_string;
  std::string cigar;
  std::string healthy_cns;
  std::string qstring;

  // augmented sam - bwa optional fields
  std::vector<std::string> bwa_fields;

  void load_non_mandatory_field(std::string field);
  // function loads a field from the sam record into the
  // corresponding position in bwa_fields


public:

  SamEntry(std::string record);
  int getFlag();
  int getAlnPos();
  int getMapq();
  int getInsSize();

  std::string getChrNum();
  std::string getMutationString();
  std::string getCigarString();
  std::string getHealthyCns();
  std::string getQString();

  bool bwaFieldContainsElem(int field_num);
  // return whether the field contains element i.e is valid, 
  // or not (is empty string);

  std::string getBwaField(int field_num);
  // returns the fiedl at index field_num from the bwa unique sam fields

  void printEntry();
  // Print the SamEntry info

  void printRecord();
};



#endif

