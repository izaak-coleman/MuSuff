#include <iostream>
#include <vector>
#include <string>


#include "string.h"
#include "SamEntry.cpp"

static const int MUT_CODE = 0;
static const int FLAG = 1;
static const int CHR = 2;
static const int ALIGNMENT = 3;
static const int QUAL = 4;
static const int CIGAR = 5;
static const int INS_SIZE = 8;
static const int AL_CNS = 9;
static const int QSTRING = 10;
static const int MISMATCHES = 18;
static const int SNV = 1;
static const int SSV = 2;
static const int LSV = 3;
static const int MUT_CNS = 4;

static const string NM = "NM";
static const string MD = "MD";
static const string AS = "AS";
static const string BC = "BC";
static const string X0 = "X0";
static const string X1 = "X1";
static const string XN = "XN";
static const string XM = "XM";
static const string XO = "XO";
static const string XG = "XG";
static const string XT = "XT";
static const string XA = "XA";
static const string XS = "XS";
static const string XF = "XF";
static const string XE = "XE";



static const int REVERSE_FLAG = 16;
static const int FORWARD_FLAG = 0;

using namespace std;


SamEntry::SamEntry(string record) {
  vector<string> fields;
  split_string(record, "\t", fields);

  this.mutation_string = fields[MUT_CODE];
  this.flag            = stoi(fields[FLAG]);
  this.chr_num         = stoi(fields[CHR]);
  this.aln_pos         = stoi(fields[ALIGNMENT]);
  this.mapq            = stoi(fields[MAPQ]);
  this.cigar           = fields[CIGAR];
  this.ins_size        = stoi(fields[INS_SIZE]);
  this.healthy_cns     = fields[AL_CNS];
  this.qstring         = fields[QSTRING];
  fields.erase(fields.begin(), fields.begin()+10); // done with mandatory fields

  for(int i=0; i < fields.size(); i++) {
    load_non_mandatory_field(fields[i]);
  }
}

void load_non_mandatory_field(string field) {
  string prefix = field.substr(0,2);

  if (prefix == NM){
    this.nm= field;
  }
  else if (prefix == MD){
    this.md= field;
  }
  else if (prefix == AS){
    this.as= field;
  }
  else if (prefix == BC){
    this.bc= field;
  }
  else if (prefix ==X0){
    this.x0= field;
  }
  else if (prefix == X1({
    this.x1= field;
  }
  else if (prefix == XN){
    this.xn= field;
  }
  else if (prefix == XM){
    this.xm= field;
  }
  else if (prefix == XO){
    this.xo= field;
  }
  else if (prefix == XG){
    this.xg= field;
  }
  else if (prefix == XT){
    this.xt= field;
  }
  else if (prefix == XA){
    this.xa=  field;
  }
  else if (prefix == XS){
    this.xs= field;
  }
  else if (prefix == XF){
    this.xf= field;
  }
  else (prefix == XE) {
    this.xe= field;
  }
}

// mandatory fields
int SamEntry::getFlag() {
  return flag;
}
int SamEntry::getChrNum() {
  return chr_num;
}
int SamEntry::getAlnPos() {
  return aln_pos;
}
int SamEntry::getMapq() {
  return mapq;
}
int SamEntry::getInsSize() {
  return ins_size;
}
string SamEntry::getMutationString() {
  return mutation_string;
}
string SamEntry::getCigarString() {
  return cigar;
}
string SamEntry::getHealthyCns() {
  return healthy_cns;
}
string SamEntry::getQString() {
  return qstring;
}

// augmented sam - bwa optional fields

bool SamEntry::validNM() {
  return !nm.empty()
}
string SamEntry::getNM() {
  return nm;
}

bool SamEntry::validMD() {
  return !md.empty();
}
string SamEntry::getMD() {
  return md;
}

bool SamEntry::validAS() {
  return !ad.empty();
}
string SamEntry::getAS() {
  return as;
}

bool SamEntry::validBC() {
  return !bc.empty();
}
string SamEntry::getBC() {
  return bc;
}

bool SamEntry::validX0() {
  return !x0.empty();
}
string SamEntry::getX0 {
  return x0;
}

bool SamEntry::validX1() {
  return !x1.empty();
}
string SamEntry::getX1() {
  return x1;
}

bool SamEntry::validXN() {
  return !xn.empty();
}
string SamEntry::getXN() {
  return xn;
}
