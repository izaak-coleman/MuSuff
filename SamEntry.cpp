#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>


#include "string.h"
#include "SamEntry.h"

using namespace std;

SamEntry::SamEntry(string record): bwa_fields(NUM_BWA_FIELDS){
  this->record = record;
  re
  vector<string> fields;
  split_string(record, "\t", fields);


  // load mandatory fields
  this->mutation_string = fields[MUT_CODE];
  this->flag            = stoi(fields[FLAG]);
  this->chr_num         = fields[CHR];
  this->aln_pos         = stoi(fields[ALIGNMENT]);
  this->mapq            = stoi(fields[MAPQ]);
  this->cigar           = fields[CIGAR];
  this->ins_size        = stoi(fields[INS_SIZE]);
  this->healthy_cns     = fields[AL_CNS];
  this->qstring         = fields[QSTRING];

  // erase mandatory fields after loading
  fields.erase(fields.begin(), fields.begin()+11); 

  // load non mandatory (bwa specific) fields
  for(int i=0; i < fields.size(); i++) {
    load_non_mandatory_field(fields[i]);
  }
}

void SamEntry::load_non_mandatory_field(string field) {
  string prefix = field.substr(0,2);

  if (prefix == PREFIX_NM){
    bwa_fields[NM] = field;
  }
  else if (prefix == PREFIX_MD){
    bwa_fields[MD] = field;
  }
  else if (prefix == PREFIX_AS){
    bwa_fields[AS] = field;
  }
  else if (prefix == PREFIX_BC){
    bwa_fields[BC] = field;
  }
  else if (prefix == PREFIX_X0){
    bwa_fields[X0] = field;
  }
  else if (prefix == PREFIX_X1){
    bwa_fields[X1] = field;
  }
  else if (prefix == PREFIX_XN){
    bwa_fields[XN] = field;
  }
  else if (prefix == PREFIX_XM){
    bwa_fields[XM] = field;
  }
  else if (prefix == PREFIX_XO){
    bwa_fields[XO] = field;
  }
  else if (prefix == PREFIX_XG){
    bwa_fields[XG] = field;
  }
  else if (prefix == PREFIX_XT){
    bwa_fields[XT] = field;
  }
  else if (prefix == PREFIX_XA){
    bwa_fields[XA] =  field;
  }
  else if (prefix == PREFIX_XS){
    bwa_fields[XS] = field;
  }
  else if (prefix == PREFIX_XF){
    bwa_fields[XF] = field;
  }
  else if (prefix == PREFIX_XE) {
    bwa_fields[XE] = field;
  }
  else {
    cout << "Error parsing sam fields" << endl;
    exit(1);
  }
}


// mandatory fields
int SamEntry::getFlag() {
  return flag;
}
string SamEntry::getChrNum() {
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

bool SamEntry::bwaFieldContainsElem(int field_num) {
  return !bwa_fields[field_num].empty();

}

string SamEntry::getBwaField(int field_num) {
  return bwa_fields[field_num];
}



void SamEntry::printEntry(){
  cout << "Flag: " << getFlag() << endl
  << "Aln pos: " << getAlnPos() << endl
  << "Mapq: " << getMapq() << endl
  << "InsSize: " << getInsSize() << endl
  << "ChrNum: " << getChrNum() << endl
  << "MutString :" << getMutationString() << endl
  << "Cigar: " << getCigarString() << endl
  << "Healthy CNS: " << getHealthyCns() << endl
  << "Qstring: " << getQString() << endl << endl << endl;

  if(bwaFieldContainsElem(NM)) {
    cout << getBwaField(NM) << endl;
  }
  if(bwaFieldContainsElem(MD)) {
    cout << getBwaField(MD) << endl;
  }
  if(bwaFieldContainsElem(AS)) {
    cout << getBwaField(AS) << endl;
  }
  if(bwaFieldContainsElem(BC)) {
    cout << getBwaField(BC) << endl;
  }
  if(bwaFieldContainsElem(X0)) {
    cout << getBwaField(X0) << endl;
  }
  if(bwaFieldContainsElem(X1)) {
    cout << getBwaField(X1) << endl;
  }
  if(bwaFieldContainsElem(XN)) {
    cout << getBwaField(XN) << endl;
  }
  if(bwaFieldContainsElem(XM)) {
    cout << getBwaField(XM) << endl;
  }
  if(bwaFieldContainsElem(XO)) {
    cout << getBwaField(XO) << endl;
  }
  if(bwaFieldContainsElem(XG)) {
    cout << getBwaField(XG) << endl;
  }
  if(bwaFieldContainsElem(XT)) {
    cout << getBwaField(XT) << endl;
  }
  if(bwaFieldContainsElem(XA)) {
    cout << getBwaField(XA) << endl;
  }
  if(bwaFieldContainsElem(XF)) {
    cout << getBwaField(XF) << endl;
  }
  if(bwaFieldContainsElem(XE)) {
    cout << getBwaField(XE) << endl;
  }
}

void SamEntry::printRecord() {
  cout << record << endl;
}

vector<single_snv> * SamEntry::generateReportableSNVs() {

  vector<string> sub_fields;
  split_string(getBwaField(XA), ":", sub_fields); // split XA into colon fields

  vector<string> alt_alignments;
  split_string(sub_fields[2], ";", alt_alignments); // alts split by ;

  for(int i = 0; i < alt_alignments.size(); i++) {
    string str_aln = alt_alignments[i];
    vector<string> aln_fields;
    split_string(str_aln, ",", aln_fields);   // separate alignment fields

    if(aln_fields[0] != "22") {
      continue;     // alt alignment is not to 22, so do not report
    }

    vector<single_snv> reported_snvs = new vector<single_snv>;
    for(int i = 0; i < SNV_pos.size(); i++) {
      single_snv aln;
      aln.chr = stoi(alt_fields[0]);
      int 

    }


  }
}

void SamEntry::skipHeader(ifstream &samfile) {
  boost::regex header("(@).*");   // regex differentiates header / non-header
  string line;
  int end_of_header_line_num = 1;
  while (getline(samfile, line)) {
    if (!boost::regex_match(line, header)) {
      break;
    }
    end_of_header_line_num++;
  }

  samfile.clear();
  samfile.seekg(0, ios::beg); // return to start
  for(; end_of_header_line_num > 1; end_of_header_line_num--) {
    getline(samfile, line);
  }
}
