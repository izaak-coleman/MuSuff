#ifndef XA_AS_LINES_H
#define XA_AS_LINES_H
using namespace std;
static const int SAM_EXT = 4;
static const int FLAG = 1;
static const int POS = 3;
static const int SEQ = 9;
static const int CIGAR = 5;
static const int CHR = 2;
static const string XA_PREFIX = "XA";
static const string REV = "16";
static const string FWD = "0";
void xaToSeparateLine(string filename, string &out_filename);
void skipSamHeader(ifstream &samfile);
bool xaPresent(vector<string> &fields);
string getXA(vector<string> &fields);
string constructAlternateSamEntry(vector<string> const &entry_fields, 
                                  string alignment);


#endif
