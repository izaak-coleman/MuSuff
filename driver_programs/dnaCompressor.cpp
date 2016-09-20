#include <iostream>
#include <string>
#include <cstring>

using namespace std;

uint8_t converter(string sequence);

int main() {

  string dna;
  cout << "Enter DNA: ";
  cin >> dna;


  for(int i=0; i < dna.size()-3; i = i + 4) {
    string chunk = dna.substr(i, 4);
    cout << (int) converter(chunk);
  }
  cout << endl;
  return 0;
}

uint8_t converter(string sequence) {
  int sum = 0;
  for(int i=0; i < sequence.size(); i++) {
    sum += sequence[i];
  }
  const char * seq = sequence.c_str();

  if(!strcmp("AAAA", seq)) {
    return 0;
  }
  if(!strcmp("AAAC", seq)) {
    return 1;
  }
  if(!strcmp("AAAG", seq)) {
    return 2;
  }
  if(!strcmp("AAAT", seq)) {
    return 3;
  }
  if(!strcmp("AACA", seq)) {
    return 4;
  }
  if(!strcmp("AACC", seq)) {
    return 5;
  }
  if(!strcmp("AACG", seq)) {
    return 6;
  }
  if(!strcmp("AACT", seq)) {
    return 7;
  }
  if(!strcmp("AAGA", seq)) {
    return 8;
  }
  if(!strcmp("AAGC", seq)) {
    return 9;
  }
  if(!strcmp("AAGG", seq)) {
    return 10;
  }
  if(!strcmp("AAGT", seq)) {
    return 11;
  }
  if(!strcmp("AATA", seq)) {
    return 12;
  }
  if(!strcmp("AATC", seq)) {
    return 13;
  }
  if(!strcmp("AATG", seq)) {
    return 14;
  }
  if(!strcmp("AATT", seq)) {
    return 15;
  }
  if(!strcmp("ACAA", seq)) {
    return 16;
  }
  if(!strcmp("ACAC", seq)) {
    return 17;
  }
  if(!strcmp("ACAG", seq)) {
    return 18;
  }
  if(!strcmp("ACAT", seq)) {
    return 19;
  }
  if(!strcmp("ACCA", seq)) {
    return 20;
  }
  if(!strcmp("ACCC", seq)) {
    return 21;
  }
  if(!strcmp("ACCG", seq)) {
    return 22;
  }
  if(!strcmp("ACCT", seq)) {
    return 23;
  }
  if(!strcmp("ACGA", seq)) {
    return 24;
  }
  if(!strcmp("ACGC", seq)) {
    return 25;
  }
  if(!strcmp("ACGG", seq)) {
    return 26;
  }
  if(!strcmp("ACGT", seq)) {
    return 27;
  }
  if(!strcmp("ACTA", seq)) {
    return 28;
  }
  if(!strcmp("ACTC", seq)) {
    return 29;
  }
  if(!strcmp("ACTG", seq)) {
    return 30;
  }
  if(!strcmp("ACTT", seq)) {
    return 31;
  }
  if(!strcmp("AGAA", seq)) {
    return 32;
  }
  if(!strcmp("AGAC", seq)) {
    return 33;
  }
  if(!strcmp("AGAG", seq)) {
    return 34;
  }
  if(!strcmp("AGAT", seq)) {
    return 35;
  }
  if(!strcmp("AGCA", seq)) {
    return 36;
  }
  if(!strcmp("AGCC", seq)) {
    return 37;
  }
  if(!strcmp("AGCG", seq)) {
    return 38;
  }
  if(!strcmp("AGCT", seq)) {
    return 39;
  }
  if(!strcmp("AGGA", seq)) {
    return 40;
  }
  if(!strcmp("AGGC", seq)) {
    return 41;
  }
  if(!strcmp("AGGG", seq)) {
    return 42;
  }
  if(!strcmp("AGGT", seq)) {
    return 43;
  }
  if(!strcmp("AGTA", seq)) {
    return 44;
  }
  if(!strcmp("AGTC", seq)) {
    return 45;
  }
  if(!strcmp("AGTG", seq)) {
    return 46;
  }
  if(!strcmp("AGTT", seq)) {
    return 47;
  }
  if(!strcmp("ATAA", seq)) {
    return 48;
  }
  if(!strcmp("ATAC", seq)) {
    return 49;
  }
  if(!strcmp("ATAG", seq)) {
    return 50;
  }
  if(!strcmp("ATAT", seq)) {
    return 51;
  }
  if(!strcmp("ATCA", seq)) {
    return 52;
  }
  if(!strcmp("ATCC", seq)) {
    return 53;
  }
  if(!strcmp("ATCG", seq)) {
    return 54;
  }
  if(!strcmp("ATCT", seq)) {
    return 55;
  }
  if(!strcmp("ATGA", seq)) {
    return 56;
  }
  if(!strcmp("ATGC", seq)) {
    return 57;
  }
  if(!strcmp("ATGG", seq)) {
    return 58;
  }
  if(!strcmp("ATGT", seq)) {
    return 59;
  }
  if(!strcmp("ATTA", seq)) {
    return 60;
  }
  if(!strcmp("ATTC", seq)) {
    return 61;
  }
  if(!strcmp("ATTG", seq)) {
    return 62;
  }
  if(!strcmp("ATTT", seq)) {
    return 63;
  }
  if(!strcmp("CAAA", seq)) {
    return 64;
  }
  if(!strcmp("CAAC", seq)) {
    return 65;
  }
  if(!strcmp("CAAG", seq)) {
    return 66;
  }
  if(!strcmp("CAAT", seq)) {
    return 67;
  }
  if(!strcmp("CACA", seq)) {
    return 68;
  }
  if(!strcmp("CACC", seq)) {
    return 69;
  }
  if(!strcmp("CACG", seq)) {
    return 70;
  }
  if(!strcmp("CACT", seq)) {
    return 71;
  }
  if(!strcmp("CAGA", seq)) {
    return 72;
  }
  if(!strcmp("CAGC", seq)) {
    return 73;
  }
  if(!strcmp("CAGG", seq)) {
    return 74;
  }
  if(!strcmp("CAGT", seq)) {
    return 75;
  }
  if(!strcmp("CATA", seq)) {
    return 76;
  }
  if(!strcmp("CATC", seq)) {
    return 77;
  }
  if(!strcmp("CATG", seq)) {
    return 78;
  }
  if(!strcmp("CATT", seq)) {
    return 79;
  }
  if(!strcmp("CCAA", seq)) {
    return 80;
  }
  if(!strcmp("CCAC", seq)) {
    return 81;
  }
  if(!strcmp("CCAG", seq)) {
    return 82;
  }
  if(!strcmp("CCAT", seq)) {
    return 83;
  }
  if(!strcmp("CCCA", seq)) {
    return 84;
  }
  if(!strcmp("CCCC", seq)) {
    return 85;
  }
  if(!strcmp("CCCG", seq)) {
    return 86;
  }
  if(!strcmp("CCCT", seq)) {
    return 87;
  }
  if(!strcmp("CCGA", seq)) {
    return 88;
  }
  if(!strcmp("CCGC", seq)) {
    return 89;
  }
  if(!strcmp("CCGG", seq)) {
    return 90;
  }
  if(!strcmp("CCGT", seq)) {
    return 91;
  }
  if(!strcmp("CCTA", seq)) {
    return 92;
  }
  if(!strcmp("CCTC", seq)) {
    return 93;
  }
  if(!strcmp("CCTG", seq)) {
    return 94;
  }
  if(!strcmp("CCTT", seq)) {
    return 95;
  }
  if(!strcmp("CGAA", seq)) {
    return 96;
  }
  if(!strcmp("CGAC", seq)) {
    return 97;
  }
  if(!strcmp("CGAG", seq)) {
    return 98;
  }
  if(!strcmp("CGAT", seq)) {
    return 99;
  }
  if(!strcmp("CGCA", seq)) {
    return 100;
  }
  if(!strcmp("CGCC", seq)) {
    return 101;
  }
  if(!strcmp("CGCG", seq)) {
    return 102;
  }
  if(!strcmp("CGCT", seq)) {
    return 103;
  }
  if(!strcmp("CGGA", seq)) {
    return 104;
  }
  if(!strcmp("CGGC", seq)) {
    return 105;
  }
  if(!strcmp("CGGG", seq)) {
    return 106;
  }
  if(!strcmp("CGGT", seq)) {
    return 107;
  }
  if(!strcmp("CGTA", seq)) {
    return 108;
  }
  if(!strcmp("CGTC", seq)) {
    return 109;
  }
  if(!strcmp("CGTG", seq)) {
    return 110;
  }
  if(!strcmp("CGTT", seq)) {
    return 111;
  }
  if(!strcmp("CTAA", seq)) {
    return 112;
  }
  if(!strcmp("CTAC", seq)) {
    return 113;
  }
  if(!strcmp("CTAG", seq)) {
    return 114;
  }
  if(!strcmp("CTAT", seq)) {
    return 115;
  }
  if(!strcmp("CTCA", seq)) {
    return 116;
  }
  if(!strcmp("CTCC", seq)) {
    return 117;
  }
  if(!strcmp("CTCG", seq)) {
    return 118;
  }
  if(!strcmp("CTCT", seq)) {
    return 119;
  }
  if(!strcmp("CTGA", seq)) {
    return 120;
  }
  if(!strcmp("CTGC", seq)) {
    return 121;
  }
  if(!strcmp("CTGG", seq)) {
    return 122;
  }
  if(!strcmp("CTGT", seq)) {
    return 123;
  }
  if(!strcmp("CTTA", seq)) {
    return 124;
  }
  if(!strcmp("CTTC", seq)) {
    return 125;
  }
  if(!strcmp("CTTG", seq)) {
    return 126;
  }
  if(!strcmp("CTTT", seq)) {
    return 127;
  }
  if(!strcmp("GAAA", seq)) {
    return 128;
  }
  if(!strcmp("GAAC", seq)) {
    return 129;
  }
  if(!strcmp("GAAG", seq)) {
    return 130;
  }
  if(!strcmp("GAAT", seq)) {
    return 131;
  }
  if(!strcmp("GACA", seq)) {
    return 132;
  }
  if(!strcmp("GACC", seq)) {
    return 133;
  }
  if(!strcmp("GACG", seq)) {
    return 134;
  }
  if(!strcmp("GACT", seq)) {
    return 135;
  }
  if(!strcmp("GAGA", seq)) {
    return 136;
  }
  if(!strcmp("GAGC", seq)) {
    return 137;
  }
  if(!strcmp("GAGG", seq)) {
    return 138;
  }
  if(!strcmp("GAGT", seq)) {
    return 139;
  }
  if(!strcmp("GATA", seq)) {
    return 140;
  }
  if(!strcmp("GATC", seq)) {
    return 141;
  }
  if(!strcmp("GATG", seq)) {
    return 142;
  }
  if(!strcmp("GATT", seq)) {
    return 143;
  }
  if(!strcmp("GCAA", seq)) {
    return 144;
  }
  if(!strcmp("GCAC", seq)) {
    return 145;
  }
  if(!strcmp("GCAG", seq)) {
    return 146;
  }
  if(!strcmp("GCAT", seq)) {
    return 147;
  }
  if(!strcmp("GCCA", seq)) {
    return 148;
  }
  if(!strcmp("GCCC", seq)) {
    return 149;
  }
  if(!strcmp("GCCG", seq)) {
    return 150;
  }
  if(!strcmp("GCCT", seq)) {
    return 151;
  }
  if(!strcmp("GCGA", seq)) {
    return 152;
  }
  if(!strcmp("GCGC", seq)) {
    return 153;
  }
  if(!strcmp("GCGG", seq)) {
    return 154;
  }
  if(!strcmp("GCGT", seq)) {
    return 155;
  }
  if(!strcmp("GCTA", seq)) {
    return 156;
  }
  if(!strcmp("GCTC", seq)) {
    return 157;
  }
  if(!strcmp("GCTG", seq)) {
    return 158;
  }
  if(!strcmp("GCTT", seq)) {
    return 159;
  }
  if(!strcmp("GGAA", seq)) {
    return 160;
  }
  if(!strcmp("GGAC", seq)) {
    return 161;
  }
  if(!strcmp("GGAG", seq)) {
    return 162;
  }
  if(!strcmp("GGAT", seq)) {
    return 163;
  }
  if(!strcmp("GGCA", seq)) {
    return 164;
  }
  if(!strcmp("GGCC", seq)) {
    return 165;
  }
  if(!strcmp("GGCG", seq)) {
    return 166;
  }
  if(!strcmp("GGCT", seq)) {
    return 167;
  }
  if(!strcmp("GGGA", seq)) {
    return 168;
  }
  if(!strcmp("GGGC", seq)) {
    return 169;
  }
  if(!strcmp("GGGG", seq)) {
    return 170;
  }
  if(!strcmp("GGGT", seq)) {
    return 171;
  }
  if(!strcmp("GGTA", seq)) {
    return 172;
  }
  if(!strcmp("GGTC", seq)) {
    return 173;
  }
  if(!strcmp("GGTG", seq)) {
    return 174;
  }
  if(!strcmp("GGTT", seq)) {
    return 175;
  }
  if(!strcmp("GTAA", seq)) {
    return 176;
  }
  if(!strcmp("GTAC", seq)) {
    return 177;
  }
  if(!strcmp("GTAG", seq)) {
    return 178;
  }
  if(!strcmp("GTAT", seq)) {
    return 179;
  }
  if(!strcmp("GTCA", seq)) {
    return 180;
  }
  if(!strcmp("GTCC", seq)) {
    return 181;
  }
  if(!strcmp("GTCG", seq)) {
    return 182;
  }
  if(!strcmp("GTCT", seq)) {
    return 183;
  }
  if(!strcmp("GTGA", seq)) {
    return 184;
  }
  if(!strcmp("GTGC", seq)) {
    return 185;
  }
  if(!strcmp("GTGG", seq)) {
    return 186;
  }
  if(!strcmp("GTGT", seq)) {
    return 187;
  }
  if(!strcmp("GTTA", seq)) {
    return 188;
  }
  if(!strcmp("GTTC", seq)) {
    return 189;
  }
  if(!strcmp("GTTG", seq)) {
    return 190;
  }
  if(!strcmp("GTTT", seq)) {
    return 191;
  }
  if(!strcmp("TAAA", seq)) {
    return 192;
  }
  if(!strcmp("TAAC", seq)) {
    return 193;
  }
  if(!strcmp("TAAG", seq)) {
    return 194;
  }
  if(!strcmp("TAAT", seq)) {
    return 195;
  }
  if(!strcmp("TACA", seq)) {
    return 196;
  }
  if(!strcmp("TACC", seq)) {
    return 197;
  }
  if(!strcmp("TACG", seq)) {
    return 198;
  }
  if(!strcmp("TACT", seq)) {
    return 199;
  }
  if(!strcmp("TAGA", seq)) {
    return 200;
  }
  if(!strcmp("TAGC", seq)) {
    return 201;
  }
  if(!strcmp("TAGG", seq)) {
    return 202;
  }
  if(!strcmp("TAGT", seq)) {
    return 203;
  }
  if(!strcmp("TATA", seq)) {
    return 204;
  }
  if(!strcmp("TATC", seq)) {
    return 205;
  }
  if(!strcmp("TATG", seq)) {
    return 206;
  }
  if(!strcmp("TATT", seq)) {
    return 207;
  }
  if(!strcmp("TCAA", seq)) {
    return 208;
  }
  if(!strcmp("TCAC", seq)) {
    return 209;
  }
  if(!strcmp("TCAG", seq)) {
    return 210;
  }
  if(!strcmp("TCAT", seq)) {
    return 211;
  }
  if(!strcmp("TCCA", seq)) {
    return 212;
  }
  if(!strcmp("TCCC", seq)) {
    return 213;
  }
  if(!strcmp("TCCG", seq)) {
    return 214;
  }
  if(!strcmp("TCCT", seq)) {
    return 215;
  }
  if(!strcmp("TCGA", seq)) {
    return 216;
  }
  if(!strcmp("TCGC", seq)) {
    return 217;
  }
  if(!strcmp("TCGG", seq)) {
    return 218;
  }
  if(!strcmp("TCGT", seq)) {
    return 219;
  }
  if(!strcmp("TCTA", seq)) {
    return 220;
  }
  if(!strcmp("TCTC", seq)) {
    return 221;
  }
  if(!strcmp("TCTG", seq)) {
    return 222;
  }
  if(!strcmp("TCTT", seq)) {
    return 223;
  }
  if(!strcmp("TGAA", seq)) {
    return 224;
  }
  if(!strcmp("TGAC", seq)) {
    return 225;
  }
  if(!strcmp("TGAG", seq)) {
    return 226;
  }
  if(!strcmp("TGAT", seq)) {
    return 227;
  }
  if(!strcmp("TGCA", seq)) {
    return 228;
  }
  if(!strcmp("TGCC", seq)) {
    return 229;
  }
  if(!strcmp("TGCG", seq)) {
    return 230;
  }
  if(!strcmp("TGCT", seq)) {
    return 231;
  }
  if(!strcmp("TGGA", seq)) {
    return 232;
  }
  if(!strcmp("TGGC", seq)) {
    return 233;
  }
  if(!strcmp("TGGG", seq)) {
    return 234;
  }
  if(!strcmp("TGGT", seq)) {
    return 235;
  }
  if(!strcmp("TGTA", seq)) {
    return 236;
  }
  if(!strcmp("TGTC", seq)) {
    return 237;
  }
  if(!strcmp("TGTG", seq)) {
    return 238;
  }
  if(!strcmp("TGTT", seq)) {
    return 239;
  }
  if(!strcmp("TTAA", seq)) {
    return 240;
  }
  if(!strcmp("TTAC", seq)) {
    return 241;
  }
  if(!strcmp("TTAG", seq)) {
    return 242;
  }
  if(!strcmp("TTAT", seq)) {
    return 243;
  }
  if(!strcmp("TTCA", seq)) {
    return 244;
  }
  if(!strcmp("TTCC", seq)) {
    return 245;
  }
  if(!strcmp("TTCG", seq)) {
    return 246;
  }
  if(!strcmp("TTCT", seq)) {
    return 247;
  }
  if(!strcmp("TTGA", seq)) {
    return 248;
  }
  if(!strcmp("TTGC", seq)) {
    return 249;
  }
  if(!strcmp("TTGG", seq)) {
    return 250;
  }
  if(!strcmp("TTGT", seq)) {
    return 251;
  }
  if(!strcmp("TTTA", seq)) {
    return 252;
  }
  if(!strcmp("TTTC", seq)) {
    return 253;
  }
  if(!strcmp("TTTG", seq)) {
    return 254;
  }
  if(!strcmp("TTTT", seq)) {
    return 255;
  }
}
