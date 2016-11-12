#ifndef SAMENTRY_H
#define SAMENTRY_H

class SamEntry {
private:

  // sam mandatory fields
  int flag;
  int chr_num;
  int aln_pos;
  int mapq;
  int ins_size;
  std::string mutation_string;
  std::string cigar;
  std::string healthy_cns;
  std::string qstring;

  // augmented sam - bwa optional fields
  std::string nm;
  std::string md;
  std::string as;
  std::string bc;
  std::string x0;
  std::string x1;
  std::string xn;
  std::string xm;
  std::string xo;
  std::string xg;
  std::string xt;
  std::string xa;
  std::string xs;
  std::string xf;
  std::string xe;


public:
  int getFlag();
  int getChrNum();
  int getAlnPos();
  int getMapq();
  int getInsSize();

  std::string getMutationString();
  std::string getCigarString();
  std::string getHealthyCns();
  std::string getQString();

  std::string getNM();
  std::string getMD();
  std::string getAS();
  std::string getBC();
  std::string getX0();
  std::string getX1();
  std::string getXN();
  std::string getXM();
  std::string getXO();
  std::string getXG();
  std::string getXT();
  std::string getXA();
  std::string getXS();
  std::string getXF();
  std::string getXE();
  
};



#endif
