#ifndef VARIANTCALL_H
#define VARIANTCALL_H 

#include <string>
#include <vector>
#include <fstream>


struct snv_aln_info {
 std::vector<int> SNV_pos;
 int flag;
 int chr;
 int position;
 std::string non_mutated_cns;
 std::string mutated_cns;
};

struct single_snv {
 int flag;
 int chr;
 int position;
 char mutation_base;
 char healthy_base;
};


struct mutation_classes{
  std::vector<int> SNV_pos;
};

struct consensus_pair {
  std::string mutated;
  std::string non_mutated;
  std::string read_freq_m;
  std::string read_freq_nm;
  int mut_offset;
  int nmut_offset;
  mutation_classes mutations;
};


void identifySNVs(std::vector<snv_aln_info> &alignments);
// iterates through alignments and calls countSNVs() to identify mutations
// handles reverse complement aligment of the healthy cns
void countSNVs(snv_aln_info &alignment);



std::string generateParseString(mutation_classes &m);
// Function generates a string of the following string (mutation string)
// [SNV:a;b;c][SSV:e;f;g][LSV:x;y;z]
// this string is stored as the header for each
// aligned healthy read. This allows identification
// of the mutation indexes directly from the SAM file

void constructSNVFastqData();
// generates the fastq file of non_mutated reads to align
// to the refrence genome. 
// a fastq element has format 

void callBWA();
// Function first calls bwa aln:
// cmd ./bwa/bwa aln -t 16 hg19.fa cns_pairs.fastq > cns_pairs.sai
// Then cals to bwa samse
// ./bwa/bwa samse hg19.fa cna_pairs.sai cns_pair.fastq > cns_pairs.sam


void parseSamFile(std::vector<snv_aln_info> &alignments, std::string filename);
void printAllAlignments(std::vector<snv_aln_info> &alignments);
void printSingleAlignment(snv_aln_info &snv);
void correctReverseCompSNV(std::vector<snv_aln_info> &alignments);
static bool compareSNVLocations(const single_snv &a, const single_snv &b);
void outputSNVToUser(std::vector<snv_aln_info> &alignments, std::string report_filename);


 void printConsensusPairs();
 // print out each mutated and non mutated string
 void printMutation(char healthy, char cancer, std::ofstream &mut_file);

#endif
