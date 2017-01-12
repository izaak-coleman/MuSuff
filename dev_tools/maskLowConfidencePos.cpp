#include <iostream>
#include <vector>


using namespace std;

static const double ALLELIC_FREQ_OF_ERROR = 0.1;


  struct consensus_pair {
    std::string mutated;
    std::string non_mutated;
    std::string read_freq_m;
    std::string read_freq_nm;
    unsigned int pair_id;
    int mut_offset;
    int nmut_offset;
    int left_ohang;
    int right_ohang;
  };

void maskLowConfidencePositions(consensus_pair &pair,
    vector< vector<int> > &healthy_base_freq,
    vector< vector<int> > &tumour_base_freq);

int main() {
  consensus_pair pair;
  pair.mutated =         "ATCCGCG";
  pair.non_mutated = "AACCTTCCGCAGCCTG";  // A vs low confidence T

  pair.left_ohang = 4;
  pair.right_ohang = 5;

  cout << "Consensus pair vals: " << endl;
  cout << "T:     " << pair.mutated << endl;
  cout << "H: " << pair.non_mutated << endl;


  vector< vector<int> > tumour_freq = {
    { -1,-1,-1,-1,4,1,1,1,1,1,0,-1,-1,-1 }, // A
    { -1,-1,-1,-1,0,4,0,0,0,0,0,-1,-1,-1 }, // T
    { -1,-1,-1,-1,0,0,4,4,0,4,0,-1,-1,-1 }, // C
    { -1,-1,-1,-1,0,0,0,0,4,0,4,-1,-1,-1 }  // G
  };

  vector< vector<int> > healthy_freq = {
    { -1,-1,4,4,0,0,1,0,0,0,0,0,3,0,0,0,0,0,-1 }, // A
    { -1,-1,0,0,0,0,9,4,0,0,0,0,1,0,0,0,4,0,-1 }, // T
    { -1,-1,0,0,4,4,1,0,4,4,0,4,0,0,4,4,0,0,-1 }, // C
    { -1,-1,0,0,0,0,1,0,0,0,4,0,0,4,0,0,0,4,-1 }  // G
  };

  maskLowConfidencePositions(pair, healthy_freq, tumour_freq);

  cout << "Consensus pair vals: " << endl;
  cout << "T:     " << pair.mutated << endl;
  cout << "H: " << pair.non_mutated << endl;


}

void maskLowConfidencePositions(consensus_pair &pair,
                                vector< vector<int> > &healthy_base_freq,
                                vector< vector<int> > &tumour_base_freq) {
  unsigned int start_h= 0, start_t= 0;

  for(int i=0; i < healthy_base_freq[0].size(); i++) {
    if(healthy_base_freq[0][i] != -1) {
      start_h = i;
      break;
    }
  }
  for(int i=0; i < tumour_base_freq[0].size(); i++) {
    if(tumour_base_freq[0][i] != -1) {
      start_t = i;
      break;
    }
  }

  // mask based on tumour cns
  for(int pos = start_t; pos < pair.mutated.size() + start_t; pos++) {
    int n_tumour_bases_above_err_freq = 0;

    if(tumour_base_freq[0][pos] == -1) {
      break;
    }

    // get total reads
    double total_reads = 0;
    for(int base=0; base < 4; base++) {
      total_reads += tumour_base_freq[base][pos];
    }

    // calc number of bases over the error frequency
    for(int base=0; base < 4; base++) {
      if(tumour_base_freq[base][pos] / total_reads > ALLELIC_FREQ_OF_ERROR) {
        n_tumour_bases_above_err_freq++;
      }
    }

    // if the number of bases with a high allelic frequency is above
    // one, then the position is of low condifence, so mask
    if (n_tumour_bases_above_err_freq > 1) {
      pair.mutated[pos - start_t] = pair.non_mutated[pos - start_t + pair.left_ohang];
    }
  }


  // mask based on healthy cns
  for(int pos = start_h + pair.left_ohang; pos < pair.mutated.size() +
      start_h + pair.left_ohang; pos++) {

    int n_healthy_bases_above_err_freq = 0;
    if(healthy_base_freq[0][pos] == -1) {
      break;
    }

    // get total reads
    double total_reads = 0;
    for(int base = 0; base < 4; base++) {
      total_reads += healthy_base_freq[base][pos];
    }

    // cals number of bases over the error frequency
    for(int base=0; base < 4; base++) {
      if(healthy_base_freq[base][pos] / total_reads > ALLELIC_FREQ_OF_ERROR) {
        n_healthy_bases_above_err_freq++;
      }
    }

    // if n bases with high allelic freq. is above one, then 
    // position is low confidence so mask
    if(n_healthy_bases_above_err_freq > 1) {
      pair.mutated[pos - pair.left_ohang - start_h] = pair.non_mutated[pos -
        start_h];
    }
  }
}
