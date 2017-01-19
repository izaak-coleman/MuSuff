import sys

def parse_consensus_pairs(filename):
    consensus_pairs = [line.rstrip('\n') for line in open(filename)]

    pairs = []

    pair_start = 0
    for i in range(len(consensus_pairs)):
        if("Mutation" in consensus_pairs[i]):
            start_block = i
            while("SNV" not in consensus_pairs[i]):	# advance to SNP line
                i = i + 1

	    i = i + 1 # move to SNV mutation line

            if i == len(consensus_pairs):
                break
            if len(consensus_pairs[i]) == 0:		# no mutation
                i = i + 1 # move to line before M
                continue

            else:
                i = i + 1
                pairs.append(consensus_pairs[start_block:i])
            

    return pairs

def extract_frequencies(cns_pairs):
    freq_list = []
    for pair in cns_pairs:
        cns_length = len(pair[1])
        SNV_list = pair[5].split(", ")
        SNV_list.pop()
        SNV_list = [(1+int(x)) for x in SNV_list]

        sym_list = []
        for snv in SNV_list:
            if(snv > (cns_length/2)):			# fold SNVs to same end
                sym_snv = (cns_length - snv)
            else:
                sym_snv = snv

            if sym_snv < 40:			# unbiased length
                freq_list.append(sym_snv)
                

    return freq_list

    
def count_frequencies(occurences):
    frequencies = [0] * 40

    for occ in occurences:
        frequencies[occ] = frequencies[occ] + 1

    return frequencies

def print_frequencies(freqs):
    pos = 0
    print "cns position, occurences"
    for freq in freqs:
        print pos, ", ", freq
        pos = pos + 1

def main():
    cns_pairs = parse_consensus_pairs(sys.argv[1])
    occurences = extract_frequencies(cns_pairs)

    frequencies = count_frequencies(occurences)

    print frequencies
    print_frequencies(frequencies)




if __name__ == "__main__":
    main()
