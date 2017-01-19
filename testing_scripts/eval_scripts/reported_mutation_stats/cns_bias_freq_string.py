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


def extract_freq_string_and_SNVs(cns_pairs):
    strip_list = []
    for pair in cns_pairs:
        # extract raw strings
        mut_freq = pair[5]
        snvs = pair[9]

        #split them
        mut_freq = mut_freq.split(", ")
        snvs = snvs.split(", ")
        mut_freq.pop()
        snvs.pop()

        # integerize them
        mut_freq = [int(x) for x in mut_freq]
        snvs     = [int(x) for x in snvs]


        strip_list.append([mut_freq, snvs])
    return strip_list


def main():

    cns_pairs = parse_consensus_pairs(sys.argv[1])
    stripped_pairs = extract_freq_string_and_SNVs(cns_pairs)
    for pair in stripped_pairs:
        print pair

if __name__ == "__main__":
    main()
   
