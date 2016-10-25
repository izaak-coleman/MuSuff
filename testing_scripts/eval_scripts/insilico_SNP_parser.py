import sys
import csv # tab file parser

def parse_snp_from_file(file_one, file_two):
    """Function parses file returning 
    tuple of (location, healthy_base, mutation_base)"""

    print "parse_snp_from_file"
    print "parsing over ", file_one, " ", file_two

    snp_list = []   # store parsed snps

    # open file as tab delimited
    with open(file_one) as snps_one:
        for line in csv.reader(snps_one, delimiter="\t"):
            if line[1] == "SNP":
                line[4] = int(line[4])
                snp_list.append(tuple(line[4:]))

    with open(file_one) as snps_two:
        for line in csv.reader(snps_two, delimiter="\t"):
            if line[1] == "SNP":
                line[4] = int(line[4])
                snp_list.append(tuple(line[4:]))


    
    return sorted(set(snp_list), key=lambda x: x[0])

def main():
    snps = parse_snp_from_file(sys.argv[1], sys.argv[2])
    for snp in snps:
        print snp


if __name__ == "__main__":
    main()
