import sys
import csv # tab file parser

def parse_snv_from_file(file_one, file_two):
    """Function parses file returning 
    tuple of (location, healthy_base, mutation_base)"""

    print "parse_snv_from_file"
    print "parsing over ", file_one, " ", file_two

    snv_list = []   # store parsed snvs

    # open file as tab delimited
    with open(file_one) as snvs_one:
        for line in csv.reader(snvs_one, delimiter="\t"):
            line[2] = int(line[2])   #convert locations to ints for sorting
            snv_list.append(tuple(line[2:]))

    with open(file_two) as snvs_two:
        for line in csv.reader(snvs_two, delimiter="\t"):
            line[2] = int(line[2]) 
            snv_list.append(tuple(line[2:]))


    snv_list.sort(key=lambda x: x[0])
    return snv_list

def main():
    snvs = parse_mutations_from_file(sys.argv[1], sys.argv[2])
    for snv in snvs:
        print snv


if __name__ == "__main__":
    main()
