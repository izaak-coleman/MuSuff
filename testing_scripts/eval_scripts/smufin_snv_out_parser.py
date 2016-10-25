import sys
import csv

def parse_smufin_snv_file(snvfile):
    print "parse_smufin_snv_from: ", snvfile
   
    snv_list = []
  
    with open(snvfile) as snvs:
        for line in csv.reader(snvs, delimiter="\t"):
            if line[0] == "Mut_ID":		#skip header
                continue

            line[3] = int(line[3])
            snv_list.append(tuple(line[3:6]))	      # ICL-SMuFin requirement
    return sorted(set(snv_list), key=lambda x: x[0])  # ICL-SMuFin-unique


def main():
    smuf_list = parse_smufin_snv_file(sys.argv[1])

    for snv in smuf_list:
        print snv


if __name__ == "__main__":
    main()
    
