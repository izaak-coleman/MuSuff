import sys 

def mem_stripper(bwa_mem_file):

    skimmed_output = []
    f = open(bwa_mem_file)
    lines = f.readlines()
    for line in lines:
        fields = line.split("\t")
        h = fields[4]
        t = fields[5]
        if h != "A" and h != "T" and h != "C" and h != "G":
            continue
        if t != "A" and t != "T" and t != "C" and t != "G":
            continue 

        skimmed_output.append(line.rstrip())
    return skimmed_output


def main():

    skimmed_output = mem_stripper(sys.argv[1])
    for line in skimmed_output:
        print line

   
    

main()
