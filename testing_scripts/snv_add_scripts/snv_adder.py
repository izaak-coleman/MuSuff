import sys

start_of_chr22 = 15888000

def add_snvs(fa_name, n_mutations, dist, after_n):
    fa_string = ""
    f = open(fa_name)
    lines = f.readlines()

    header =  lines.pop(0) # remove header
    for line in lines:
        fa_string += line.rstrip()


    line_len = len(lines[0])-1
  
    #convert fa_strin to muatable list
    fa_string = list(fa_string)

    # at every i * dist interval add a mutation
    print "location, healthy, tumour"
    for i in range(n_mutations):
        mutation = ""
        if (i % 2) == 0:
            mutation = transversion_mutation(fa_string[after_n + i*dist])
        else:
            mutation = transition_mutation(fa_string[after_n + i*dist])
 
        print (i*dist) + start_of_chr22,",", fa_string[after_n + i*dist],",", mutation
        fa_string[after_n + i*dist] = mutation
 
 
    print header.rstrip()
   
    # reslice the singe string into the length of lines, and add crg return

    # reconvert fa_string to list
    fa_string = ''.join(fa_string)
    for i in range(len(lines)):
        print fa_string[i*(line_len) : (i+1)*line_len]


def transversion_mutation(base):
    return {
        "A":"T",
        "T":"A",
        "C":"G",
        "G":"C",
        "N":"N"
    }[base]

def transition_mutation(base):
    return {
        "A":"G",
        "G":"A",
        "T":"C",
        "C":"T",
        "N":"N"
    }[base]



def main():
    if (sys.argv) != 4:
        print "useage: <exe> <fasta> <n mutations> <dist> <number_after_n>"

    add_snvs(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))



if __name__ == "__main__":
    main()
