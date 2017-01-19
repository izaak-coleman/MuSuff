import sys

import insilico_SNV_parser as gen_snv_psr 
import insilico_SNP_parser as snp_psr
import smufin_snv_out_parser as smuf_psr

def print_list(l):
    for i in l:
        print i

def unique_sorted(l, field):
   return sorted(set(l), key=lambda x: x[field])

def add_inclusion_check(l):
    inc_list = []
    for elem in l:
        inc_list.append([False, elem])
    return inc_list


def exact_mutation_match(a, b):
    for i in range(len(a)):	#if any field not identical
        if a[i] != b[i]:
            return False
    return True

def location_match(a, b):
    return a[0] == b[0]


def sort_location_hits(non_exact_matches, genuine_snvs, genuine_snps):
    snv_location_hits = []
    snp_location_hits = []
    complete_misses = []

    for non_exact_match in non_exact_matches:
        snv_location_hit = False
        snp_location_hit = False

        for gen_snv in genuine_snvs:
            if location_match(gen_snv, non_exact_match):
                snv_location_hits.append(non_exact_match)
                snv_location_hit = True
                break
        if snv_location_hit:
            continue


        for gen_snp in genuine_snps:
            if location_match(gen_snp, non_exact_match):
                snp_location_hits.append(non_exact_match)
                snp_location_hit = True
                break
        if snp_location_hit:
            continue


        complete_misses.append(non_exact_match)

    return snv_location_hits, snp_location_hits, complete_misses


def sort_reported_hits(reported_snvs, genuine_snvs, genuine_snps):
    exact_snv_hits = []
    exact_snp_hits = []
    non_exact_matches = []

    for reported_snv in reported_snvs:
        snv_hit = False
        snp_hit = False

        for gen_snv in genuine_snvs:
            if exact_mutation_match(gen_snv, reported_snv):
                exact_snv_hits.append(reported_snv)
                snv_hit = True
                break

        if snv_hit:		# move to next reported snv
           continue

        for gen_snp in genuine_snps:
            if exact_mutation_match(gen_snp, reported_snv):
                exact_snp_hits.append(reported_snv)
                snp_hit = True
                break

        if snp_hit:
           continue


        # if reached here in inner loop, then no hits
        non_exact_matches.append(reported_snv)


    return exact_snv_hits, exact_snp_hits, non_exact_matches


def sort_snvs():

    print "sort_snvs from report_sorter"


    # parse data
    snv_list = gen_snv_psr.parse_snv_from_file(sys.argv[1], sys.argv[2])
    snp_list = snp_psr.parse_snp_from_file(sys.argv[3], sys.argv[4])
    reported_list = smuf_psr.parse_smufin_snv_file(sys.argv[5])
    print "Unique list len: ", len(reported_list)

    # identify exact snp, snv hits in smuf_list
    snv_matches, snp_matches, non_matches = sort_reported_hits(reported_list, snv_list, snp_list)


    # identify reported snvs that match a snv or snp location but report the wrong base
    snv_location_hits, snp_location_hits, complete_misses = sort_location_hits(non_matches, snv_list, snp_list)

    return snv_matches, snp_matches, non_matches, snv_location_hits, snp_location_hits, complete_misses, snv_list, snp_list, reported_list

if __name__ == "__main__":
    sort_snvs()
