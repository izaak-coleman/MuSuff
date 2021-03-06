import sys
import report_sorter

def print_list(l):
    for i in l:
        print i


def ambiguous_report(snv_matches, location_matches):

    count = 0

    ambiguous_locs = []

    for snv in snv_matches:
        for loc_match in location_matches:
            if report_sorter.location_match(snv, loc_match):
                count = count + 1
                ambiguous_locs.append(loc_match)
    return count, ambiguous_locs


def near_hit(complete_misses, hit_list, error):
    down_stream = []
    up_stream = []

    for rep_snv in complete_misses:
        for gen_hit in hit_list:
            dist = rep_snv[0] - gen_hit[0] # calc distance between two
            if abs(dist) == error:
                if rep_snv[1] == gen_hit[1] and rep_snv[2] == gen_hit[2]: # bases equiv
                    #determine whether it was an up or down stream error
                    if dist < 0:
                        down_stream.append([rep_snv, gen_hit])
                    else:
                        up_stream.append([rep_snv, gen_hit])

    return down_stream, up_stream

def reportedSNVAtDistance(reportedSNVs, genuineList, distance):
    downStream, upStream = [], []
    for reportedSNV in reportedSNVs:
        for genuineHit in genuineList:
          if abs(reportedSNV[0] - genuineHit[0]) == distance:
              if reportedSNV[0] - genuineHit[0] < 0:
                downStream.append(reportedSNV)
              else:
                upStream.append(reportedSNV)

    return [downStream, upStream]

def main():

    chr22_len = 51304566

    if(len(sys.argv) != 6):
        print "usage: <gen_snv1> <gen_snv2> <snp_1> <snp_2> <output_file>"

    # statistics buffer
    results = []




    snv_matches, snp_matches, non_matches, snv_location_hits, snp_location_hits, complete_misses, snv_list, snp_list, reported_list = report_sorter.sort_snvs()

    location_matches = snv_location_hits + snp_location_hits
    # count the number of ambigous snv matches
    ambiguous_snv_matches, ambig_snv_list = ambiguous_report(snv_matches, location_matches)
    ambiguous_snp_matches, ambig_snp_list = ambiguous_report(snp_matches, location_matches)


    #calculate any complete misses that are a distance of X from some
    # genuine mutation, and share the same mutation structure as the 
    # genuine mutation
    # SNV comparison...
    dwn_err, up_err = near_hit(complete_misses, snv_list, 1)
    dwn_err2, up_err2 = near_hit(complete_misses, snv_list, 2)
    dwn_err3, up_err3 = near_hit(complete_misses, snv_list, 3)
    dwn_err4, up_err4 = near_hit(complete_misses, snv_list, 4)
    dwn_err5, up_err5 = near_hit(complete_misses, snv_list, 5)

    # SNP comparison...
    locAndStructureHit = []
    for i in range (1, 6):
      locAndStructureHit += near_hit(complete_misses, snp_list, i)


    # generate list of complete_misses that are a distance of X form
    # some mutation (they do not have to share the same mutation structure)
    # SNP comparison...
    completeMissToSNP = []
    for i in range(1, 6):
      completeMissToSNP += reportedSNVAtDistance(complete_misses, snp_list, i)

    # SNV comparison
    completeMissToSNV = []
    for i in range(1, 6):
      completeMissToSNV += reportedSNVAtDistance(complete_misses, snv_list, i)

    # record results in results_table
    results.append("exact_snv_matches, " + str(len(snv_matches)))
    results.append("exact_snp_matches, " + str(len(snp_matches)))
    results.append("non_exact_matches, " + str(len(non_matches)))
    results.append("complete_misses, " + str(len(complete_misses)))
    results.append("ambiguous_snv_matches, " + str(ambiguous_snv_matches))
    results.append("ambiguous_snp_matches, " + str(ambiguous_snp_matches))
    results.append("location_hits_with_snv, " + str(len(snv_location_hits)))
    results.append("location_hits_with_snp, " + str(len(snp_location_hits)))
    results.append("minus_one_snv_match, " + str(len(dwn_err)))
    results.append("plus_one_snv_match, " + str(len(up_err)))

    results.append("loc and struct hit: minus_2_snv_match, " + str(len(dwn_err2)))
    results.append("loc and struct hit: plus_2_snv_match, " + str(len(up_err2)))
    results.append("loc and struct hit: minus_3_snv_match, " + str(len(dwn_err3)))
    results.append("loc and struct hit: plus_3_snv_match, " + str(len(up_err3)))
    results.append("loc and struct hit: minus_4_snv_match, " + str(len(dwn_err4)))
    results.append("loc and struct hit: plus_4_snv_match, " + str(len(up_err4)))
    results.append("loc and struct hit: minus_5_snv_match, " + str(len(dwn_err5)))
    results.append("loc and struct hit: plus_5_snv_match, " + str(len(up_err5)))

    results.append("loc and struct hit: minus_1_snp_match, " + str(len(locAndStructureHit[0])))
    results.append("loc and struct hit: plus_1_snp_match, " + str(len(locAndStructureHit[1])))
    results.append("loc and struct hit: minus_2_snp_match, " + str(len(locAndStructureHit[2])))
    results.append("loc and struct hit: plus_2_snp_match, " + str(len(locAndStructureHit[3])))
    results.append("loc and struct hit: minus_3_snp_match, " + str(len(locAndStructureHit[4])))
    results.append("loc and struct hit: plus_3_snp_match, " + str(len(locAndStructureHit[5])))
    results.append("loc and struct hit: minus_4_snp_match, " + str(len(locAndStructureHit[6])))
    results.append("loc and struct hit: plus_4_snp_match, " + str(len(locAndStructureHit[7])))
    results.append("loc and struct hit: minus_5_snp_match, " + str(len(locAndStructureHit[8])))
    results.append("loc and struct hit: plus_5_snp_match, " + str(len(locAndStructureHit[9])))

    results.append("loc hit: minus_1_snp_match, " + str(len(completeMissToSNP[0])))
    results.append("loc hit: plus_1_snp_match, " + str(len(completeMissToSNP[1])))
    results.append("loc hit: minus_2_snp_match, " + str(len(completeMissToSNP[2])))
    results.append("loc hit: plus_2_snp_match, " + str(len(completeMissToSNP[3])))
    results.append("loc hit: minus_3_snp_match, " + str(len(completeMissToSNP[4])))
    results.append("loc hit: plus_3_snp_match, " + str(len(completeMissToSNP[5])))
    results.append("loc hit: minus_4_snp_match, " + str(len(completeMissToSNP[6])))
    results.append("loc hit: plus_4_snp_match, " + str(len(completeMissToSNP[7])))
    results.append("loc hit: minus_5_snp_match, " + str(len(completeMissToSNP[8])))
    results.append("loc hit: plus_5_snp_match, " + str(len(completeMissToSNP[9])))
            
            
    results.append("loc hit: minus_1_snv_match, " + str(len(completeMissToSNV[0])))
    results.append("loc hit: plus_1_snv_match, " + str(len(completeMissToSNV[1])))
    results.append("loc hit: minus_2_snv_match, " + str(len(completeMissToSNV[2])))
    results.append("loc hit: plus_2_snv_match, " + str(len(completeMissToSNV[3])))
    results.append("loc hit: minus_3_snv_match, " + str(len(completeMissToSNV[4])))
    results.append("loc hit: plus_3_snv_match, " + str(len(completeMissToSNV[5])))
    results.append("loc hit: minus_4_snv_match, " + str(len(completeMissToSNV[6])))
    results.append("loc hit: plus_4_snv_match, " + str(len(completeMissToSNV[7])))
    results.append("loc hit: minus_5_snv_match, " + str(len(completeMissToSNV[8])))
    results.append("loc hit: plus_5_snv_match, " + str(len(completeMissToSNV[9])))


    # calculate sensitivity and selectivity
    sensitivity = float(len(snv_matches)) / float(len(snv_list))

    selectivity = float((chr22_len - len(location_matches) - len(complete_misses) - len(snp_matches))) / float(chr22_len)
    selectivity_snp_filter = float((chr22_len - len(location_matches) - len(complete_misses))) / float(chr22_len)
   
    sel = float((chr22_len - (len(snv_matches) + len(snp_matches) + len(complete_misses))) / 
                (float(chr22_len) - (len(snv_matches) + len(snp_matches))))
    results.append("sensitivity, " + str(sensitivity))
    results.append("selectivity, " + str(selectivity))
    results.append("selectivity_snp_filter, " + str(selectivity_snp_filter))
    results.append("sel_actual, " + str(sel))


    print "snv matches..."
    print_list(snv_matches)
    print "snp matches..."
    print_list(snp_matches)
    print "non_exact_matches..."
    print_list(non_matches)
    print "snv_location_hits..."
    print_list(snv_location_hits)
    print "snp_location_hits..."
    print_list(snp_location_hits)
    print "complete_misses..."
    print_list(complete_misses)
    print "ambiguous snv matches..."
    print_list(ambig_snv_list)
    print "ambiguous snp matches..."
    print_list(ambig_snp_list)

    print "minus_one_snv_match..."
    print_list(dwn_err)
    print "plus_one_snv_match..."
    print_list(up_err)

    print "minus_2_snv_match..."
    print_list(dwn_err2)
    print "plus_2_snv_match..."
    print_list(up_err2)


    print "minus_3_snv_match..."
    print_list(dwn_err3)
    print "plus_3_snv_match..."
    print_list(up_err3)

    print "minus_4_snv_match..."
    print_list(dwn_err4)
    print "plus_4_snv_match..."
    print_list(up_err4)

    print "minus_5_snv_match..."
    print_list(dwn_err5)
    print "plus_5_snv_match..."
    print_list(up_err5)
    print "plus_1 close proximity SNV hit"
    print_list(completeMissToSNP[1])

    # print results to user
    print_list(results)


if __name__ == "__main__":
    main()





