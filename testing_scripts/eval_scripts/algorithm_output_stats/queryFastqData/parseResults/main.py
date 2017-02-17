# -*- coding: utf-8 -*-
#!/usr/bin/env python

import sys
import AnalyseReadsCoveringFN as fn

def main():
  if len(sys.argv) != 3:
    print "Usage: <exe> <reads_covering_false_negs> <tuple_point_file>"
    sys.exit()

  analysis = fn.AnalyseReadsCoveringFN(sys.argv[1])
  print "Less than 4"
  analysis.countThan(lambda x:x < 4, "Count < 4: ")

  print "\n\n"
  print "Greater than 80"
  analysis.countThan(lambda x:x > 80, "Count > 80: ")
  print "\n\n"
  print "Greater than 70"
  analysis.countThan(lambda x: x > 79, "Count > 79: ")
  print "\n\n"
  print "\n\n"
  
  analysis.quantifyLostReads(sys.argv[2])

main()
