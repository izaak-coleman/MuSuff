import sys
import AnalyseReadsCoveringFN as fn

def main():
  if len(sys.argv) != 2:
    print "Usage: <exe> <fn_data>"
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

main()
