import sys
import AnalyseReadsCoveringFN as fn

def main():
  if len(sys.argv) != 2:
    print "Usage: <exe> <fn_data>"
    sys.exit()

  analysis = fn.AnalyseReadsCoveringFN(sys.argv[1])
  analysis.printKeyVals()


main()
