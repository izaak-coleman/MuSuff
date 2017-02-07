#!/usr/bin/python
import sys

def main():
  if len(sys.argv) != 2:
    print "Usage: <exe> <tuple_file>"
    sys.exit()

  tuples = []
  with open(sys.argv[1], 'r') as fsock:
    tuples = [line.strip() for line in fsock]
    tuples = [(int(line[1:line.find(',')]), line[-2:-1]) for line in tuples] 
    tuples = list(set(tuples))
 
  tuples = sorted(tuples, key=lambda x: (x[1], x[0]))
  for t in tuples:
    print("(%d,%s)" % (t[0], t[1]))

main()
