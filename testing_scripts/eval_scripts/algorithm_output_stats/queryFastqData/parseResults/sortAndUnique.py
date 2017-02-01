import sys
from operator import itemgetter, attrgetter, methodcaller

def f7(seq):
  seen = set()
  seen_add = seen.add
  return [x for x in seq if not (x in seen or seen_add(x))]


def main():
  hTuples = []
  cTuples = []
  with open(sys.argv[1], 'r') as fsock:
    tuples = [line.strip() for line in fsock]
    hTuples = [(int(line[1:line.find(',')]), line[-2:-1]) for line in tuples if line[-2:-1] == "H"]
    cTuples = [(int(line[1:line.find(',')]), line[-2:-1]) for line in tuples if line[-2:-1] == "T"]


 
  hTup = list(set(hTuples))
  cTup = list(set(cTuples))
  for t in sorted(hTup, key=itemgetter(0)):
    print("(%d,%s)" % (t[0], t[1]))
  for t in sorted(cTup, key=itemgetter(0)):
    print("(%d,%s)" % (t[0], t[1]))

main()


