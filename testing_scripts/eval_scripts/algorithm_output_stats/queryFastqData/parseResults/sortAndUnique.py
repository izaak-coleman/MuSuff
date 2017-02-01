import sys
from operator import itemgetter, attrgetter, methodcaller

def f7(seq):
  seen = set()
  seen_add = seen.add
  return [x for x in seq if not (x in seen or seen_add(x))]


def main():
  tuples = []
  with open(sys.argv[1], 'r') as fsock:
    tuples = [line.strip() for line in fsock]
    tuples = [(int(line[1:line.find(',')]), line[-2:-1]) for line in tuples] 
    tuples = list(set(tuples))
 
  tuples = sorted(tuples, key=itemgetter(1,0))
  for t in tuples:
    print("(%d,%s)" % (t[0], t[1]))

main()


