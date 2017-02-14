import sys
import fileinput
import itertools

class FalseNegativeEntry(dict):
  HEADER_IDX = 0     # ONLY CHANGE IF OUTPUT FORMAT CHANGES
  SNIPPET_IDX = 1
  COORD_IDX = 2
  READ_IDX = 4

  def __init__(self, entryData):
    self.loadEntryData(entryData)

  def loadEntryData(self, entryData):
    entryData = [line.strip() for line in entryData]
    header, hBase, cBase, tissueType = self.extractHeaderInfo(entryData[self.HEADER_IDX])
    self["header"] = header
    self["hBase"] = hBase
    self["cBase"] = cBase
    self["snippet"] = entryData[self.SNIPPET_IDX]
    self["coordinate"] = \
    int(entryData[self.COORD_IDX][entryData[self.COORD_IDX].find(":") + 2:])

    reads = [tuple(read.split(" :: ")) for read in entryData[self.READ_IDX:]]
    reads = [(read, int(idx), nm, fr, int(offset), int(rel)) 
    for (read, idx, nm, fr, offset, rel) in reads]

    if tissueType == "H":
      self["hReads"] = reads
    else:
      self["cReads"] = reads

  def extractHeaderInfo(self, line):
    header = line[line.find(">") : line.find("-") -1]
    mutationPair = line[line.find("-") + 2 : line.find("-") + 5]
    tissueType = line[len(line) - 1]
    return header, mutationPair[0], mutationPair[2], tissueType

  def updateInfo(self, entryData):
    entryData = [line.strip() for line in entryData]
    _, _, _, tissueType = self.extractHeaderInfo(entryData[self.HEADER_IDX])
    reads = [tuple(read.split(" :: ")) for read in entryData[self.READ_IDX:]]
    reads = [(read, int(idx), nm, fr, int(offset), int(rel)) 
    for (read, idx, nm, fr, offset, rel) in reads]
    if tissueType == "H":
      self["hReads"] = reads
    else:
      self["cReads"] = reads


class AnalyseReadsCoveringFN:
  COORD_IDX = 2
  def __init__(self, fname):
    self.falseNegData = {}
    self.parseReadsFile(fname)


  def parseReadsFile(self, fname):
    fileBuf = [line for line in open(fname, "r")]
    groupedData = [list(x) for a, x in itertools.groupby(fileBuf, lambda x:
        x.startswith("Header"))]
    entries = [groupedData[i]+groupedData[i+1] for i in range(0,
        len(groupedData), 2)]

    for entry in entries:
      coordinate = int(entry[self.COORD_IDX].strip()[entry[self.COORD_IDX].find(":") + 2:])
      if coordinate not in self.falseNegData.keys():
        self.falseNegData[coordinate] = FalseNegativeEntry(entry)
      else:
        self.falseNegData[coordinate].updateInfo(entry)


  def printEntry(self):
    print self.falseNegData[16075001]


  def countThan(self, fun, resultMessage):
    count = 0
    for (k, v) in self.falseNegData.items():
      totalReads = len(v["hReads"]) + len(v["cReads"])
      if fun(totalReads):
        print "Coordi : ", str(k)
        print "H count: ", len(v["hReads"])
        print "T count: ", len(v["cReads"])
        count = count + 1
    print resultMessage, count


  def computeStats(self, reads, tcond, rcond):
    results = {}

    def percentage(remain, total):  # quick perc cal
      try:
        return (remain / float(total)) * 100
      except ZeroDivisionError:
        return 0 

    results["Total"] = len([idx for (_, idx, found) in reads if tcond(found)])
    results["Rem"]  = len([idx for (_, idx, found) in reads 
        if rcond(idx) and tcond(found)])
    results["Perc"] = percentage(results["Rem"], results["Total"])
    return results


  def quantifyLostReads(self, filename):
    statFile = filename + ".stat"
    readFile = filename + ".stat.reads"
    remainingReads = []
    results = []
    with open(filename, "r") as fileHandle:
      remainingReads = [(int(line[1:line.find(",")]), line[-3:-2]) for line in fileHandle]
      remainingReads = set(remainingReads)
    
    # compute stats
    for k, v in self.falseNegData.items():
       result = {}
       hSet = lambda x: (x, "H") in remainingReads
       cSet = lambda x: (x, "T") in remainingReads
       nSet = lambda x: x == "N"
       mSet = lambda x: x == "M"
       v["HealthyTot"] = self.computeStats(v["hReads"], lambda x: True, hSet)
       v["HealthyM"]   = self.computeStats(v["hReads"], mSet, hSet) 
       v["HealthyN"]   = self.computeStats(v["hReads"], nSet, hSet)
       v["CancerTot"]  = self.computeStats(v["cReads"], lambda x: True, cSet)
       v["CancerN"]    = self.computeStats(v["cReads"], nSet, cSet)
       v["CancerM"]    = self.computeStats(v["cReads"], mSet, cSet) 


    # print stats
    statFileHandle = open(statFile, "w")
    for k, v in self.falseNegData.items():
      self.printStats(statFileHandle, v)

    # print stats and reads
    readFileHandle = open(readFile, "w")
    for k, v in self.falseNegData.items():
      self.printStats(readFileHandle, v, v["hReads"], v["cReads"])


  def printReads(self, reads, snippet, fsock):
    # first sort the reasd into reverse and forward
    fwdReads, revReads = self.binarySplit(reads, lambda x: x[3] == 'F')
    revSnippet = self.reverseComplement(snippet)
    alignedFwd, snippet  = self.alignReads(fwdReads, snippet)
    alignedRev, revSnippet  = self.alignReads(revReads, revSnippet)
    self.printMismatches(alignedFwd, fwdReads, snippet, fsock)
    self.printMismatches(alignedRev, revReads, revSnippet, fsock)


  def alignReads(self, reads, snippet, fsock):
    # Find the maximum offset read
    maxOffset = max([(offset - relativeTo) for (_,_,_,_, offset, relativeTo) in reads])
    aligned = ["-" * (maxOffset + (relativeTo - offset)) + read
      for (read,_,_,_, offset, relativeTo) in reads]
    return aligned, ("-"*maxOffset + snippet)


  def binarySplit(self, reads, cond):
    return [t for t in reads if cond(t)], [f for f in reads if not cond(f)]
    

  def printMismatches(self, aligned, readsData, snippet, fsock):
    fsock.write("Snippet:\n")
    fsock.write(snippet + "\n")
    for read in aligned:
      # determine mismatch positions (excluding gap differences)
      mismatch_idx = [i for i,x in enumerate(zip(read, snippet),0) if x[0] != x[1] and
      x[0] != '-' and x[1] != '-']


      # determine the length of matching characters between mismatches
      gaps = [(l[i] - l[i-1] - 1) for i in range(0, len(mismatch_idx))]
      gaps.insert(0, diffs[0])

      # print the read
      fsock.write(read + "\n")  # print the read
      # print a pointer to every mismatch below it in the file
      for gap in gaps:
        fsock.write(' '*gap[0] + '^')
      fsock.write('\n')


  def reverseComplement(self, s):
    complement = ""
    for ch in s:
      complement = complement + {'A':'T', 'T':'A', 'C':'G', 'G':'C'}[ch]
    return complement[::-1]


  def printStats(self, fileHandle, result, hReads = None, cReads = None, nmData
      = False):
    fileHandle.write("Cood: %d, %s -> %s\n" % (result["coordinate"],
          result["hBase"], result["cBase"]))

    fileHandle.write("Healthy Total Stats (Total, Remain, Perc)\n")
    fileHandle.write("%d, %d, %.2f\n" % (
      result["HealthyTot"]["Total"],
      result["HealthyTot"]["Rem"], 
      result["HealthyTot"]["Perc"])
    )

    # switch to print the NM data
    if nmData:
      fileHandle.write("Healhty N Stats (Total, Remain, Perc)\n")
      fileHandle.write("%d, %d, %2.f\n" % (
        result["HealthyN"]["Total"], 
        result["HealthyN"]["Rem"],
        result["HealthyN"]["Perc"])
      )

      fileHandle.write("Healthy M Stats (Total, Remain, Perc)\n")
      fileHandle.write("%d, %d, %.2f\n" % (
        result["HealthyM"]["Total"], 
        result["HealthyM"]["Rem"],
        result["HealthyM"]["Perc"])
      )

    if hReads != None:
      printReads(hReads, result["snippet"], fileHandle)
    
    fileHandle.write("Cancer Total Stats (Total, Remain, Perc)\n")
    fileHandle.write("%d, %d, %.2f\n" % (
      result["CancerTot"]["Total"],
      result["CancerTot"]["Rem"], 
      result["CancerTot"]["Perc"])
    )

    # switch to print the NM data
    if nmData:
      fileHandle.write("Healhty N Stats (Total, Remain, Perc)\n")
      fileHandle.write("%d, %d, %2.f\n" % (
        result["CancerN"]["Total"], 
        result["CancerN"]["Rem"],
        result["CancerN"]["Perc"])
      )

      fileHandle.write("Cancer M Stats (Total, Remain, Perc)\n")
      fileHandle.write("%d, %d, %.2f\n" % (
        result["CancerM"]["Total"], 
        result["CancerM"]["Rem"],
        result["CancerM"]["Perc"])
      )

    if cReads != None:
      printReads(cReads, result["snippet"], fileHandle)

    fileHandle.write("\n\n")
