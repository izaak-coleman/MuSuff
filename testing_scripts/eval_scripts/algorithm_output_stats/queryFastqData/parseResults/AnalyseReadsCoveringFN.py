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
    reads = [(read, int(read_id), found) for (read, read_id, found) in reads]

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
    reads = [(read, int(read_id), found) for (read, read_id, found) in reads]
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
      if fun(len(v["hReads"])) or fun(len(v["cReads"])):
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
       result["HealthyTot"] = self.computeStats(v["hReads"], lambda x: True, hSet)
       result["HealthyM"]   = self.computeStats(v["hReads"], mSet, hSet) 
       result["HealthyN"]   = self.computeStats(v["hReads"], nSet, hSet)
       result["CancerTot"]  = self.computeStats(v["cReads"], lambda x: True, cSet)
       result["CancerN"]    = self.computeStats(v["cReads"], nSet, cSet)
       result["CancerM"]    = self.computeStats(v["cReads"], mSet, cSet) 
       result["hBase"] = v["hBase"]
       result["cBase"] = v["cBase"]

       result["hReads"] = "\n".join(["%s :: %d :: %s" % (read, idx, found)
        for (read, idx, found) in v["hReads"]])

       result["cReads"] = "\n".join(["%s :: %d :: %s" % (read, idx, found)
        for (read, idx, found) in v["cReads"]])

       result["coordinate"] = k

       results.append(result)


    # print stats
    statFileHandle = open(statFile, "w")
    for result in results:
      self.printStats(statFileHandle, result)

    # print stats and reads
    readFileHandle = open(readFile, "w")
    for result in results:
      self.printStats(readFileHandle, result, result["hReads"], result["cReads"])



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
      fileHandle.write(result["hReads"])
    
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
      fileHandle.write(result["cReads"])
    fileHandle.write("\n\n")
