import sys
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
    reads = [(read[:read.find("$") + 1], read[read.find("::")+3:])
              for read in entryData[self.READ_IDX:]]
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
    reads = [(read[:read.find("$") + 1], read[read.find("::")+3:])
              for read in entryData[self.READ_IDX:]]
    if tissueType == "H":
      self["hReads"] = reads
    else:
      self["cReads"] = reads



class AnalyseReadsCoveringFN:
  COORD_IDX = 2
  def __init__(self, fname):
    self.falseNegData = {}
    self.parseFile(fname)


  def parseFile(self, fname):
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

  def quantifyLostReads(self, filename):
    remainingReads = []
    results = []
    with open(filename, "r") as fileHandle:
      remainingReads = [int(line.strip()) for line in fileHandle]
    for k, v in self.falseNegData.items():
       result = {}
       result["remainingHealthy"] = "\n".join(["%s :: %d :: %r" % 
           (read, idx, (idx in remainingReads)) for (read, idx) in v["hReads"])
       result["remainingCancer"] = "\n".join(["%s :: %d :: %r" % 
           (read, idx, (idx in remainingReads)) for (read, idx) in v["cReads"])

       result["cPerc"] = len([idx for (_, idx) in v["cReads"] if 
         (idx in remainingReads)]) / float(len(v["cReads"]))

       result["hPerc"] = len([idx for (_, idx) in v["hReads"] if 
         (idx in remainingReads)]) / float(len(v["hReads"]))

       result["tPerc"] = (cPercent * len(v["cReads"]) + 
                         hPercent * len(v["hReads"])) / 
                         float(len(v["cReads"]) + len(v["hReads"]))
       result["coordinate"] = k

       results.append(result)
