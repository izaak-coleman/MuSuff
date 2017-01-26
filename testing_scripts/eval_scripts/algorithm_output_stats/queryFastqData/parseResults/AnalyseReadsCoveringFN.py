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
    reads = [read[:read.find("$") + 1] for read in entryData[self.READ_IDX:]]
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
    reads = [read[:read.find("$") + 1] for read in entryData[self.READ_IDX:]]
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

  def printKeyVals(self):
    k, v  = self.falseNegData.items()[0]
    print k
    print v
