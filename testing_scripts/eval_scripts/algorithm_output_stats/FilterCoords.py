import sys
class FilterCoords:
  """Class reads an SNV locations files and reported SNV files from ICSMuFin.
    It parses the files, and reports properies of the reported SNVs.
    For example, it can reported the coords of all the SNVs that were not 
    reported in the required form by GenomeCoordinatesExtractor.py:
    (header, loc, _, _). Output is to stdout.  
      header is specified at command line, and corresponds to the 
    SNV_locations header. """

  def __init__(self, header):
    self.header = header
    self.SNV_loc_coords = []
    self.reported_SNV_coords = []
    self.IDX_CALIB = -1

  def getCoordsFromSNV_locations(self):
    print("Entering coordinates from SNV_locations file...")
    fileName = raw_input("Enter co-ord file, or done to finish: ")
    while fileName != "done" and fileName != "Done":
      ifileHandle = open(fileName, 'r')
      for fields in ifileHandle:
        fields = fields.strip()
        fields = fields.split("\t")
        self.SNV_loc_coords.append( \
            (self.header, int(fields[2]) + self.IDX_CALIB, fields[3], fields[4]) \
        )
      fileName = raw_input("Enter co-ord file, or done to finish: ")
    self.SNV_loc_coords = sorted(self.SNV_loc_coords, key=lambda x: x[1])
    print("Finished entering coordinates from SNV_locations file...\n\n")

  def getReportedSNVs(self):
    print("Entering coordinates from algorithm output file...")
    fileName = raw_input("Enter co-ord file, or done to finish: ")
    while fileName != "done" and fileName != "Done":
      ifileHandle = open(fileName, 'r')
      for fields in ifileHandle:
        fields = fields.strip()
        fields = fields[1:-1]       # slice off tuple brackets
        fields = fields.split(",")
        fields[2] = fields[2][2:-1]  # slice off the ' ' quotes
        fields[3] = fields[3][2:-1]  # slice off the ' ' quotes
        self.reported_SNV_coords.append( \
          (fields[0], int(fields[1]) + self.IDX_CALIB, fields[2], fields[3], ) \
        )
      fileName = raw_input("Enter co-ord file, or done to finish: ")
    self.reported_SNV_coords = sorted(self.reported_SNV_coords, key=lambda x: x[1])
    print("reported SNV len: %d", len(self.reported_SNV_coords))
    print("Finished entering coordinates from algorithm output file...")

  def printUnreportedCoordinates(self):
    allCoords = list(set(self.SNV_loc_coords)) + self.reported_SNV_coords
    unReportedCoords = \
      [coord for coord in allCoords if allCoords.count(elem) == 1]
    for coord in unReportedCoords:
      print coord

  def printList(self, l):
    for e in l:
      print e 
    


if __name__ == "__main__":
  if len(sys.argv) != 2:
    print("Usage: <exe> <SNV_locations>")
    sys.exit()
  obj = FilterCoords(sys.argv[1])
  obj.getCoordsFromSNV_locations()
  obj.getReportedSNVs()
  print("SNV_loc")
  obj.printList(list(set(obj.SNV_loc_coords)))
  print("Rep")
  obj.printList(obj.reported_SNV_coords)
