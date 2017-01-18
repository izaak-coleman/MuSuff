import sys
import pprint
class GenomeCoordinateExtractor:
  """GenomeCoordinatesExtractor takes a fastq file and a list of coordinates,
  and prints out the DNA sequence from -flanking_dist to +flanking_dist
  either side of the coordinate from the fastq file"""

  def __init__(self, fasta_filename, coordinate_filename, flanking_dist):
    self.flanking_dist = flanking_dist
    self.extractFastaData(fasta_filename)
    self.extractCoordinates(coordinate_filename)
    print "hello"

  def extractFastaData(self, fasta_filename):
    fasta_handle = open(fasta_filename, 'r')
    self.fasta_data = {}
    header, sequence = "", ""
    for line in fasta_handle:
      if line[0] == ">":      # start new entry
        if header != "":
          self.fasta_data[header] = (self.nonNStartPos(sequence), sequence)
        header = line.strip()
        sequence = ""
      else:
        sequence = sequence + line.strip()
    self.fasta_data[header] = (self.nonNStartPos(sequence), sequence)
    fasta_handle.close()

  def nonNStartPos(self, sequence):
    n = 0
    for char in sequence:
      if char == "N":
        n = n + 1
      else:
        return n

  def extractCoordinates(self, coordinate_filename):
    coordinate_handle = open(coordinate_filename, 'r')
    self.coordinates = []
    for line in coordinate_handle:
      coordinate = line[1:line.find(",")]
      self.coordinates.append(int(coordinate))
    coordinate_handle.close()

  def print_fasta_data(self):
    pprint.pprint(self.fasta_data)
  
  def print_coordinates(self):
    pprint.pprint(self.coordinates)

  def print_flanking_dist(self):
    print self.flanking_dist

if __name__ == "__main__":
  if len(sys.argv) != 4:
    print "Usage <exe> <fasta_file> <coordinate_file> <flanking_dist>"

  ge = GenomeCoordinateExtractor(sys.argv[1], sys.argv[2], sys.argv[3])
  ge.print_fasta_data()
  ge.print_coordinates()
  ge.print_flanking_dist()
