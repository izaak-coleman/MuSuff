import sys
import pprint
class GenomeCoordinateExtractor:
  """GenomeCoordinatesExtractor takes a fastq file and a list of coordinates,
  and prints out the DNA sequence from -flanking_dist to +flanking_dist
  either side of the coordinate from the fastq file"""
  def __init__(self, fasta_filename, coordinate_filename, flanking_dist, 
      out_filename):
    self.flanking_dist = int(flanking_dist)
    self.out_filename = out_filename
    self.fasta_data = {}
    self.coordinates = []
    self.coord_seq_pairs = []

    self.extractFastaData(fasta_filename)
    self.extractCoordinates(coordinate_filename)

  def extractFastaData(self, fasta_filename):
    """ Loads each fasta sequence in the file into dictionary self.fastq_data.
      The sequence header is used as the key, whilst values are of the
      for (endNPos, sequence), where endNPos is the index of the first leftmost
      non-N character in sequence, and sequence is the fasta string."""
    fasta_handle = open(fasta_filename, 'r')
    header, sequence = "", ""
    for line in fasta_handle:
      if line[0] == ">":      # new sequence, new dict entry
        if header != "":      # skip loop initialization val
          self.fasta_data[header] = (self.nonNStartPos(sequence), sequence)
        header = line.strip()
        sequence = ""
      else:
        sequence = sequence + line.strip()
    self.fasta_data[header] = (self.nonNStartPos(sequence), sequence) # add last
    fasta_handle.close()

  def nonNStartPos(self, sequence):
    """ Returns the index of the first leftmost non-N char in sequence."""
    n = 0
    for char in sequence:
      if char == "N":
        n = n + 1
      else:
        return n

  def extractCoordinates(self, coordinate_filename):
    """ Takes a file with a list of coodinates of form (header,loc, __, __) where
      __ can be anything and loc is an integer value representing the
      coordinate of the mutation. Note, format MUST remain comma delim, 
      or modify line slicing. """
    coordinate_handle = open(coordinate_filename, 'r')
    SWITCH_INDEX = -1                   # bwa 1-based, pystring 0-based
    for line in coordinate_handle:
      header = line[1:line.find(",")]
      line = line[line.find(",")+1:]
      coordinate = line[:line.find(",")]
      self.coordinates.append( (header, int(coordinate) + SWITCH_INDEX) ) 
    coordinate_handle.close()

  def extractGenomicSequences(self):
    """Function loops through coordinates and extracts the (nEndPos, seq) value
      associated with the coordinate header (key). The function then
      extracts the section of DNA around the coordinate that is +/- the 
      flanking distance of the coordinate and loads data into a list."""
    for (header, coordinate) in self.coordinates:
      (n_pos, sequence) = self.fasta_data[header] # retrieve the loc, seq tuple

      left_most = coordinate - self.flanking_dist
      right_most = coordinate + self.flanking_dist + 1
      if left_most < 0:                     # boundary check
        left_most = 0
      if right_most > len(sequence):
        right_most = len(sequence)

      subseq = sequence[left_most:coordinate]
      subseq = subseq + sequence[coordinate:right_most]
      self.coord_seq_pairs.append( (header, coordinate, subseq) )

  def printCoordSeqPairs(self):
    out_handle = open(self.out_filename, 'w')
    for (header, coordinate, subseq) in self.coord_seq_pairs:
      metaData = header + ": " + str(coordinate)
      out_handle.write(metaData)
      out_handle.write("\n")
      out_handle.write(subseq)
      out_handle.write("\n")
      out_handle.write("\n")
    out_handle.close()

  def print_fasta_data(self):
    pprint.pprint(self.fasta_data)
  
  def print_coordinates(self):
    pprint.pprint(self.coordinates)

  def print_flanking_dist(self):
    print self.flanking_dist

if __name__ == "__main__":
  if len(sys.argv) != 5:
    print "Usage <exe> <fasta_file> <coordinate_file> <flanking_dist> <oFileName>"

  ge = GenomeCoordinateExtractor(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
  ge.extractGenomicSequences()
  ge.printCoordSeqPairs()
