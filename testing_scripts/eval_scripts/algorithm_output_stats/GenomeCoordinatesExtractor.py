class GenomeCoordinatesExtractor:
  """GenomeCoordinatesExtractor takes a fastq file and a list of coordinates,
  and prints out the DNA sequence from -flanking_dist to +flanking_dist
  either side of the coordinate from the fastq file"""

  def __init__(self, fasta_filename, coordinate_filename, flanking_dist):
    self.flanking_dist = flanking_dist
    extractFastaData(fasta_filename)
    extractCoordinate(coordinate_filename)

  def extractFastaData(self, fasta_filename):
    fasta_handle = open(fasta_filename, 'r')
    self.fasta_data = {}
  
    header, sequence = "", ""
    for line in fasta_handle:
      if line[0] == ">":      # start new entry
        if header != "":
          self.fasta_data[header] = (nonNStartPos(sequence), sequence)
        header = line.strip()
        sequence = ""
      else:
        sequence = sequence + line.strip()
    self.fasta_data[header] = (nonNStartPos(sequence), sequence)
    fasta_handle.close()
