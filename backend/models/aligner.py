from Bio.Align.Applications import ClustalOmegaCommandline
from random import choice
import string
import os
from models.file_name_generator import FileNameGenerator

class Aligner:
  def __init__(self, path=None, sequence_name=None, sequence_1=None, pdb_key=None, sequence_2=None):
    self.__path = path
    self.__sequence_name = sequence_name
    self.__pdb_key = pdb_key
    self.__sequence_1 = sequence_1
    self.__sequence_2 = sequence_2

  def file_align(self):
    fasta_path = self.__in_file(lambda file: self.write_fasta(file))
    aligment_file_path = self.random_file_name()

    clustalomega_cline = ClustalOmegaCommandline(infile=fasta_path, outfile=aligment_file_path, verbose=True, auto=True)
    os.system(str(clustalomega_cline))

    return aligment_file_path

  def __in_file(self, write_block):
    file_path = self.random_file_name()
    file = open(file_path, 'w')
    write_block(file)
    file.close()
    return file_path

  def write_fasta(self, file):
      file.write(">P1;%s" % self.__sequence_name) 
      file.write("\n")
      file.write(self.__sequence_1)
      file.write("\n")
      file.write("\n")
      file.write(">P1;%s" % self.__pdb_key) 
      file.write("\n")
      file.write(self.__sequence_2)

  def random_file_name(self):
    return FileNameGenerator().random(extension='fasta', path=self.__path)
