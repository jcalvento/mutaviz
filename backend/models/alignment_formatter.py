from modeller import *


class AlignmentFormatter:
    def __init__(self, in_file=None, out_file=None):
        self.__in_file = in_file
        self.__out_file = out_file

    def to_pir(self):
        e = environ()

        a = alignment(e, file=self.__in_file, align_codes='all', alignment_format='FASTA')
        a.write(file=self.__out_file, alignment_format='PIR')
        return self.__out_file
