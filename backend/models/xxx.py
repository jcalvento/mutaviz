import os

from Bio.Blast import NCBIWWW, NCBIXML
from Bio.PDB import PDBList
from models.synth import Synthesizer
from functools import reduce

from models.aligner import Aligner
from models.alignment_formatter import AlignmentFormatter
from models.file_name_generator import FileNameGenerator


class XXX:
    def __init__(self, seq_string):
        self.__seq_string = seq_string
        self.__protein_chain = None
        self.__matching_sequence = None
        self.__pdb_key = None

    def synthesize(self):
        self.__protein_chain = Synthesizer.accepting(Synthesizer.ADN, self.__seq_string[0:]).run()
        return self.__protein_chain

    def blast(self):
        scan_result = NCBIWWW.qblast(
            "blastp", "pdb", self.__protein_chain, word_size=2, threshold=200000, matrix_name="BLOSUM62", gapcosts="11 1"
        )
        blast_records = NCBIXML.read(scan_result)
        most_similar_structure = reduce(lambda result, output: self.__most_similar_between(result, output), blast_records.alignments)
        self.__matching_sequence = most_similar_structure.hsps[0].sbjct

    def fetch_pdb(self):
        self.__pdb_key = self.__matching_sequence.accession.split("_")[0]
        # me traigo el pdb del primer resultado
        pdb_file_path = PDBList().retrieve_pdb_file(self.__pdb_key, pdir='backend/atom_files', file_format="pdb")

        # no renombro para que tenga el nombre que espera modeller (xxxx.pdb)
        new_pdb_path = 'backend/atom_files/%s.pdb' % self.__pdb_key
        os.rename(pdb_file_path, new_pdb_path)

        return new_pdb_path

    def align(self):
        align_file_path = Aligner(
            path="backend/alignments",
            sequence_name="gilada",
            sequence_1=self.__protein_chain,
            pdb_key=self.__pdb_key,
            sequence_2=self.__matching_sequence
        ).file_align()

        pir_file_path = FileNameGenerator().random(extension='pir', path='backend/alignments')
        AlignmentFormatter(align_file_path, pir_file_path).to_pir()

    def __most_similar_between(self, alignment, another_alignment):
        if self.__identity_percentage(alignment) > self.__identity_percentage(another_alignment):
            return alignment
        return another_alignment

    def __identity_percentage(self, alignment):
        hsp = alignment.hsps[0]
        return round(hsp.identities / hsp.align_length, 2) == 1

# exact_protein = next((alignment for alignment in blast_records.alignments if is_same_protein(alignment)), False)