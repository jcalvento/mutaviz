from functools import reduce

from Bio.Blast import NCBIWWW, NCBIXML

from backend.models.aligner import Aligner
from backend.models.alignment_formatter import AlignmentFormatter
from backend.models.file_name_generator import FileNameGenerator
from backend.models.modeller import Modeller
from backend.models.synth import Synthesizer


class Mutaviz:
    def __init__(self, seq_string, mutations, sequence_name):
        self.__mutations = mutations
        self.__seq_string = seq_string
        self.__sequence_name = sequence_name
        self.__protein_chain = None
        self.__matching_sequence = None
        self.__mutated_sequence = None
        self.__original_pdb_filename = None

    def process(self):
        self.__protein_chain = self.synthesize(self.__seq_string)
        self.blast()
        if self.__is_same_protein():
            self.mutate()
            mutated_protein = self.synthesize(self.__mutated_sequence)
            alignment_file = self.align(mutated_protein)
            model_filename = Modeller().execute(
                alignment_file=alignment_file, pdb_id=self.__pdb_key(), sequence=self.__sequence_name
            )['name']
            original_pdb_filename = self.__original_pdb_filename

    def synthesize(self, sequence):
        return Synthesizer.accepting(Synthesizer.ADN, sequence[0:]).run()

    def mutate(self):
        self.__mutated_sequence = self.__seq_string

        for (position, mutation) in self.__mutations.items():
            changing_seq = list(self.__mutated_sequence)
            changing_seq[position] = mutation
            self.__mutated_sequence = ''.join(changing_seq)

        return self.__mutated_sequence

    def blast(self):
        scan_result = NCBIWWW.qblast(
            "blastp", "pdb", self.__protein_chain, word_size=2, threshold=200000, matrix_name="BLOSUM62", gapcosts="11 1"
        )
        blast_records = NCBIXML.read(scan_result)
        most_similar_structure = reduce(lambda result, output: self.__most_similar_between(result, output), blast_records.alignments)
        self.__matching_sequence = most_similar_structure.hsps[0].sbjct

    def __pdb_key(self):
        return self.__matching_sequence.accession.split("_")[0]

    def align(self, protein):
        aligner = Aligner(path="backend/alignments", sequence_name=self.__sequence_name, sequence_1=protein,
                          pdb_key=self.__pdb_key(), sequence_2=self.protein_chain)
        align_file_path = aligner.file_align()
        self.__original_pdb_filename = aligner.pdb_file_path

        pir_file_path = FileNameGenerator().random(extension='pir', path='backend/alignments')
        return AlignmentFormatter(align_file_path, pir_file_path).to_pir()

    @property
    def protein_chain(self):
        return self.__protein_chain

    @property
    def original_sequence(self):
        return self.__seq_string

    def __most_similar_between(self, alignment, another_alignment):
        if self.__identity_percentage(alignment) > self.__identity_percentage(another_alignment):
            return alignment
        return another_alignment

    def __identity_percentage(self, alignment):
        hsp = alignment.hsps[0]
        return round(hsp.identities / hsp.align_length, 2)

    def __is_same_protein(self):
        # Hsp_identity / Hsp_align - len
        return self.__identity_percentage(self.__matching_sequence) == 1
