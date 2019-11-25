from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.PDB import *
from Bio.Align.Applications import ClustalOmegaCommandline
import os
from models.synth import Synthesizer
from models.aligner import Aligner
from models.alignment_formatter import AlignmentFormatter
from models.file_name_generator import FileNameGenerator
from modeller import *

def read_seq(input_file):
    with open(input_file, "r") as f:
        lines = f.read().splitlines(True)
    seq = [line for line in lines if not line.startswith(">")]
    seq = "".join(seq)
    seq = seq.replace("\n", "")
    seq = seq.replace("\r", "")
    return seq

def is_same_protein(alignment):
    # Hsp_identity / Hsp_align - len
    hsp = alignment.hsps[0]
    return round(hsp.identities / hsp.align_length, 2) == 1


if __name__ == "__main__":
    seq_string = read_seq("backend/rat_albumin_dna.fasta")
    sequence_data = Synthesizer.accepting(Synthesizer.ADN, seq_string[0:]).run()
    print(sequence_data)
    result_sequence = NCBIWWW.qblast(
        "blastp", "pdb", sequence_data, word_size=2, threshold=200000, matrix_name="BLOSUM62", gapcosts="11 1"
    )
    # Parametrizar word size, matrix_name
    blast_records = NCBIXML.read(result_sequence)
    # exact_protein = next((alignment for alignment in blast_records.alignments if is_same_protein(alignment)), False)
    exact_protein = blast_records.alignments[0]
    if exact_protein:
        print("Esta es igual")
        # mutamos el adn
        # sintetizamos la proteina
        # alineamos la proteina sintetizada con la proteina original
        # pasamos el alineamiento a pir
        # modelamos con el pdb que tira el resultado del blast que coincide 100%
        # elegimos el mejor modelo (es el pdb con "DOPE score" mas bajo)
        # le mostramos al usuario los dos pdb's el original y el modelo

    # sino
    #   le damos a elegir al usuario con que pdb de los primeros 5 resultados del blast quiere modelar
    #   alineamos la secuencia que tenemos con la del resultado que el usuario eligio
    #   pasamos el alineamiento a pir
    #   modelamos con el pdb del resultado que eligio el usuario
    #   elegimos el mejor modelo (es el pdb con "DOPE score" mas bajo)
    #   mutamos el adn
    #   sintetizamos la proteina
    #   alineamos la proteina sintetizada con la proteina original
    #   pasamos el alineamiento a pir
    #   modelamos con el pdb resultado del modelado anterior
    #   elegimos el mejor modelo (es el pdb con "DOPE score" mas bajo)
    #   le mostramos al usuario los dos pdb's el modelado de la proteina original y el modelado de la mutante

    pdb_key = exact_protein.accession.split("_")[0]
    print(pdb_key)
    
    # me traigo el pdb del primer resultado
    pdb_file_path = PDBList().retrieve_pdb_file(pdb_key, pdir='backend/atom_files', file_format="pdb")
    coso = 'backend/atom_files/%s.pdb' % pdb_key

    # no renombro para que tenga el nombre que espera modeller (xxxx.pdb)
    os.rename(pdb_file_path, coso)

    structure = PDBParser().get_structure(pdb_key, coso)
    for chain in structure.get_chains():
        print(chain.get_id())

    # busco la secuencia problema en el primer resultado del blast
    my_sequence = blast_records.alignments[0].hsps[0].query

    # buco la secuencia que machea en el primer resultado del blast
    matching_sequence = blast_records.alignments[0].hsps[0].sbjct

    align_file_path = Aligner(
        path="backend/alignments",
        sequence_name="gilada",
        sequence_1=my_sequence,
        pdb_key=pdb_key,
        sequence_2=matching_sequence
    ).file_align()

    pir_file_path = FileNameGenerator().random(extension='pir', path='backend/alignments')
    AlignmentFormatter(align_file_path, pir_file_path).to_pir()

