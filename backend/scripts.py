from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.PDB import *

from backend.models.synth import Synthesizer


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
    seq_string = read_seq("serum_albumin_dna.fasta")
    sequence_data = Synthesizer.accepting(Synthesizer.ADN, seq_string[110:]).run()
    print(sequence_data)
    result_sequence = NCBIWWW.qblast(
        "blastp", "pdb", sequence_data, word_size=2, threshold=200000, matrix_name="BLOSUM62", gapcosts="11 1"
    )
    # Parametrizar word size, matrix_name
    blast_records = NCBIXML.read(result_sequence)
    exact_protein = next(alignment for alignment in blast_records.alignments if is_same_protein(alignment))
    if exact_protein:
        print("Esta es igual")
    pdb_key = exact_protein.accession.split("_")[0]
    print(pdb_key)
    pdb_file_path = PDBList().retrieve_pdb_file(pdb_key, file_format="pdb")
    structure = PDBParser().get_structure(pdb_key, pdb_file_path)
    for chain in structure.get_chains():
        print(chain.get_id())


# Si es igual
#   mutamos el adn
#   obtenemos la proteina mutada
#   traemos la estructura pdb
#   con foldx modelamos la mutada

# Si no es igual
#   agarramos la mas parecida (identidad > y evalue <)
#   modelamos la original
#   con foldx modelamos la mutada
