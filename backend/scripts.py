from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

from backend.synth import Synthesizer


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
    seq_string = read_seq("serum_albumin.fasta")
    sequence_data = Synthesizer(seq_string[110:]).run()
    print(sequence_data)
    result_sequence = NCBIWWW.qblast(
        "blastp", "pdb", sequence_data, word_size=2, threshold=200000, matrix_name="BLOSUM62", gapcosts="11 1"
    )
    # Parametrizar word size, matrix_name
    print(result_sequence)
    blast_records = NCBIXML.read(result_sequence)
    exact_protein = next(alignment for alignment in blast_records.alignments if is_same_protein(alignment))
    if exact_protein:
        print("Esta es igual")


# Si es igual
#   mutamos el adn
#   obtenemos la proteina mutada
#   traemos la estructura pdb
#   con foldx modelamos la mutada

# Si no es igual
#   agarramos la mas parecida (identidad > y evalue <)
#   modelamos la original
#   con foldx modelamos la mutada
