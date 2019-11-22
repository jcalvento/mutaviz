

class Synthesizer:
    STOP = '_'

    def __init__(self, seq):
        self.__adn_table = self.adn_table()
        self.__seq = seq

    def run(self):
        protein = ""
        for i in range(0, len(self.__seq), 3):
            codon = self.__seq[i:i + 3]
            if self.is_stop_codon(codon):
                break
            protein += self.__adn_table[codon]
        return protein

    def adn_table(self):
        return {
            'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
            'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
            'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
            'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
            'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
            'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
            'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
            'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
            'TAC': 'Y', 'TAT': 'Y', 'TGG': 'W', 'TGC': 'C',
            'TGT': 'C', 'TAA': self.STOP, 'TAG': self.STOP, 'TGA': self.STOP
        }

    def is_stop_codon(self, codon):
        return self.__adn_table[codon] == self.STOP

# ADN:
    # TGA
    # TAA
    # TAG
    #
    # ARN:
    # UAG
    # UAA
    # UGA
