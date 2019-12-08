from abc import abstractmethod


class Synthesizer:
    STOP = '_'
    ADN = 'ADN'
    ARN = 'ARN'
    PROTEIN = 'PROTEIN'

    @classmethod
    def accepting(cls, synth_type, seq):
        return next(
            (synth for synth in cls.__subclasses__() if synth.has_type(synth_type)),
            AdnSynthesizer
        )(seq)

    @classmethod
    def has_type(cls, synth_type):
        return synth_type == cls.synth_type()

    @classmethod
    def synth_type(cls):
        raise NotImplementedError

    def __init__(self, seq):
        self.__seq = seq

    def run(self):
        protein = ""
        for i in range(0, len(self.__seq), 3):
            codon = self.__seq[i:i + 3]
            if self.__is_stop_codon(codon):
                break
            protein += self.table()[codon]
        return protein

    def __is_stop_codon(self, codon):
        return self.table()[codon] == self.STOP

    @abstractmethod
    def table(self):
        raise NotImplementedError

    @property
    def seq(self):
        return self.__seq


class AdnSynthesizer(Synthesizer):
    @classmethod
    def synth_type(cls):
        return cls.ADN

    def table(self):
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


class ArnSynthesizer(Synthesizer):
    @classmethod
    def synth_type(cls):
        return cls.ARN

    def table(self):
        return {
            'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
            'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
            'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
            'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
            'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
            'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
            'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
            'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
            'UAC': 'Y', 'UAU': 'Y', 'UGG': 'W', 'UGC': 'C',
            'UGU': 'C', 'UAA': self.STOP, 'UAG': self.STOP, 'UGA': self.STOP
        }


class ProteinSynthesizer(Synthesizer):
    @classmethod
    def synth_type(cls):
        return cls.PROTEIN

    def table(self):
        pass

    def run(self):
        return self.seq