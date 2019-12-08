import os

from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.PDB import PDBList, PDBParser

from backend.models.file_name_generator import FileNameGenerator


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

        clustalomega_cline = ClustalOmegaCommandline(
            infile=fasta_path, outfile=aligment_file_path, verbose=True, auto=True
        )
        os.system(str(clustalomega_cline))

        return aligment_file_path

    def __in_file(self, write_block):
        file_path = self.random_file_name()
        file = open(file_path, 'w')
        write_block(file)
        file.close()
        return file_path

    def write_fasta(self, file):
        start, end = self.__structure_info()
        file.write(">P1;%s" % self.__sequence_name)
        file.write("\n")
        file.write("sequence:%s:1:%s::::0.00:0.00" % (self.__sequence_name, len(self.__sequence_1) + 1))
        file.write("\n")
        file.write(self.__sequence_1)
        file.write("\n")
        file.write(">P1;%s" % self.__pdb_key)
        file.write("\n")
        file.write('structureX:{key}:{start}:{end}::::0.00:0.00'.format(key=self.__pdb_key, start=start, end=end))
        file.write("\n")
        file.write(self.__sequence_2)

    def random_file_name(self):
        return FileNameGenerator().random(extension='fasta', path=self.__path)

    def __structure_info(self):
        new_file_name = self.__fetch_pdb()

        structure = PDBParser().get_structure(self.__pdb_key, new_file_name)
        chain = list(structure.get_chains())[0]  # que cadena tenemos que usar?
        residues = list(chain.get_residues())
        atoms = list(filter(lambda residue: not residue.get_id()[0].strip(), residues))
        first_atom_residue = atoms[0].parent.get_id()[1]
        last_atom_residue = atoms[-1].parent.get_id()[1]

        return {'chain': chain.get_id(), 'start': first_atom_residue, 'end': last_atom_residue}

    def __fetch_pdb(self):
        pdb_file_path = PDBList().retrieve_pdb_file(self.__pdb_key, pdir='atom_files', file_format="pdb")
        new_file_name = 'atom_files/%s.pdb' % self.__pdb_key
        # renombro para que tenga el nombre que espera modeller (xxxx.pdb)
        os.rename(pdb_file_path, new_file_name)
        self.__pdb_file_path = new_file_name
        return self.__pdb_file_path

    @property
    def pdb_file_path(self):
        return self.__pdb_file_path
