from modeller import *
from modeller.automodel import *

# target = 'result'

# in_file = "backend/human_vs_rat_albumin.fasta"
# out_file = "backend/%s.fasta" % target

# clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=True)
# print(clustalomega_cline)
# os.system(str(clustalomega_cline))

# # FASTA to PIR
# e = environ()

# a = alignment(e, file=out_file, align_codes='all', alignment_format='FASTA')
# a.write(file='backend/%s.pir'%target, alignment_format='PIR')

# file = PDBList().retrieve_pdb_file('4BKE', pdir='backend/atom_files', file_format="pdb")
# os.rename(file, 'backend/atom_files/4BKE.pdb')

# Homology modeling with multiple templates
from backend.models.modeller import Modeller

if __name__ == "__main__":
    # log.verbose()    # request verbose output
    # env = environ()  # create a new MODELLER environment to build this model in
    #
    # # directories for input atom files
    # env.io.atom_files_directory = 'atom_files'
    # env.io.hetatm = True
    #
    # a = automodel(
    #   env,
    #   alnfile = 'aln_rat_3v03.pir', # alignment filename
    #   knowns = ('3V03'),     # codes of the templates
    #   sequence='NM_134326',
    #   assess_methods=(assess.DOPE),
    #
    # )
    # # code of the target
    # a.starting_model= 1                 # index of the first model
    # a.ending_model  = 5                 # index of the last model
    #                                    # (determines how many models to calculate)
    # a.make()                            # do the actual homology modeling
    #
    # from functools import reduce
    #
    #
    # def funcion(x, y):
    #   if x['DOPE score'] > y['DOPE score']:
    #     return y
    #   return x
    #
    #
    # reduce(lambda result, output: funcion(result, output), a.chains[0].seq.outputs)
    print(Modeller().execute(alignment_file='aln_rat_3v03.pir', pdb_id='3V03', sequence='NM_134326'))
