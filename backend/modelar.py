from modeller import *
from modeller.automodel import *    # Load the automodel class
from Bio.Align.Applications import ClustalOmegaCommandline
import os

target = 'result'

in_file = "human_vs_rat_albumin.fasta"
out_file = "%s.fasta" % target

clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=True)
print(clustalomega_cline)
os.system(str(clustalomega_cline))

# FASTA to PIR
e = environ()

a = alignment(e, file='%s.fasta'%target, alignment_format='FASTA')
a.write(file='%s.pir'%target, alignment_format='PIR')


# Homology modeling with multiple templates
log.verbose()    # request verbose output
env = environ()  # create a new MODELLER environment to build this model in
 
# directories for input atom files
env.io.atom_files_directory = './:../atom_files'
env.io.hetatm = True


a = automodel(
  env,
  alnfile = '%s.pir'%target, # alignment filename
  knowns = '4BKE',     # codes of the templates
  sequence = 'coso')               # code of the target
a.starting_model= 1                 # index of the first model
a.ending_model  = 100                 # index of the last model
                                   # (determines how many models to calculate)
a.make()                            # do the actual homology modeling
