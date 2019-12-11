# def to_pir(self):
#     e = environ()
#     a = alignment(e, file=self.__in_file, align_codes='all', alignment_format='FASTA')
#     a.write(file=self.__out_file, alignment_format='PIR')
#     with open(self.__out_file, "r") as f:
#         file_lines = f.read().splitlines()
#     sequences = list(filter(lambda line: not (line.startswith("sequence") or line.startswith(">")), file_lines))
#     return self.__in_file(lambda file: self.__write_pir(file, sequences))
# 
# def __write_pir(self, file, sequences):
#     structure_info = self.__structure_info
#     sequence_name = structure_info["seq_name"]
#     file.write(">P1;%s" % sequence_name)
#     file.write("\n")
#     file.write("sequence:%s:1:%s::::0.00:0.00" % (sequence_name, len(sequences[0]) + 1))
#     file.write("\n")
#     file.write("\n")
#     file.write(sequences[0])
#     file.write("\n")
#     file.write("\n")
#     pdb_key = structure_info["pdb_key"]
#     file.write(">P1;%s" % pdb_key)
#     file.write("\n")
#     file.write(
#         f"structureX:{pdb_key}:{structure_info['start']}:{structure_info['chain']}:{structure_info['end']}:{structure_info['chain']}::::0.00:0.00")
#     file.write("\n")
#     file.write(sequences[1])
# 
# def __in_file(self, write_block):
#     file_path = self.random_file_name()
#     file = open(file_path, 'w')
#     write_block(file)
#     file.close()
#     return file_path
# 
# def random_file_name(self):
#     return FileNameGenerator().random(extension='pir', path=self.__path)