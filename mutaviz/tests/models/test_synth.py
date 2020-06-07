from mutaviz.models.synth import Synthesizer, InvalidSequenceLengthException


class TestSynth:
    def test_given_an_adn_seq_synthesizes_the_seq_returning_the_correct_protein(self):
        albumin = 'ATGGATGCACACAAGAGTGAGGTTGCTCATCGGTTTAAAGATTTGGGAGAAGAAAATTTCAAAGCCTTGGTGTTGATTGCCTTTGCTCAGTATCTTC' \
                  'AGCAGTGTCCATTTGAAGATCATGTAAAATTAGTGAATGAAGTAACTGAATTTGCAAAAACATGTGTTGCTGATGAGTCAGCTGAAAATTGTGACAA' \
                  'ATCACTTCATACCCTTTTTGGAGACAAATTATGCACAGTTGCAACTCTTCGTGAAACCTATGGTGAAATGGCTGACTGCTGTGCAAAACAAGAACCT' \
                  'GAGAGAAATGAATGCTTCTTGCAACACAAAGATGACAATCCAAATCTCCCCCGATTGGTGAGACCAGAGGTTGATGTGATGTGCACTGCTTTTCATG' \
                  'ACAATGAAGAGACATTTTTGAAAAAATACTTATATGAAATTGCCAGAAGACATCCTTACTTTTATGCCCCGGAACTCCTTTTCTTTGCTAAAAGGTA' \
                  'TAAAGCTGCTTTTACAGAATGTTGCCAAGCTGCTGATAAAGCAGCCTGCCTGTTGCCAAAGCTCGATGAACTTCGGGATGAAGGGAAGGCTTCGTCT' \
                  'GCCAAACAGAGACTCAAGTGTGCCAGTCTCCAAAAATTTGGAGAAAGAGCTTTCAAAGCATGGGCAGTAGCTCGCCTGAGCCAGAGATTTCCCAAAG' \
                  'CTGAGTTTGCAGAAGTTTCCAAGTTAGTGACAGATCTTACCAAAGTCCACACGGAATGCTGCCATGGAGATCTGCTTGAATGTGCTGATGACAGGGC' \
                  'GGACCTTGCCAAGTATATCTGTGAAAATCAAGATTCGATCTCCAGTAAACTGAAGGAATGCTGTGAAAAACCTCTGTTGGAAAAATCCCACTGCATT' \
                  'GCCGAAGTGGAAAATGATGAGATGCCTGCTGACTTGCCTTCATTAGCGGCTGATTTTGTTGAAAGTAAGGATGTTTGCAAAAACTATGCTGAGGCAA' \
                  'AGGATGTCTTCTTGGGCATGTTTTTGTATGAATATGCAAGAAGGCATCCTGATTACTCTGTCGTACTGCTGCTGAGACTTGCCAAGACATATGAAAC' \
                  'CACTCTAGAGAAGTGCTGTGCCGCTGCAGATCCTCATGAATGCTATGCCAAAGTGTTCGATGAATTTAAACCTCTTATGGAAGAGCCTCAGAATTTA' \
                  'ATCAAACAAAATTGTGAGCTTTTTGAGCAGCTTGGAGAGTACAAATTCCAGAATGCGCTATTAGTTCGTTACACCAAGAAAGTACCCCAAGTGTCAA' \
                  'CTCCAACTCTTGTAGAGGTCTCAAGAAACCTAGGAAAAGTGGGCAGCAAATGTTGTAAACATCCTGAAGCAAAAAGAATGCCCTGTGCAGAAGACTA' \
                  'TCTATCCGTGGTCCTGAACCAGTTATGTGTGTTGCATGAGAAAACGCCAGTAAGTGACAGAGTCACCAAATGCTGCACAGAATCCTTGGTGAACAGG' \
                  'CGACCATGCTTTTCAGCTCTGGAAGTCGATGAAACATACGTTCCCAAAGAGTTTAATGCTGAAACATTCACCTTCCATGCAGATATATGCACACTTT' \
                  'CTGAGAAGGAGAGACAAATCAAGAAACAAACTGCACTTGTTGAGCTTGTGAAACACAAGCCCAAGGCAACAAAAGAGCAACTGAAAGCTGTTATGGA' \
                  'TGATTTCGCAGCTTTTGTAGAGAAGTGCTGCAAGGCTGACGATAAGGAAACCTGCTTTGCCGAGGAGGGTAAAAAACTTGTTGCTGCAAGTCAAGCT' \
                  'GCCTTAGGCTTATAA'

        seq = Synthesizer.accepting(Synthesizer.DNA, albumin).run()

        assert seq == 'MDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCA' \
                      'KQEPERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLFFAKRYKAAFTECCQAADKAACLLPKLDEL' \
                      'RDEGKASSAKQRLKCASLQKFGERAFKAWAVARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLKEC' \
                      'CEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYARRHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYA' \
                      'KVFDEFKPLMEEPQNLIKQNCELFEQLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVVLNQLCVLH' \
                      'EKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTLSEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEK' \
                      'CCKADDKETCFAEEGKKLVAASQAALGL'

    def test_given_an_arn_seq_synthesizes_the_seq_returning_the_correct_protein(self):
        uvi1d = 'AUGGAUGCACACAAGAGUGAGGUUGCUCAUCGGUUUAAAGAUUUGGGAGAAGAAAAUUUCAAAGCCUUGGUGUUGAUUGCCUUUGCUCAGUAUCUUC' \
                 'AGCAGUGUCCAUUUGAAGAUCAUGUAAAAUUAGUGAAUGAAGUAACUGAAUUUGCAAAAACAUGUGUUGCUGAUGAGUCAGCUGAAAAUUGUGACAA' \
                 'AUCACUUCAUACCCUUUUUGGAGACAAAUUAUGCACAGUUGCAACUCUUCGUGAAACCUAUGGUGAAAUGGCUGACUGCUGUGCAAAACAAGAACCU' \
                 'GAGAGAAAUGAAUGCUUCUUGCAACACAAAGAUGACAAUCCAAAUCUCCCCCGAUUGGUGAGACCAGAGGUUGAUGUGAUGUGCACUGCUUUUCAUG' \
                 'ACAAUGAAGAGACAUUUUUGAAAAAAUACUUAUAUGAAAUUGCCAGAAGACAUCCUUACUUUUAUGCCCCGGAACUCCUUUUCUUUGCUAAAAGGUA' \
                 'UAAAGCUGCUUUUACAGAAUGUUGCCAAGCUGCUGAUAAAGCAGCCUGCCUGUUGCCAAAGCUCGAUGAACUUCGGGAUGAAGGGAAGGCUUCGUCU' \
                 'GCCAAACAGAGACUCAAGUGUGCCAGUCUCCAAAAAUUUGGAGAAAGAGCUUUCAAAGCAUGGGCAGUAGCUCGCCUGAGCCAGAGAUUUCCCAAAG' \
                 'CUGAGUUUGCAGAAGUUUCCAAGUUAGUGACAGAUCUUACCAAAGUCCACACGGAAUGCUGCCAUGGAGAUCUGCUUGAAUGUGCUGAUGACAGGGC' \
                 'GGACCUUGCCAAGUAUAUCUGUGAAAAUCAAGAUUCGAUCUCCAGUAAACUGAAGGAAUGCUGUGAAAAACCUCUGUUGGAAAAAUCCCACUGCAUU' \
                 'GCCGAAGUGGAAAAUGAUGAGAUGCCUGCUGACUUGCCUUCAUUAGCGGCUGAUUUUGUUGAAAGUAAGGAUGUUUGCAAAAACUAUGCUGAGGCAA' \
                 'AGGAUGUCUUCUUGGGCAUGUUUUUGUAUGAAUAUGCAAGAAGGCAUCCUGAUUACUCUGUCGUACUGCUGCUGAGACUUGCCAAGACAUAUGAAAC' \
                 'CACUCUAGAGAAGUGCUGUGCCGCUGCAGAUCCUCAUGAAUGCUAUGCCAAAGUGUUCGAUGAAUUUAAACCUCUUAUGGAAGAGCCUCAGAAUUUA' \
                 'AUCAAACAAAAUUGUGAGCUUUUUGAGCAGCUUGGAGAGUACAAAUUCCAGAAUGCGCUAUUAGUUCGUUACACCAAGAAAGUACCCCAAGUGUCAA' \
                 'CUCCAACUCUUGUAGAGGUCUCAAGAAACCUAGGAAAAGUGGGCAGCAAAUGUUGUAAACAUCCUGAAGCAAAAAGAAUGCCCUGUGCAGAAGACUA' \
                 'UCUAUCCGUGGUCCUGAACCAGUUAUGUGUGUUGCAUGAGAAAACGCCAGUAAGUGACAGAGUCACCAAAUGCUGCACAGAAUCCUUGGUGAACAGG' \
                 'CGACCAUGCUUUUCAGCUCUGGAAGUCGAUGAAACAUACGUUCCCAAAGAGUUUAAUGCUGAAACAUUCACCUUCCAUGCAGAUAUAUGCACACUUU' \
                 'CUGAGAAGGAGAGACAAAUCAAGAAACAAACUGCACUUGUUGAGCUUGUGAAACACAAGCCCAAGGCAACAAAAGAGCAACUGAAAGCUGUUAUGGA' \
                 'UGAUUUCGCAGCUUUUGUAGAGAAGUGCUGCAAGGCUGACGAUAAGGAAACCUGCUUUGCCGAGGAGGGUAAAAAACUUGUUGCUGCAAGUCAAGCU' \
                 'GCCUUAGGCUUAUAA'

        seq = Synthesizer.accepting(Synthesizer.RNA, uvi1d).run()

        assert 'MDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPER' \
               'NECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLFFAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQR' \
               'LKCASLQKFGERAFKAWAVARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLKECCEKPLLEKSHCIAEVEND' \
               'EMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYARRHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLMEEPQNLIKQNCELF' \
               'EQLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVVLNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVD' \
               'ETYVPKEFNAETFTFHADICTLSEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLVAASQAALGL' == seq

    def test_given_an_already_synthesized_protein_it_returns_the_same(self):
        protein = 'MDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCA'

        seq = Synthesizer.accepting(Synthesizer.PROTEIN, protein).run()

        assert seq == protein

    def test_fails_when_the_given_sequence_length_is_not_multiple_of_3(self):
        try:
            Synthesizer.accepting(Synthesizer.DNA, "AA").run()
        except InvalidSequenceLengthException as e:
            assert e.message == "Sequence length has to multiple a multiple of 3"

    def test_fails_when_the_given_sequence_is_shorter_than_60_chars(self):
        try:
            Synthesizer.accepting(Synthesizer.DNA, "AAA").run()
        except InvalidSequenceLengthException as e:
            assert e.message == "Sequence length not valid. Length has to be greater than 60"
