#!/usr/bin/env python3

# Usage: python3 fasta_grader.py file.fasta [MAX_SCORE]

import sys
from itertools import product
from Bio.Seq import Seq
from Bio.Data import IUPACData

IUPACData.extended_protein_letters += "*"  # Adds stop codon

def fasta_parser(filename_to_grade):
    """
    Parses a FASTA formatted file and returns a dict with {sequence_name:
    sequence, ...}
    """
    infile = open(filename_to_grade, "r")
    sequence_data = {}
    count = 0
    _score_mod = 0
    alphabet_in_use = [IUPACData.ambiguous_dna_letters,
                       IUPACData.unambiguous_dna_letters,
                       IUPACData.unambiguous_dna_letters,
                       IUPACData.unambiguous_dna_letters,
                       IUPACData.unambiguous_dna_letters,
                       IUPACData.unambiguous_dna_letters,
                       IUPACData.unambiguous_dna_letters,
                       IUPACData.extended_protein_letters]
    for lines in infile:
        if lines.startswith(">"):
            seqname = lines.strip()[1:]
            sequence_data[seqname] = Seq("")
            count += 1
        elif lines.startswith("\n"):
            if count >= 8:
                break
        else:
            lines = lines.upper()
            if "U" in lines and count < 8:
                _score_mod += -1
                sys.stderr.write("Found 'U'. This is not RNA. -1\n")
                lines = lines.replace("U", "T")
            corrected = "".join([x for x in lines if x in alphabet_in_use[count
                                 - 1]])
            if lines.strip() != corrected:
                _score_mod += -2
                sys.stderr.write("Wrong characters found and removed. -2\n")
            sequence_data[seqname] += corrected.strip().upper()

    return sequence_data, _score_mod


def auto_grader(sequence_data, _score=100, seq_length=21):
    """
    Grades the homework, take a sequence dictionary as input, returns a score
    Default starting score is 100, every error will subtract from this
    Score should start lower for resubmissions, or late submissions
    """
    if any([len(x) != seq_length for x in list(sequence_data.values())[:-1]]):
        _score += -10
        sys.stderr.write("Not all sequences are %s bp long!\n" % seq_length)


    def _correct_encoding(seq, data_type):
        """
        Checks if all letters are valid DNA codes
        Takes a sequence as input and returns the score difference
        """
        if any([x not in data_type for x in seq]):
            sys.stderr.write("Some bases are not correctly encoded:\n%s\n" % seq)
            minus = -10
        else:
            minus = 0

        return minus


    def _extend_ambiguous_dna(seq):
        """
        return list of all possible sequences given an ambiguous DNA
        input
        Taken from https://stackoverflow.com/a/27552377/3091595
        """
        ambigs = IUPACData.ambiguous_dna_values
        return list(map("".join, product(*map(ambigs.get, seq))))


    seqlist = list(sequence_data.keys())

    # Ambiguous sequence
    ambig_seq = sequence_data[seqlist[0]]
    # Are all letters valid DNA?
    _score += _correct_encoding(ambig_seq, IUPACData.ambiguous_dna_values)
    # Are there 4 ambiguities?
    ambig_bases = sum([ambig_seq.count(x) for x in
                       IUPACData.ambiguous_dna_letters[4:]])
    ambig_penalty = -abs(ambig_bases - 4) * 3
    _score += ambig_penalty
    if ambig_penalty != 0:
        sys.stderr.write("Sequence %s does not have 4 ambiguities!\n" % ambig_seq)
    # Generate all possible disambiguations:
    disambiguations = _extend_ambiguous_dna(ambig_seq)

    # Disambiguation
    for disamb in seqlist[1:4]:
        _score += _correct_encoding(sequence_data[disamb],  # Valid DNA?
                                   IUPACData.unambiguous_dna_letters)
        if sequence_data[disamb] not in disambiguations:
            _score += -10
            sys.stderr.write("Sequence %s is not a disambiguation of %s.\n" %
                             (sequence_data[disamb], ambig_seq))

    # Reverse
    rev_seq = sequence_data[seqlist[4]]
    _score += _correct_encoding(rev_seq,  # Valid DNA?
                               IUPACData.unambiguous_dna_letters)
    if rev_seq[::-1] not in [sequence_data[x] for x in seqlist[1:4]]:
        _score += -10
        sys.stderr.write("Sequence %s is not the reverse of any of the "
                         "disambiguations.\n" % rev_seq)

    # Complement
    comp_seq = sequence_data[seqlist[5]]
    _score += _correct_encoding(comp_seq,  # Valid DNA?
                               IUPACData.unambiguous_dna_letters)
    if comp_seq.complement() not in [sequence_data[x] for x in seqlist[1:4]]:
        _score += -10
        sys.stderr.write("Sequence %s is not the complement of any of the "
                         "disambiguations.\n" % comp_seq)

    # Rev Comp
    revcomp_seq = sequence_data[seqlist[6]]
    _score += _correct_encoding(revcomp_seq,  # Valid DNA?
                               IUPACData.unambiguous_dna_letters)
    if revcomp_seq.reverse_complement() not in [sequence_data[x] for x in
                                                seqlist[1:4]]:
        _score += -10
        sys.stderr.write("Sequence %s is not the reverse-complement of any of the "
                         "disambiguations.\n" % revcomp_seq)

    # Translation
    translation = sequence_data[seqlist[7]]
    _score += _correct_encoding(translation,  # Valid Protein?
                               IUPACData.extended_protein_letters)
    if translation not in [sequence_data[x].translate() for x in seqlist[1:4]]:
        _score += -10
        sys.stderr.write("Sequence %s is not the translation of any of the "
                         "disambiguations.\n" % translation)
    #print([sequence_data[x].translate() for x in seqlist[1:4]])
    return _score


if __name__ == "__main__":
    try:
        SEQUENCES, score_mod = fasta_parser(sys.argv[1])
    except:
        sys.exit("This is not spec compliant! It has to be manually graded...")
    if len(sys.argv) == 2:
        SCORE = 100 + score_mod
    else:
        SCORE = int(sys.argv[2]) + score_mod
    FINAL_SCORE = auto_grader(SEQUENCES, SCORE)
    print("Final score is:\t%s\t%s/%s" % (sys.argv[1], FINAL_SCORE, sys.argv[2]))
