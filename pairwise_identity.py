import argparse
import copy
import multiprocessing

from Bio import SeqIO
from Bio.Alphabet import generic_protein
from lib.align import muscle_wrapper


def main():
    """
    Pairwise identity is a script which takes a protein sequence file that's full of protein sequences and calculates
    the percent identity of all possible pairwise combinations (basically the cartesian product) of the proteins
    """

    sh_parse = argparse.ArgumentParser(description="Cartesian product pairwise identity of a protein sequence file")

    sh_parse.add_argument("-i", "--input", dest="infile", help="Input FILE", metavar="FILE", required=True)
    sh_parse.add_argument("-o", "--output", dest="outfile", help="Output FILE", metavar="FILE", required=True)
    sh_parse.add_argument("-t", "--type", dest="type", help="Input file TYPE", metavar="TYPE", default="fasta")
    sh_parse.add_argument("-c", "--cpu", dest="cores", help="Number of CORES to use", metavar="CORES", default=2,
                          type=int)

    sh_args = sh_parse.parse_args()

    pairwise_id(sh_args.infile, sh_args.outfile, in_type=sh_args.type, cores=sh_args.cores)


def pairwise_id(infile, outfile, in_type="fasta", cores=1):
    """
    pairwise_id is the worker function which implements the pairwise identity script

    Required Argument:

    :param infile: str
        Path to the input sequence file
    :param outfile: str
        Path to the output file

    Keyword Argument:

    :param in_type: str
        Input file type (defaults to "fasta", can use anything BioPython takes)
    :param cores: int
        Number of separate cores to use
    """

    # Read the sequence file into a dict of SeqRecords keyed by ID
    sequence_list = {}
    with open(infile, mode="rU") as in_fh:
        for seq_record in SeqIO.parse(in_fh, format=in_type):
            sequence_list[seq_record.id] = seq_record

    # Generator which yields all the possible unique pairwise combinations (not permutations)
    def _cart_generator(seq_list):
        right_list = copy.copy(seq_list)
        for left_id in seq_list:
            left_side_seq = right_list.pop(left_id)
            for right_id in right_list:
                right_side_seq = seq_list[right_id]
                yield left_side_seq, right_side_seq

    mp = multiprocessing.Pool(processes=cores)

    # Write a TSV file with the pairwise identity and similarity scores
    with open(outfile, mode="w+") as out_fh:
        print("Protein\tProtein\tIdentity\tSimilarity", file=out_fh)
        for left_id, right_id, identity, similar in mp.imap_unordered(pairwise_muscle_and_identity,
                                                                      _cart_generator(sequence_list)):
            print("{}\t{}\t{}\t{}".format(left_id, right_id, identity, similar), file=out_fh)


def pairwise_muscle_and_identity(input):
    """
    pairwise_muscle_and_identity calls the MUSCLE wrapper. Is a separete function to facilitate multiprocessing

    :param input: SeqRecord, SeqRecord
        The two sequences to align and determine identity for
    :return: str, str, float, float
        Returns the first sequence ID string, the second sequence ID string, the identity ratio, and the similarity
        ratio
    """

    left_seq, right_seq = input
    align, identity_ratio, similar_ratio = muscle_wrapper([left_seq, right_seq], alphabet=generic_protein)

    return left_seq.id, right_seq.id, identity_ratio, similar_ratio


if __name__ == '__main__':
    main()
