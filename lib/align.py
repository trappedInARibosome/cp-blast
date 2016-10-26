import tempfile
import subprocess
import os

from Bio.Alphabet import generic_dna
from Bio import AlignIO
from Bio import SeqIO

def muscle_wrapper(list_of_seqrecords, alphabet=generic_dna):

    """
    muscle_wrapper uses the command line MUSCLE application, with the default options for alignment, to align a list of
    SeqRecord objects. Implemented with tempfiles instead of pipelining because of reasons I can't remember right now.

    The identity and similarity ratios are calculated so that they're the LARGEST possible number because that seemed
    more important for IP considerations. In some cases this may not be the most appropriate metric

    Edgar, R.C. (2004). MUSCLE: a multiple sequence alignment method with reduced time and space complexity.
    BMC bioinformatics 5, 113.

    Required Arguments:

    :param list_of_seqrecords:  [Bio.SeqRecord]
        Sequence record objects in a list that will be used to build an alignment.

    Keyword Arguments:

    :param alphabet:    Bio.Alphabet
        BioPython alphabet object to use when reading in the alignment
        Should have just implemented detection but didnt

    :return aligned: Bio.Align.MultipleSeqAlignment
        Multiple sequence alignment object containing the MUSCLE alignment

    :return identity_ratio: float
        Identity ratio of the MSA. Calculated by counting up the "*" characters in the ClustalW output.

    :return similar_ratio: float
        Similarity ratio of the MSA. Calculated by counting up the ":" and the "*" characters in the ClustalW output.
        Not generally informative in the case of DNA sequence.
    """

    # Create and write a temporary fasta file with the sequences
    seq_file = tempfile.mkstemp(suffix=".fasta")
    with open(seq_file[0], mode="w") as seq_fh:
        SeqIO.write(list_of_seqrecords, seq_fh, "fasta")

    # Create a temporary file to hold MUSCLE output and pass the MUSCLE command to the shell.
    muscle_file = tempfile.mkstemp(suffix=".muscle")
    muscle_cmd = ["muscle", "-in", seq_file[1], "-out", muscle_file[1], "-clw"]
    subprocess.call(muscle_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    ident = 0
    similar = 0
    align_chars = ["*",":","."]

    with open(muscle_file[0], mode="rU") as muscle_fh:

        #Read MUSCLE alignment into a MultipleSeqAlign object

        try:
            aligned = AlignIO.read(muscle_fh, "clustal", alphabet=alphabet)
        except ValueError:
            print("Multiple alignments in file. Skipping.")
            aligned = None

        #Read through the MUSCLE alignment

        muscle_fh.seek(0)
        for line in muscle_fh:
            try:

                #Pass lines that don't have degree of conservation symbols
                if line.strip()[0] not in align_chars:
                    continue

                #Count up the identity and conserved residue symbols
                else:
                    ident += line.count("*")
                    similar += line.count(":")

            except IndexError:
                continue

    identity_ratio = ident / min(map(len,aligned))
    similar_ratio = (ident + similar) / min(map(len,aligned))

    # Delete the temporary files
    os.remove(seq_file[1])
    os.remove(muscle_file[1])

    return aligned, identity_ratio, similar_ratio
