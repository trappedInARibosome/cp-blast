import io
import subprocess

from Bio import AlignIO
from Bio import SeqIO
from Bio.Alphabet import generic_dna


def muscle_wrapper(list_of_seqrecords, alphabet=generic_dna):
    """
    muscle_wrapper uses the command line MUSCLE application, with the default options for alignment, to align a list of
    SeqRecord objects.

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
    seq_file = io.StringIO()
    SeqIO.write(list_of_seqrecords, seq_file, "fasta")

    # Create a temporary file to hold MUSCLE output and pass the MUSCLE command to the shell.
    muscle_cmd = ["muscle", "-clw"]
    muscle_proc = subprocess.Popen(muscle_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                                   stdin=subprocess.PIPE, universal_newlines=True)
    muscle_output = muscle_proc.communicate(input=seq_file.getvalue())[0]

    ident = 0
    similar = 0
    align_chars = ["*", ":", "."]

    # Get the MultipleSeqAlign object through StringIO
    try:
        aligned = AlignIO.read(io.StringIO(muscle_output), "clustal", alphabet=alphabet)
    except ValueError:
        print("Multiple alignments in file. Skipping.")
        aligned = None

    for line in muscle_output.splitlines():
        try:
            # Pass lines that don't have degree of conservation symbols
            if line.strip()[0] not in align_chars:
                continue

            # Count up the identity and conserved residue symbols
            else:
                ident += line.count("*")
                similar += line.count(":")

        except IndexError:
            continue

    identity_ratio = ident / min(map(len, aligned))
    similar_ratio = (ident + similar) / min(map(len, aligned))
    return aligned, identity_ratio, similar_ratio
