from Bio.Alphabet import generic_dna
from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord


class ORF_Scanner:

    """
    The ORF Scanner class takes a nucleotide sequence and finds Open Reading Frames (ORFs) according to a handful of
    simple rules. Done as a class to facilitate multiprocessing.
    """

    start = ["ATG"]
    stops = ["TAA", "TAG", "TGA"]

    def __init__(self, min_aa_length=100, max_overlap=None):
        """
        Instantiation method for the ORF_Scanner class

        Keyword Arguments:

        :param min_aa_length: int
            Defines the minimum required length in amino acids to constitute a potential ORF
        :param max_overlap: int
            Defines the largest range that ORFs are allowed to overlap on the ends
        """

        self.min_aa_length = min_aa_length
        self.max_overlap = max_overlap

    def orf_scan(self, window_seq):
        """
        The orf_scan method takes a window of nucleotides and finds ORFs, returning them as a list of ORFs that consist
        of a Seq object, a start position, a stop position, and an offset number

        Required Argument:

        :param window_seq: Bio.Seq / Bio.SeqRecord
            Takes a Bio.Seq object (nucleotide) or Bio.SeqRecord

        Returns:

        :return: [(Bio.Seq, int, int, int)]
            Returns a list of ORF records, which have the orf nucleotide sequence as a Seq object, an integer start
            position (relative to the window start), an integer stop position (relative to the window end), and an
            offset integer
        """

        if isinstance(window_seq, SeqRecord):
            window_seq = window_seq.seq

        try:
            window_seq = Seq(window_seq, generic_dna)
        except TypeError:
            pass

        open_orf = [0, 0, 0]
        possible_starts = [[], [], []]
        orfs = []

        # Feed the sequence into a generator which walks one character at a time
        for codon, pos, frame in self._slider(window_seq):

            # If a start codon has been reached
            if codon in self.start:

                # Append the start location to all of the possible start positions in this frame
                possible_starts[frame].append(pos)

                # Mark the start of a potential ORF if there is nothing else in-frame
                if open_orf[frame] == 0:
                    open_orf[frame] = pos

            # If a stop codon has been reached and there's a potential ORF in-frame
            if codon in self.stops and open_orf[frame] != 0:

                # Determine if the longest possible in-frame ORF can meet the minimum length requirement
                if int(pos - open_orf[frame]) >= (self.min_aa_length * 3):

                    # Add the ORF start and stop locations to the list
                    orfs.append((open_orf[frame], pos + 3))

                    # Go through the list of start positions for each frame and keep the longest one that meets the
                    # overlap criteria
                    for frame_check in range(2):
                        for maybe_start in possible_starts[frame_check]:
                            if self.max_overlap is None or pos - maybe_start <= self.max_overlap:
                                open_orf[frame_check] = maybe_start
                                break

                # Close the ORF in this frame
                possible_starts[frame] = []
                open_orf[frame] = 0

        # Return sequences, start position, stop position, and the orf number
        seq_list = []
        for i, (start, stop) in enumerate(orfs):
            seq_list.append((window_seq[start:stop], start, stop, i))
        return seq_list

    def _slider(self, seq, step=3):
        """
        _slider is an internal generator which takes the window sequence and slides a window forward

        Required Arguments:

        :param seq: A Seq object or string containing the nucleotide sequence

        :param step:
        :return:
        """

        frame = -1
        string = str(seq)

        for i_pos in range(len(seq) - step):

            if frame > 1:
                frame = 0
            else:
                frame += 1

            yield string[i_pos:i_pos + step], i_pos, frame
