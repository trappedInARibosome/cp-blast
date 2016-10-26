from Bio.Alphabet import generic_protein
from lib.align import muscle_wrapper
from lib.blast import blast_query
from xml.etree.ElementTree import ParseError

class IdentityScreener:
    """
    IdentityScreener is a class that will calculate maximum identity of a SeqRecord object compared to a protein
    database, using BLAST and MUSCLE

    Implemented as a class because I find it easier to use multiprocessing with a class
    """

    def __init__(self, database, cores=2):
        """
        Instantiates the IdentityScreener class object

        :param database: str
            Path to a BLAST database to search against

        :param cores: int
            Number of cores to use per BLAST thread

        """
        self.database = database
        self.cores = cores

    def ident_screen(self, input_tuple):
        """
        ident_screen is a wrapper for the ident_query function that's a helper to unpack a tuple so I don't have to use
        starmap when I'm working with multiprocess.

        It just passes most of the input data back anyway. Don't worry bout it.

        :param input_tuple:

        """

        seq_record, offset, bin_num = input_tuple
        max_ident = self.ident_query(seq_record)

        return seq_record, max_ident, offset, bin_num


    def ident_query(self, seq_record):
        """
        ident_query takes a SeqRecord and determines the maximum identity percentage against a blast db

        Required Arguments:

        :param seq_record: Bio.SeqRecord
            Sequence record to use as a database query
        :param comp_database: str
            Database path to search

        Keyword Arguments:

        :param cores: int
            Number of cores to use for blast

        Returns:

        :return: float
            Ratio of identical amino acids to the total query sequence length
        """

        try:
            bq = blast_query(seq_record, self.database, cores=self.cores)
        except ParseError:
            return 0

        max_ident = 0
        for hit in bq.iterhits():
            ident_percent = calculate_identity(hit)
            if ident_percent > max_ident:
                max_ident = ident_percent

        return max_ident


def calculate_identity(hit):
    """
    Calculates identity of a query and hit pair by aligning with MUSCLE

    Required Arguments:

    :param hit: Hit
        Hit object containing a blast query hit

    Returns:

    :return: float
        Ratio of identical amino acids to the total query sequence length

    """

    for hsp in hit:
        query_seqrecord = hsp.query
        hit_seqrecord = hsp.hit

        _, identity_ratio, _ = muscle_wrapper([query_seqrecord, hit_seqrecord], alphabet=generic_protein)

        return identity_ratio