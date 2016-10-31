import argparse
import multiprocessing
import os
import subprocess
import tempfile
from xml.etree.ElementTree import ParseError

from Bio import SearchIO
from Bio import SeqIO
from Bio.Alphabet import generic_protein
from Bio.Data.CodonTable import TranslationError
from lib.align import muscle_wrapper
from lib.blast import blast_record_grab
from lib.identity import IdentityScreener


def main():
    """
    The Direct Blast function takes an input query file and tblastn, blastp or psi-blasts it against a target protein
    database. It then filters results based on size, aligns the hit against the query, and calculates percent identity
    and similarity to the query.

    It also, if set, will check the hits against a second database and calculate the percent identity to that database
    Useful for determining similarity to existing IP
    """

    sh_parse = argparse.ArgumentParser(description="BLAST a protein database")

    sh_parse.add_argument("-f", "--file", dest="infile", help="Input sequence FILE", metavar="FILE", required=True)
    sh_parse.add_argument("-t", "--type", dest="intype", help="Input sequence file TYPE (default is FASTA)",
                          metavar="TYPE", default="fasta")
    sh_parse.add_argument("-o", "--out", dest="outfile", help="Output sequence FILE", metavar="FILE", default=None)
    sh_parse.add_argument("-d", "--db", dest="database", help="BLAST database to search", metavar="BLAST", default="nr")
    sh_parse.add_argument("--minaa", dest="minaa",
                          help="Minimum SIZE of final candidate protein in amino acids to be considered",
                          metavar="SIZE", default=0, type=int)
    sh_parse.add_argument("--cpu", dest="cores", help="Number of CORES to use per thread",
                          metavar="CORES", default=4, type=int)
    sh_parse.add_argument("--patchk", dest="patchk", help="DATABASE to screen for patentability", metavar="DATABASE",
                          default=None)
    sh_parse.add_argument("--psi", dest="psi", help="Run the PSI-BLAST algorithm", default=False, const=True,
                          action="store_const")
    sh_parse.add_argument("--bitscore", dest="bitscore", help="Minimum required HSP BITSCORE", default=50, type=int,
                          metavar="BITSCORE")

    sh_args = sh_parse.parse_args()

    direct_blast(sh_args.infile, sh_args.database, outfile=sh_args.outfile, in_type=sh_args.intype, cores=sh_args.cores,
                 patent_db=sh_args.patchk, min_aa_size=sh_args.minaa, psi_blast=sh_args.psi,
                 min_bitscore=sh_args.bitscore)


def direct_blast(infile, database, outfile=None, in_type="fasta", cores=1, patent_db=None, min_aa_size=0,
                 psi_blast=False, min_bitscore=50):
    """
    direct_blast is the worker function which takes an input filename and blasts it against a target database, then
    calculates identity against the query and against some other database

    Positional Arguments:

    :param infile:  str
        Input file path
    :param database:    str
        Input database path. Goes directly into command line blast, requires the database name as well.

    Keyword Arguments:

    :param outfile: str
        Output file path
    :param in_type: str
        Input file type. Takes anything SeqIO can take.
    :param cores: int
        Passed to blast command line as num_threads, and number of pool processes to spawn to analyze hits
    :param patent_db: str
        Screening database path. Goes directly into command line blast, requires the database name as well.
    :param min_aa_size: int
        Minimum size of the protein in amino acids
    :param psi_blast: bool
        Flag to use psi-blast instead of blastp
    :param min_bitscore: int
        The minimum required bitscore of a BLAST high scoring pair. HSPs with lower scores will be filtered

    """

    # Use a generic name for the outfile if it hasn't been explicitly set
    if outfile is None:
        outfile = infile + ".db.out"

    # Make sure the databases exist
    # Easier then catching it way down the line
    if not (os.path.isfile(database + ".phd") or os.path.isfile(database + ".00.phd")):
        raise FileNotFoundError("BLAST database {} not located".format(database))
    if patent_db is not None and not (os.path.isfile(patent_db + ".phd") or os.path.isfile(patent_db + ".00.phd")):
        raise FileNotFoundError("BLAST database {} not located".format(database))

    # Set the BLAST hsp minimums
    hsp_filter = lambda hsp: hsp.aln_span > 50 and hsp.bitscore > min_bitscore

    # Open and read in the query file as SeqRecords
    with open(outfile, mode="w") as out_fh, open(infile, mode="rU") as in_fh:
        for query_sequence in SeqIO.parse(in_fh, format=in_type):

            # Count the DNA bases to do a lazy job determining if this is a protein or a DNA sequence
            dna = 0
            for base in ["A", "T", "G", "C", "N"]:
                dna += str(query_sequence.seq).upper().count(base)

            blast_out_file = tempfile.mkstemp(suffix=".blast.xml")

            # Decide which blast command line arguments to use for the query sequence
            # Also set sequence alphabet

            # If the sequence looks like protein
            if dna / len(query_sequence) < 0.95:
                query_sequence.seq.alphabet = generic_protein

            # If it's probably DNA
            else:
                # Translate the sequence if it looks like DNA and use the translated sequence for downstream
                try:
                    query_sequence.seq = query_sequence.seq.translate(cds=True)
                    query_sequence.seq.alphabet = generic_protein
                except TranslationError as trans_err:
                    print("Input Sequence {} Not a CDS: {}".format(infile, trans_err.args))
                    try:
                        query_sequence.seq = query_sequence.seq.translate()
                        query_sequence.seq.alphabet = generic_protein
                    except TranslationError:
                        print("Translation Error")
                        exit(0)

                # Write the translated protein to a file to use as the query
                blast_query_temp = tempfile.mkstemp(suffix=".fasta")
                with open(blast_query_temp[0], mode="w") as blast_temp_fh:
                    SeqIO.write(query_sequence, blast_temp_fh, format="fasta")
                infile = blast_query_temp[1]

            # Run psiblast if the psi flag is set otherwise blastp
            if psi_blast:
                blast_cmd = ["psiblast", "-db", database, "-query", infile, "-outfmt", "5", "-out",
                             blast_out_file[1],
                             "-num_threads", str(cores), "-max_target_seqs", str(5000)]
            else:
                blast_cmd = ["blastp", "-db", database, "-query", infile, "-outfmt", "5", "-out", blast_out_file[1],
                             "-num_threads", str(cores), "-task", "blastp", "-max_target_seqs", str(5000)]

            subprocess.call(blast_cmd)

            # Read in the blast output file as a QueryResult
            with open(blast_out_file[0], mode="rU") as blast_fh:
                try:
                    query = SearchIO.read(blast_fh, format='blast-xml')
                except ParseError:
                    print("BLAST Command Failed")
                    exit(0)

            # Preprocess the BLAST query result with the hsp filter object
            filter_query = query.hsp_filter(hsp_filter)
            print("{} BLAST results [{} filtered]".format(len(filter_query), len(query) - len(filter_query)))

            # Pass the control arguments to the hit processor and then multiprocess hits through mp.Pool
            heavy_hitter = HitProcess(database, query_sequence, min_aa=min_aa_size, patentdb=patent_db)
            blast_process_runner = multiprocessing.Pool(processes=cores).imap_unordered(heavy_hitter.process_hits,
                                                                                        (hit for hit in filter_query))

            print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format("Query ID",
                                                      "Hit ID",
                                                      "Hit Description",
                                                      "Hit Identity",
                                                      "Hit Similarity",
                                                      "Patent DB Identity",
                                                      "Hit Sequence"),
                  file=out_fh)

            # Iterate through the processing results and print them
            for hit_id, hit_seq, hit_ident, hit_simil, patent_ident in blast_process_runner:
                if hit_ident is not None:
                    print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(query_sequence.id,
                                                              hit_id,
                                                              hit_seq.description,
                                                              hit_ident,
                                                              hit_simil,
                                                              patent_ident,
                                                              str(hit_seq.seq)),
                          file=out_fh)

            os.remove(blast_out_file[1])


class HitProcess:
    """
    HitProcess is the analysis module for BLAST hits. It's written as a class to facilitate multiprocessing.
    """

    def __init__(self, db, query, min_aa=None, patentdb=None, alphabet=generic_protein):
        """
        Instantiation method for the HitProcess class.

        Required Arguments:

        :param db: str
            Path to the primary database that the query was searched against

        :param query: SeqRecord
            SeqRecord containing the query sequence

        Keyword Arguments:

        :param min_aa: int
            If not none, any hits smaller then this will be discarded

        :param patentdb: str
            Path to the secondary database to search with primary hits and determine identity percent

        :param alphabet: Bio.Alphabet
            Alphabet to use. Pretty sure that this class doesnt work with DNA but I didn't hard-code protein
        """

        self.database = db
        self.query = query
        self.min_aa = min_aa
        self.patentdb = patentdb
        self.alphabet = alphabet

    def process_hits(self, hit):
        """
        Process_hits retrieves the hit sequence from the blast database with blastdbcmd and calculates the identity and
        similarity. If set in the instantiation method, it also searches a second database with the hit sequence and
        calculates identity from that.

        :param hit: Bio.SearchIO.hit
            BLAST search hit object to be processed

        :return hit.id: str
            Database ID of the hit object that was passed in

        :return db_seq: Bio.SeqRecord
            SeqRecord object that contains the entire hit from the BLAST database

        :return ident:  float
            Identity ratio between the query passed in at class instantiation and the hit

        :return simil:  float
            Similarity ratio between the query passed in at class instantiation and the hit

        :return patent_identity:    float
            Identity ratio between the primary database hit, and the highest result from searching a secondary database
        """

        db_seq = blast_record_grab(hit.id, self.database, alphabet=self.alphabet)

        # If min_aa is set, return Nones for the hit sequence ratios if it's too small
        if self.min_aa is not None and len(db_seq) < self.min_aa:
            return hit, db_seq, None, None, None

        # If patentdb is set, search the secondary database with the primary hit and get the maximum identity
        # Otherwise return 0
        patent_identity = 0
        if self.patentdb is not None:
            papers_please = IdentityScreener(self.patentdb)
            patent_identity = papers_please.ident_query(db_seq)

        _, ident, simil = muscle_wrapper([self.query, db_seq], alphabet=self.alphabet)

        return hit.id, db_seq, ident, simil, patent_identity


if __name__ == '__main__':
    main()
