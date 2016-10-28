import argparse
import os
import subprocess
import tempfile
from multiprocessing import Pool
from xml.etree.ElementTree import ParseError

from Bio import SearchIO
from Bio import SeqIO
from Bio.Alphabet import generic_protein
from Bio.Data.CodonTable import TranslationError
from Bio.SeqRecord import SeqRecord
from lib.align import muscle_wrapper
from lib.blast import blast_record_grab
from lib.identity import IdentityScreener
from lib.orf import ORF_Scanner


def main():
    """

    Shift_blast is a BLAST wrapper that takes a seed query and blasts it against a genomic BLAST database. It then
    identifies nearby ORFs and filters them (currently based principally on size), allowing for the identification of
    pathway elements in physically adjacent operons

    It takes either a nucleotide seed sequence that gets blastn or a protein sequence that gets tblastn against a
    nucleotide database

    """

    sh_parse = argparse.ArgumentParser(
        description="Shift-BLAST, retreiving hits from a seed sequence and finding nearby ORFs")

    sh_parse.add_argument("-f", "--file", dest="infile", help="Input sequence FILE", metavar="FILE", required=True)
    sh_parse.add_argument("-t", "--type", dest="intype", help="Input sequence file TYPE (default is FASTA)",
                          metavar="TYPE", default="fasta")
    sh_parse.add_argument("-o", "--out", dest="outfile", help="Output sequence FILE", metavar="FILE", default=None)
    sh_parse.add_argument("-d", "--db", dest="database", help="BLAST database to search", metavar="BLAST", default="nr")
    sh_parse.add_argument("-w", "--window", dest="window", help="Window SIZE around BLAST hits in first round",
                          metavar="SIZE", default=10000, type=int)
    sh_parse.add_argument("--offset", dest="offset", help="Number of ORFs to consider on each side of hit",
                          metavar="SIZE", default=2, type=int)
    sh_parse.add_argument("--finalaa", dest="finalaa",
                          help="Minimum SIZE of final candidate protein in amino acids to be considered",
                          metavar="SIZE", default=100, type=int)
    sh_parse.add_argument("--cpu", dest="cores", help="Number of CORES to use per thread",
                          metavar="CORES", default=4, type=int)
    sh_parse.add_argument("--patchk", dest="patchk", help="DATABASE to screen for patentability", metavar="DATABASE",
                          default=None)
    sh_parse.add_argument("--ident", dest="ident_db", help="DATABASE to screen for gene id & description",
                          metavar="DATABASE", default=None)

    sh_args = sh_parse.parse_args()

    shift_blast(sh_args.infile, sh_args.database, infile_type=sh_args.intype, window=sh_args.window,
                outfile_name=sh_args.outfile, cores=sh_args.cores, ident_db=sh_args.ident_db,
                patent_db=sh_args.patchk, final_aa=sh_args.finalaa, max_offset=sh_args.offset)


def shift_blast(infile_name, database, infile_type="fasta", window=10000, outfile_name=None,
                cores=4, min_aa=100, patent_db=None, final_aa=None, max_offset=2, ident_db=None):
    """
    shift_blast is the worker function which takes a seed sequence and blasts it against a nucleotide database, then
    finds ORFs near the blast hit. It then filters those hits and compares them to an additional BLAST protein database

    Required Arguments:

    :param infile_name: str
        Path to the input file
    :param database:    str
        Path to the BLAST nucleotide database

    Keyword Arguments:

    :param infile_type: str
        The input file type
    :param window:  int
        Window size in base pairs
    :param outfile_name:    str
        Output file path
    :param cores:   int
        Number of cores to use
    :param min_aa:  int
        Minimum size of the orf in amino acids to be considered
    :param patent_db:   str
        Path to the secondary protein database to search for identity scoring
    :param final_aa:    int
        Filter the final protein to size if
    :param max_offset:  int
        How many ORFs to consider (+/- from the seed sequence hit)
    :param ident_db:    str
        Path to a protein database to recover any associated protein-level data.
        This is SUPER SLOW right now

    :return:
    """

    # Set the hit filter for hsp. Fiddling with this may alter signal/noise
    hit_filter = lambda hsp: hsp.aln_span > 20 and hsp.bitscore > 40

    # Retreive the query sequence from the input file
    with open(infile_name, mode="rU") as in_fh:
        try:
            query_sequence = SeqIO.read(in_fh, format=infile_type)
        except ValueError:
            print("Only one query sequence at a time")
            exit(0)

    # Make a best guess at the input sequence DNA vs protein
    dna = 0
    seq = str(query_sequence).upper()

    for base in ["A", "T", "G", "C", "N"]:
        dna += seq.count(base)

    blast_out_file = tempfile.mkstemp(suffix=".blast.xml")

    # If this looks like protein, run tblastn
    if dna / len(seq) < 0.95:
        blast_cmd = ["tblastn", "-db", database, "-query", infile_name, "-outfmt", "5", "-out", blast_out_file[1],
                     "-num_threads", str(cores), "-max_target_seqs", str(2000)]
        query_sequence.alphabet = generic_protein
    # Otherwise run blastn
    else:
        blast_cmd = ["blastn", "-db", database, "-query", infile_name, "-outfmt", "5", "-out", blast_out_file[1],
                     "-num_threads", str(cores), "-task", "blastn", "-max_target_seqs", str(2000)]
        try:
            query_sequence.seq = query_sequence.seq.translate()
            query_sequence.alphabete = generic_protein
        except TranslationError:
            pass

    subprocess.call(blast_cmd)

    # Process the BLAST result with SearchIO
    with open(blast_out_file[0], mode="rU") as blast_fh:
        try:
            query = SearchIO.read(blast_fh, format='blast-xml')
        except ParseError:
            print(" ".join(blast_cmd) + "BLAST COMMAND FAILED... EXITING")
            exit(1)

    # Filter the results by HSP stats
    filter_query = query.hsp_filter(hit_filter)
    print("{} BLAST results [{} filtered]".format(len(filter_query), len(query) - len(filter_query)))

    # Open a process pool and create a generator for BLAST hits
    process_pool = Pool(processes=cores, maxtasksperchild=100)

    def _hit_generator(bq):
        hsp_num = 0
        for hit in bq.iterhits():
            for hsp in hit:
                hsp_num += 1
                hsp.hsp_number = hsp_num
                yield hsp

    # Multiprocess the BLAST hit processing through the BlastProcessor class
    blast_processor = BlastProcessor(database, query_sequence, window=window,
                                     orf_scanner=ORF_Scanner(min_aa_length=min_aa))
    blast_process_runner = process_pool.imap_unordered(blast_processor.blast_process, _hit_generator(filter_query))

    orfs = {}
    seen_orfs = {}

    # Take the BlastProcessor results and iterate through them
    for orf_tuple in blast_process_runner:
        if orf_tuple is not None:
            for orf_seq, hit_num, orf_num in orf_tuple:

                # Skip if the ORF is too far offset from the seed query
                if abs(orf_num) <= max_offset:

                    # Skip if the ORF doesnt meet the minimum length requirement
                    if final_aa is not None and final_aa * 3 > len(orf_seq):
                        continue

                    # Skip if the ORF has been seen before
                    try:
                        if seen_orfs[orf_seq.id]:
                            continue
                    except KeyError:
                        seen_orfs[orf_seq.id] = True

                    # Translate ORF to protein sequence
                    try:
                        orf_seq.seq = orf_seq.seq.translate(cds=True)
                        orf_seq.seq.alphabet = generic_protein
                    except TranslationError:
                        print("Translation error [{}]: {}...{}".format(orf_seq.id, orf_seq[:10], orf_seq[:-10]))
                        continue

                    # Throw out garbage alignment ORFs that are full of Xs
                    if str(orf_seq.seq).upper().count("X") / len(orf_seq) > 0.1:
                        print("Translation QC [{}]: {}...{}".format(orf_seq.id, orf_seq[:10], orf_seq[:-10]))
                        continue

                    # Save the ORF sequence
                    try:
                        orfs[orf_num].append(orf_seq)
                    except KeyError:
                        orfs[orf_num] = [orf_seq]

    print("BLAST Analysis: {} ORFs detected in {} offsets".format(len(seen_orfs), len(orfs)))

    # Create a generator with all the potential hits identified by the BLAST result processor
    def _orf_generator(orf_dict):
        for orf_idx in orf_dict:
            for orf in orf_dict[orf_idx]:
                yield orf

    # Process the potential hits with the HitProcess class to filter and determine identity
    heavy_hitter = HitProcess(database, reference=query_sequence, min_aa=final_aa, patentdb=patent_db, id_db=ident_db)
    analysis_runner = process_pool.imap_unordered(heavy_hitter.process_hits,
                                                  _orf_generator(orfs))

    with open(outfile_name, mode="w+") as out_fh:
        print("{}\t{}\t{}\t{}\t{}".format("Hit ID",
                                          "Hit Description",
                                          "Hit Location",
                                          "Patent DB Identity",
                                          "Hit Sequence"),
              file=out_fh)

        # Iterate through the retained hits and print them to the output file
        for hit_seq, hit_ident, hit_simil, patent_ident in analysis_runner:
            if hit_ident is not None:
                print("{}\t{}\t{}\t{}\t{}".format(hit_seq.id,
                                                  hit_seq.description,
                                                  hit_seq.name,
                                                  patent_ident,
                                                  str(hit_seq.seq)),
                      file=out_fh)


class BlastProcessor:
    """
    The BlastProcessor class takes a BLAST hsp (high scoring pair), retrieves the genomic nucleotide data from nearby,
    identifies the ORFs in that opened genomic window, and returns them
    """

    def __init__(self, db, query, window=10000, orf_scanner=None):
        """
        Instantiation of the BlastProcessor class

        Required Arguments:

        :param db:  str
            Path to the primary genomic BLAST database
        :param query:   Bio.SeqRecord
            The seed query that's been translated to protein

        Keyword Arguments:

        :param window:  int
            Integer size of the nucleotide window to open from the genomic hit
        :param orf_scanner: ORF_Scanner object
            ORF_scanner class passed in with the various options set at instantiation. Use defaults if not passed.
        """

        self.query = query
        self.query_nuclen = 3 * len(self.query)
        self.window = window
        self.db = db

        if orf_scanner is None:
            self.orf_scanner = ORF_Scanner()
        else:
            self.orf_scanner = orf_scanner

    def blast_process(self, hsp):
        """
        blast_process takes a BLAST search high-scoring pair, opens a window from the genomic DNA retrieved from the
        BLAST database, and passes it to the ORF scanner

        Required Arguement:

        :param hsp: Bio.Blast.Record.HSP
            HSP object from the BLAST query results

        Returns:

        :return: [(Bio.Seq), int, int, int)]
            Returns a list of ORF tuples that consist of a Bio.Seq object, a start integer, a stop integer, and
        """

        # Store the hsp in the class object
        self.hsp = hsp
        self.hsp_num = hsp.hsp_number

        # Open a window from the genomic nucleotide data
        try:
            window_seq, hit_start, hit_end, hit_frame = self._open_window()
        except ValueError:
            return None

        self.window_len = len(window_seq)
        self.hit_start = hit_start
        self.hit_end = hit_end

        # Pass the genomic nucleotide window to the ORF scanner
        orf_list = self.orf_scanner.orf_scan(window_seq)

        # Process the identified ORFs
        process_orfs = self.orf_process(orf_list)

        return process_orfs

    def _open_window(self):
        """
        _open_window takes the class hsp object, gets the hit_id from it, and use that ID to get the nucleotide window
        from the BLAST database. Uses the HSP set in the class object - don't call this method yourself

        :return return_seq: Bio.SeqRecord
            The sequence record of the opened window of genomic DNA
        :return hit_start:  int
            The integer start of the window in the genomic DNA
        :return hit_end:    int
            The integer end of the window in the genomic DNA
        :return hit_frame:  int
            The strand of the hit returned within the genomic DNA record (1 for sense, -1 for sense)
        """

        # Get the genomic DNA sequence
        record_seq = blast_record_grab(self.hsp.hit_id, self.db)
        if record_seq is None:
            return None

        # Set values for the start and stop of the window
        hit_frame = int(self.hsp.hit_frame)
        hit_start = self.hsp.hit_start - int(self.window / 2)
        hit_end = self.hsp.hit_end + int(self.window / 2)

        # Make sure that the start and stop is within the genomic DNA record
        if hit_start < 0:
            hit_start = 0
        elif hit_end > len(record_seq.seq):
            hit_end = len(record_seq.seq)

        # Slice the sequence using the start and end
        return_seq = record_seq.seq[hit_start:hit_end]

        # Reverse and complement to get the genomic window sense to the hit
        if hit_frame > 0:
            pass
        elif hit_frame < 0:
            return_seq = return_seq.reverse_complement()
        else:
            raise AttributeError("Hit frame is 0")

        return return_seq, hit_start, hit_end, hit_frame

    def orf_process(self, orf_list):
        """
        orf_process takes the list of total ORFs from the ORF scanner and filters based on maximum offset to the
        seed sequence

        Required Arguments:

        :param orf_list: Seq, int, int, int

        Returns:

        :return:
        """

        find_query = (self.window_len, 0, 0)
        left_sash = self.window / 2
        right_sash = self.window_len - left_sash

        # Iterate through the list of ORFs and try to figure out which one is most likely the query that seeded
        for orf, start, stop, count in orf_list:

            # Find the center of the specific ORF
            orf_midpoint = (start + stop) / 2

            # Figure out how far away from the center of the window the center of the ORF is
            offset = max([left_sash - orf_midpoint, orf_midpoint - right_sash, 0])

            # If this ORF is closest to the center of the window, save the offset number
            if offset < find_query[0]:
                find_query = (offset, count, len(orf))
            elif offset == find_query[0] and abs(len(orf) - self.query_nuclen) < abs(find_query[0] - self.query_nuclen):
                find_query = (offset, count, len(orf))
            else:
                pass

        return_orfs = []

        # Replace the 1...N orf count ID with an offset centered around 0, which should be the query, and return
        for orf, start, stop, count in orf_list:
            new_id = "{}[Window={}:{}][ORF={}:{}]".format(self.hsp.hit_id, self.hit_start, self.hit_end, start, stop)
            return_orfs.append((SeqRecord(orf, id=self.hsp.hit_id, name=new_id), self.hsp_num, count - find_query[1]))
        return return_orfs


class HitProcess:
    """
    The HitProcess class takes the hits identified by BLAST and scans a second database for identity, and a third
    database for protein ID & description
    """

    def __init__(self, db, reference=None, min_aa=100, patentdb=None, alphabet=generic_protein, id_db=None, cores=2):
        """
        Instantiates the HitProcess class

        Required Arguments:

        :param db: str
            Path to the primary nucleotide genomic BLAST database

        Keyword Arguments:

        :param reference: Bio.SeqRecord
            The query sequence to reference against
        :param min_aa: int
            Minimum protein length, smaller proteins will be discarded
        :param patentdb: str
            Path to a protein BLAST database to screen for identity %
        :param alphabet: Bio.Alphabet
            Only protein alphabets work. Never implemented nucleotide
        :param id_db: str
            Path to a protein BLAST database to recover protein ID & description from
        :param cores: int
            Number of cores to use per HitProcess class instance
        """

        self.database = db
        self.reference = reference
        self.min_aa = min_aa
        self.patent_db = patentdb
        self.alphabet = alphabet
        self.id_db = id_db
        self.cores = cores

    def process_hits(self, hit_sequence):
        """
        process_hits is the worker function for the HitProcess class. It takes a hit sequence and BLASTs it at the
        database

        :param hit_sequence: Bio.SeqRecord
            The protein sequence to screen against percent identity and protein description databases
        :return: Bio.SeqRecord, float, float, float
            Returns the hit_sequence, the identity ratio to the query, the similarity ratio to the query, and the
            identity ratio to the secondary screening patent_db
        """

        # Screen protein database patent_db to find the highest % identity to the candidate protein
        patent_identity = 0
        if self.patent_db is not None:
            papers_please = IdentityScreener(self.patent_db)
            patent_identity = papers_please.ident_query(hit_sequence)

        _, ident, simil = muscle_wrapper([self.reference, hit_sequence], alphabet=self.alphabet)

        # Screen protein database id_db to find the description of the closest protein
        if self.id_db is not None:
            found_id = self.get_id(hit_sequence)
            if found_id is not None:
                hit_sequence.description = found_id

        # If there's no id_db set, just set the description to be the seqid from the genomic database hit
        else:
            accession = str(hit_sequence.id).strip().split(sep="|")[-2]
            blastdbcmd_cmd = ["blastdbcmd", "-db", self.database, "-dbtype", "nucl", "-entry", accession, "-outfmt",
                              "\"%a %t\""]
            blast_call = subprocess.Popen(blastdbcmd_cmd, stdout=subprocess.PIPE)
            hit_sequence.description = blast_call.communicate()[0].decode().strip()

        return hit_sequence, ident, simil, patent_identity

    def get_id(self, hit):
        """
        get_id takes a hit SeqRecord and blasts it against a protein database

        Required Arguments:

        :param hit: SeqRecord
            SeqRecord of an unidentified protein sequence retreived from a genomic nucleotide database and translated

        Returns:

        :return top_hit_id: str
            The blast sequence description of the highest-scoring hit from the BLAST search

        """

        # Write the candidate protein sequence into a fasta tempfile
        id_file = tempfile.mkstemp(suffix=".fasta")
        with open(id_file[0], mode="w+") as id_fh:
            SeqIO.write(hit, id_fh, format="fasta")

        # Blast it
        blast_id_file = tempfile.mkstemp(suffix=".blast.xml")
        blast_cmd = ["blastp", "-db", self.id_db, "-query", id_file[1], "-outfmt", "5", "-out", blast_id_file[1],
                     "-num_threads", str(self.cores), "-task", "blastp", "-max_target_seqs", str(2000)]

        subprocess.call(blast_cmd)
        os.remove(id_file[1])

        # Read in the BLAST results
        blast_fh = open(blast_id_file[0], mode="rU")

        try:
            query = SearchIO.read(blast_fh, format='blast-xml')
        except ParseError:
            return None
        finally:
            blast_fh.close()
            os.remove(blast_id_file[1])

        # Return the highest scoring BLAST hit as a hit_id string
        top_hit_id = query.hit_map(lambda hit: hit[:-1]).hit_keys()[0]
        return top_hit_id


if __name__ == '__main__':
    main()
