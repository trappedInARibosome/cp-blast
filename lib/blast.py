import os
import io
import subprocess
import tempfile
import time
from xml.etree.ElementTree import ParseError

from Bio import SearchIO
from Bio import SeqIO
from Bio.Alphabet import ProteinAlphabet
from Bio.Alphabet import generic_dna
from Bio.Alphabet.IUPAC import ExtendedIUPACProtein
from Bio.Blast.Applications import NcbiblastnCommandline, NcbiblastpCommandline
from Bio.Seq import Seq

def blast_query(sequence, database, delete=True, cores=2, max_target_seqs=500):
    """
    blast_query is a wrapper for NCBI command line blastn or blastp

    Required Arguments:

    :param sequence: Bio.Seq
        A Seq object which contains the query sequence. Will also take strings containing sequence data
    :param database: str
        A string that identifies the path of the blast database to search

    Keyword Arguments:

    :param delete: bool
        If set to TRUE, will clean up the temporary files generated by the script.


    Returns:

    :return query: Bio.SearchIO
        A SearchIO object created by reading in the blast output XML file

    """

    if isinstance(sequence, list):
        raise TypeError("blast_query takes single sequences, not lists of sequences")

    if isinstance(sequence, str):
        sequence = Seq(sequence, generic_dna)

    # Make a temp file and flush the sequence to it

    temp_file = tempfile.mkstemp()
    temp_fh = open(temp_file[0], mode="w+")
    SeqIO.write(sequence, temp_fh, "fasta")
    temp_fh.close()

    out_name = temp_file[1] + ".blast"

    if is_protein(seq=sequence):
        blast_command = NcbiblastpCommandline(query=temp_file[1], db=database, out=out_name, outfmt="5", task="blastp",
                                              num_threads=cores, max_target_seqs=max_target_seqs)
    else:
        blast_command = NcbiblastnCommandline(query=temp_file[1], db=database, out=out_name, outfmt="5", task="blastn",
                                              num_threads=cores, max_target_seqs=max_target_seqs)

    start_time = time.time()
    subprocess.call(str(blast_command), stdout=subprocess.DEVNULL, shell=True)

    run_time = int(100 * (time.time() - start_time)) / 100
    print("BLAST run completed on database {} in {}s".format(database, run_time))

    with open(out_name, mode="r") as blast_fh:
        try:
            query = SearchIO.read(blast_fh, format='blast-xml')
        except ParseError:
            print(str(blast_command))
            os.remove(out_name)
            raise

    if delete:
        os.remove(out_name)
        os.remove(temp_file[1])

    return query


def blast_record_grab(accession, database, alphabet=generic_dna):
    """
    blast_record_grab uses the command line blastdbcmd to obtain sequence data from the blast database

    Required Arguments:

    :param accession:   str
        String to use to search through a blast database for a record
    :param database:    str
        Path to the BLAST database

    Returns:

    :return db_seq: Bio.SeqRecord
        SeqRecord object created from the BLAST database entry
    """

    # Pull the sequence of the hit out of the primary BLAST database using blastdbcmd and put it into a string
    blastdbcmd_cmd = ["blastdbcmd", "-db",database, "-entry", accession]
    blast_call = subprocess.Popen(blastdbcmd_cmd, stdout=subprocess.PIPE)
    blast_fasta = blast_call.communicate()[0].decode().strip()

    # Turn the string with the hit sequence into a SeqRecord object
    try:
        db_seq = SeqIO.read(io.StringIO(blast_fasta), format="fasta", alphabet=alphabet)
    except ValueError:
        return None

    return db_seq

def is_protein(seq=None, alphabet=None):
    """
    Is_protein returns true if the sequence or alphabet correspond to a protein. Necessary because there are a bunch
    of protein alphabets in BioPython. Have to pass in one or the other or a ValueError is raised

    Keyword Arguements:

    :param seq: Seq
        Sequence object to be checked
    :param alphabet: Alphabet
        Alphabet object to be checked

    Returns:

    :return: bool
        Returns True if the sequence or alphabet has a protein alphabet

    """

    if seq is None and alphabet is None:
        raise ValueError

    if seq is not None:
        try:
            alphabet = seq.alphabet
        except AttributeError:
            pass

        try:
            alphabet = seq.seq.alphabet
        except AttributeError:
            pass

    if isinstance(alphabet, ExtendedIUPACProtein):
        return True
    elif isinstance(alphabet, ProteinAlphabet):
        return True
    else:
        return False