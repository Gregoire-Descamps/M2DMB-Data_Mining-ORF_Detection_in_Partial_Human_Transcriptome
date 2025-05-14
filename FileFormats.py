from __future__ import annotations
import io
import os
from typing import Union
import datetime
import logging
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import ORFinfder


# logging config
logging.basicConfig(level=logging.DEBUG)

class GffFile:
    """
        A class representing an GFF File

        Args:

         path: The output file path

        Methods:

        - set_header: -> Set the header information
        -  add_entry -> append a new line to the
        - write_orf_file -> write a gff file from a list of orf objects
        """

    def __init__(self, path):
        self.path = path
        self.__header: bool = False

        if os.path.exists(path):
            logging.warning(f'File {path} already exists and have been overwritten by GffFile object!')

        file = open(path, "w")
        file.close()

    def __open_file(self):
        return open(self.path, "a")

    def set_header(self, desc: str = "", provider: str = "", contact: str = "",
                   date: Union[str, datetime.date] = "Today"):
        """
        A method to overwrite the file and generate a header for the gff file

        !!! Using this method will erase all the file content and replace the header !!!

        Args:
            desc: the description string
            provider: the provider name
            contact: contact info string
            date: The date of generation
        """
        if isinstance(date, str) and date == "Today":
            date = datetime.date.today()
        elif not (isinstance(date, datetime.date)):
            raise TypeError('date provided should be a datetime.date type or "Today"')

        with open(self.path, "w") as file:
            logging.warning(f'File {file.name} have been completely overwritten to modify the header!')
            file.write(f'#gff-version 3\n'
                       f"#description: {desc}\n"
                       f"#provider: {provider}\n"
                       f"#contact: {contact}\n"
                       f"#format: gff3\n"
                       f"#date: {date}\n")
            file.close()
        self.__header = True

    def add_entry(self,
                  seqid: str,
                  source: str,
                  seq_type: str,
                  start: int,
                  end: int,
                  score: str = ".",
                  strand: str = ".",
                  phase: Union[str, int] = ".",
                  attributes: Union[dict, str] = ".",
                  open_file: io.TextIOWrapper = None) -> object:
        """
            Append a new entry at the end of the file

            :param source: The source sequence id
            :type source: str
            :param seq_type:  type of sequence (gff categories)
            :type seq_type: str
            :param start: start position of the sequence
            :type start: int
            :param end: end position of the sequence
            :type end: int
            :param score: Sequence score
            :type score: str
            :param strand:  strand of the sequence
            :type strand: str
            :param phase:  phase of the sequence
            :type phase: str
            :param attributes: Dictionary of sequence attributes (metadata)
            :type attributes: dict

            seqid (object):

        """
        # converting dictionary into attribute string
        if isinstance(attributes, dict):
            new_attributes = ""
            for key, value in attributes.items():
                new_attributes += f';{key}={value}'

            attributes = new_attributes[1:]
        if not open_file:
            with self.__open_file() as file:
                file.write(f'{seqid}\t{source}\t{seq_type}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{attributes}\n')

            file.close()
        else:
            open_file.write(f'{seqid}\t{source}\t{seq_type}\t{start}\t{end}\t{score}\t\t{strand}\t{phase}\t{attributes}\n')

    def write_orf_file(self,
                       orfEntries: Union[list[ORFinfder.ORF], ORFinfder.ORF] = None) -> object:
        """
        A method to write a gff file from a list of ORF objects

        Args:
            orfEntries: A list of ORF object to add in the gff file

        Returns:
            object:
        """

        if isinstance(orfEntries, ORFinfder.ORF):
            orfEntries= [orfEntries]

        self.set_header(desc="Potential CDS entries extracted from human transcriptome assembly",
                              provider="Gregoire Descamps",
                              contact="gregÔire.desc@mps@thisisnotarealemail.com",
                              date="Today")
        with self.__open_file() as file:
            for orf in orfEntries:
                self.add_entry(*orf.to_gff_tuple(), open_file = file)
        file.close()



class FastaFile:
    """
        A class representing a Fasta File

        Args:

         path: The output file path

        Methods:

        -  add_entry -> append a new line to the
        - write_orf_file -> write a gff file from a list of orf objects
        """

    def __init__(self, path):
        self.path = path


    def __open_file(self):
        return open(self.path, "a")

    def add_entry(self,
                  seqid: str,
                  desc: str,
                  seq: str):

        """
            Append a new entry at the end of the file

            :param seqid: The sequence id
            :param desc:  Entry description
            :param seq: The sequence
            """

        sequence = SeqRecord(seq,id = seqid, description= desc)
        SeqIO.write(sequence, self.__open_file(), 'fasta')

    def write_orf_file(self,
                       orfEntries: Union[list[ORFinfder.ORF], ORFinfder.ORF] = None) -> object:
        """
        A method to write a fasta file from a list of ORF objects

        Args:
            orfEntries: A list of ORF object to add in the fasta file

        """
        seq_list =[]

        if isinstance(orfEntries, ORFinfder.ORF):
            orfEntries = [orfEntries]

        for orf in orfEntries:
            seq_list.append(SeqRecord(orf.seq, id=orf.id, description=f"{orf.src}:{orf.start_pos}-{orf.end_pos}({orf.frame[0]})"))


        SeqIO.write(seq_list, self.path, "fasta")



