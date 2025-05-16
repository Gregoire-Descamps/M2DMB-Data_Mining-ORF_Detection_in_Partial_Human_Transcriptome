from __future__ import annotations
import io
import os
from typing import Union
import datetime
import logging
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import ORFinder


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

    def open_file(self):
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
                  open_file: io.TextIOWrapper = None):
        """
            Append a new entry at the end of the file

            :param seqid: The sequence id
            :type seqid: str
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
            :param open_file: A file object
            :type open_file: str

            seqid (object):



        """

        # converting dictionary into attribute string
        if isinstance(attributes, dict):
            new_attributes = ""
            for key, value in attributes.items():
                new_attributes += f';{key}={value}'

            attributes = new_attributes[1:]

        entry_string = f'{seqid}\t{source}\t{seq_type}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{attributes}\n'

        if not open_file:
            with self.open_file() as file:
                file.write(entry_string)

            file.close()
        else:
            open_file.write(entry_string)

    def write_orf_file(self, orfEntries: Union[dict[ORFinder.ORF], ORFinder.ORF]):
        """
        A method to write a gff file from a list of ORF objects

        Args:
            orfEntries: A list of ORF object to add in the gff file

        Returns:
            object:
        """

        if isinstance(orfEntries, ORFinder.ORF):
            orfEntries= {orfEntries.id: orfEntries}



        with self.open_file() as file:
            for orf in orfEntries.values():
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


    def open_file(self):
        return open(self.path, "a")

    def add_entry(self,
                  seqid: str,
                  desc: str,
                  seq: str,
                  open_file: io.TextIOWrapper = None):

        """
            Append a new entry at the end of the file

            :param seqid: The sequence id
            :param desc:  Entry description
            :param seq: The sequence
            """

        sequence = SeqRecord(seq,id = seqid, description= desc)


        if not open_file:
            with self.open_file() as file:
                SeqIO.write(sequence, file, 'fasta')
            file.close()

        else:
            SeqIO.write(sequence, open_file, 'fasta')



    def write_orf_file(self, orfEntries: Union[dict[ORFinder.ORF], ORFinder.ORF]):
        """
        A method to write a fasta file from a list of ORF objects

        Args:
            orfEntries: A list of ORF object to add in the fasta file

        """
        seq_list =[]

        if isinstance(orfEntries, ORFinder.ORF):
            orfEntries = {orfEntries.id: orfEntries}

        for orf in orfEntries.values():
            seq_list.append(SeqRecord(orf.seq, id=orf.id, description=f"{orf.seq_id}:{orf.start_pos}-{orf.end_pos}({orf.frame[0]})"))


        SeqIO.write(seq_list, self.path, "fasta")


class BlastTsvFile:
    """
            A class representing a blast result tsv file

            Args:

             path: The file path

            Methods:

            - parse: -> Set the header information
            -  add_entry -> append a new line to the
            - write_orf_file -> write a gff file from a list of orf objects
            """
    def __init__(self, path):
        if not os.path.isfile(path):
            raise FileNotFoundError(f"File not found: {path}")
        self.path = path
        self.comment_lines = []
        self.__line = None

    def parse(self):
        """
        A method that yields each BLAST result line
         return hits as dict.
        Skips comment lines starting with '#'.
        """
        with open(self.path, 'r') as f:
            for self.__line in f:
                if self.__line.startswith("#"):
                    self.comment_lines.append(self.__line.strip())
                    if len(self.comment_lines)>5:
                        self.comment_lines = [self.comment_lines[-1]]
                    if self.__line.startswith("# Query"):
                        src = self.__line.split(":")[1].split(" ")[2]
                    continue

                yield src, self._parse_line(self.__line.strip())
                self.comment_lines =[]

    def _parse_line(self, line):
        """
        Parses a single line of BLAST outfmt 7 output.
        Returns a dictionary of relevant fields.
        """
        fields = [
            "query_id", "subject_id", "percent_identity", "alignment_length",
            "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end",
            "evalue", "bit_score"
        ]
        parts = line.split('\t')
        return dict(zip(fields, parts))



    def top_hit_extract(self, outpath: str):
        """
        Yields the top hit per source (cdna) on best (lowest) e-value and store it in a new file.
        Assumes there is only one CDS per source (cdna)
        """
        best_hit = None
        best_comment = None
        cdna = None

        with open(outpath,"wt") as fout:
            for src, hit in self.parse():
                if src == cdna:
                    #  TEST
                    # print(f" type of q_start : {type(float(hit['q_start']))}")
                    # print(f" value of q_start : {hit['q_start']}")
                    # print(f" q_start % 3 : {float(hit['q_start'])%3}")
                    # print(f"{float(hit['q_start']) < float(hit['q_end'])}")
                    # print(f" type of evalue : {type(float(hit['evalue']))}")
                    # print(f" value of evalue : {float(hit['evalue'])}")

                    # END TEST

                    # checking for right frame
                    if (float(hit["q_start"])<float(hit["q_end"])) and float(hit["q_start"])%3 == 1.0:
                        if not best_hit or float(best_hit["evalue"]) > float(hit["evalue"]):
                            best_hit = hit
                            best_comment = self.comment_lines
                    # else:
                    #     fout.write("WRONG\t"+hit["query_id"]+"\n")
                else:
                    # append best hit to gff output
                    # print(f"New source, from {cdna} to {src}\n{self.__line}Best hit\t{best_hit}\n")
                    if best_hit:
                        fout.write("\n".join(best_comment)+"\n")
                        fout.write(f"{best_hit['query_id']}\t"
                                   f"{best_hit['subject_id']}\t"
                                   f"{best_hit['percent_identity']}\t"
                                   f"{best_hit['alignment_length']}\t"
                                   f"{best_hit['mismatches']}\t"
                                   f"{best_hit['gap_opens']}\t"
                                   f"{best_hit['q_start']}\t"
                                   f"{best_hit['q_end']}\t"
                                   f"{best_hit['s_start']}\t"
                                   f"{best_hit['s_end']}\t"
                                   f"{best_hit['evalue']}\t"
                                   f"{best_hit['bit_score']}\n")
                    cdna = src
                    best_hit = hit
                    best_comment = self.comment_lines
        fout.close()
        return BlastTsvFile(outpath)
