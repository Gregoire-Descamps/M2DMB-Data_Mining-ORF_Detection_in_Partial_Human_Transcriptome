from __future__ import annotations
import io
import FileFormats
from typing import Union, Self
import logging
import os

from Bio import SeqIO, SeqRecord
from Bio.Seq import Seq, MutableSeq
from Bio.Data.CodonTable import standard_dna_table
from Bio.Application import _Option
from Bio.Blast.Applications import NcbiblastxCommandline


from FileFormats import GffFile, FastaFile

# logging config
logging.basicConfig(level=logging.DEBUG)


class ORFExtractor:

    def __init__(self, file: str):
        """
        A class that extract ORFs from fasta file
        :param file: path to a fasta file
        """
        # Initialize the dictionaries containing the ORFs
        self.file_path: str = file
        self.seqIO: SeqIO.Iterable = self._load_seqio_fasta(file, "fasta")
        self.ORF_dict: dict[str,ORF] = {}
        self.__seq_orf_list: list[ORF] = []
        self.non_orf_seq = 0

        if os.path.exists("output"):
            logging.warning("output directory already existing, file(s) inside might be replaced")
        else:
            os.mkdir("output")

        self.output = False
        self.gff_file: FileFormats.GffFile = None
        self.fasta_file: FileFormats.FastaFile = None
        self.gff_handler: io.TextIOWrapper = None
        self.fasta_handler: io.TextIOWrapper = None
        self.__last_orf_id: int = 0

        # self.blastresults: FileFormats.BlastTsvFile = None

        # for test purpose (skip blastx)
        self.blastresults: FileFormats.BlastTsvFile = FileFormats.BlastTsvFile("/analysis/output/blast_results.tsv")

    def _load_seqio_fasta(self,
                          file: str,
                          file_format: str):
        """
        Load the fasta file in a SeqIO object
        :param file: path to the fasta file to load
        :param file_format: type of file format to load (default : fasta)
        :return: seqIO object
        """

        obj: Seq = SeqIO.parse(file, file_format)
        return obj

    def __extract_all_codons(self,
                             seq: Seq,
                             codons: Union[list[str], str],
                             strand: str,
                             start_pos: int = None,
                             end_pos: int = None):

        """
        Extract all codons positions from a given sequence.
        The position indexes are related to the direction of the sequence (3' -> 5') and range from 0 to N-1

        :param seq: The sequence to search in, it must be a Bio.Seq object
        :param codons: Either a string or a list of string to search for
        :param strand: used to build the dictionary key and distinguish positive and negative strand positions
        :param start_pos: Optional, the lowest initial position to search codon
        :param end_pos: Optional, the highest position to search codon
        :return: dictionary of position list by ORF
            (if no position has been found the ORF doesn't appear in the dictionary keys)

        """

        pos_dict = {
            "+1": [],
            "+2": [],
            "+3": [],
            "-1": [],
            "-2": [],
            "-3": [],
        }

        # reformat str single codon as list
        if isinstance(codons, str):
            codons = [codons]

        # find all codons
        for codon in codons:  # type: str
            found = True
            last_pos = start_pos
            while found:
                pos: int = seq.find(codon, last_pos, end_pos)
                # if a position is found
                if pos != -1:
                    pos_dict[strand + str(pos % 3 + 1)].append(pos)
                    last_pos = pos + 1
                else:
                    found = False

        # remove empty ORFs pos lists and sort non-empty lists in asc order
        del_keys = []
        for orf_key, orf in pos_dict.items():
            if not orf:
                del_keys.append(orf_key)
            else:
                orf.sort()
        for key in del_keys:
            del pos_dict[key]

        return pos_dict

    def __clean_orf_pos_dict(self,
                             dict1: dict[str, list],
                             dict2: dict[str, list]):

        """
        remove non-paired key between two dictionaries
            return the 2 cleared dicts

            if an ORF contains only start or only stop codons, it can't translate

        :param dict1: first codon position dictionary to compare
        :param dict2: second codon position dictionary to compare
        :return: Cleared dict1, Cleared dict2
        """

        remove_list1 = []
        remove_list2 = []

        for key in dict1.keys():
            if not (key in dict2.keys()):
                remove_list1.append(key)

        for key in dict2.keys():
            if not (key in dict1.keys()):
                remove_list2.append(key)

        [dict1.pop(key) for key in remove_list1]
        [dict2.pop(key) for key in remove_list2]

        return dict1, dict2

    def __orf_extract(self,
                      frame: str,
                      start_pos_list: list[int],
                      stop_pos_list: list[int],
                      seq: SeqRecord.SeqRecord,
                      rev_seq: SeqRecord.SeqRecord,
                      trigger: int = 0):
        """
        Extract all ORF of a given frame from the potential start and stop positions
        Position indexes should range from 0 to N-1
        All eligible ORF get stored in ORF objects before removing overlapping ORF

        :param frame: the frame to search in
        :param start_pos_list: the list of start codons positions found for this frame
        :param stop_pos_list: the list of start codons positions found for this frame
        :param seq: the sequence of the positive strand (in 3' -> 5' orientation)
        :param rev_seq: the sequence of the negative strand (in 3' -> 5' orientation)
        :param trigger: the minimum length of an ORF to be added in the dict
        :return: None
        """

        while start_pos_list and stop_pos_list:

            # alt_start is not used in this project but can help identify nested ORFs
            alt_start = []
            init_pos = start_pos_list.pop(0)

            # remove preceding stop codons
            while stop_pos_list and stop_pos_list[0] < init_pos:
                stop_pos_list.pop(0)

            if stop_pos_list:
                # correct end pos to the last nucleotide of the seqence
                end_pos = stop_pos_list[0] + 2
                # only process orf longer than trigger
                if end_pos - init_pos > trigger:
                    # extract alternative starts in the ORF
                    while start_pos_list and start_pos_list[0] < end_pos:
                        alt_start.append(start_pos_list.pop(0))
                    if frame[0] == "+":
                        orf_seq = seq.seq[init_pos:end_pos + 1]
                        self.__seq_orf_list.append(ORF(seq=Seq(orf_seq), seq_id=seq.id, start_pos=init_pos + 1,
                                                       end_pos=end_pos + 1, frame=frame, parent=seq))

                        # convert alt_start positions
                        alt_start = [x + 1 for x in alt_start]
                        self.__seq_orf_list[-1].add_alt_start(alt_start)

                    else:
                        orf_seq = rev_seq.seq[init_pos:end_pos + 1]
                        self.__seq_orf_list.append(
                            ORF(seq=Seq(orf_seq), seq_id=seq.id, start_pos=len(seq) - end_pos,
                                end_pos=len(seq) - init_pos, frame=frame, parent=seq))

                        # self.__add_orf_to_dict(orf_seq, seq.id, len(seq) - end_pos, len(seq) - init_pos, frame)

                        # convert alt_start positions
                        alt_start = [len(seq) - x for x in alt_start]
                        self.__seq_orf_list[-1].add_alt_start(alt_start)

    def __aggregate_overlapping_orf(self):
        """
        Aggregate all ORF of the sequence to remove overlapping ORF
        :return: None .
        """
        while len(self.__seq_orf_list) > 1:
            non_overlapping_orf: list[ORF] = []
            max_orf_length = 0
            longest_orf = None
            index_to_remove = None
            # select longest ORF
            index = 0
            for orf in self.__seq_orf_list:
                if orf.length > max_orf_length:
                    longest_orf = orf
                    index_to_remove = index
                index += 1

            # remove the longest orf from list
            self.__seq_orf_list.pop(index_to_remove)

            max_end_pos = longest_orf.end_pos
            max_start_pos = longest_orf.start_pos

            # putting aside non overlapping ORFs, register overlapping ORfs in the longest ORF
            for orf in self.__seq_orf_list:
                if (orf.frame == longest_orf.frame and
                        (max_start_pos <= orf.start_pos <= max_end_pos
                         or max_start_pos <= orf.end_pos <= max_end_pos)):
                    # removing orf sequence to reduce object size before aggregating it
                    orf.seq = None
                    longest_orf.add_nested_orf(orf)
                else:
                    non_overlapping_orf.append(orf)

            # export orf entry to dict and output files
            self.entry_export(longest_orf)


            self.__seq_orf_list = non_overlapping_orf

        if len(self.__seq_orf_list) == 1:
            self.entry_export(self.__seq_orf_list.pop(0))

        if len(self.__seq_orf_list) > 0:
            raise Exception("remaining sequences in the list at the end of aggregation, this shouldn't be the case")


    def extract(self):
        """
        Method to extract all ORF from the fasta file
        All information get stored in the ORF_dict of the object
        Create output fasta and gff files
        :return: self.ORf_dict
        """
        self.gff_file = GffFile("output/ORF_extract.gff")
        self.fasta_file = FastaFile("output/ORF_extract.fasta")

        self.gff_file.set_header(desc="Potential CDS entries extracted from human transcriptome assembly",
                                 provider="Gregoire Descamps",
                                 contact="gregÔire.desc@mps@thisisnotarealemail.com",
                                 date="Today")

        self.gff_handler = self.gff_file.open_file()
        self.fasta_handler = self.fasta_file.open_file()

        max_seq_length = 0
        for seq in self.seqIO:  # type: SeqRecord.SeqRecord
            self.gff_file.add_entry(seqid=seq.id, source="manual", seq_type="transcript", start=1, end=len(seq.seq), score=".", strand=".", phase=".", attributes={"ID": seq.id})
            rev_seq: SeqRecord = seq.reverse_complement()
            self.__seq_orf_list = []
            if len(seq.seq) > max_seq_length:
                max_seq_name = seq.id
                max_seq_length = len(seq.seq)
            start_codons_pos = {}
            stop_codons_pos = {}

            # extract start codons
            for frame, pos_list in self.__extract_all_codons(seq.seq, "ATG", "+").items():
                start_codons_pos[frame] = pos_list
            for frame, pos_list in self.__extract_all_codons(rev_seq.seq, "ATG", "-").items():
                start_codons_pos[frame] = pos_list

            # extract stop codons
            for frame, pos_list in self.__extract_all_codons(seq.seq, standard_dna_table.stop_codons, "+").items():
                stop_codons_pos[frame] = pos_list
            for frame, pos_list in self.__extract_all_codons(rev_seq.seq, standard_dna_table.stop_codons, "-").items():
                stop_codons_pos[frame] = pos_list

            # cleaning codon pos dicts
            start_codons_pos, stop_codons_pos = self.__clean_orf_pos_dict(start_codons_pos, stop_codons_pos)

            # extract ORF:
            for frame in stop_codons_pos.keys():
                self.__orf_extract(frame, start_codons_pos[frame], stop_codons_pos[frame], seq, rev_seq)
            # aggregate and filter  overlapping ORFs
            self.__aggregate_overlapping_orf()

            if self.ORF_dict["ORF_"+str(self.__last_orf_id)].seq_id != seq.id:
                self.non_orf_seq += 1

        # closing files
        self.fasta_handler.close()
        self.gff_handler.close()

        self.output = True

        return self.ORF_dict

    def entry_export(self, orf: ORF):
        """
        Append gff and multifasta files with the extracted ORF

        Args:
            orf: an ORF object
        """
        # set orf ID
        self.__last_orf_id += 1
        orf.set_id("ORF_" + str(self.__last_orf_id))

        # write entry to gff and fasta file and append orf to dict
        self.gff_file.add_entry(*orf.to_gff_tuple(), open_file=self.gff_handler)
        self.fasta_file.add_entry(*orf.to_fasta_tuple(), open_file=self.fasta_handler)
        self.ORF_dict[orf.id] = orf

    def blastx(self, database: str, outputfile: str, **kwargs):
        """
        Run a blastx against the source fasta file using Biopython NcbitblastxCommandline

        Args:
            database: path to the blast database (should be a directory and not a file)
            outputfile: results output path
            **kwargs: arguments used to run Blastx see https://www.ncbi.nlm.nih.gov/books/NBK279684/#appendices.Options_for_the_commandline_a for the Ncbi command line arguments
                    and https://biopython.org/docs/1.76/api/Bio.Blast.Applications.html for the biopython documentation

        Returns: class:tuple (STDOUT,STDERR)

        """
        if not self.output:
            raise Exception("Unable to find an output file use extract method to generate the corresponding fasta file")

        # add the taxids argument as extra parameter
        extra_parameters = []
        if "taxids" in kwargs.keys() or "-taxids" in kwargs.keys():

            extra_parameters.append(_Option(
                ["-taxids", "taxids"],
                "Filtering algorithm for soft masking (integer).\n\n"
                "Filtering algorithm ID to apply to BLAST database as soft masking. "
                "Incompatible with: db_hard_mask, subject, subject_loc",
                equate=False,
            ))

        # add parameters to the biopython blastx object
        bio_blast = NcbiblastxCommandline(query=self.fasta_file.path, db = database,
                                          out=outputfile, num_threads=os.cpu_count()-2)
        bio_blast.parameters = bio_blast.parameters + extra_parameters

        for key, value in kwargs.items():
            bio_blast.set_parameter(key, value)

        bio_blast()
        self.blastresults = FileFormats.BlastTsvFile(outputfile)
        return self.blastresults


    def orf_validate(self):
        self.blastresults.top_hit_extract()

        for src, hit in self.blastresults.parse():
            try:
                self.ORF_dict[hit["query_id"]].to_CDS()
                # add UniProtKB/Swiss-Prot Dbxref
                self.ORF_dict[hit["query_id"]].add_attr({"Dbxref":"UniProtKB/Swiss-Prot:"+hit["subject_id"]})
            except:
                print(src, hit)
                raise

        self._gff_update()

    def _gff_update(self):
        tempfile = FileFormats.GffFile("output/temp.gff")
        tempfile.set_header(desc="Potential CDS entries extracted from human transcriptome assembly",
                            provider="Gregoire Descamps",
                            contact="gregÔire.desc@mps@thisisnotarealemail.com",
                            date="Today")

        for entry in self.gff_file.parse():
            id = entry[8]["ID"]
            if id.startswith("cdna_"):
                tempfile.add_entry(*entry)
            else:
                dict_entry=self.ORF_dict[id]
                if dict_entry.type == "CDS":
                    tempfile.add_entry(*dict_entry.utr_3.to_gff_tuple())
                    tempfile.add_entry(*dict_entry.utr_5.to_gff_tuple())

                # keep also putative ORFs
                tempfile.add_entry(*dict_entry.to_gff_tuple())

        gff_path = self.gff_file.path
        os.remove(gff_path)
        self.gff_file = tempfile.rename(gff_path.split("/")[-1])
        return self





class ORF:
    """
    A class representing an ORF

    Attributes:

    - :class:`Seq` seq -> The sequence of the ORF
    - :class:`str` seq_id -> the contig source id
    - :class:`int` length -> length of the ORF sequence
    - :class:`int` start_pos -> position of the first ORF's base on the source contig
    - :class:`int` end_pos -> position of the last ORF's base on the source contig
    - :class:`str` frame -> reading frame of the ORF
    - :class:`list[int]` _alt_start -> position of alternative start codons within the ORF
    - :class:`list[ORF]` _nested -> list of nested or overlapping ORFs in other reading frames

    Methods:

    - add_alt_start -> Add alternative start codon(s) position in the ORF object
    - add_nested_orf -> Add nested ORF objects in the ORF object
    - to_dict_tuple -> Return OFR information as tuple for dictionary entry
    - to_gff_tuple -> Return OFR information as tuple for gff file entry
    """

    def __init__(self,
                 seq: Seq,
                 seq_id: str,
                 start_pos: int,
                 end_pos: int,
                 src: str = "manual",
                 score: str = ".",
                 frame: str = ".",
                 id: str = ".",
                 parent: SeqRecord.SeqRecord = None,
                 attributes: dict = None,
                 seq_type: str = "ORF"):

        self.seq = seq
        self.seq_id = seq_id
        self.src = src
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.score = score
        self.frame = frame
        self.length = end_pos - start_pos + 1
        self._alt_start = []
        self._nested: list[ORF] = []
        self.id = id
        self.parent = parent

        # Append attributes
        if isinstance(attributes, dict):
            self.attributes = attributes
        else:
            self.attributes = {}

        # Add ID attribute
        if self.id != "":
            self.attributes["ID"] = self.id

        # Add parent attribute
        if parent is not None:
            self.attributes["parent"] = parent.id

        self.type = seq_type

    def __len__(self):
        return len(self.seq)

    def set_id(self, id: str):
        self.id = id
        self.attributes["ID"] = id

    def add_attr(self, attr: Union[dict, tuple]):
        """
        Add attributes to the sequence, can be directly as dict or a tuple of (key, value)
        Args:
            attr: Attribute to add
        """
        if isinstance(attr, tuple) and len(attr) == 2:
            self.attributes[attr[0]] = attr[1]

        elif isinstance(attr, dict):
            for key, value in attr.items():
                self.attributes[key] = value

        else:
            raise Exception(f"The attribute doesn't have the right format, it should be a dict or a (key, value) pair."
                            f"{type(attr)} was provided")

    def add_alt_start(self, pos: Union[int, list[int]]):
        if isinstance(pos, list) and all(type(x) is int for x in pos):
            for position in pos:
                self._alt_start.append(position)
        elif isinstance(pos, int):
            self._alt_start.append(pos)
        else:
            raise TypeError(f"pos must be either int or list[int], not {type(pos)}")

    def add_nested_orf(self, orf: Self):
        if not isinstance(orf, ORF):
            raise TypeError(f"orf must be an ORF object, not {type(orf)}")
        self._nested.append(orf)

    def to_dict_tuple(self):
        """
        DEPRECATED  a method to provide ORF properties as tuple

        Returns: (ORF.id, ORF.seq, ORF.seq_id, ORF.start_pos, ORF.end_pos, ORF.frame)

        """
        return self.id, self.seq, self.seq_id, self.start_pos, self.end_pos, self.frame

    def to_gff_tuple(self):
        """
        :return: a tuple to fill a GffFile.add_entry method
        """

        if self.frame == ".":
            phase = "."
        else:
            phase = str(int(self.frame[-1]) - 1)

        return (self.seq_id,
                self.src,
                self.type,
                self.start_pos,
                self.end_pos,
                self.score,
                self.frame[0],
                phase,
                self.attributes)
    def to_fasta_tuple(self):
        """
       :return: a tuple to fill a GffFile.add_entry method
       """
        desc = f"{self.seq_id}:{self.start_pos}-{self.end_pos}({self.frame[0]})"
        return self.id, desc, self.seq

    def to_CDS(self):
        self.type = "CDS"
        self.set_id("CDS_"+self.id.split("_")[1])
        self.src = "BLASTX 2.16.0+"

        if self.frame[0] == "+":
            start3 = 1
            end3 = self.start_pos - 1
            start5 = self.end_pos + 1
            end5 = len(self.parent.seq)
            seq3 = Seq(self.parent.seq[start3:end3 + 1])
            seq5 = Seq(self.parent.seq[start5:end5 + 1])


        else:
            start3 = len(self.parent.seq)
            end3 = self.end_pos + 1
            start5 = self.start_pos - 1
            end5 = 1
            seq3 = MutableSeq(self.parent.seq[end3:start3 + 1]).reverse_complement(inplace=True)
            seq5 = MutableSeq(self.parent.seq[end5:start5 + 1]).reverse_complement(inplace=True)

        self.utr_3 = ORF(seq=seq3, seq_id=self.seq_id, seq_type="UTR3", start_pos=start3, end_pos=end3,
                         frame=self.frame, id="UTR3a", parent=self.parent, attributes={})
        self.utr_5 = ORF(seq=seq5, seq_id=self.seq_id, seq_type="UTR5", start_pos=start5, end_pos=end5,
                         frame=self.frame, id="UTR5a", parent=self.parent, attributes={})

