import FileFormats
from typing import Union, Self
import logging

from Bio import SeqIO, SeqRecord
from Bio.Seq import Seq
from Bio.Data.CodonTable import standard_dna_table

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
        self.ORF_dict = {}
        self.ORF_list: list[ORF] = []
        self.non_orf_seq = 0

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
                end_pos = stop_pos_list[0] + 2
                # only process orf longer than trigger
                if end_pos - init_pos > trigger:
                    while start_pos_list and start_pos_list[0] < end_pos:
                        alt_start.append(start_pos_list[0])
                        start_pos_list.pop(0)
                    if frame[0] == "+":
                        orf_seq = seq.seq[init_pos:end_pos + 1]
                        self.__seq_orf_list.append(ORF(Seq(orf_seq), seq.id, init_pos + 1, end_pos + 1, frame))
                        # self.__add_orf_to_dict(orf_seq, seq.id, init_pos + 1, end_pos + 1, frame)

                        # convert alt_start positions
                        alt_start = [x + 1 for x in alt_start]
                        self.__seq_orf_list[-1].add_alt_start(alt_start)

                    else:
                        orf_seq = rev_seq.seq[init_pos:end_pos + 3]
                        self.__seq_orf_list.append(
                            ORF(Seq(orf_seq), seq.id, len(seq) - end_pos, len(seq) - init_pos, frame))
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
            self.ORF_list.append(longest_orf)
            self.__seq_orf_list = non_overlapping_orf

        if len(self.__seq_orf_list) == 1:
            self.ORF_list.append(self.__seq_orf_list.pop(0))

        if len(self.__seq_orf_list) > 0:
            raise Exception("remaining sequences in the list at the end of aggregation, this shouldn't be the case")

    def __add_orf_to_dict(self,
                          id: str,
                          seq: str,
                          src: str,
                          start_pos: int,
                          stop_pos: int,
                          frame: str):
        """
        Add ORF information into the ORF dict
        :param seq: The ID string of the ORF
        :param seq: The string of the complete sequence of the ORF
        :param src: The id of the contig or the source sequence the ORF is found in
        :param start_pos: The first base position of the ORF, should be in the positive strand direction,
            (3' -> 5' or 5' -> 3' for negative strand) and the indexes should range from 1 to N.
        :param stop_pos: The last base position of the ORF, should be in the positive strand direction,
            (3' -> 5' or 5' -> 3' for negative strand) and the indexes should range from 1 to N.
        :param frame:
        :return:
        """
        self.ORF_dict[id] = (seq, stop_pos - start_pos + 1, src, (start_pos, stop_pos), frame)

    def extract(self):
        """
        Method to extract all ORF from the fasta file
        All information get stored in the ORF_dict of the object
        :return: self.ORf_dict
        """
        max_seq_length = 0
        for seq in self.seqIO:  # type: SeqRecord.SeqRecord
            rev_seq: SeqRecord = seq.reverse_complement()
            self.__seq_orf_list: list[ORF] = []
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

            if self.ORF_list[-1].src != seq.id:
                self.non_orf_seq += 1

        # add all saved ORF to dict
        i = 1
        for orf in self.ORF_list:
            orf.id = "ORF_" + str(i)
            self.__add_orf_to_dict(*orf.to_dict_tuple())
            i += 1

        return self.ORF_dict

    def result_export(self,
                     output_name: str):
        """
        Create gff and multifasta files with all the extracted ORF

        Args:
            output_name: Name for the generated files
        """
        gffFile = GffFile(output_name+".gff")
        fastaFile = FastaFile(output_name + ".fasta")

        gffFile.write_orf_file(self.ORF_list)
        fastaFile.write_orf_file(self.ORF_list)


class ORF:
    """
    A class representing an ORF

    Attributes:

    - :class:`Seq` seq -> The sequence of the ORF
    - :class:`str` src -> the contig source id
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
                 src: str,
                 start_pos: int,
                 end_pos: int,
                 frame: str):

        self.seq = seq
        self.src = src
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.frame = frame
        self.length = end_pos - start_pos + 1
        self._alt_start = []
        self._nested: list[ORF] = []
        self.id = None

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
        return self.id, self.seq, self.src, self.start_pos, self.end_pos, self.frame

    def to_gff_tuple(self):
        """
        :return: a tuple to fill a GffFile.add_entry method
        """
        return (self.src,
                "custom_ORF_finder",
                "ORF",
                self.start_pos,
                self.end_pos,
                ".",
                self.frame[0],
                int(self.frame[-1]) - 1,
                {"ID": self.id})




