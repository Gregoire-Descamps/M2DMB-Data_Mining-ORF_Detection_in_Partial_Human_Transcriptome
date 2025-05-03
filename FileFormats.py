import io
import os
from typing import Union
import datetime
import logging
from ORFinfder import ORF


# logging config
logging.basicConfig(level=logging.DEBUG)

class GffFile:
    """
        A class representing an GFF File

        Attributes:

        - :class:`str` path -> The output file path


        Methods:

        -  set_header -> Set the header information
        - add_entry -> append a new line to the file
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

        # converting dictionary into attribute string
        if isinstance(attributes, dict):
            new_attributes = ""
            for key, value in attributes.items():
                new_attributes += f';{key}={value}'

            attributes = new_attributes[1:]
        if not open_file:
            with self.__open_file() as file:
                file.write(f'{seqid}\t{source}\t{seq_type}\t{start}\t{end}\t{score}\t\t{strand}\t{phase}\t{attributes}\n')

            file.close()
        else:
            open_file.write(f'{seqid}\t{source}\t{seq_type}\t{start}\t{end}\t{score}\t\t{strand}\t{phase}\t{attributes}\n')

    def write_orf_file(self,
              orfEntries : Union[list[ORF], str] = None):

        if isinstance(orfEntries, str):
            orfEntries= [orfEntries]

        self.set_header(desc="Potential CDS entries extracted from human transcriptome assembly",
                              provider="Gregoire Descamps",
                              contact="gregÔire.desc@mps@thisisnotarealemail.com",
                              date="Today")
        with self.__open_file() as file:
            for orf in orfEntries:
                self.add_entry(*orf.to_gff_tuple(), open_file = file)
        file.close()
