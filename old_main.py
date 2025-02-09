from typing import Union
from Bio import SeqIO, SeqRecord
from Bio.Seq import Seq
from Bio.Data.CodonTable import standard_dna_table
fastaFile = "Homo_sapiens_cdna_assembed.fasta"

# print(f'{standard_dna_table.start_codons}')
# print(f'{standard_dna_table.stop_codons}')
# print(standard_dna_table)


class ORFExtractor:

    # A class that extract ORFs from fasta file

    def __init__(self, file):

        # Initialize the dictionaries containing the ORFs

        self.open: dict[str, dict[int, list[ORF]]] = \
                    {"Plus": {
                             1: [],
                             2: [],
                             3: []},
                     "Minus": {
                             1: [],
                             2: [],
                             3: []}
                     }
        self.closed: dict[str, dict[int, list[ORF]]] = {
            "Plus": {
                    1: [],
                    2: [],
                    3: []
            },
            "Minus": {
                    1: [],
                    2: [],
                    3: []
            }
        }
        self.file: str = file
        self.seqIO: Seq = self._LoadSeqIOFast(self.file, "fasta")
        self.ORF_dict = {}

    def _LoadSeqIOFast(self,
                       file: str,
                       format: str):

        # Load the fasta file in a SeqIO object
        #     file : (str) path to the fasta file to load
        #     format : (str) type of file format to load
        #     obj : (seqIO) the seqIO object returned

        # load a fastq file into an iterable SeqIO object
        obj: Seq = SeqIO.parse(file, format)
        return obj

    def __n_strand_append(self,
                          codon: Seq,
                          frame: int):

        # append codon in all opened ORF for negative strand
        #     codon : (seqIO) codon to insert
        #     frame : (int) number of the frame the codon originate
        #     orf : (list) list of ORF in the dictionary

        for orf in self.open["Minus"][frame]:  # type: ORF
            orf.appendCodon(codon)

    def __p_strand_append(self,
                          codon: Seq,
                          frame: int):

        # append codon in all opened ORF for positive strand
        #     codon : (SeqRecord) codon to insert
        #     frame : (int) number of the frame the codon originate
        #     orf : (list) list of ORF in the dictionary

        for orf in self.open["Plus"][frame]:  # type: ORF
            orf.appendCodon(codon)

    def __n_strand_open(self,
                        codon: Seq,
                        frame: int,
                        pos: int,
                        src_id: str):

        # Open a new ORF for negative strand
        #     codon : (Seq) start codon
        #     frame : (int) number of the frame the codon originate
        #     pos : (int) position of the start codon
        #     src_id : (str) id of the source contig

        dest = self.open["Plus"][frame]
        if len(dest) == 0:
            dest.append(ORF(id=f"ORF_{src_id}_N_{pos}",startPos=pos, strand="plus",codon=codon))
        else:
            dest.append(ORF(id=f"ORF_{src_id}_N_{pos}",startPos=pos,codon=codon,strand="plus", parent=dest[-1]))

    def __p_strand_open(self,
                        codon: Seq,
                        frame: int,
                        pos: int,
                        src_id: str):

        # Open a new ORF for positive strand
        #     codon : (Seq) start codon
        #     frame : (int) number of the frame the codon originate
        #     pos : (int) position of the start codon
        #     src_id : (str) id of the source contig

        dest = self.open["Plus"][frame]
        if len(dest) == 0:
            dest.append(ORF(id=f"ORF_{src_id}_P_{pos}",startPos=pos, codon=codon, strand="minus"))
        else:
            dest.append(ORF(id=f"ORF_{src_id}_P_{pos}",startPos=pos,codon=codon, strand="minus", parent=dest[-1]))

    def __p_strand_close(self, frame, pos):

        # Close all ORF for positive strand
        #     codon : (Seq) start codon
        #     frame : (int) number of the frame the codon originate
        #     pos : (int) position of the start codon
        #     src_id : (str) id of the source contig

        dest_open: list[ORF] = self.open["Plus"][frame]
        dest_close = self.closed["Plus"][frame]

        for orf in dest_open: #type: ORF
            orf.closeSeq(pos)
            dest_close.append(orf)

        dest_open.clear()
        # print(f'num of orf+{frame} remaining : {len(self.open["Plus"][frame])}')

    def __n_strand_close(self, frame, pos):

        # Close all ORF for negative strand
        #     codon : (Seq) start codon
        #     frame : (int) number of the frame the codon originate
        #     pos : (int) position of the start codon
        #     src_id : (str) id of the source contig

        dest_open: list[ORF] = self.open["Minus"][frame]
        dest_close = self.closed["Minus"][frame]

        for orf in dest_open:  # type: ORF
            orf.closeSeq(pos)
            dest_close.append(orf)

        dest_open.clear()
        # print(f'num of orf-{frame} remaining : {len(self.open["Minus"][frame])}')

    def extract(self):

        # Extract ORFs from the file in the object
        #     seq : (SeqRecord sequence entry from the fasta file
        #     rev_seq : (SeqRecord) reverse complement of seq
        #     pos : (int) position of the first base of the codon in the sequence (in transcription direction)
        #     codon : (Seq) seq of the codon in the "forward" sequence
        #     rev_codon : (Seq) seq of the codon in reverse complement sequence
        #     read_frm : (int) number of the frame (1, 2, 3)

        for seq in self.seqIO:  # type: SeqRecord.SeqRecord
            pos = 1
            # make complementary reverse of the sequence
            rev_seq: SeqRecord = seq.reverse_complement()

            # parsing sequence (in fwd and reverse)
            while pos < len(seq)-3:  # type: int
                codon:Seq = Seq(seq[pos-1:pos+2].seq)
                rev_codon: Seq = Seq(rev_seq[pos-1:pos+2].seq)

                # number of the frame (1, 2 or 3)
                read_frm = pos % 3 + 1

                # Insert corresponding codon in all opened ORFs
                self.__n_strand_append(rev_codon, read_frm)
                self.__p_strand_append(codon, read_frm)

                # print(codon)
                # open new ORF if Start codon
                if codon == 'ATG': # in standard_dna_table.start_codons:
                    self.__p_strand_open(codon, read_frm, pos, seq.id)
                    # print(f"found new ORF! {seq.id}(ORF{read_frm}) : {codon}")
                if rev_codon == 'ATG' : #in standard_dna_table.start_codons:
                    self.__n_strand_open(rev_codon, read_frm, len(seq)-pos+1, seq.id)
                    # print(f"found new ORF! {seq.id}(ORF{read_frm}) : {rev_codon}")

                # close all ORF if stop codon found
                if codon in standard_dna_table.stop_codons:
                    self.__p_strand_close(read_frm, pos)
                    # print(f"closed ORF! {seq.id}(ORF{read_frm}) : {codon}")
                if rev_codon in standard_dna_table.stop_codons:
                    self.__n_strand_close(read_frm, len(seq) - pos + 1)
                    # print(f"closed ORF! {seq.id}(ORF{read_frm}) : {rev_codon}")


                pos +=1

            # emptying opened ORFs
            for strand in self.open.items():
                for frame in strand[1].items():
                    if len(frame[1])>0:
                        print(f'max Unclosed ORF length for frame{strand[0]}{frame[0]} : {max(len(orf.codSeq) for orf in frame[1])}')
                    frame[1].clear()

            print(f'Sequence {seq.id} parsed! '
                  f'{len(self.closed["Minus"][1]) + len(self.closed["Minus"][2]) + len(self.closed["Minus"][3]) + len(self.closed["Plus"][1]) + len(self.closed["Plus"][2]) + len(self.closed["Plus"][3])} ORF found! '
                  f'{len(self.open["Minus"][1]) + len(self.open["Minus"][2]) + len(self.open["Minus"][3]) + len(self.open["Plus"][1]) + len(self.open["Plus"][2]) + len(self.open["Plus"][3])} opened ORFs')


class ORF:
    def __init__(self, id, startPos, codon: SeqRecord, strand: str, parent= None):
        self.id = id
        self.strand = strand
        self.codSeq = codon
        self.startPos = startPos
        self.endPos = None
        if isinstance(parent, (ORF, type(None))):
            self.parent = parent
        else:
            raise TypeError("parent must be ORF object or None")

    def appendCodon(self, codon):
        self.codSeq += codon

    def closeSeq(self, end_pos):
        self.endPos = end_pos
        self.seq = self._codToSeq()


    def _codToSeq(self):
        seq=""
        return [seq + str(x) for x in self.codSeq]

fasta = ORFExtractor(fastaFile)
fasta.smartExtract()