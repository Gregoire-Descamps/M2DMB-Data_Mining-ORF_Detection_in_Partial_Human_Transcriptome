import logging

import matplotlib
from matplotlib import pyplot as plt
from ORFinder import ORFExtractor
import FileFormats
matplotlib.set_loglevel("WARNING")

fastaFileName = "Homo_sapiens_cdna_assembed.fasta"
outputFileName = "/analysis/output/ORF_extracted"
# outputFileName = "output/ORF_extracted"
blastoutput = "/analysis/output/blast_results.tsv"
blastdb = "/db/swissprot/swissprot"


fasta_to_ORF = ORFExtractor(fastaFileName)
ORF_Dict = fasta_to_ORF.extract()
# fasta_to_ORF.blastx(database=blastdb, outputfile=blastoutput, taxids=9606, evalue=0.05, outfmt=7)
fasta_to_ORF.orf_validate()


# parser = FileFormats.BlastTsvFile(blastoutput)
#
# parser.top_hit_extract("output/test_top_hits2.tsv")

max_length_dict = {}
for key, orf in ORF_Dict.items():

    if orf.seq_id in max_length_dict.keys():
        if len(orf) > max_length_dict[orf.seq_id]:
            max_length_dict[orf.seq_id] = len(orf)
    else:
        max_length_dict[orf.seq_id] = len(orf)
global_length_dict = {}
for key, orf in ORF_Dict.items():
    global_length_dict[key] = len(orf)

print(f"Contig(s) missing ORF : {fasta_to_ORF.non_orf_seq}")
print(f'maximum ORF length : {max(max_length_dict, key=max_length_dict.get)} {max(max_length_dict.values())}')
# print(sorted(ORF_Dict.items(), key=lambda i: i[1][5])[-10:])
print(sorted(max_length_dict.items(), key=lambda i: i[1])[-5:])
print(f'number orf ORF : {len(ORF_Dict.keys())}')
plt.hist(x=global_length_dict.values(), bins=range(1000))
plt.savefig("/analysis/output/orf_length_distribution")
