from matplotlib import pyplot as plt
from ORFinfder import ORFExtractor


fastaFileName = "Homo_sapiens_cdna_assembed.fasta"
outputFileName = "/analysis/output/ORF_extracted"

fasta_to_ORF = ORFExtractor(fastaFileName)
ORF_Dict = fasta_to_ORF.extract()
fasta_to_ORF.result_export(outputFileName)

max_length_dict = {}
for key, tup in ORF_Dict.items():

    if tup[2] in max_length_dict.keys():
        if tup[1] > max_length_dict[tup[2]]:
            max_length_dict[tup[2]] = tup[1]
            max_ORF = tup
    else:
        max_length_dict[tup[2]] = tup[1]
global_length_dict = {}
for key, tup in ORF_Dict.items():
    global_length_dict[key] = tup[1]

print(f"Contig(s) missing ORF : {fasta_to_ORF.non_orf_seq}")
# print(f'maximum ORF length : {max(max_length_dict, key=max_length_dict.get)} {max(max_length_dict.values())}')
# print(sorted(ORF_Dict.items(), key=lambda i: i[1][5])[-10:])
print(sorted(max_length_dict.items(), key=lambda i: i[1])[-5:])
print(f'number orf ORF : {len(ORF_Dict.keys())}')
plt.hist(x=global_length_dict.values(), bins=range(1000))
plt.savefig("/analysis/output/orf_length_distribution")

# command for blastx : blastx -db /db/swissprot -query /analysis/output/ORF_extracted.fasta -taxids 9606 -evalue 0.05 -outfmt 7 -out /analysis/output/blast_results.tsv


