import logging
from ORFinder import ORFExtractor
import FileFormats


# define file variables
fastaFileName = "Homo_sapiens_cdna_assembed.fasta"
outputFileName = "/analysis/output/ORF_extracted"
blastoutput = "/analysis/output/blast_results.tsv"
blastdb = "/db/swissprot/swissprot"

# Extract ORF from fasta source file
fasta_to_ORF = ORFExtractor(fastaFileName)
ORF_Dict = fasta_to_ORF.extract()

# run BLASTx against putative ORF extracted
fasta_to_ORF.blastx(database=blastdb, outputfile=blastoutput, taxids=9606, evalue=0.05, outfmt=7, max_target_seqs=1 )

# Validate ORF found and annotate source sequences
fasta_to_ORF.orf_validate()

# Assess the quality of the extraction
fasta_to_ORF.orf_assessment()

