NRS pipeline identifies novel segments of a primary assembly not aligning to the reference genome
and compares the variant calling outputs of two pipelines; one using the reference and the other using the reference modified with novel segments
Detailed description of the project and other open-source tools utilized in the project can be found at: 
https://www.biorxiv.org/content/10.1101/2023.10.23.563520v1

Scripts in this project in running order: \
unaligned_contigs.R \
lift_database.R \
annotation lifts: lift_RNA.R,lift_Refseq.R, lift_annovar.R \
LostGained.R \

unaligned_contig.R:
uses MUMmer alignment output to extract aligned sequences. 
MUMmer tools: https://mummer.sourceforge.net/manual/
input is generated using show-coords tool of MUMmmer
show-coords -T -l path/delta.filter > path/delta.filter.tab
running R
source("path/unaligned_contig.R")
F("path/delta.filter.tab")
outputs: output_alignment.csv and output_alignment_summary.csv  

For reference editing with NRS segments please refer to: https://systemsbio.ucsd.edu/perEditor/

lift_database.R:
Needed to perform variant calling using edited genome reference with GATK.
Lifts over the SNP database in GATK bundle.
input files: dbsnp file (header removed), list of segments inserted into the genome reference.
source("path/lift_database.R")
lift_database("path/dbsnp_138.hg38.sort.filtered.no_header.vcf","path/insertions_all_chr.txt”)

Three scripts to lift-over databases used by ANNOVAR
ANNOVAR: https://annovar.openbioinformatics.org/en/latest/
lift_RNA.R:
lift_RNA("path/hg38_refGeneMrna.fa", "path/insertions_all_chr.xlsx")
lift_RefSeq.R:
lift_RefSeq("path/hg38_refGene.txt", "path/insertions_all_chr.xlsx", "path/RefSeq_changes.xlsx")
RefSeq_changes.xlsx describes the genes with insertions. The file used in project is available in the folder.
lift_annovar.R:
lift_annovar("path/hg38_avsnp150.txt", "path/insertions_all_chr.xlsx")

LostGained.R:
Compared the the reference and reference-NRS annotated variant data.
input files: chromosome based annotation files for reference and reference-NRS, insertion files. add chromosome number as an argument
find_gained_lost_snps("path/chr10_hg38_multianno.txt", "path/chr10_hg38NRS_multianno.txt", "path/insertions_all_chr.xlsx", "10")

Authors: Cigdem Koroglu, Serdar Altok



