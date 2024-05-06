NRS pipeline identifies novel segments of a primary assembly not aligning to the reference genome
and compares the variant calling outputs of two pipelines; one using the reference and the other using the reference modified with novel segments

Detailed description of the project and other open-source tools utilized in the project can be found at: 
https://www.biorxiv.org/content/10.1101/2023.10.23.563520v1

Scripts containing the functions in this project in running order: 
- unaligned_contigs.R 
- lift_database.R 
- annotation lifts: lift_RNA.R, lift_Refseq.R, lift_annovar.R 
- LostGained.R

Files used:
- insertions file: contains locations of novel segments inserted in the reference
- dbsnp file: dbsnp file from GATK bundle
- ANNOVAR database files: hg38_refGeneMrna.fa, hg38_refGene.txt, hg38_avsnp150.txt

unaligned_contig.R:
- Uses MUMmer alignment output to extract unaligned sequences. 
- MUMmer tools: https://mummer.sourceforge.net/manual/ 
- Input is generated using show-coords tool of MUMmer: \
&emsp;show-coords -T -l path/delta.filter > path/delta.filter.tab
- and then running the following R function: \
&emsp;find_unaligned_parts(mummer_output_file)     
- Output files: \
&emsp;output_unalignment.csv and output_unalignment_summary.csv  

For reference editing with NRS segments please refer to: https://systemsbio.ucsd.edu/perEditor/

lift_database.R:
- Needed to perform variant calling using edited genome reference with GATK.
- Lifts over the SNP database in GATK bundle.
- Sample call: \
&emsp; lift_database(dbsnp_file, insertions_file)

Three functions to lift-over databases used by ANNOVAR (https://annovar.openbioinformatics.org/en/latest/): \
- lift_RNA.R: \
&emsp;lift_RNA(hg38_refGeneMrna, insertions_file)
- lift_RefSeq.R: \
&emsp;lift_RefSeq(hg38_refGene, insertions_file)
- lift_annovar.R: \
&emsp;lift_annovar(hg38_avsnp150, insertions_file)

LostGained.R: 
- Compared the reference and the reference-NRS annotated variant data.
- Input files: \
&emsp; Insertion file, chromosome based annotation files for reference and reference-NRS 
- find_gained_lost_snps(hg38_multianno_file, hg38NRS_multianno_file, insertions_file, chromosome_number)

Authors: Cigdem Koroglu, Serdar Altok



