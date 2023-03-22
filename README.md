# KCNN1_ewing_sarcoma

## Chipseq data

ChIPseq data are from GSE176400, GSE129155, and GSE94278 from the Gene Expression
Omnibus database. GSE176400 contains ChIPseq from 18 Ewing's sarcoma cell lines targeting FLI1 or
ERG based on cell line fusion; histone marks H3K4me3, H3K27me3, and H3K27ac; and an input
serving as a control. GSE129155 contains ChIPseq targeting the H3K27ac mark of the ASP14 cell line
without sh induction targeting EWS-FLI1 (day 0) and after 7 days of sh induction, and ChIPseq
targeting FLI1 on the same cell line starting at the end of sh induction (day 7), and then on days 9, 10,
11, 14, and 17; in order to observe the effects of gradual reappearance of EWS-FLI1. These data are
supplemented with an input from day 0 that served as a control. GSE94278 contains Bigwig files of
normalized ChIPseq values of Mesenchymal Stem Cells (MSCs) transfected with a control plasmid,
FLI1 plasmid, and EWS-FLI1 plasmid.
A few corrupted reads (nucleotide sequence length different from the phred score) were
removed using Awk 4.0.2 (https://doi.org/10.1002/spe.4380090403). ChIPseq data were all filtered
with fastp 0.23.2 (https://doi.org/10.1093/bioinformatics/bty560) by removing reads with more than
40% of bases with a phred score below 30 and removing PCR adapters if detected. Reads were aligned
with BWA-MEM 0.7.17 (https://doi.org/10.1093/bioinformatics/btp324) against the GRCh38.p14
reference genome "analysis set" version from NCBI. After conversion to BAM format and sorting of
reads according to their coordinates with samtools 1.16.1
(https://pubmed.ncbi.nlm.nih.gov/33590861/), peaks were detected with Macs 2.2.7.1
(https://doi.org/10.1186/gb-2008-9-9-r137) using the controls associated with each lineage. Sequencing
depths of the alignment files to construct the graphs were retrieved through python scripts using Pysam
0.20.0 (https://github.com/pysam-developers/pysam). Macs2 sequencing depths and peaks coordinates
plots were constructed with ggplot2 3.4.0 (https://doi.org/10.1007/978-0-387-98141-3), and transcript
representations with gggenes 0.4.1 (https://github.com/wilkox/gggenes) from the RefSeq GTF file
associated with the reference genome used. These graphs were assembled with Cowplot 1.1.1
(https://github.com/wilkelab/cowplot). As for the number of microsatellite repeats at the EWS-FLI1 or EWS-
ERG binding site, the consensus sequence of the region was obtained with bcftools mpileup 1.16
(https://pubmed.ncbi.nlm.nih.gov/33590861/) by keeping only mutations present at least 5 times, then
the counting was done by TandemRepeatFinder 4.09.1 (https://academic.oup.com/nar/article/27/2/573/1061099) 
by strongly penalizing the mismatches and indels.

