# MappingReadsTargetRegion
Automatic downloading NGS reads from SRA, mapping to the reference and extracting the reads mapped to the targets.
The purpose of this script is to focuse the reads mapped at the desired regions, instead of whole genome, to discover potential polymorphisms against the reference. This is specially helpful for identifying functional site in the narrowed QTL interval.

The rational is to balance the time and computing resources (memory, disk space et al) because of the large size of sra fand fastq files, after downloading sra from NCBI SRA. The whole fastq file will be split into small chunks, then each chunk will be mapped indepently.
After all chunks are finished, the bam file from each chunk will be merged into a one file and indexed. 

Prerequisites
1. bwa https://github.com/lh3/bwa
2. Btrim64 and paired_end_trim.pl from http://graphics.med.yale.edu/trim/
3. indexed reference by bwa 

Only required input file
---------------------------
bwa_Group	Strain	SRR

1	80M	SRR5271055

2	Hegari	SRR5271058

3	IS3620C	SRR5271059

4	BTx623	SRR5271056

---------------------------
The final combined bam file, for example, will be '80M.bam' for the 1st accession, 'Hegari.bam' for the 2nd accession.
