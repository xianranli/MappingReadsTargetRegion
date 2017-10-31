# MappingReadsTargetRegion
Automatic downloading NGS reads from SRA, mapping to the reference and extracting the reads mapped to the targets.
The purpose of this script is to focuse the reads mapped at the desired regions, instead of whole genome, to discover potential polymorphisms against the reference. This is specially helpful for identifying functional site in the narrowed QTL interval.

The rational is to balance the time and computing resources (memory, disk space et al) because of the large size of sra fand fastq files, after downloading sra from NCBI SRA. The whole fastq file will be split into small chunks, then each chunk will be mapped indepently.
After all chunks are finished, the bam file from each chunk will be merged into a one file and indexed. 
