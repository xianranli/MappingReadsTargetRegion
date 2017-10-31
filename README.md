# MappingReadsTargetRegion
Automatically downloading NGS reads from SRA, mapping to the reference and extracting the reads mapped to the targets.
The purpose of this script is to focus the reads mapped at the desired regions, instead of whole genome, to discover potential polymorphisms between among the interested strains and reference. This is helpful for identifying functional site in the narrowed QTL interval.
The rationale behind this script is a) the wet-bench oriented people may have limited knowledge to run NGS, 2) their compute ring resources is limited (memory, disk space), 3) but the time is not a big concern.
because of the large size of sra and fastq files, after downloading sra from NCBI SRA, this script balance the time and computing resources by split the whole fastq file into small chunks, then each chunk will be mapped independently, which requires less memories, the reads mapped to target region will be extracted with awk then stored as bam file.
After all chunks are finished, the bam file from each chunk will be merged into one file and indexed.
Currently only support in Unix/Linux environment (with wget and awk built-in the system)


Prerequisites
1. bwa https://github.com/lh3/bwa
2. Btrim64 and paired_end_trim.pl from http://graphics.med.yale.edu/trim/
3. indexed reference by bwa 

Required input file
---------------------------
This example file try to map reads from 4 strains
bwa_Group	Strain	SRR

1	80M	SRR5271055

2	Hegari	SRR5271058

3	IS3620C	SRR5271059

4	BTx623	SRR5271056

---------------------------

Usage:
perl bp_run_bwa_mem_TP.pl project_dir srr_list reference sam_header

For example, assuming the script is stored at /home/xxx/Test/scripts/.

a). The reference is stored in /home/xxx/refs/Test_ref_bwa;

b). The project directory is /home/xxx/Test/;

c). The srr_list file is /home/xxx/Test/srr_list;

d). The SAM header file is /home/xxx/Test/sam_header


perl /home/xxx/Test/scripts/bp_run_bwa_mem_TP.pl /home/xxx/Test/ /home/xxx/Test/srr_list /home/xxx/refs/Test_ref_bwa /home/xxx/Test/sam_header

The final combined bam file, for example, will be stored  under the project fold, '/home/xxx/Test/80M.bam' for the 1st strain, '/home/xxx/Test/Hegari.bam' for the 2nd strain.

