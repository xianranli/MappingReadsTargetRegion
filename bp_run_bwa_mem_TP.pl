#!/usr/local/bin/perl -w

use strict;

my $pre_dir =  $ARGV[0]; ##'/home/lixr/Sb/qHT7/';
my $srr_list_file = $ARGV[1]; ##'/home/lixr/Sb/sra_info';
my $genome_ref_4_bwa = $ARGV[2];; ### '/home/lixr/Sb/Refs/Sb_V3_all_bwa';
## sam header file. futhure work: this file can be generated within pipline to reduce one argument passing to perl
my $sam_header_file = $ARGV[3]; ##'/home/lixr/Sb/Refs/Sb_SAM_header';

my $bwa = 'bwa';
my $btrim  = 'Btrim64';
my $PE_pl = 'paired_end_trim.pl';
my $split_lines = 8000000; #the number can be modified
my $split_filter = '--filter='."'".'gzip > $FILE.gz'."'";

### example of one region
my $awk_std = "'".'$3~/Chr07/&&$4>56430000&&$4<56475000'."'";
### example of multiple regions: Ch04:55000000-65000000 and Chr8: 50000000-60000000
#my $awk_std = "'".'$3~/Chr04/&&$4>55000000&&$4<65000000||$3~/Chr08/&&$4>50000000&&$4<60000000'."'";            

for (my $bwa_group = 1; $bwa_group <= 4; $bwa_group ++) {
	my ($sra_arrayref, $sra_hashref) = SRA_list_DIR($bwa_group, $srr_list_file);
	my @srs = @$sra_arrayref; #keys %$sra_hashref;
	for (my $i = 0; $i <= $#srs; $i ++) { ## 0 - 3; 4 - 6;
		my $sample = $srs[$i];
		my $sample_dir = $pre_dir.$sample.'/'; 
		mkdir $sample_dir unless -e $sample; 
		my $merged_bam = $pre_dir.$sample.'.bam';
		next if -e $merged_bam;
		
		my @remote_files = @{ $$sra_hashref{$sample}};
		my $srs = $remote_files[0]; ##shift (@remote_files); 
		my ($pre_srs, $x) = $srs =~ /(SRR\d\d\d)(\d+)/;
		my $remote_ftp_pre = 'ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/'.$pre_srs.'/'; 
		my $k; 
		foreach my $run (@remote_files) {
			$k ++;
#			next if $k > 1;
			my $strain = $sample.'_'.$k;
			my $run_dir = $sample_dir.$run.'/';
			mkdir $run_dir unless -e $run_dir;
			my $ori_run_sra_file = $run_dir.$run.'.sra';
			my $ln_run_sra_file = $sample_dir.$strain.'.sra';		
			my $left_file = $sample_dir.$strain.'_1.fastq';
			my $right_file = $sample_dir.$strain.'_2.fastq';
			my $left_gz_file = $left_file.'.gz';
			my $right_gz_file = $right_file.'.gz';
			my $left_chunk_prex = $sample_dir.$strain.'_1_';
			my $right_chunk_prex = $sample_dir.$strain.'_2_';
			my $run_quality_file = $sample_dir.$strain.'_fastq_sample';
			
			my $remote_ftp_file = $remote_ftp_pre.$srs.'/'.$run.'.sra';
			if ($k > 0) {
				system("wget -q -nH $remote_ftp_file -P $run_dir") unless -e $ori_run_sra_file;
				system("ln -s  $ori_run_sra_file $ln_run_sra_file") unless -e $ln_run_sra_file;
				system("fastq-dump -I -split-files --gzip $ln_run_sra_file -O $sample_dir") unless -e $left_file && -e $right_file;
				system("gzip -cd $left_gz_file | head -n 4000 > $run_quality_file") unless -e $run_quality_file;
				system("rm $ln_run_sra_file") ;
				system("rm $ori_run_sra_file -f");
				system("gunzip -c $left_gz_file | split -l $split_lines -d -a 3 $split_filter - $left_chunk_prex");
				system("rm $left_gz_file");
				system("gunzip -c $right_gz_file | split -l $split_lines -d -a 3 $split_filter - $right_chunk_prex");
				system("rm $right_gz_file");
			}
			my $quality_code = Quality_Code($run_quality_file);
			for (my $i = 0; $i < 1000; $i ++) {
#				next unless $i == 10;
				my $suffix = sprintf "%03d", $i;
				my $left_chunk_file = $left_chunk_prex.$suffix.'.gz';
				next unless -e $left_chunk_file;
#### preparing necessary temp files during the process. All temp files will be deleted to reduce the disk usage 				
				my $right_chunk_file = $right_chunk_prex.$suffix.'.gz';
				my $trim_left_chunk = $left_chunk_file.'.fq_trm';
				my $trim_left_chunk_s = $trim_left_chunk.'_sum';
				my $trim_right_chunk = $right_chunk_file.'.fq_trm';
				my $trim_right_chunk_s = $trim_right_chunk.'_sum';
				
				my $trim_PE_1_chunk = $trim_left_chunk.'.pe';
				my $trim_PE_1_chunk_s = $sample_dir.$strain.'_1_'.$suffix.'.pe_sum';
				
				my $trim_PE_2_chunk = $trim_right_chunk.'.pe';
				my $trim_PE_2_chunk_s = $sample_dir.$strain.'_2_'.$suffix.'.pe_sum';
				
				my $trim_SE_1_chunk = $trim_left_chunk.'.se';
				my $trim_SE_2_chunk = $trim_right_chunk.'.se';
				my $trim_SE_chunk = $sample_dir.$strain.'_0_'.$suffix.'.se';
			
				my $PE_target_bam = $sample_dir.$strain.'_pe_Target.bam_'.$suffix;
				my $SE_target_bam = $sample_dir.$strain.'_se_Target.bam_'.$suffix;
  	
				my $trim_PE_chunk_sam = $sample_dir.$strain.'_pe_all.sam_'.$suffix;
				my $trim_SE_chunk_sam = $sample_dir.$strain.'_se_all.sam_'.$suffix;
				
				my $trim_PE_chunk_uni_sam = $sample_dir.$strain.'_uni_pe.sam_'.$suffix;
				my $trim_SE_chunk_uni_sam = $sample_dir.$strain.'_uni_se.sam_'.$suffix;
			
				my $trim_PE_chunk_uni_bam = $sample_dir.$strain.'_pe_uni.bam_'.$suffix;
				my $trim_SE_chunk_uni_bam = $sample_dir.$strain.'_se_uni.bam_'.$suffix;
				my $chunk_uni_bam = $sample_dir.$strain.'_uni.bam_'.$suffix;
				
				my $chunk_uni_srt_bam = $sample_dir.$strain.'_Target_'.$suffix.'.bam';
				next if -e $chunk_uni_srt_bam;
				unlink $trim_PE_chunk_sam if -e $trim_PE_chunk_sam;
				unlink $trim_SE_chunk_sam if -e $trim_SE_chunk_sam;
				
				if ($quality_code eq 'I') {
					system("Btrim64 -i -q -t $left_chunk_file -Z -l 40 -o $trim_left_chunk -s $trim_left_chunk_s");
					system("rm $left_chunk_file");
					system("Btrim64 -i -q -t $right_chunk_file -Z -l 40 -o $trim_right_chunk -s $trim_right_chunk_s");
					system("rm $right_chunk_file");
					}
					elsif ($quality_code eq 'S') {
						system("Btrim64 -q -t $left_chunk_file -Z -l 40 -o $trim_left_chunk -s $trim_left_chunk_s");
						system("rm $left_chunk_file");
						system("Btrim64 -q -t $right_chunk_file -Z -l 40 -o $trim_right_chunk -s $trim_right_chunk_s");
						system("rm $right_chunk_file");
						}
				
				system("$PE_pl $trim_left_chunk_s $trim_right_chunk_s $trim_left_chunk $trim_right_chunk");
				system("rm $trim_left_chunk_s $trim_right_chunk_s $trim_left_chunk $trim_right_chunk");
			
				system("cat $trim_SE_1_chunk $trim_SE_2_chunk > $trim_SE_chunk");
				system("rm $trim_SE_1_chunk $trim_SE_2_chunk");
				
				system("bwa mem -t 2 -v 0 $genome_ref_4_bwa $trim_PE_1_chunk $trim_PE_2_chunk | awk $awk_std >> $trim_PE_chunk_sam");
				system("cat $sam_header_file $trim_PE_chunk_sam | samtools view - -q 10 -Sb | samtools sort - -o $PE_target_bam")	;
				system("rm $trim_PE_1_chunk $trim_PE_2_chunk $trim_PE_chunk_sam");
  		
				system("bwa mem -t 2 -v 0 $genome_ref_4_bwa $trim_SE_chunk | awk $awk_std >> $trim_SE_chunk_sam");
				system("cat $sam_header_file $trim_SE_chunk_sam | samtools view - -q 10 -Sb | samtools sort - -o $SE_target_bam")	;
				system("rm $trim_SE_chunk $trim_SE_chunk_sam");
  		
				system("samtools merge $chunk_uni_bam $PE_target_bam $SE_target_bam");
  		
				system("rm $PE_target_bam $SE_target_bam");
				
				system("samtools sort $chunk_uni_bam -o $chunk_uni_srt_bam");
				system("rm $chunk_uni_bam");
				
				}
			
			}
			my $target_bams = $sample_dir.'*Target*bam';	
			system("samtools merge $merged_bam $target_bams -f");	
			system("samtools index $merged_bam");
		}
	}

###########################################
sub SRA_list_DIR {
	my ($bwa_group, $f) = @_;
	my (%hash, @array);
	open (F, $f) || die;
	while (<F>) {
		chomp;
		my @t = split /\t/;
		next if $t[0] eq 'bwa_Group';
#		  next unless $#t == 1;
		next unless $t[0] == $bwa_group;  
		push @array, $t[1];
		  @{ $hash{$t[1]} } = @t[2..$#t];  
		}
	close F;
	return (\@array, \%hash);	
	}

sub Quality_Code {
	my ($f) = @_;
	my ($j, $flag, %hash, $code) ;
	open (F, $f) || die;
	while (<F>) {
		chomp;
		my $line = $_;
		$flag = 0 if $line =~ /^\+/;
		$flag ++; 
		$j ++;
		if ($flag) {
			my @t = split //, $line;
			foreach my $t(@t) {
				my $n = ord($t);
				  $hash{$n} ++;
				}
			}
		last if $j > 1000	
		}
	close F;	
	if (exists $hash{90}) { $code = 'I'} elsif (exists $hash{48}) {$code = 'S'};
	return $code;		
	}


