#!/usr/local/bin/perl -w

### SRS378444_2.fastq (SRR999007) is Sanger code quality
### SRS1615359 is P898012; SRS1615291 is Tx430; SRS1523163 is Tx623

###
##	ID	               Region_V313	          Alias
##	Sobic.001G152901	Chr01:12171177..12377291	Sh1
##	Sobic.002G076600	Chr02:7875937..8085221	Tan2
##	Sobic.004G280800	Chr04:62215396..62418779	Tan1
##	Sobic.006G004400	Chr06:597459..800101	Ma6/Ghd7
##	Sobic.006G047700	Chr06:33448755..33651604	FT12
##	Sobic.006G057866	Chr06:40204883..40416799 	PRR37
##	unknown          	Chr07:56232000..56528079	qHT
##	Sobic.007G163800	Chr07:59721905..59929910	Dw3
##	Sobic.009G229800	Chr09:56938653..57141166	Dw1
##	Sobic.010G045100	Chr10:3399965..3602278	FT
##	Sobic.009G259100	Chr09:59163972..59367819	CCT
##	Sobic.009G257300	Chr09	59055166..59260131	ELF3


# for Ma1 and FT
#my $target_chr = "'".'@SQ\|Chr06\|Chr10'."'";
#my $target_region = 'Chr06:40205117-40416799 Chr10:3399965-3602278'; ##  Ma1(Sobic.006G057866): Chr06:40305117..40316799; FT: Chr10:3499965..3502278 forward 
###my $awk_std = "'".'$3~/Chr06/&&$4>40205117&&$4<40416799||$3~/Chr10/&&$4>3399965&&$4<3602278'."'";            
##                      sh1                                 Tan2                                     Tan1                                      Ma6/Ghd7                       FT12
#my $awk_std1 = '$3~/Chr01/&&$4>12171177&&$4<12377291||$3~/Chr02/&&$4>7875937&&$4<8085221||$3~/Chr04/&&$4>62215396&&$4<62418779||$3~/Chr06/&&$4>597459&&$4<800101||$3~/Chr06/&&$4>33448755&&$4<33651604||';
####                    PRR37                                qHT7                                Dw3                                      Dw1                                FT
#my $awk_std2 = '$3~/Chr06/&&$4>40204883&&$4<40416799||$3~/Chr07/&&$4>56232000&&$4<56528079||$3~/Chr07/&&$4>59721905&&$4<59929910||$3~/Chr09/&&$4>56938653&&$4<57141166||$3~/Chr10/&&$4>3399965&&$4<3602278';
#my $awk_std = "'".$awk_std1.$awk_std2."'";

## for ch4 and 8
##my $target_chr = "'".'@SQ\|Chr04\|Chr08'."'";
##my $target_region = 'Chr04:40205117-40416799 Chr10:3399965-3602278'; ##  Ma1(Sobic.006G057866): Chr06:40305117..40316799; FT: Chr10:3499965..3502278 forward 
#my $awk_std = "'".'$3~/Chr04/&&$4>55000000&&$4<65000000||$3~/Chr08/&&$4>50000000&&$4<60000000'."'";            
### for ch8 CO population 
#my $awk_std = "'".'$3~/Chr08/&&$4>58000000&&$4<60000000'."'";            

### for CCT in chr9
#my $awk_std = "'".'$3~/Chr09/&&$4>57000000&&$4<60000000'."'";

### for ch6
#my $awk_std = "'".'$3~/Chr07/'."'";


use strict;
use Net::FTP;
my $bwa = 'bwa';
my $btrim  = 'Btrim64';
my $PE_pl = 'paired_end_trim.pl';
#my $filter_pl = '/home/lixr/PollenSeq/scripts/filter_unique_reads_bwa.pl';
my $genome_ref_4_bwa = '/home/lixr/Sb/Refs/Sb_V3_all_bwa';
#my $genome_ref_4_samtools = '/work/agron_1/Sb/Refs/Sb_V2_all.fa';
my $pre_dir = '/home/lixr/Sb/qHT7/';

###SRR_list_by_SRX();
##SRR_list_by_SRS();
my $split_lines = 8000000; #8000000; ###8_000_000; ## total of BW350-1_1 is 353_607_940. 45 files 
my $split_filter = '--filter='."'".'gzip > $FILE.gz'."'";

my $awk_std = "'".'$3~/Chr07/&&$4>56430000&&$4<56475000'."'";

my $srr_list_file = '/home/lixr/Sb/qHT7_otherpops';
#my $srr_list_file = '/home/lixr/Sb/IS3620C_wSRRs';
#my $bwa_group = 3;
for (my $bwa_group = 1; $bwa_group <= 4; $bwa_group ++) {
	next unless $bwa_group == 4;
	my ($sra_arrayref, $sra_hashref) = SRA_list_DIR($bwa_group, $srr_list_file);
	my @srs = @$sra_arrayref; #keys %$sra_hashref;
	my $sam_header_file = '/home/lixr/Sb/Refs/Sb_SAM_header';
	my $local_SRA_chunk_dir = '/home/lixr/Sb/SRAs/';
	for (my $i = 0; $i <= $#srs; $i ++) { ## 0 - 3; 4 - 6;
		my $sample = $srs[$i];
		my $sample_dir = $pre_dir.$sample.'/'; 
		mkdir $sample_dir unless -e $sample; 
		my $merged_bam = $pre_dir.$sample.'.bam';
		next if -e $merged_bam;
		
		my @remote_files = @{ $$sra_hashref{$sample}};
		my $srs = $remote_files[0]; ##shift (@remote_files); 
#		print $srs."\n";
		my ($pre_srs, $x) = $srs =~ /(SRR\d\d\d)(\d+)/;
		my $remote_ftp_pre = 'ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/'.$pre_srs.'/'; 
	
#	my @remote_files = qw/SRR2759161/;
		my $k; 
		foreach my $run (@remote_files) {
			$k ++;
#			next if $k > 1;
			my $strain = $sample.'_'.$k;
			my $run_dir = $sample_dir.$run.'/';
			mkdir $run_dir unless -e $run_dir;
			my $ori_run_sra_file = $run_dir.$run.'.sra';
			my $ln_run_sra_file = $sample_dir.$strain.'.sra';
#			
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
#				system("split -l $split_lines -d -a 3 $split_filter $left_file $left_chunk_prex");
				system("rm $left_gz_file");
				system("gunzip -c $right_gz_file | split -l $split_lines -d -a 3 $split_filter - $right_chunk_prex");
#				system("split -l $split_lines -d -a 3 $split_filter $right_file $right_chunk_prex");
				system("rm $right_gz_file");
			}
			my $quality_code = Quality_Code($run_quality_file);
			for (my $i = 50; $i < 60; $i ++) {
#				next unless $i == 10;
#				print $i." ";
				my $suffix = sprintf "%03d", $i;
				my $left_chunk_file = $left_chunk_prex.$suffix.'.gz';
#				my $left_chunk_file = $local_SRA_chunk_dir.$sample.'/'.$strain.'_1_'.$suffix.'.gz';
				next unless -e $left_chunk_file;
				my $right_chunk_file = $right_chunk_prex.$suffix.'.gz';
#				my $right_chunk_file = $local_SRA_chunk_dir.$sample.'/'.$strain.'_2_'.$suffix.'.gz';
				
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
#				unlink $left_chunk_file if -e  $chunk_uni_srt_bam;
#				unlink $right_chunk_file if -e  $chunk_uni_srt_bam;
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


sub SRR_list_by_SRX {
#	my ($bwa_group) = @_;
#	my (%hash, @array);
	open (F, '/home/lixr/Sb/NC_SRAs') || die;
	open (O, '>/home/lixr/Sb/NC_SRAs_wSRRs') || die;
	my $ftp_site = 'ftp-trace.ncbi.nih.gov';
	my $ftp = Net::FTP->new($ftp_site);
     $ftp->login();
	while (<F>) {
		chomp;
		my $line = $_;
		my @t = split /\t/, $line;
		next if $t[0] eq 'bwa_Group';
##		  next unless $#t == 1;
#		next unless $t[-1] == $bwa_group;  
#		push @array, $t[1];
	my $srx = $t[2];
		my ($pre_srs, $x) = $srx =~ /(SRX\d\d\d)(\d+)/;
		my $remote_ftp_dir = '/sra/sra-instant/reads/ByExp/sra/SRX/'.$pre_srs.'/'.$srx;
		  $ftp->cwd($remote_ftp_dir);
		my @remote_files = $ftp->ls();
		print O $line;
		print O "\t".$_ foreach (@remote_files);
		print O "\n";
		}
	close F;
#	return (\@array, \%hash);	
	   $ftp->quit;

	}
	
sub SRR_list_by_SRS {
#	my ($bwa_group) = @_;
#	my (%hash, @array);
	open (F, '/home/lixr/Sb/SbNAM_founders_1') || die;
	open (O, '>/home/lixr/Sb/SbNAM_founders_wSRRs') || die;
	my $ftp_site = 'ftp-trace.ncbi.nih.gov';
	my $ftp = Net::FTP->new($ftp_site);
     $ftp->login();
	while (<F>) {
		chomp;
		my $line = $_;
		my @t = split /\t/, $line;
		next if $t[0] eq 'bwa_Group';
##		  next unless $#t == 1;
#		next unless $t[-1] == $bwa_group;  
#		push @array, $t[1];
	my $srx = $t[1];
		my ($pre_srs, $x) = $srx =~ /(SRS\d\d\d)(\d+)/;
		my $remote_ftp_dir = '/sra/sra-instant/reads/ByStudy/sra/SRS/'.$pre_srs.'/'.$srx;
		  $ftp->cwd($remote_ftp_dir);
		my @remote_files = $ftp->ls();
		print O $line;
		print O "\t".$_ foreach (@remote_files);
		print O "\n";
		}
	close F;
#	return (\@array, \%hash);	
	   $ftp->quit;

	}	