#!/usr/bin/perl

###### Print USAGE if no argument is given
#my $usage = "\nUSAGE: Extract_concensus_based_on_filter_and_unfilter_vcf.pl <input_reference_fasta_file> <input_unfilter_vcf_file> <input_filter_vcf_file>

###### Get the arguments:  
my $input_reference_fasta_file = shift; # or die $usage;
my $input_unfilter_vcf_file = shift; # or die $usage;
my $input_filter_vcf_file = shift;
my $input_bam_file = shift;
my $input_filetr_bam_file = shift;
my $output_concensus = shift;

my $ST_start_pos = 367;
my $ST_end_pos = 718;

my $depth_cutoff = 3;

open (CONFIG, "config.txt") or die;
local $/;
my $config = <CONFIG>;
my ($samtools_path) = $config =~ /samtools_path\s*=\s*(.*?)\n/;
local $/ = "\n";

my $ST_length = $ST_end_pos-$ST_start_pos;
my $output_fasta = $output_concensus.".fasta";
my $output_ST = $output_concensus.".ST_MompS_res.txt";
open (RES, ">$output_ST") or die;
open (FASTA, ">$output_fasta") or die;

#check for regunes without reads mapping coverage at the non-filter bam fila
#my $pnd '$samtools_path/samtools depth -a -r Paris_mompS_R:367-718 $input_bam_file | grep "0$" | wc -l';
my $any_cov_nonfiletr = qx($samtools_path/samtools depth -a -r Paris_mompS_R:367-718 $input_bam_file);
my $num_of_zero_nonfiletr = qx($samtools_path/samtools depth -a -r Paris_mompS_R:367-718 $input_bam_file | grep '\\s0\$' | wc -l );
my $num_of_one_nonfiletr = qx($samtools_path/samtools depth -a -r Paris_mompS_R:367-718 $input_bam_file | grep '\\s1\$' | wc -l );
print "$samtools_path/samtools depth -a -r Paris_mompS_R:367-718 $input_bam_file | grep '\\s1\$' | wc -l\n";
my $num_of_two_nonfiletr = qx($samtools_path/samtools depth -a -r Paris_mompS_R:367-718 $input_bam_file | grep '\\s2\$' | wc -l );
print "depth: 0 - $num_of_zero_nonfiletr, 1 - $num_of_one_nonfiletr, 2 - $num_of_two_nonfiletr\n";
if($num_of_zero_nonfiletr>0 or $num_of_one_nonfiletr>0 or $num_of_two_nonfiletr>0 or $any_cov_nonfiletr eq ""){
	print RES "!!!Too low coverage at non filter file at file $input_bam_file !!!!!!!";
}

#check for regunes without reads mapping coverage at the filter bam fila
my $flage_coverage = 0;

#read the reference
my $ref_string = "";
open (REF, $input_reference_fasta_file) or die;
while (my $line = <REF>){
	chomp $line;
	if($line =~ /^>/){
		next;
	}else{
		$ref_string = $ref_string.$line;
	}
}
my $befor_ref = substr($ref_string,$ST_start_pos-1 ,$ST_length+1);
print "befor:\n$befor_ref\n";
#read the unfilter vcf
my $unfiletr_hash;
my $unfiletr_hash_ALT;
my $filetr_hash_ALT;
my $filter_hash;
open (UNFILT_VCF, $input_unfilter_vcf_file) or die;
while (my $line = <UNFILT_VCF>){
	chomp $line;
	if($line =~ /^#/){
		next;
	}else{
		my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, $info) = split("\t",$line);
		if ($POS> $ST_start_pos and  $POS < $ST_end_pos){
			my ($zigot) = $INFO =~ /AB=(.*?);/;
		
		#homozigot mutation - take it
			if($zigot == 0){
				my $alt_length = length($ALT);
				print substr($ref_string,$POS-1,$alt_length), " have to replace by $ALT\n";
				substr($ref_string,$POS-1,$alt_length) = $ALT;
				print substr($ref_string,$POS-1,$alt_length), "\n";
				print "$POS\t$REF\t$ALT\tAB = $zigot\t$alt_length\n";
			}
			else{#an hetrozigot mutation found
				my $depth_in_hetro_pos_term = qx($samtools_path/samtools depth -a -r Paris_mompS_R:$POS-$POS $input_filetr_bam_file);
				my ($depth_in_hetro_pos) = $depth_in_hetro_pos_term =~ /\d+\s+(\d+)/;
				if($depth_in_hetro_pos<$depth_cutoff){ #first, check the coverage of the filter (board) bam file.
					$flage_coverage = 1;
					print RES "\n The coverage at point $POS is $depth_in_hetro_pos\n";
				}
				print RES "Avidamce to hetrozigot mutation: posotin $POS, zigot: $zigot\n";
				print "hetro - $POS\n";
				
				if ($ALT =~ /,/){  
					$unfiletr_hash->{$POS} = 0;
					$unfiletr_hash_ALT->{$POS} = ",";
					print "There is , in the ALT\n";
				}else{
				$unfiletr_hash->{$POS} = $zigot;
				$unfiletr_hash_ALT->{$POS} = $ALT;
				}
			}
		}
	}
}

if($flage_coverage == 1){ #first, check the coverage of the filter (board) bam file.

	print RES "\n!!!Too low coverage at filter bam file! cannot trust the results!!!!!\n\n";
}

my @hetro_pos = sort keys %{$unfiletr_hash};
if (@hetro_pos != NULL){
	print RES "there are hetrozigot positions!\n";
	open (FILT_VCF, $input_filter_vcf_file) or die;
	while (my $line = <FILT_VCF>){
		chomp $line;
		if($line =~ /^#/){
			next;
		}else{
			my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, $info) = split("\t",$line);
			my ($zigot) = $INFO =~ /AB=(.*?);/;
			$filter_hash->{$POS} = $zigot;
			$filetr_hash_ALT->{$POS} = $ALT;
		}
	}
	my $alt_length;
	foreach my $i (@hetro_pos){
		print "exists $i\n";
		if (exists $filter_hash->{$i}){
		print RES "zigot befor clean: $unfiletr_hash->{$i}, zigot after clean: $filter_hash->{$i}\n";
			if ($unfiletr_hash_ALT->{$i} eq "," and exists $filter_hash->{$i}){ #if there are hetrozigot SNP and both are diferent from the reference, it will take the SNP as it appear in the filter file
				$alt_length = length($filetr_hash_ALT->{$i});
				print substr($ref_string,$i-1,$alt_length), " have to replace by $filetr_hash_ALT->{$i}\n";
				substr($ref_string,$i-1,$alt_length) = $filetr_hash_ALT->{$i};
				print substr($ref_string,$i-1,$alt_length), "\n";
				print "$i\t$filetr_hash_ALT->{$i}\tAB = $filter_hash->{$i} > $unfiletr_hash->{$i}\n";
			}
			elsif ($filter_hash->{$i} > $unfiletr_hash->{$i} or $filter_hash->{$i} == 0){
				$alt_length = length($unfiletr_hash_ALT->{$i});
				print substr($ref_string,$i-1,$alt_length), " have to replace by $unfiletr_hash_ALT->{$i}\n";
				substr($ref_string,$i-1,$alt_length) = $unfiletr_hash_ALT->{$i};
				print substr($ref_string,$i-1,$alt_length), "\n";
				print "$i\t$unfiletr_hash_ALT->{$i}\tAB = $filter_hash->{$i} > $unfiletr_hash->{$i}\n";
			}
		}
	}
	
}



my $after_ref = substr($ref_string,$ST_start_pos-1 ,$ST_length+1);
print RES ">momps\n$after_ref\n";
print FASTA ">momps\n$after_ref\n";


