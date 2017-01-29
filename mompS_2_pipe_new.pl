#!/usr/bin/perl

my $fastq1 = shift;
my $fastq2 = shift;
my $prefix = shift;
my $output_Path = shift;

my $fasta_of_ref_momps_2_regen = "Ref_Paris_mompS_2.fasta";
my $log = $output_Path.$prefix.".log";
open (RES, ">>$log") or die;


my $output_sam = $prefix.".sam";
my $output_sam_header_fix = $prefix."_headerFix.sam";
my $output_bam = $prefix.".bam";
my $output_dup = $prefix.".dedup.bam";
my $output_proper_paireds = $prefix.".properPaired.sam";
my $output_board_resds = $prefix.".board.sam";
my $output_consensus = $prefix.".consensus";
my $vcf_unfilter = $prefix.".unfilter.vcf";
my $vcf_filter = $prefix.".filter.vcf";

open (CONFIG, "config.txt") or die;
local $/;
my $config = <CONFIG>;
my ($bwa_path) = $config =~ /bwa_path\s*=\s*(.*?)\n/;
my ($samtools_path) = $config =~ /samtools_path\s*=\s*(.*?)\n/;
my ($picard_path) = $config =~ /picard\s*=\s*(.*?)\n/;
my ($freebayse_path) = $config =~ /freebayse_path\s*=\s*(.*?)\n/;


#Mapping fastq paired file to mompS_2_regen by bwa
`$bwa_path/bwa mem Ref_Paris_mompS_2.fasta -R '\@RG\tID:group1\tSM:sample81\tPL:illumina\tLB:Technion\tPU:unit1'  $fastq1 $fastq2 > $output_Path$output_sam`;
#reheater sam file
my $grep_io = qx(grep -v '\@PG' $output_Path$output_sam > $output_Path$output_sam_header_fix);

#Sam to sort bam
`$picard_path SortSam I=$output_Path$output_sam_header_fix O=$output_Path$output_bam SORT_ORDER=coordinate`;
#remove duplicales
`$picard_path MarkDuplicates M=matrix_ncbi.txt I=$output_Path$output_bam O=$output_Path$output_dup`;
#index
`$picard_path BuildBamIndex I=$output_Path$output_dup`;

#filter proper paireds
`$samtools_path/samtools view -h -f 0x2 $output_Path$output_dup > $output_Path$output_proper_paireds`;

#select only boards reads
`perl Filter_mompS2_boards_reads.pl $output_Path$output_proper_paireds $output_Path$output_board_resds`;

#Sam to sort bam, index
`$picard_path SortSam INPUT=$output_Path$output_board_resds OUTPUT=$output_Path$output_board_resds".bam" SORT_ORDER=coordinate`;
`$picard_path BuildBamIndex I=$output_Path$output_board_resds".bam"`;

#freebayes vcf for filtered and unfiltered bams
`$freebayse_path/freebayes -f $fasta_of_ref_momps_2_regen $output_Path$output_dup > $output_Path$vcf_unfilter`;
`$freebayse_path/freebayes -f $fasta_of_ref_momps_2_regen $output_Path$output_board_resds".bam" > $output_Path$vcf_filter`;

#find the concensus
`perl Extract_concensus_based_on_filter_and_unfilter_vcf.pl $fasta_of_ref_momps_2_regen $output_Path$vcf_unfilter $output_Path$vcf_filter $output_Path$output_dup $output_Path$output_board_resds".bam" $output_Path$output_consensus`;
print RES "perl Extract_concensus_based_on_filter_and_unfilter_vcf.pl $fasta_of_ref_momps_2_regen $output_Path$vcf_unfilter $output_Path$vcf_filter $output_Path$output_dup $output_Path$output_board_resds .bam $output_Path$output_consensus\n";
#find the momps ST (fq to fasta, makeblastdb, blast mompS representor, extract seq due to blast res, find the fit ST)
`perl Find_mompS_ST_from_mompS_regen_fq.pl $output_Path$output_consensus`;
print RES "perl Find_mompS_ST_from_mompS_regen_fq.pl $output_Path$output_consensus\n";