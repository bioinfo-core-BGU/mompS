#!/usr/bin/perl
use Getopt::Long;

GetOptions
	(
	"f=s"=>\$reads_F,# forwad fastq file name
	"r=s"=>\$reads_R,# reverse fastq file name
	"p=s"=>\$prefix,# prefix
	"o=s"=>\$outputpath,# output path
	"a=s"=>\$assembly_fasta # asembly genome 
	);
#This script combine the mompS pipe and the MLST of the rest 6 genes

my $output_blast = $outputpath.$prefix.".BLAST_res.txt";
my $output_MLST = $outputpath.$prefix.".MLST_res.txt";
my $log = $outputpath.$prefix.".log";


open (CONFIG, "config.txt") or die;
local $/;
my $config = <CONFIG>;
my ($blast_path) = $config =~ /blast_path\s*=\s*(.*?)\n/;
open (RES, ">$log") or die;
#call the mompS pipe
`perl mompS_2_pipe_new.pl $reads_F $reads_R $prefix $outputpath`;
open (mompS,$outputpath.$prefix.".consensus.ST_MompS_res.txt") or die;
local $/;
my $mompS_res = <mompS>;
my ($mompS_ST) = $mompS_res =~ /mompS ST number: (\d+)/;
print RES "$mompS_ST\n";

#make blastDB for assembly_fasta
`$blast_path/makeblastdb -in $assembly_fasta -dbtype nucl`;
print RES "$blast_path/makeblastdb -in $assembly_fasta -dbtype nucl\n";
#do tha blast of representors MLST without mompS againse the assembly_fasta
`$blast_path/blastn -query MLST_Representative_Legionella.fasta -db $assembly_fasta -outfmt "6 qseqid sallseqid qlen slen qstart qend sstart send length evalue bitscore score pident qframe" -out $output_blast`;
print RES "$blast_path/blastn -query MLST_Representative_Legionella.fasta -db $assembly_fasta -outfmt 6 qseqid sallseqid qlen slen qstart qend sstart send length evalue bitscore score pident qframe -out $output_blast\n";
#call the MLST pipe
`perl Extract_MLST_Legionella_new_mompS.pl $assembly_fasta $output_blast $mompS_ST $output_MLST`;
print RES "perl Extract_MLST_Legionella_new_mompS.pl $assembly_fasta $output_blast $mompS_ST $output_MLST\n";
