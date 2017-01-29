#!/usr/bin/perl

#my $input_fq_file = shift;
my $input_file = shift;

my $fasta_file = $input_file.".fasta";
my $outoutBlastFile = $fasta_file.".BlastMompS_res.txt";
my $outout_ST_File = $input_file.".ST_MompS_res.txt";
my $outout_log = $fasta_file.".log";
my $query = "MLST_mompS_Representative_Legionella.fasta";
my $mompS_MLST_schem = "schema/mompS.fas";

open (CONFIG, "config.txt") or die;
local $/;
my $config = <CONFIG>;
my ($blast_path) = $config =~ /blast_path\s*=\s*(.*?)\n/;
local $/ = "\n";
my $no_coverage_flag = qx(grep "Too low coverage" $outout_ST_File);
print "$no_coverage_flag\n";

if ($no_coverage_flag eq ""){

	`$blast_path/makeblastdb -in $fasta_file -dbtype nucl`;  # make blastDB

	#`blastn -query $query -db $fasta_file -evalue 0.0001 -outfmt 6 > $outoutFile`; # blast aginest mompS representator (1)

	`$blast_path/blastn -query $query -db $fasta_file -outfmt "6 qseqid sallseqid qlen slen qstart qend sstart send length evalue bitscore score pident qframe" -out $outoutBlastFile`;

	`perl Extract_MLST_mompS_Legionella.pl $fasta_file $outoutBlastFile $outout_log`;

	open (ST_RES, $outout_log) or die;
	my $title = <ST_RES>;
	my $seq = <ST_RES>;
	chomp $seq;
	my $mompS_type = <ST_RES>;
	chomp $mompS_type;
	my $mompS = "none";
	close (ST_RES);
	if ($mompS_type =~ /mompS ST number: \d+/){
			($mompS) = $mompS_type =~ /mompS ST number: (\d+)/;
	}
	open (BLAST_RES, ">>$outout_ST_File") or die;
	print "mompS - $mompS\n";
	if ($mompS eq "none"){
		my $MompSfasta = $outout_ST_File.".fasta";
		my $MompS_BLAST_RES = $outout_ST_File.".BLASTvsMompS.txt";
		open (mompSFasta, ">$MompSfasta") or die;
		print mompSFasta ">my_momS\n$seq\n";
		close (mompSFasta);
		`$blast_path/blastn -query $MompSfasta -db $mompS_MLST_schem -outfmt "6 qseqid sallseqid qlen slen qstart qend sstart send length evalue bitscore score pident qframe" -out $MompS_BLAST_RES`;
		
		#my_momS	14	352	352	1	352	1	352	352	0.0	 640	346	99.43	1
		open (BLASTmompS_RES, $MompS_BLAST_RES) or die;
		
		my $first_blast_res_line = <BLASTmompS_RES>;
		my $second_blast_res_line = <BLASTmompS_RES>;
		my $third_blast_res_line = <BLASTmompS_RES>;
		my ($qseqid,$sallseqid,$qlen,$slen,$qstart,$qend,$sstart,$send,$length,$evalue,$bitscore,$score,$pident,$qframe) = split("\t",$first_blast_res_line);
		print BLAST_RES "1st most close ST:\t$sallseqid\t$pident\n";
		my ($qseqid,$sallseqid,$qlen,$slen,$qstart,$qend,$sstart,$send,$length,$evalue,$bitscore,$score,$pident,$qframe) = split("\t",$second_blast_res_line);
		print BLAST_RES "2st most close ST:\t$sallseqid\t$pident\n";
		my ($qseqid,$sallseqid,$qlen,$slen,$qstart,$qend,$sstart,$send,$length,$evalue,$bitscore,$score,$pident,$qframe) = split("\t",$third_blast_res_line);
		print BLAST_RES "3st most close ST:\t$sallseqid\t$pident\n";	
		
	}else{
		print BLAST_RES "$mompS_type\n";
	}
}
