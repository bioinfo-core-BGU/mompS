#!/usr/bin/perl

###### Print USAGE if no argument is given
#my $usage = "\nUSAGE: Extract_SBT_seq_BaseOnBlastRes.pl <input_fasta_file> <input_blast_file>

 
###### Get the arguments:  
my $input_fasta_file = shift; # or die $usage;
my $input_blast_file = shift; # or die $usage;
my $output_log_file = shift;


####open SBT files 

open (mompS, "schema/mompS.fas") or die;



 sub reverse_complement {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}

#read mompS allel file
sub read_SBT_gene {
	my $hash;
	my $hendl = "mompS";
	while (my $first_line = <$hendl>){
		chomp $line;
		if ($first_line =~ /^>/){
			my ($num) = $first_line =~ /^>(.*)\n/;
			my $second_line = <$hendl>;
			chomp $second_line;
			$hash->{$hendl}{$second_line} = $num;
			#print "from read_SBT_gene\n$hendl\t$second_line\t$num\n";
		}
	}
	return $hash;
}


######read the SBT files
my $SBT_hash_map = read_SBT_gene();


#asd_1	4_S7_L001_R1_001_(paired)_contig_38	473	73398	1	473	9258	8786	473	0.0	 835	452	98.52	1
open(CONTIG,$input_fasta_file) or die;
open(BLAST,$input_blast_file) or die;
open(LOG, ">$output_log_file") or die;

#Read blast results
my $blastHash;
my $length_gene_hash;
while (my $line = <BLAST>){
	chomp $line;
	my ($qseqid,$sallseqid,$qlen,$slen,$qstart,$qend,$sstart,$send,$length,$evalue,$bitscore,$score,$pident,$qframe) = split("\t",$line);
	
	print "$sallseqid\n";
	if($qseqid =~ /\w+_\d+/){($qseqid) = $qseqid =~ /(\w+)_\d+/;} #for multi size SBT genes as neuA gene
	my $geneHeader = $qseqid."  Contig = ".$sallseqid." location on contig = ".$sstart."_".$send."_pident = ".$pident." eval = ".$evalue;
	$length_gene_hash->{$qseqid} = $length;
	if($sstart>$send){#complement
			$sstart = $sstart+ ($qstart-1);
			$send = $send - ($qlen-$qend);
			my $pos = $sstart."_".$send;
			$blastHash->{$sallseqid}{$pos} = $geneHeader;
		}else{
			$sstart = $sstart-($qstart-1);
			$send = $send+$qlen-$qend;
			my $pos = $sstart."_".$send;
			$blastHash->{$sallseqid}{$pos} = $geneHeader;
		}
}

local $/ = '>';

my $sub_seq;
my $SBT_hash;
my $junk = <CONTIG>; # Discard the ">" at the begining of the file

 while (my $line = <CONTIG>) {
 my $temp = "";
     chomp $line;
	 my $seq;
	
	my ($header, @seqLines) = split /\n/, $line;
	#	Join the individual sequence lines into one single sequence and store in a scalar variable.
	my $seq = join('',@seqLines); # Concatenates all elements of the @seqLines array into a single string.
	if ($header =~ /.*?\s+/){
		($header) = $header =~ /(.*?)\s+/;
	}
	
	 if (exists $blastHash->{$header}){
		my @positions = sort keys %{$blastHash->{$header}};
		foreach my $i (@positions){
			
			#print LOG ">$i\t";
			my ($start,$end) = split ("_",$i);
			my ($gene_name_line) = $blastHash->{$header}{$i};
			my ($gene_neme) = $gene_name_line =~ /(.*?)\s/; 
			if ($start>$end) {
				$temp = $start; 
				$start = $end-1; 
				$end = $temp;
				$sub_seq = substr $seq, $start, $end - $start;
				$sub_seq = reverse_complement($sub_seq);
			}else{
				
				$sub_seq = substr $seq, $start - 1, $end - $start+1;
			}
			print LOG ">$blastHash->{$header}{$i}\n";
			print LOG "$sub_seq\n";
			#print LOG "$seq\n";
			$SBT_hash->{$gene_neme} = $sub_seq;
			#print "Gene_neme\t$gene_neme\t$sub_seq\n";
		}
	 }
 }
 
 #print ">$input_blast_file\n";
 #print LOG ">Concat_seq\n";
 my @keys = sort keys %{$SBT_hash};
 my $profile = "";
 foreach my $key (@keys){
	my $print_length = $SBT_hash->{$key};
	#print LOG "$print_length";
	$profile = $SBT_hash_map->{$key}{$print_length};
 }
 #print "\n";
 print LOG "mompS ST number: $profile\n";

