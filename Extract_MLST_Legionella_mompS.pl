#!/usr/bin/perl

###### Print USAGE if no argument is given
#my $usage = "\nUSAGE: Extract_SBT_seq_BaseOnBlastRes.pl <input_fasta_file> <input_blast_file>

 
###### Get the arguments:  
my $input_fasta_file = shift; # or die $usage;
my $input_blast_file = shift; # or die $usage;
my $mompS_ST = shift;
my $output_log_file = shift;


####open SBT files 
open (SBT, "schema/profiles.txt") or die;
open (asd, "schema/asd.fas") or die;
open (flaA, "schema/flaA.fas") or die;
open (mip, "schema/mip.fas") or die;
open (neuA, "schema/neuA.fas") or die;
open (pilE, "schema/pilE.fas") or die;
open (proA, "schema/proA.fas") or die;



 sub reverse_complement {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}

sub read_SBT_gene {
	my $hash;
	my @hendls = ("flaA","pilE","asd","mip","proA","neuA");
	foreach my $hendl (@hendls){
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
	}
	return $hash;
}

sub read_ST_profile_map {
	my $hash;
	my $hendl = shift;
	while(my $line = <$hendl>){
		chomp $line;
		$line =~ s/[\r\n]+//g;
		if($line =~ /^\d/){
			my($ST,$flaA,$pilE,$asd,$mip,$mompS,$proA,$neuA) = split(/\t/,$line);
			my $profile = "mompS".$mompS."asd".$asd."flaA".$flaA."mip".$mip."neuA".$neuA."pilE".$pilE."proA".$proA; #sorted
			$hash->{$profile} = $ST;
		}
	}
	return $hash;
}

######read the SBT files
my $SBT_hash_map = read_SBT_gene();


###### read ST profile
my $ST_hash = read_ST_profile_map("SBT");

#asd_1	4_S7_L001_R1_001_(paired)_contig_38	473	73398	1	473	9258	8786	473	0.0	 835	452	98.52	1
open(CONTIG,$input_fasta_file) or die;
open(BLAST,$input_blast_file) or die;
open(OUT, ">$output_log_file") or die;

my $blastHash;
my $length_gene_hash;
while (my $line = <BLAST>){
	chomp $line;
	my ($qseqid,$sallseqid,$qlen,$slen,$qstart,$qend,$sstart,$send,$length,$evalue,$bitscore,$score,$pident,$qframe) = split("\t",$line);
	
	print "$sallseqid\n";
	if($qseqid =~ /\w+_\d+/){($qseqid) = $qseqid =~ /(\w+)_\d+/;} #for multi size SBT genes as neuA gene
	my $geneHeader = $qseqid."  Contig = ".$sallseqid." location on contig = ".$sstart."_".$send."_pident = ".$pident." eval = ".$evalue;
	if (exists $length_gene_hash->{$qseqid}){
		my $exist_length = $length_gene_hash->{$qseqid};
		if ($length>$exist_length){
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
	}else{
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
		print "I'm here\n";
		my @positions = sort keys %{$blastHash->{$header}};
		foreach my $i (@positions){
			
			#print OUT ">$i\t";
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
			print OUT ">$blastHash->{$header}{$i}\n";
			print OUT "$sub_seq\n";
			print OUT "$seq\n";
			$SBT_hash->{$gene_neme} = $sub_seq;
		}
	 }
 }
 
  print OUT ">Concat_seq\n";
 my @keys = sort keys %{$SBT_hash};
 my $profile = "mompS".$mompS_ST;
 foreach my $key (@keys){
	my $print_length = $SBT_hash->{$key};
	print OUT "$print_length";
	$profile = $profile.$key.$SBT_hash_map->{$key}{$print_length};
 }
 #print "\n";
 print OUT "\n\n";
 print OUT "Profile = $profile\n";
 print OUT "ST = $ST_hash->{$profile}\n";
