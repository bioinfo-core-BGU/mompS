#!/usr/local/bin/perl
use strict ;  
use warnings ;

my $sam_file = shift;
my $sam_out = shift;

my $border_start = 15;
my $border_end = 972;

my $border_reads_hash;

open (BAM, $sam_file) or die;
open (RES, ">$sam_out") or die;
while (my $line = <BAM>){
	chomp $line;
	my @lineArr = split("\t",$line);
	my $readID = $lineArr[0];
	my $readStart = $lineArr[3];
	my ($read_length) = $lineArr[5] =~ /(\d+)M/; # only maches
	
	my $end_border = $border_end-$read_length;
	
	if ($readID =~ "1114:12980:7947"){
		print "this is the read:$readID, read start:$readStart read_length:$read_length, border: $border_start-$end_border\n";
	}
	if ($readStart < $border_start or $readStart > $end_border){
			
		$border_reads_hash->{$readID}++;
	}
	
}
close (BAM);
open (BAM, $sam_file) or die;
while (my $line = <BAM>){
	chomp $line;
	if ($line =~ \\@\){
		print RES $line;
	}else{
		my @lineArr = split("\t",$line);
		my $readID = $lineArr[0];
		if (exists $border_reads_hash->{$readID}){
			print RES "$line\n";
		}
	}
}

