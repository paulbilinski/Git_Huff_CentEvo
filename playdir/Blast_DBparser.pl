#!/usr/bin/perl
#Script that will extract all hits from blast

my $input = $ARGV[0];
open(input, "<$input");

while (<input>){
	chomp;
	my ($queryid, $subjectid, $identity, $alignment_length, $mismatches, $gapopens, $qstart, $qend, $sstart, $send, $evalue, $bitscore) = split("\t", $_);
	if ($alignment_length > 30){
		print "$subjectid\n";
	}
}
