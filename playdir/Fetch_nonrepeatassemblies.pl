#!/usr/bin/perl
#Script that takes the names of the repeats that don't hit the repeat database and extracts them from the original file

my $input = $ARGV[0];
my $originals = $ARGV[1];
open(input, "<$input");
open(orig, "<$originals");

my @allseq = <orig>;
my @nonrepeats = <input>;

foreach (@nonrepeats) {
	my $search_for = $_;
	my ( $index ) = grep {$allseq[$_] eq $search_for} 0..$#allseq; 
	print "$allseq[$index]$allseq[$index+1]";
#		my ($index) = grep { $IDs[$tag] eq $_ } 0..$#IDs;
}

