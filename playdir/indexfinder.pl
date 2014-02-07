#!/usr/bin/perl

my @array = qw( your array here );

my @search_for = qw(stuff here);

foreach (@search_for) {
	my $temp = $_;
	my( $index )= grep { $array[$_] eq $temp } 0..$#array;
	print "$index\n";
}
