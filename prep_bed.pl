#!/usr/bin/perl
use strict;
use warnings;

if ($#ARGV != 2){
	die "use the program as perl $0 <Depth_sum> <Desired minimal length> <Desired minimal average depth>\n";
}

my $length_min=$ARGV[1];
my $depth_min=$ARGV[2];

open (INFILE,$ARGV[0]) || die $!;
readline INFILE;
while (<INFILE>){
	chomp;
	my @list = split(/\t/,$_);
	my ($length,$depth) = ($list[3],$list[4]);
#	print "$length\n";
#	print "$depth\n";
	if ($length >= $length_min && $depth >= $depth_min){
		#print join("\t",$list[0],$list[1],$list[2])."\n";
		print join("\t",@list)."\n";
	}
}
close INFILE;
