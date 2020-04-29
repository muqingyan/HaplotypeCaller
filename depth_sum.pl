#!/usr/bin/perl
use warnings;
use strict;
use POSIX;

if ($#ARGV != 1){
	die "Use it as 'perl $0 <Depth_file> <Result_file>'\n";
	
}
open(INFILE, $ARGV[0]) or die "Could not open the file: $!\n";
open(OUTFILE, ">$ARGV[1]");

my @new_list;
my $chr=1;
my $start=0;
my $length=1;
my $end=0;
my $depth_sum=0;
my $average_depth=0;
my $string="";
print OUTFILE ("chr\tstart\tend\tlength\taverage_depth\n");
while(<INFILE>){
	chomp;
	my @list = split(/\t/,$_);
	my ($new_chr,$new_pos,$new_depth)=@list;
	if ($new_chr eq $chr)
	{
		if($new_pos - $end eq 1)
		{
			$length+=1;
			$end=$new_pos;
			$depth_sum+=$new_depth;
			$average_depth=$depth_sum/($end-$start+1);
			$average_depth=ceil($average_depth);
		}
		else
		{
			print OUTFILE ("$chr\t$start\t$end\t$length\t$average_depth\n");
			$start=$new_pos;
			$end=$new_pos;
			$depth_sum=$new_depth;
			$length=1;
		}
	}
	else
	{
		$string=join("\t",$chr,$start,$end,$length,$average_depth);
		print OUTFILE ("$string\n");
		$chr=$new_chr;
		$start=$new_pos;
		$end=$new_pos;
		$depth_sum=$new_depth;
		$length=1;
	}
}
print OUTFILE ("$chr\t$start\t$end\t$length\t$average_depth\n");
close (INFILE);
close (OUTFILE);
