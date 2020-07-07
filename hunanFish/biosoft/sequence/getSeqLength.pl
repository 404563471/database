#!/usr/bin/perl -w
use warnings;
use strict;
###############################################################
#get length of files in fasta format 
my $reads_number;
my $base_number;
my $s=$/;
$/="\n>";
my @len;
while(<>){
	chomp;
	my $id;
	my $sequence;
	my @seq;
	if(/(.*?)\n(.*)/ms){
		 $id=$1;
		 $sequence=$2;
		 @seq=split "",$sequence;
     }
	my $len=$#seq+1;
	push @len,$len;
	$base_number+=$len;
	$reads_number++;
#	print  "$id\t$len\n";
}
$/=$s;
my %hash;
$hash{$_}++ for @len;
print "$_\t".$hash{$_}."\n" for (sort {$a<=>$b} keys %hash);
print  "all reads number is : $reads_number\n";
print  "all base number is : $base_number\n";
