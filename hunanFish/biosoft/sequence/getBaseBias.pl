#!/usr/bin/perl -w
use warnings;
use strict;
###############################################################
#get the first base bias of all reads.
my $s=$/;
my %hash;
$/="\n>";
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
	push @{$hash{$len}},$seq[0];
}
$/=$s;
print "sRNA length\tA\tC\tG\tT\n";
#print "sRNA length\tA\tC\tG\tT\tsRNA length\tA\tC\tG\tT\n";
foreach my $key(sort keys %hash){
	print "$key\t";
	my %statistcs;
	my @percent;
	foreach my $item (@{$hash{$key}})
	{
		   $statistcs{$item}++;
	} 
	foreach my $item (sort keys(%statistcs))
	 {
		   print "$statistcs{$item}\t";
		   #push @percent,$statistcs{$item};
	 }
    #print "$key\t";
	#percent (\@percent);
    print "\n";
}
sub percent{
     my $list=shift @_;
	 my  $sum;
	 my @percent;
	 foreach (@$list){
		 $sum+=$_;
		 }
     foreach (@$list){
		 my $per;
		 $per=$_/$sum;
		 printf ("%.2f",$per*100);
		 print "%\t";
	 }
}
