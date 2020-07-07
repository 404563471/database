#!/usr/bin/perl -w
use warnings;
use strict;
###############################################################
#get bias of every position bases.
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
		$sequence =~ s/\n//;
        @seq=split "",$sequence;
     }
	my $i;
#	print "@seq\t$#seq\n";
	for($i=0;$i<=$#seq;$i++){
        push @{$hash{$i}},$seq[$i];
	}
}
$/=$s;
print "sRNA position\tA\tC\tG\tT\n";
#print "sRNA position\tA\tC\tG\tT\tsRNA position\tA\tC\tG\tT\n";
foreach my $key(sort {$a<=>$b} keys %hash){
	print $key+1;
	my %statistcs=("A"=>0,"T"=>0,"C"=>0,"G"=>0);
	my @percent;
	foreach my $item (@{$hash{$key}}){
         $statistcs{$item}++;
	} 
	my @item=keys %statistcs;
	my @value=values %statistcs;
#	print "@value\n";
	foreach my $item (sort keys(%statistcs)){
		 print "\t$statistcs{$item}"; 
    }
	#print $key+1,"\t";
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
