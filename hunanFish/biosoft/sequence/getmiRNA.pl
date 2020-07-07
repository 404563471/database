#!usr/bin/perl
$fileA=$ARGV[0];
$fileB=$ARGV[1];
open IN,"<$fileA";
open SEQ,"<$fileB";
my %hash;
my $s=$/;
$/="\n>";
while(<SEQ>){
	s/>//g;
	if($_=~/(.*?)\n(.*)/g){
    	$hash{$1}=$2;
	   }
	}
$/=$s;


while(<IN>){
	my $line1=$_;
    my $line2=<IN>;
    $line1=~s/\n//g;
	$line2=~s/\n//g;
	if($line2=~/mir|let/g){
		foreach my $key(sort keys %hash){
			if($line1 eq $key){
				print ">$line1\n$hash{$key}\n";
				}
			}
		}
	}
