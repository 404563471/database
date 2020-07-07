#!usr/bin/perl
print "count sRNA number from Rfam blast out:\n";
my $rRNA;
my $scRNA;
my $snRNA;
my $snoRNA;
my $tRNA;
my $other;
my $total;
#open O,"reads_md_collapse.fa" or die "can't open the reads file.$!\n";
$fileA=$ARGV[0];
$fileB=$ARGV[1];
open O,"<$fileA";
open B,"<$fileB";
my $s=$/;
$/="\n>";
while(<O>){
   $total++;
}
close O;
$/=$s;
print "rRNA\tscRNA\tsnRNA\tsnoRNA\ttRNA\tother\ttotal\n";
while(<B>){
	chomp;
	if($_=~/>.*rRNA/six){
		$rRNA++;
	}elsif($_=~/>.*?scRNA/six){
	    $scRNA++;		
	}elsif($_=~/>.*?(snR|sn\d)/six){
		$snRNA++;
	}elsif($_=~/>.*?(sno|Afu|ceN|DdR)/six){
		$snoRNA++;
	}elsif($_=~/>.*?tRNA/six){
		$tRNA++;
	}
}
$other=$total-$rRNA-$scRNA-$snRNA-$snoRNA-$tRNA;
print "$rRNA\t$scRNA\t$snRNA\t$snoRNA\t$tRNA\t$other\t$total\n";
print "--------------------------------\n";
