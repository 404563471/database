#!/usr/bin/perl
use Cwd;
my %num;
my $curdir=getcwd();
my $path=$1 if($curdir=~/(.*)\/.*/g);
&find_fileindir("$path");
sub find_fileindir(){
	  local($dir) = @_;
	  my $id=$1 if($dir=~/.*\/(.*)/);
	  opendir(DIR,"$dir"|| die "can't open this $dir");
      local @files =readdir(DIR);
	  closedir(DIR);
	  for $file (@files){
		  next if($file=~m/\.$/ || $file =~m/\.\.$/);
			 if ($file =~/reads.fa/i){
				  my $count;
				  #print "$dir\/$file \n";
				  open O,"$dir/$file" or die "cannot open the file!$!\n";
				  while(<O>){
					  $count++ if($_=~/>/g)
				  }
                  $num{$id}=$count;
				 #push @num,$count;
				  close O;
			}
			elsif(-d "$dir/$file"){
				find_fileindir("$dir/$file" );
			}
	}
}
#foreach my $key(sort keys %num){
#   print "$key\t$num{$key}\n";	
#}


my @files=glob("*.txt");
my @f;
my @id;
open O,">miRNA_express.out" or die "cannot open the out file!$!\n";
for my $i ( 0 .. $#files ) {
	if($files[$i]=~/.*_(.*).txt/g){
	   push @id,$1;
    }
    open my $fh, $files[$i];
    $f[$i] = $i ? { map  split, <$fh> } : [ grep chomp, <$fh> ];
}
my $s=join "\t",@id;
#print O join("\t", 'miRNAs', @files[1..$#files]), $/;
print O "$s\n";
for my $n ( @{ shift @f } ) {
	print O  join("\t", $n, map $$_{$n} // '0', @f), $/;
}
close O;

my $n;
my @id;
open A,"miRNA_express.out" or die "can't open miRNA_express.out!$!\n";
open B,">miRNA_express_normalize.out" or die "can't open miRNA_express_normalize.out!$!\n";
while(<A>){
 chomp;
 $n++;
 s/\s/\t/g;
 if($n == 1){
	 print B "$_\n";
	 @id=split "\t",$_;
	}else{
	my @data=split "\t",$_;
#	print "@data\n";
	print B "$data[0]";
	for($i=1;$i<=$#data;$i++){
	#	print "\t$data[$i]**";
		print B "\t";
		foreach my $key(sort keys %num){
		   if($id[$i] eq $key){
			  #print "$id[$i]***$key**$data[$i]***$num{$key}\t";
			  print B ($data[$i]/$num{$key}*10**6);
		   }	
		}
	}
	#print "\n";
	print B "\n";
 } 
}
close A;
close B;
