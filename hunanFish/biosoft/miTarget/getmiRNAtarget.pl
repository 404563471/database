#!usr/bin/perl
#(1) sRNA与靶基因间的错配不得超过4个(G-U配对认为0.5个错配) 
#(2) 在miRNA/靶基因复合体中不得有超过2处发生相邻位点的错配 
#(3) 在miRNA/靶基因复合体中，从miRNA的5’端起第2-12个位点不得有相邻位点都发生错配 
#(4) miRNA/靶基因复合体的第10-11个位点不得发生错配 
#(5) 在miRNA/靶基因复合体中，从miRNA的5’端起第1-12个位点不得有超过2.5个错配 
open O,">$ARGV[1]/$ARGV[2].miRNA_target_alinged.out" or die "cannot open the align file.$!\n";
open T,">$ARGV[1]/$ARGV[2].miRNA_target_count.out" or die "cannot open the count file.$!\n";
open N,">$ARGV[1]/$ARGV[2].miRNA_target_number.out" or die "cannot open the number file.$!\n";
open W,">$ARGV[1]/$ARGV[2].miRNA_network.txt" or die "cannot open the network file.$!\n";
my %hash;
my $s=$/;
$/="Performing Scan";
print O "mirna_name\ttarget_id\tmirna_alignment\talingment\ttarget_alignment\tscore\tenergy\tmirna_start\tmirna_end\ttarget_start\ttarget_end\talign_length\tmirna_conservation\ttarget_conservation\n";
open(FILE,$ARGV[0]);
while(<FILE>){
  if($_ =~/No\sHits\sFound/mg){
	  next;
  }else{
	   if($_ =~ />/mg){	
	    my $line=$_;
		my @query;$query;@ref;$ref;$query_name;$seq_name;$message;$match;
		if($line=~/Query.*?\d.*?\s(.*?)\s.*?\n(.*?)\n/mg){
			$query=$1;
			@query=split "",$query;
			$match=$2;
			$match=~s/\s{16}//;
		}
		if($line=~/Ref.*?\d.*?\s(.*?)\s.*/mg){
			$ref=$1;
			@ref=split "",$ref;
		}
		if($line=~/>(.*?)\t(.*?)\t(.*?)\n/mg){
			$query_name=$1;
			$seq_name=$2;
			$message=$3;
		}
		my @data;
#统计匹配信息，完全匹配标为1，G-U为0.5，错配或缺失为0
		for(my $i==0;$i<=$#query;$i++){
			#$data[$i]=join "_",$query[$i],$ref[$i];
			if($query[$i]=~/a/i && $ref[$i]=~/t/i){
				push @data,1;
			}elsif($query[$i]=~/u/i && $ref[$i]=~/a/i){
				push @data,1;
			}elsif($query[$i]=~/g/i && $ref[$i]=~/c/i){
			    push @data,1;
			}elsif($query[$i]=~/c/i && $ref[$i]=~/g/i){
			    push @data,1;
			}elsif($query[$i]=~/u/i && $ref[$i]=~/g/i){
			    push @data,0.5;	
			}else{
			    push @data,0;	
			}
		} 
		my $wrong_pairs=0;  #整个复合体中的错配数目,要求不能超过4
		my $a_wr_num=0;     #相邻错配的pairs对数目，要求不超过2
		my $a_wr_12=0;      #miRNA 5'端2-12位点相邻错配的数目，要求为0
		my $wrong_pairs_12=0; #miRNA 5'端1-12位点的错配数目，要求不超过2.5
		for(my $i=0;$i<=$#data;$i++){
			if($data[$i] == 0){
				$wrong_pairs++;
			}elsif($data[$i] == 0.5){
				$wrong_pairs+=0.5;
			}
		}
		for(my $i=0;$i<$#data;$i++){
		    if($data[$i] == 0 && $data[$i+1] == 0){
				$a_wr_num++;
			}
		}
		for(my $i=($#data-1);$i>=($#data-10);$i--){
			$a_wr_12++ if($data[$i]==0 && $data[$i-1]==0);
		}
        for(my $i=$#data;$i>=($#data-11);$i--){
			if($data[$i] == 0){
				$wrong_pairs_12++;	
			}elsif($data[$i] == 0.5){
			    $wrong_pairs_12+=0.5;	
			}
		}
		if($wrong_pairs <= 4 && $a_wr_num <= 2 && $a_wr_12 == 0 && $wrong_pairs_12 <= 2.5 && $data[$#data-9] == 1 && $data[$#data-10] == 1){
			print O "$query_name\t$seq_name\t$query\t$match\t$ref\t$message\n";
			#print "$query_name\t$seq_name\n@query\n@ref\n@data\n";
			#print "wrong_pairs:$wrong_pairs\ta_wr_num:$a_wr_num\twrong_pairs_12:$wrong_pairs_12\n";
			push @{$hash{$query_name}},$seq_name;
  		} 
    }
}
}
close(FILE);
my $key_num;
my %target;
my $target_num;
foreach my $key(sort keys %hash){
	my %line;#每个miRNA的target列表
	my @line;
	$key_num++;
	foreach(@{$hash{$key}}){
		$target{$_}=1;
		$line{$_}=1;
	}
	@line=keys %line;
	print W "$key\t$_\n" foreach(@line);
	$target_num+=$#line+1;
	print N "$key\t".($#line+1)."\n";
	print T "$key\t".($#line+1)."\n@line\n";
	print T "--------------------\n";
}
my @target=keys %target;
print T "miRNA Number:$key_num\nuniq target Number:".($#target+1)."\n";
#print T "miRNA Number:$key_num\nsum target Number:".$target_num."\n";
close O;
close T;
close N;
close W;
