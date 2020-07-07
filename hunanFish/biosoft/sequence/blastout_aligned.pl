#!/usr/bin/perl

use warnings;
use strict;

use Getopt::Std;

my $usage ="$0 blastoutputfile mis_max_num > outfile";

open IN, "< $ARGV[0]";
my $mis_max= $ARGV[1];

my $query;
my $query_lng;
my $subject;
my $subject_lng;
my $bits;
my $expect;
my $identities;
my $gaps;
my $frame;
my $q_st;
my $q_ed;
my $s_st;
my $s_ed;




while (<IN>) {
    
  X:if (/^Query\s*=\s*(\S+)/) {
      #resetting variables to catch parsing errors
     $query_lng="";
     
     $query=$1;
 }
    
    if(/\s+\(([\d\,]+)\s+letters\)/) {
	$query_lng=$1;
    }
    

    #Requiring that query_lng is defined is to avoid parsing parts of the query descriptor
    #containing '>' as being the subject and going into the do while loop. Is currently a test.
    if (/^>(\S+)/ && $query_lng) {
	#resetting variables to catch parsing errors
	$subject_lng="";
	
	$subject=$1;
	
	
     	do {
     	    $_ = <IN>;
     	    tr/,//d;
	    
	        if (/\s+Length\s*=\s*(\d+)/) {$subject_lng=$1;}
	    
  	    } while (!/\s+Length\s*=\s*(\d+)/);
    }
    
    
    if(/^\s{0,2}Score/){
	if (/^\s+Score\s*=\s*(\S+)\s*.*Expect\S*\s*=\s*(\S+)$/) {
	    $bits=$1;
	    $expect=$2;
	}else{
	    print STDERR "problem parsing score and expect line\n";
	    exit;
	}

	$_ = <IN>;
	if (!/Identities\s*=\s*(\S+)/) {
	    print STDERR "problem parsing identities line\n";
	    exit;
	}
	else { $identities=$1; } #?
	
	if (/Gaps\s*=\s*(\d+)\//) { $gaps = $1; } #?
	else                      { $gaps = 0;  }
	
	
	$_ = <IN>;
	if (/\s+Frame\s*=\s*(\S+)\s*\n/){ #?
	    $frame=$1;
	}elsif(/\s*Strand\s*=\s*(\w+\s*\/\s*\w+)/){
	    $frame=$1;
	}else{
	  print STDERR "error reading the strand or frame\n";
      } 
	
	# read positions
	$q_st = "";
	$q_ed = "";
	$s_st = "";
	$s_ed = "";
	do {
	    $_ = <IN>;
	    if (/^Query:\s+(\d+).*\s(\d+)\s*$/)  {
		if ($q_st eq "")                 { $q_st = $1; }
		if ($1 < $q_st)                  { $q_st = $1; }
		if ($2 < $q_st)                  { $q_st = $2; }
		if ($q_ed eq "" or $1 > $q_ed)   { $q_ed = $1; }
		if ($2 > $q_ed)                  { $q_ed = $2; }
	    }
	    if (/^Sbjct:\s+(\d+).*\s(\d+)\s*$/)  {
		if ($s_st eq "")	         { $s_st = $1; }
		if ($1 < $s_st)		         { $s_st = $1; }
		if ($2 < $s_st)		         { $s_st = $2; }
		if ($s_ed eq "" or $1 > $s_ed)   { $s_ed = $1; }
		if ($2 > $s_ed)		         { $s_ed = $2; }
	    }
	} while ((/^\s*$/ or /^\s{8,}(\w|\s|\+|\||\*)+$/ or /^Query:/ or /^Sbjct:/) and !/^>/);
	
	my $strand=findstrand($frame);
	if($strand eq "-"){($s_st,$s_ed)=reversepositions($subject_lng,$s_st,$s_ed);}

	
	my @v = split(/\//, $identities);
	my $pid = $v[0]/$v[1];

	my $mis = $v[1] - $v[0];
	my $mis_total = mismatch($query_lng, $q_st, $q_ed, $mis);

	if ($mis_total <= $mis_max and $frame =~ /Plus \/ Plus/) {
		print "$query\n";
		print "\>$subject\n";
		do{
			$_ = <IN>;
		} until (/^Query\s*=\s*(\S+)/ or eof);
		
	}
	
	unless($query and $query_lng and $q_st and $q_ed and $subject and $subject_lng and $s_st and $s_ed and $expect and $pid and $bits){
	    print STDERR "problem parsing entry\n";
	    exit;
	}
	
	#resetting variables to catch parsing errors
	$bits="";
	$expect="";
	$identities="";
	$gaps="";
	$frame="";
	$q_st="";
	$q_ed="";
	$s_st="";
	$s_ed="";

	goto X;
	
    }
}




sub findstrand{

    #A subroutine to find the strand, parsing different blast formats
    my($other)=@_;

    my $strand="+";

    if($other=~/-/){
	$strand="-";
    }

    if($other=~/minus/i){
	$strand="-";
    }

    return($strand);
}



sub mismatch{
	my ($query_lng, $q_st, $q_ed, $mis) = @_;
	my $mis_total = $mis + ($q_st - 1) + ($query_lng - $q_ed);
	return ($mis_total);
}




sub reversepositions{

    #A subroutine to find positions relative to the minus strand
    my($length,$begin,$end)=@_;

    my $new_end=$length-$begin+1;
    my $new_beg=$length-$end+1;

    return($new_beg,$new_end);
}
