#!/usr/bin/perl -w

open(FILE,"test.fq");
open(OUT,">a_out.txt");
while(<FILE>){
	chomp;
	print OUT $_.'\n';
}	
close(FILE);close(OUT);
