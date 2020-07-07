#! /bin/bash

#$1 == input file path
#$2 == output file path

com1=./bin/clustalw2/clustalw2
com2=./bin/FastTree
com3=./bin/tree.R

#func_judge(){
#        source /etc/init.d/functions
#	if [ $? -ne 0 ]; then
#	    action "${1##*/} failed" /bin/false
#	else
#	    action "${1##*/} success" /bin/true
#	fi
#}

if [ ! -d $2 ]; then
    mkdir $2
fi

$com1 -INFILE=$1/multest.fasta  -OUTPUT=FASTA  -OUTFILE=$2/multest.out -SEQNOS=ON

#func_judge $com1

$com2 -nt -gtr $2/multest.out > $2/test.tree

#func_judge $com2

Rscript $com3 $2/test.tree $2/tree.pdf 2> $2/file.error

#func_judge $com3
