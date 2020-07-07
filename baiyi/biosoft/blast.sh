#!/bin/bash

#example:./blast/bin/blastn -query {$sequence} -db $rootpath/blast/database/mitochondrial-genome -out {$outfile} -outfmt=0 -evalue {$faval} -num_descriptions {$num}
#$1 == program
#$2 == fasta
#$3 == database
#$4 == resultpath
#$5 == evalue 
#$6 == number_result

# useage example : ./blast.sh blastn example/test.fasta mitochondrial-genome out/out1 1 2

blast_path=/biosoft/blast/blast/bin/
database_path=/biosoft/blast/database/

$blast_path$1 -query $2 -db $database_path$3 -out $4 -outfmt=0 -evalue $5 -num_descriptions $6 
