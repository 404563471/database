
#predict target seq using miRNAs
miranda miRNA.fa $UTR -out mirna_target.out

perl $program/getmiRNAtarget.pl mirna_target.out
