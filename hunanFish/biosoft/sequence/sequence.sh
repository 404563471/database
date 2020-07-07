sample="test.fq"   #in fastaq format
spece="hsa"            #the abbsract of the sepece.
adaptor="AGATCGGAAGAGCACACGTCT"             #3' adaptor
id=`echo $sample | sed s/.fq//`
pic_path="picture"     #the path of picture

tools="/share/home/lije/bin"

#if [ ! -d "out"];then
     mkdir ../out
#fi
cd ../out
mkdir $pic_path
#quality of data
ln -s ../sequence/$sample ./
fastqc $sample
#remove adorpter
cutadapt -a $adaptor --discard-untrimmed -m 10 $sample -o reads_mv_Adaptor.fastq
#IlluQC.pl -se reads_mv_Adaptor.fastq 5 5 -t 2 -z t
#AmbiguityFiltering.pl -i IlluQC_Filtered_files/*filtered -c 0 -o reads.fastq
#fastq to fasta
perl ../FastqToFasta.pl -i reads_mv_Adaptor.fastq -o reads.fa
#get the length of clean data
perl ../getSeqLength.pl reads.fa > length_$id.txt
echo "aaaa".$id."aaa".$pic_path."aa"
Rscript ../reads_len.R $id $pic_path


#get unique reads 
perl ../collapse_reads_md.pl reads.fa $spece -a > reads_md.fa

#blast with Rfam database
blastall -p blastn -d ../Rfam.fa -i reads_md.fa -o reads_md_Rfam.out -a 12 -m 8 -e 0.01 -F F
perl ../getBlastDB.pl -c reads_md.fa -b  reads_md_Rfam.out -id $id -iden 100 -cov 1.0 -o Rfam_$id
#perl ../blastout_aligned.pl reads_md_Rfam.out 3 > blast_Rfam.txt
#perl ../getRfamNumber.pl reads_md.fa blast_Rfam.txt >Rfam.txt
#cat Rfam.txt >>count_number.txt
#sed -n "2,3p" Rfam.txt > Rfam_R.txt
#Rscript ../Rfam.R $id $pic_path
#get miRNA sequence 
#perl ../getmiRNA.pl blast_Rfam.txt reads_md.fa > miRNA_predict.fa

#get the first base bias of the clean data
perl ../getBaseBias.pl miRNA_predict_$id.fa > FirstBaseBias_$id.txt
#get position base bias of all clean data
perl ../getBaseBiasAll.pl miRNA_predict_$id.fa > AllBaseBias_$id.txt
Rscript ../getBasebias.R $id $pic_path

