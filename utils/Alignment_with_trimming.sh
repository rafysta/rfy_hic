#!/bin/bash

if hash module 2>/dev/null; then
	module load bowtie2
	module load samtools
fi

DIR_LIB=$(dirname $0)
PROGRAM_CIGAR=${DIR_LIB}/CigarFilter.pl
cd ${DIR_DATA}

echo "#=========================="
echo "# alignment of ${FILE_fastq}"
echo "#=========================="
date
echo
FILE_first=$(echo ${FILE_fastq} | cut -f1 -d',')
let READ_LENGTH=$(head -n 2 ${FILE_first} | tail -n 1 | wc -m)
let MAX2_TRIM=$READ_LENGTH-25
let MAX_TRIM=$READ_LENGTH-20

echo "alignment of ${READ_LENGTH}bp read" 
bowtie2 -x ${BOWTIE2_INDEX} -U ${FILE_fastq} -q -p 12 --no-unal --un ${OUT}_unaligned.fastq > ${OUT}.sam
mv ${OUT}_unaligned.fastq ${OUT}_tmp.fastq
echo

for LEN in $(seq 5 5 $MAX2_TRIM)
do
	echo "trimming $LEN bp" 
	bowtie2 -x ${BOWTIE2_INDEX} -U ${OUT}_tmp.fastq -q -3 $LEN -p 12 --no-hd --no-unal --un ${OUT}_unaligned.fastq >> ${OUT}.sam
	mv ${OUT}_unaligned.fastq ${OUT}_tmp.fastq
	echo
done

echo "trimming $MAX_TRIM bp" 
bowtie2 -x ${BOWTIE2_INDEX} -U ${OUT}_tmp.fastq -q -3 ${MAX_TRIM} -p 12 --no-hd  >> ${OUT}.sam
rm ${OUT}_tmp.fastq
echo

echo "alignment of partial alignment reads" 
perl ${PROGRAM_CIGAR} -s ${OUT}.sam -o ${OUT}_fastqList.txt > ${OUT}_tmp.sam
mv ${OUT}_tmp.sam ${OUT}.sam
LIST_NUM=$(cat ${OUT}_fastqList.txt | wc -l)
for i in $(seq 1 ${LIST_NUM})
do
	TRIM_OPTION=$(head -n $i ${OUT}_fastqList.txt | tail -n 1 | cut -f1)
	FASTQ_REANALYZE=$(head -n $i ${OUT}_fastqList.txt | tail -n 1 | cut -f2)
	echo "Command : bowtie2 -x ${BOWTIE2_INDEX} -U ${FASTQ_REANALYZE} -q ${TRIM_OPTION} -p 12 --no-hd" 
	bowtie2 -x ${BOWTIE2_INDEX} -U ${FASTQ_REANALYZE} -q ${TRIM_OPTION} -p 12 --no-hd >> ${OUT}.sam
	rm ${FASTQ_REANALYZE}
	echo
done
rm ${OUT}_fastqList.txt

echo "convert to bam file" 
samtools view -bS ${OUT}.sam > ${OUT}.bam
echo "sort bam file" 
samtools sort -n ${OUT}.bam -o ${OUT}_sort.bam -T tmpBamSort_${OUT}
echo "convert again to sam file" 
samtools view ${OUT}_sort.bam > ${OUT}.sam 
rm ${OUT}_sort.bam ${OUT}.bam

echo "finished" 
date 
echo
echo

