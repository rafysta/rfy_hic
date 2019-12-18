#!/bin/bash

if hash module 2>/dev/null; then
	module load bowtie2
	module load samtools
fi

DIR_LIB=$(dirname $0)
PROGRAM_CIGAR=${DIR_LIB}/CigarFilter.pl
VERBOSE=${VERBOSE:-FALSE}
DIR_tmp=$(mktemp -d /tmp/tmp_${OUT}.XXXXX)
FILE_log=${DIR_DATA}/${OUT}_bowtie2.log
FILE_tmp=${DIR_tmp}/${OUT}_bowtie2_tmp.log


cd ${DIR_DATA}

echo "Length Total % NoAlign % Unique % Multiple %" | tr ' ' '\t' > $FILE_log
function getLog(){
	[ $VERBOSE = "TRUE" ] && cat $FILE_tmp
	cat $FILE_tmp | tr -d '()' | awk -v OFS='\t' 'NR>1&&NR<6{print $1,$2}' | xargs | tr ' ' '\t' >> $FILE_log && rm $FILE_tmp
}

FILE_first=$(echo ${FILE_fastq} | cut -f1 -d',')
let READ_LENGTH=$(head -n 2 ${FILE_first} | tail -n 1 | wc -m)
let MAX2_TRIM=$READ_LENGTH-25
let MAX_TRIM=$READ_LENGTH-20

[ $VERBOSE = "TRUE" ] && echo "Align entire ${READ_LENGTH}bp read"
echo -n "${READ_LENGTH}bp	" >> $FILE_log
bowtie2 -x ${BOWTIE2_INDEX} -U ${FILE_fastq} -q -p 12 --no-unal --un ${DIR_tmp}/${OUT}_unaligned.fastq > ${DIR_tmp}/${OUT}.sam 2> $FILE_tmp && getLog
mv ${DIR_tmp}/${OUT}_unaligned.fastq ${DIR_tmp}/${OUT}_tmp.fastq

### Change to temp directory
cd ${DIR_tmp}

for LEN in $(seq 5 5 $MAX2_TRIM)
do
	[ $VERBOSE = "TRUE" ] && echo && echo "Trim ${LEN}bp and align"
	echo -n "Trim ${LEN}bp	" >> $FILE_log
	bowtie2 -x ${BOWTIE2_INDEX} -U ${OUT}_tmp.fastq -q -3 $LEN -p 12 --no-hd --no-unal --un ${OUT}_unaligned.fastq >> ${OUT}.sam 2> $FILE_tmp && getLog
	mv ${OUT}_unaligned.fastq ${OUT}_tmp.fastq
done

[ $VERBOSE = "TRUE" ] && echo && echo "Trim ${MAX_TRIM}bp and align"
echo -n "Trim ${MAX_TRIM}bp	" >> $FILE_log
bowtie2 -x ${BOWTIE2_INDEX} -U ${OUT}_tmp.fastq -q -3 ${MAX_TRIM} -p 12 --no-hd  >> ${OUT}.sam 2> $FILE_tmp && getLog
rm ${OUT}_tmp.fastq

perl ${PROGRAM_CIGAR} -s ${OUT}.sam -o ${OUT}_fastqList.txt > ${OUT}_tmp.sam
mv ${OUT}_tmp.sam ${OUT}.sam
LIST_NUM=$(cat ${OUT}_fastqList.txt | wc -l)
for i in $(seq 1 ${LIST_NUM})
do
	TRIM_OPTION=$(head -n $i ${OUT}_fastqList.txt | tail -n 1 | cut -f1)
	FASTQ_REANALYZE=$(head -n $i ${OUT}_fastqList.txt | tail -n 1 | cut -f2)
	# echo "Command : bowtie2 -x ${BOWTIE2_INDEX} -U ${FASTQ_REANALYZE} -q ${TRIM_OPTION} -p 12 --no-hd"
	[ $VERBOSE = "TRUE" ] && echo && echo "Bowtie2 option: ${TRIM_OPTION} and align"
	echo -n "Option:${TRIM_OPTION}	" >> $FILE_log
	bowtie2 -x ${BOWTIE2_INDEX} -U ${FASTQ_REANALYZE} -q ${TRIM_OPTION} -p 12 --no-hd >> ${OUT}.sam 2> $FILE_tmp && getLog
	rm ${FASTQ_REANALYZE}
done
rm ${OUT}_fastqList.txt

samtools view -bS ${OUT}.sam > ${OUT}.bam

### read check
TOTAL_READ=$(cat $FILE_log | awk 'NR==2{print $2}')
TOTAL_BAM=$(samtools view -c ${OUT}.bam)
if [ $TOTAL_READ -ne $TOTAL_BAM ]; then
	echo "Total reads: $TOTAL_READ and Total Bam output: $TOTAL_BAM not matched"
	exit 1
fi

samtools sort -n ${OUT}.bam -o ${OUT}_sort.bam -T tmpBamSort_${OUT}
samtools view ${OUT}_sort.bam > ${OUT}.sam 
rm ${OUT}_sort.bam ${OUT}.bam

### Transfer to data directory
mv ${OUT}.sam ${DIR_DATA}/${OUT}.sam

rm -r ${DIR_tmp}
exit 0
