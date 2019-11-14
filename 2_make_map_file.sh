#!/bin/bash
# Making map file


get_usage(){
	cat <<EOF

Usage : $0 [OPTION]

Description
	-h, --help
		show help

	-v, --version
		show version

	-d, --directory [data directory]
		directory name of analysis file locate

	-n, --name [sample name]
		sample name

	--f1 [fastq file(s)]
		Comma-separated list of fastq to be aligned
	
	--f2 [fastq file(2)]
		Comma-separated list of fastq to be aligned

	-x, --ref [ex. hg19]
		organism name

	-r, --restriction [HindIII|MboI|MboI-HinfI]
		name for restriction

	-o, --log [log directory]
		log file directory

	-m, --mapq [mapq threshold (default:30)]
		threshold mapQ to make map

	-q, --fastqc
		TRUE for doing fastqc analysis or FALSE for not doing (default FALSE)

EOF

}

get_version(){
	echo "sh ${0} version 1.0"
}

SHORT=hvd:n:x:r:o:m:q:
LONG=help,version,directory:,name:,f1:,f2:,ref:,restriction:,log:,mapq:,fastqc:
PARSED=`getopt --options $SHORT --longoptions $LONG --name "$0" -- "$@"`
if [[ $? -ne 0 ]]; then
	exit 2
fi
eval set -- "$PARSED"

while true; do
	case "$1" in
		-h|--help)
			get_usage
			exit 1
			;;
		-v|--version)
			get_version
			exit 1
			;;
		-d|--directory)
			DIR_DATA="$2"
			shift 2
			;;
		-n|--name)
			NAME="$2"
			shift 2
			;;
		--f1)
			FILE_fastq1="$2"
			shift 2
			;;
		--f2)
			FILE_fastq2="$2"
			shift 2
			;;
		-x|--ref)
			REF="$2"
			shift 2
			;;
		-r|--restriction)
			RESTRICTION="$2"
			shift 2
			;;
		-o|--log)
			DIR_LOG="$2"
			shift 2
			;;
		-m|--mapq)
			MAPQ_THRESHOLD="$2"
			shift 2
			;;
		-q|--fastqc)
			FLAG_fastqc="$2"
			shift 2
			;;
		--)
			shift
			break
			;;
		*)
			echo "Programming error"
			exit 3
			;;
	esac
done

DIR_LIB=$(dirname $0)
TIME_STAMP=$(date +"%Y-%m-%d")

[ ! -n "${NAME}" ] && echo "Please specify NAME" && exit 1
[ ! -n "${REF}" ] && echo "Please specify ref" && exit 1
[ ! -n "${RESTRICTION}" ] && echo "Please specify restriction" && exit 1
# [ ! -n "${DIR_LOG}" ] && echo "Please specify log directory" && exit 1
[ ! -n "${DIR_DATA}" ] && echo "Please specify data directory" && exit 1
[ ! -n "${FILE_fastq1}" ] && echo "Please specify fastq file1" && exit 1
[ ! -n "${FILE_fastq2}" ] && echo "Please specify fastq file2" && exit 1
MAPQ_THRESHOLD=${MAPQ_THRESHOLD:-30}
FLAG_fastqc=${FLAG_fastqc:-FALSE}


cd ${DIR_DATA}

#-----------------------------------------------
# Load setting
#-----------------------------------------------
source ${DIR_LIB}/utils/load_setting.sh -x $REF -r $RESTRICTION


#-----------------------------------------------
# Alignment
#-----------------------------------------------
export BOWTIE2_INDEX DIR_DATA
export FILE_fastq=${FILE_fastq1} OUT=${NAME}_1
sh ${DIR_LIB}/utils/Alignment_with_trimming.sh
export FILE_fastq=${FILE_fastq2} OUT=${NAME}_2
sh ${DIR_LIB}/utils/Alignment_with_trimming.sh


#-----------------------------------------------
# fastqc
#-----------------------------------------------
if hash module 2>/dev/null; then
	module load fastqc
fi
if [ "$FLAG_fastqc" = "TRUE" ]; then
	[ ! -e "${DIR_DATA}/fastqc" ] && mkdir "${DIR_DATA}/fastqc"
	fastqc -o fastqc/ --nogroup -t 12 $(echo ${FILE_fastq1} | tr ',' ' ')
	fastqc -o fastqc/ --nogroup -t 12 $(echo ${FILE_fastq2} | tr ',' ' ')
fi


#-----------------------------------------------
# assign to nearest restriction enzyme site
#-----------------------------------------------
# "-" direction => shift align position + read length 
# Repeat or Unique was categorized based on the tag "XS:i"
# Location side column is L for "+" direction, R for "-" direction
perl ${DIR_LIB}/utils/Assign_nearest_enzymeSites.pl -a ${NAME}_1.sam -b ${NAME}_2.sam -o ${NAME}.map -e ${FILE_enzyme_def} -d ${FILE_enzyme_index} > ${NAME}_alignment.log


#-----------------------------------------------
# Extract unique map file and register to database
#-----------------------------------------------
# swap left and right to make left side always smaller than right
# if restriction sites were missing, remove those pairs
# entire chromosome were roughly split to 100 and split map file based 
perl ${DIR_LIB}/utils/Split_MapFile.pl -i ${NAME}.map -l ${CHROM_LENGTH} -o ${NAME}_list.txt
cat ${NAME}_list.txt | xargs -P12 -I@ sh -c "sort -k2,2 -k3,3n -k9,9 -k10,10n @ | awk -v OFS='\t' '{print \$0,\$2,\$3,\$9,\$10}' | uniq -f 15 | cut -f1-15 > @_sorted && mv @_sorted @"
echo "id, chr1, position1, direction1, mapQ1, restNum1, restLoc1, uniq1, chr2, position2, direction2, mapQ2, restNum2, restLoc2, uniq2" | tr ',' '\t' > ${NAME}.map
cat $(cat ${NAME}_list.txt) >> ${NAME}.map
rm $(cat ${NAME}_list.txt) ${NAME}_list.txt
Rscript --vanilla --slave ${DIR_LIB}/utils/file2database_large.R -i ${NAME}.map --db ${NAME}.db --table map


#-----------------------------------------------
# Summarize read filtering
#-----------------------------------------------
export DIR_DATA SAMPLE=${NAME}
sh ${DIR_LIB}/utils/Count_reads.sh > ${NAME}_read_filtering.log



