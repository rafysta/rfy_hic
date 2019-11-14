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
	echo "${0} version 1.0"
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



source ${DIR_LIB}/utils/load_setting.sh -x $REF -r $RESTRICTION


# [ "$FLAG_fastqc" = "TRUE" ] && [ ! -e "${DIR_DATA}/fastqc" ] && mkdir "${DIR_DATA}/fastqc"


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
# [ "$FLAG_fastqc" = "TRUE" ] &&  [ ! -e ${DIR_DATA}/fastqc/${NAME}_fastqc ] && sbatch -n 12 --job-name=fastqc_${NAME}_1 $(sq --node) -o "${DIR_LOG}/${TIME_STAMP}_fastqc_${NAME}_1.log" --open-mode append --wrap="cd ${DIR_DATA}; /applications/fastqc/current/fastqc -o fastqc/ --nogroup -t 12 ${NAME}_1.fastq" && sbatch -n 12 --job-name=fastqc_${NAME}_2 -o "${DIR_LOG}/${TIME_STAMP}_fastqc_${NAME}_2.log" --open-mode append --wrap="cd ${DIR_DATA}; /applications/fastqc/current/fastqc -o fastqc/ --nogroup -t 12 ${NAME}_2.fastq"

