#!/bin/bash
# Merge map files from technical replicates (remove PCR duplicate)


get_usage(){
	cat <<EOF

Usage : $0 [OPTION]

Description
	-h, --help
		show help

	-v, --version
		show version

	-i, --in [map files]
		map files. Separated with ,. Map files could be gziped

	-d, --directory [data directory]
		directory name of analysis file locate

	-n, --name [merged sample name]
		sample name after merged

EOF

}

get_version(){
	echo "sh ${0} version 1.0"
}

SHORT=hvi:d:n:
LONG=help,version,in:,directory:,name:
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
		-i|--in)
			FILE_IN="$2"
			shift 2
			;;
		-d|--directory)
			DIR_DATA="$2"
			shift 2
			;;
		-n|--name)
			NAME="$2"
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
[ ! -n "${DIR_DATA}" ] && echo "Please specify data directory" && exit 1
[ ! -n "${FILE_IN}" ] && echo "Please specify map files for technical replicates for merging" && exit 1

cd ${DIR_DATA}

#-----------------------------------------------
# Merge map files
#-----------------------------------------------
NUM_MAP_file=$(echo $FILE_IN | tr ',' ' ' | xargs -n1 | wc -l)
if [ $NUM_MAP_file -eq 1 ]; then
	zcat $FILE_IN > ${NAME}.map
else
	DIR_tmp=$(mktemp -d /tmp/tmp_${NAME}.XXXXX)
	perl ${DIR_LIB}/utils/MergeMaps.pl -i "$FILE_IN" -o "${DIR_DATA}/${NAME}.map" -t ${DIR_tmp}
	[ $? -ne 0 ] && echo "Merge failed" && exit 1
	rm -r ${DIR_tmp}
fi


#-----------------------------------------------
# Register to database
#-----------------------------------------------
Rscript --vanilla --slave ${DIR_LIB}/utils/file2database_large.R -i ${NAME}.map --db ${NAME}.db --table map
gzip ${NAME}.map

#-----------------------------------------------
# Summarize read filtering
#-----------------------------------------------
export DIR_DATA SAMPLE=${NAME}
sh ${DIR_LIB}/utils/Count_reads.sh > ${NAME}_read_filtering.log

#-----------------------------------------------
# DNA amount estimation
#-----------------------------------------------
perl ${DIR_LIB}/utils/Count_DNA_amount.pl -i ${NAME}.db -o ${NAME}_DNA_amount.bed
