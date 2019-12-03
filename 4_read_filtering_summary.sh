#!/bin/bash
# output read filtering summary

get_usage(){
	cat <<EOF

Usage : $0 [OPTION] [target sample names. separated by space]

Description
	-h, --help
		show help

	-v, --version
		show version

	-d, --directory [directory of data]
		log file directory
	
	-o, --out [output file name]
		output file
EOF

}

get_version(){
	echo "${0} version 1.0"
}

SHORT=hvd:o:
LONG=help,version,directory:,out:
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
		-o|--out)
			FILE_OUT="$2"
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

[ ! -n "${DIR_DATA}" ] && echo "Please specify log file directory" && exit 1
[ ! -n "${FILE_OUT}" ] && echo "Please specify output file name" && exit 1

cd ${DIR_DATA}

tmpfile=$(mktemp tmp.XXXXXX)
tmpfile2=$(mktemp tmp.XXXXXX)

FLAG=0
for NAME in "$@"
do
	FILE_map=${NAME}_alignment.log
	FILE_count=${NAME}_read_filtering.log

	if [ $FLAG -eq 0 ]
	then
		echo "SAMPLE name" > $FILE_OUT
		tail -n 4 $FILE_map | cut -f1 -d: >> $FILE_OUT
		cut -f1 -d: ${FILE_count} >> $FILE_OUT
		FLAG=1
	fi
	
	echo $NAME > $tmpfile
	tail -n 4 $FILE_map | cut -f2 >> $tmpfile
	cut -f2 -d: ${FILE_count} | sed 's/^ //' >> $tmpfile
	paste $FILE_OUT $tmpfile > $tmpfile2
	mv $tmpfile2 $FILE_OUT
done

rm $tmpfile
