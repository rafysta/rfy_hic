#!/bin/bash
# Generate matrix

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
	
	-x, --ref [ex. hg19]
		organism name

	--include [including chromsome list]
		if only specific chromosome should be calculated, specify the list. Separated by ,.

	--exclude [exclude chromsome list]
		list of excluding chromosomes. Separated by ,. Ex. chrM,chrY

	-r, --resolution [ex. 500kb]
		map resolution

	-t, --intra [TRUE/FALSE]
		onlyt intra chromosome (TRUE) or all (FALSE). Default : TRUE
	
	-e, --normalization [TRUE/FALSE]
		Do normalization. Default : TRUE
	
	-c, --raw [TRUE/FALSE]
		Generate raw matrices. If already exists, overwrite. Default : TRUE
	
	--use_blacklist [TRUE/FALSE]
		use fragment blacklist to eliminate strange ones. Default : TRUE
EOF

}

get_version(){
	echo "${0} version 1.0"
}

SHORT=hvd:n:x:r:t:e:c:
LONG=help,version,directory:,name:,ref:,include:,exclude:,resolution:,intra:,normalization:,raw:,use_blacklist:
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
		-x|--ref)
			REF="$2"
			shift 2
			;;
		--include)
			CHR_include="$2"
			shift 2
			;;
		--exclude)
			CHR_exclude="$2"
			shift 2
			;;
		-r|--resolution)
			RESOLUTION_string="$2"
			RESOLUTION=${RESOLUTION_string/kb/000}
			RESOLUTION=${RESOLUTION/bp/}
			shift 2
			;;
		-t|--intra)
			# if only intra chromosome TRUE, otherwise FALSE (default TRUE)
			FLAG_INTRA="$2"
			shift 2
			;;
		-e|--normalization)
			FLAG_NORM="$2"
			shift 2
			;;
		-c|--raw)
			# whether making raw matrices or not
			FLAG_RAW="$2"
			shift 2
			;;
		--use_blacklist)
			FLAG_blacklist="$2"
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

[ ! -n "${NAME}" ] && echo "Please specify NAME" && exit 1
[ ! -n "${RESOLUTION}" ] && echo "Please specify resolution" && exit 1
[ ! -n "${DIR_DATA}" ] && echo "Please specify data directory" && exit 1
[ ! -n "${REF}" ] && echo "Please specify ref" && exit 1
FLAG_INTRA=${FLAG_INTRA:-TRUE}
FLAG_NORM=${FLAG_NORM:-TRUE}
FLAG_RAW=${FLAG_RAW:-TRUE}
FLAG_blacklist=${FLAG_blacklist:-TRUE}
CHR_include=${CHR_include:-NA}
CHR_exclude=${CHR_exclude:-NA}
[ "$FLAG_blacklist" = "TRUE" ] && [ ! -e ${DIR_DATA}/${NAME}_bad_fragment.txt ] && echo "bad fragment list not exists" && exit 1


TIME_STAMP=$(date +"%Y-%m-%d")
DIR_LIB=$(dirname $0)


#-----------------------------------------------
# Load setting
#-----------------------------------------------
source ${DIR_LIB}/utils/load_setting.sh -x $REF -r NA


#-----------------------------------------------
# Load chromosome length
#-----------------------------------------------
CHR_TABLE=$(Rscript --vanilla --slave ${DIR_LIB}/utils/Chromosome_length.R --in $FILE_CHROME_LENGTH --include $CHR_include --exclude $CHR_exclude)
CHRs=($(echo $CHR_TABLE | xargs -n1 | awk 'NR==1' | tr ',' ' '))
LENGTH=($(echo $CHR_TABLE | xargs -n1 | awk 'NR==2' | tr ',' ' '))

#==============================================================
# matrix用のフォルダの作成
#==============================================================
if [ ! -e ${DIR_DATA}/${NAME} ]; then
	mkdir ${DIR_DATA}/${NAME}
fi
if [ ! -e ${DIR_DATA}/${NAME}/${RESOLUTION_string} ]; then
	mkdir ${DIR_DATA}/${NAME}/${RESOLUTION_string}
fi
if [ ! -e ${DIR_DATA}/${NAME}/${RESOLUTION_string}/Raw ]; then
	mkdir ${DIR_DATA}/${NAME}/${RESOLUTION_string}/Raw
fi
if [ ! -e ${DIR_DATA}/${NAME}/${RESOLUTION_string}/ICE ]; then
	mkdir ${DIR_DATA}/${NAME}/${RESOLUTION_string}/ICE
fi
if [ ! -e ${DIR_DATA}/${NAME}/${RESOLUTION_string}/ICE2 ]; then
	mkdir ${DIR_DATA}/${NAME}/${RESOLUTION_string}/ICE2
fi
if [ ! -e ${DIR_DATA}/${NAME}/${RESOLUTION_string}/InterBin ]; then
	mkdir ${DIR_DATA}/${NAME}/${RESOLUTION_string}/InterBin
fi


#==============================================================
# Raw matrixの作成
#==============================================================
if [ $FLAG_RAW = "TRUE" ]; then
	cd ${DIR_DATA};
	if [ $FLAG_INTRA = "TRUE" ]; then
		PRO_RAW_matrix=${DIR_LIB}/utils/Make_association_from_fragmentdb_onlyIntraChr.pl
	else
		PRO_RAW_matrix=${DIR_LIB}/utils/Make_association_from_fragmentdb_allChromosome.pl
	fi
	if [ "$FLAG_blacklist" = "TRUE" ] && [ -e ${DIR_DATA}/${NAME}_bad_fragment.txt ]; then
		perl $PRO_RAW_matrix -i ${NAME}_fragment.db -o ${NAME}/${RESOLUTION_string}/Raw/  -r ${RESOLUTION} -b ${NAME}_bad_fragment.txt
	else
		perl $PRO_RAW_matrix -i ${NAME}_fragment.db -o ${NAME}/${RESOLUTION_string}/Raw/  -r ${RESOLUTION}
	fi

	if [ $FLAG_INTRA = "TRUE" ]; then
		cd ${DIR_DATA}/${NAME}/${RESOLUTION_string}/Raw
		for i in $(seq 1 ${#CHRs[@]})
		do
			let index=i-1
			CHR=${CHRs[index]}
			Rscript --slave --vanilla ${DIR_LIB}/utils/Convert_matrix_to_object.R -i ${CHR}.matrix
		done
	else
		cd ${DIR_DATA}/${NAME}/${RESOLUTION_string}/Raw
		Rscript --slave --vanilla ${DIR_LIB}/utils/Convert_matrix_to_object.R -i ALL.matrix
	fi
fi


#==============================================================
# binごとのinter-chromosomeのデータを計算
#==============================================================
if [ $FLAG_INTRA = "TRUE" ]; then
	[ ! -e ${DIR_DATA}/${NAME}/${RESOLUTION_string}/InterBin/${CHRs[0]}.txt ] && cd ${DIR_DATA} && perl ${DIR_LIB}/utils/Make_association_from_fragmentdb_interChromosome_perBin.pl -i ${NAME}_fragment.db -o ${NAME}/${RESOLUTION_string}/InterBin/  -r ${RESOLUTION} -b ${NAME}_bad_fragment.txt
fi

#==============================================================
# ICE normalization
#==============================================================
if [ $FLAG_NORM = "TRUE" ]; then
	cd ${DIR_DATA}/${NAME}/${RESOLUTION_string}
	if [ $FLAG_INTRA = "TRUE" ]; then
		for i in $(seq 1 ${#CHRs[@]})
		do
			let index=i-1
			CHR=${CHRs[index]}
			[ ! -e ${DIR_DATA}/${NAME}/${RESOLUTION_string}/ICE/${CHR}.rds ] && Rscript --vanilla --slave ${DIR_LIB}/utils/Bias_normalization.R -i Raw/${CHR}.matrix -o ICE/${CHR}.matrix --inter InterBin/${CHR}.txt --times 30 && Rscript --slave --vanilla ${DIR_LIB}/utils/Convert_matrix_to_object.R -i ICE/${CHR}.matrix
			[ ! -e ${DIR_DATA}/${NAME}/${RESOLUTION_string}/ICE2/${CHR}.rds ] && Rscript --vanilla --slave ${DIR_LIB}/utils/Bias_normalization_ICE2.R -i Raw/${CHR}.matrix -o ICE2/${CHR}.matrix --inter InterBin/${CHR}.txt --times 30 && Rscript --slave --vanilla ${DIR_LIB}/utils/Convert_matrix_to_object.R -i ICE2/${CHR}.matrix
		done
	else
		Rscript --vanilla --slave ${DIR_LIB}/utils/Bias_normalization.R -i Raw/ALL.matrix -o ICE/ALL.matrix --times 30 && Rscript --slave --vanilla ${DIR_LIB}/utils/Convert_matrix_to_object.R -i ICE/ALL.matrix
		Rscript --vanilla --slave ${DIR_LIB}/utils/Bias_normalization_ICE2.R -i Raw/ALL.matrix -o ICE2/ALL.matrix --times 30 && Rscript --slave --vanilla ${DIR_LIB}/utils/Convert_matrix_to_object.R -i ICE2/ALL.matrix
	fi
fi

