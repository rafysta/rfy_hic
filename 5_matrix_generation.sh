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

	-o, --out [output directory]
		default: same as data directory

	-n, --name [sample name]
		sample name
	
	-x, --ref [ex. hg19]
		organism name

	--include [including chromsome list]
		if only specific chromosome should be calculated, specify the list. Separated by ,.

	--exclude [exclude chromsome list]
		list of excluding chromosomes. Separated by ,. Ex. chrM,chrY

	-r, --resolution [ex. 500kb. 2f for fragment resolution]
		map resolution.

	-f, --fragment [TRUE/FALSE]
		fragment resolution analysis mode. resolution value were recognized as fragment number. (default : false)

	--intra [TRUE/FALSE]
		onlyt intra chromosome (TRUE) or all (FALSE). Default : TRUE
	
	-e, --normalization [TRUE/FALSE]
		Do normalization. Default : TRUE

	--normalization_old [TRUE/FALSE]
		Do old ICE normalization. Default : FALSE
	
	-c, --raw [TRUE/FALSE]
		Generate raw matrices. If already exists, overwrite. Default : TRUE
	
	--use_blacklist [TRUE/FALSE]
		use fragment blacklist to eliminate strange ones. Default : TRUE

	-t, --threshold [threshold. default: 10000]
		threshold to remove different direction reads to remove potential self ligation. (default 10kb)
		We use 2kb for 3 restriction enzyme Hi-C

	--max_distance [maximum distance]
		maximum distance of output.
EOF

}

get_version(){
	echo "${0} version 1.0"
}

SHORT=hvd:o:n:x:r:f:e:c:t:
LONG=help,version,directory:,out:,name:,ref:,include:,exclude:,resolution:,fragment:,intra:,normalization:,normalization_old:,raw:,use_blacklist:,threshold:,max_distance:
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
			DIR_OUT="$2"
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
			RESOLUTION=${RESOLUTION_string/Mb/000000}
			RESOLUTION=${RESOLUTION/kb/000}
			RESOLUTION=${RESOLUTION/bp/}
			RESOLUTION=${RESOLUTION/f/}
			shift 2
			;;
		-f|--fragment)
			FLAG_fragment="$2"
			shift 2
			;;
		--intra)
			# if only intra chromosome TRUE, otherwise FALSE (default TRUE)
			FLAG_INTRA="$2"
			shift 2
			;;
		-e|--normalization)
			FLAG_NORM="$2"
			shift 2
			;;
		--normalization_old)
			FLAG_NORM_OLD="$2"
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
		-t|--threshold)
			THRESHOLD_SELF="$2"
			shift 2
			;;
		--max_distance)
			FLAG_dataframe="TRUE"
			THRESHOLD_MAX_DISTANCE="$2"
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
FLAG_NORM_OLD=${FLAG_NORM_OLD:-FALSE}
FLAG_RAW=${FLAG_RAW:-TRUE}
FLAG_blacklist=${FLAG_blacklist:-TRUE}
CHR_include=${CHR_include:-NA}
CHR_exclude=${CHR_exclude:-NA}
[ "$FLAG_blacklist" = "TRUE" ] && [ ! -e ${DIR_DATA}/${NAME}_bad_fragment.txt ] && echo "bad fragment list not exists" && exit 1
FLAG_dataframe=${FLAG_dataframe:-FALSE}
[ "$FLAG_dataframe" = "TRUE" ] && [ ! -n "$THRESHOLD_MAX_DISTANCE" ] && echo "If output as dataframe format, please specify maximum distance" && exit 1
THRESHOLD_SELF=${THRESHOLD_SELF:-10000}
DIR_OUT=${DIR_OUT:-"${DIR_DATA}/${NAME}/${RESOLUTION_string}"}
FLAG_fragment=${FLAG_fragment:-"FALSE"}

[ ! -e ${DIR_OUT} ] && mkdir -p ${DIR_OUT}


DIR_LIB=$(dirname $0)

### load module
# module purge
# module load intel 2> /dev/null
# module load R/4.0.2
module load perl/5.24.2

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
CHRs_list=$(echo ${CHRs[@]} | tr ' ' ',')




#==============================================================
# 最大距離のしきい値がある場合は、dataframe formatで出力する
#==============================================================
if [ $FLAG_dataframe = "TRUE" ]; then
	if [ ! -e ${DIR_OUT}/Raw ]; then
		mkdir -p ${DIR_OUT}/Raw
	fi
	if [ ! -e ${DIR_OUT}/ICE2 ]; then
		mkdir -p ${DIR_OUT}/ICE2
	fi

	if [ "$FLAG_fragment" = "TRUE" ]; then
		PROGRAM_ASSO=${DIR_LIB}/utils/Make_association_from_fragmentdb_onlyIntraChr_lessThanDistance_fragResolution.pl
	else
		PROGRAM_ASSO=${DIR_LIB}/utils/Make_association_from_fragmentdb_onlyIntraChr_lessThanDistance.pl
	fi

	cd ${DIR_DATA};
	for i in $(seq 1 ${#CHRs[@]})
	do
		let index=i-1
		CHR=${CHRs[index]}
		if [ "$FLAG_blacklist" = "TRUE" ] && [ -e ${DIR_DATA}/${NAME}_bad_fragment.txt ]; then
			perl ${PROGRAM_ASSO} -i ${NAME}_fragment.db -o ${DIR_OUT}/Raw/${CHR}_$THRESHOLD_MAX_DISTANCE.txt -c ${CHR} -r ${RESOLUTION} -b ${DIR_DATA}/${NAME}_bad_fragment.txt -m $THRESHOLD_MAX_DISTANCE -t $THRESHOLD_SELF
		else
			perl ${PROGRAM_ASSO} -i ${NAME}_fragment.db -o ${DIR_OUT}/Raw/${CHR}_$THRESHOLD_MAX_DISTANCE.txt -c ${CHR} -r ${RESOLUTION} -m $THRESHOLD_MAX_DISTANCE -t $THRESHOLD_SELF
		fi

		Rscript --vanilla --slave ${DIR_LIB}/utils/Bias_normalization_ICE2_distanceRestrict.R -i ${DIR_OUT}/Raw/${CHR}_$THRESHOLD_MAX_DISTANCE.txt -o ${DIR_OUT}/ICE2/${CHR}_$THRESHOLD_MAX_DISTANCE.txt --times 30

		gzip ${DIR_OUT}/Raw/${CHR}_$THRESHOLD_MAX_DISTANCE.txt
		gzip ${DIR_OUT}/ICE2/${CHR}_$THRESHOLD_MAX_DISTANCE.txt
	done
	exit
fi


#==============================================================
# Raw matrixの作成
#==============================================================
if [ $FLAG_RAW = "TRUE" ]; then
	if [ ! -e ${DIR_OUT}/Raw ]; then
		mkdir -p ${DIR_OUT}/Raw
	fi
	cd ${DIR_DATA};
	if [ $FLAG_INTRA = "TRUE" ]; then
		PRO_RAW_matrix=${DIR_LIB}/utils/Make_association_from_fragmentdb_onlyIntraChr.pl
	else
		PRO_RAW_matrix=${DIR_LIB}/utils/Make_association_from_fragmentdb_allChromosome.pl
	fi
	if [ "$FLAG_blacklist" = "TRUE" ] && [ -e ${DIR_DATA}/${NAME}_bad_fragment.txt ]; then
		perl $PRO_RAW_matrix -i ${NAME}_fragment.db -o ${DIR_OUT}/Raw/  -r ${RESOLUTION} -b ${DIR_DATA}/${NAME}_bad_fragment.txt -c $CHRs_list -t ${THRESHOLD_SELF}
	else
		perl $PRO_RAW_matrix -i ${NAME}_fragment.db -o ${DIR_OUT}/Raw/  -r ${RESOLUTION} -c $CHRs_list -t ${THRESHOLD_SELF}
	fi

	if [ $FLAG_INTRA = "TRUE" ]; then
		cd ${DIR_OUT}/Raw
		for i in $(seq 1 ${#CHRs[@]})
		do
			let index=i-1
			CHR=${CHRs[index]}
			Rscript --slave --vanilla ${DIR_LIB}/utils/Convert_matrix_to_object.R -i ${CHR}.matrix
		done
	else
		cd ${DIR_OUT}/Raw
		Rscript --slave --vanilla ${DIR_LIB}/utils/Convert_matrix_to_object.R -i ALL.matrix
	fi
fi


#==============================================================
# binごとのinter-chromosomeのデータを計算
#==============================================================
if [ $FLAG_INTRA = "TRUE" ] && [ $FLAG_NORM = "TRUE" ]; then
	if [ ! -e ${DIR_OUT}/InterBin ]; then
		mkdir -p ${DIR_OUT}/InterBin
	fi
	[ ! -e ${DIR_OUT}/InterBin/${CHRs[0]}.txt ] && cd ${DIR_DATA} && perl ${DIR_LIB}/utils/Make_association_from_fragmentdb_interChromosome_perBin.pl -i ${NAME}_fragment.db -o ${DIR_OUT}/InterBin/ -r ${RESOLUTION} -b ${NAME}_bad_fragment.txt
fi


#==============================================================
# ICE normalization
#==============================================================
if [ $FLAG_NORM = "TRUE" ]; then
	if [ ! -e ${DIR_OUT}/ICE2 ]; then
		mkdir -p ${DIR_OUT}/ICE2
	fi
	cd ${DIR_OUT}
	if [ $FLAG_INTRA = "TRUE" ]; then
		for i in $(seq 1 ${#CHRs[@]})
		do
			let index=i-1
			CHR=${CHRs[index]}
			[ ! -e ${DIR_OUT}/ICE2/${CHR}.rds ] && Rscript --vanilla --slave ${DIR_LIB}/utils/Bias_normalization_ICE2.R -i ${DIR_OUT}/Raw/${CHR}.matrix -o ${DIR_OUT}/ICE2/${CHR}.matrix --inter ${DIR_OUT}/InterBin/${CHR}.txt --times 30 && Rscript --slave --vanilla ${DIR_LIB}/utils/Convert_matrix_to_object.R -i ${DIR_OUT}/ICE2/${CHR}.matrix
		done
	else
		Rscript --vanilla --slave ${DIR_LIB}/utils/Bias_normalization_ICE2.R -i ${DIR_OUT}/Raw/ALL.matrix -o ${DIR_OUT}/ICE2/ALL.matrix --times 30 && Rscript --slave --vanilla ${DIR_LIB}/utils/Convert_matrix_to_object.R -i ${DIR_OUT}/ICE2/ALL.matrix
	fi
fi


if [ $FLAG_NORM_OLD = "TRUE" ]; then
	if [ ! -e ${DIR_OUT}/ICE ]; then
		mkdir ${DIR_OUT}/ICE
	fi
	cd ${DIR_OUT}
	if [ $FLAG_INTRA = "TRUE" ]; then
		for i in $(seq 1 ${#CHRs[@]})
		do
			let index=i-1
			CHR=${CHRs[index]}
			[ ! -e ${DIR_OUT}/ICE/${CHR}.rds ] && Rscript --vanilla --slave ${DIR_LIB}/utils/Bias_normalization.R -i ${DIR_OUT}/Raw/${CHR}.matrix -o ${DIR_OUT}/ICE/${CHR}.matrix --inter ${DIR_OUT}/InterBin/${CHR}.txt --times 30 && Rscript --slave --vanilla ${DIR_LIB}/utils/Convert_matrix_to_object.R -i ${DIR_OUT}/ICE/${CHR}.matrix
		done
	else
		Rscript --vanilla --slave ${DIR_LIB}/utils/Bias_normalization.R -i ${DIR_OUT}/Raw/ALL.matrix -o ${DIR_OUT}/ICE/ALL.matrix --times 30 && Rscript --slave --vanilla ${DIR_LIB}/utils/Convert_matrix_to_object.R -i ${DIR_OUT}/ICE/ALL.matrix
	fi
fi
