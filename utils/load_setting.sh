#!/bin/bash
# Loading setting

get_usage(){
	cat <<EOF

Usage : $0 [OPTION]

Description
	-h, --help
		show help

	-v, --version
		show version

	-x, --ref [reference seq]
		reference seq name

	-r, --restriction [restriction enzyme]
		restriction enzyme
EOF

}

get_version(){
	echo "${0} version 1.0"
}

SHORT=hvx:r:
LONG=help,version,ref:,restriction:
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
		-x|--ref)
			REF="$2"
			shift 2
			;;
		-r|--restriction)
			RESTRICTION="$2"
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

[ ! -n "${REF}" ] && echo "Please specify reference" && exit 1
[ ! -n "${RESTRICTION}" ] && echo "Please specify restriction" && exit 1

case $REF in
	pombe)	BOWTIE2_INDEX=${HOME}/Genome/data/pombe/2018/pombe
			CHROM_LENGTH=12571820
			FILE_CHROME_LENGTH=${HOME}/Genome/data/pombe/2018/LENGTH.txt
			case $RESTRICTION in 
				MboI)	FILE_enzyme_index=${HOME}/Genome/data/pombe/2018/Sectioning_MboI.txt
						FILE_enzyme_def=${HOME}/Genome/data/pombe/2018/MboI_sites.txt ;;
				MboI-Hinf1)	FILE_enzyme_index=${HOME}/Genome/data/pombe/2018/Sectioning_MboI-Hinf1.txt
							FILE_enzyme_def=${HOME}/Genome/data/pombe/2018/MboI-Hinf1_sites.txt ;;
				MboI-Hinf1-MluCl)	FILE_enzyme_index=${HOME}/Genome/data/pombe/2018/Sectioning_MboI-Hinf1-MluCl.txt
									FILE_enzyme_def=${HOME}/Genome/data/pombe/2018/MboI-Hinf1-MluCl_sites.txt ;;
				NA) ;;
				*)	echo "$RESTRICTION is not registered for $ORGANISM"
					exit ;;
			esac
			;;
	OR74A)	BOWTIE2_INDEX=${HOME}/Genome/data/neurospora_crassa/OR74A/Bowtie2/or74a
			CHROM_LENGTH=40463072
			FILE_CHROME_LENGTH=${HOME}/Genome/data/neurospora_crassa/OR74A/LENGTH_mainChromosome.txt
			case $RESTRICTION in 
				DpnII)	FILE_enzyme_index=${HOME}/Genome/data/neurospora_crassa/OR74A/Sectioning_DpnII.txt
						FILE_enzyme_def=${HOME}/Genome/data/neurospora_crassa/OR74A/DpnII_sites.txt ;;
				HindIII)	FILE_enzyme_index=${HOME}/Genome/data/neurospora_crassa/OR74A/Sectioning_HindIII.txt
						FILE_enzyme_def=${HOME}/Genome/data/neurospora_crassa/OR74A/HindIII_sites.txt ;;
				NA) ;;
				*)	echo "$RESTRICTION is not registered for $ORGANISM"
					exit ;;
			esac
			;;
	hg19)	BOWTIE2_INDEX=${HOME}/Genome/data/human/hg19/Bowtie2/hg19
			CHROM_LENGTH=3095677412
			FILE_CHROME_LENGTH=${HOME}/Genome/data/human/hg19/LENGTH.txt
			case $RESTRICTION in 
				HindIII)	FILE_enzyme_index=${HOME}/Genome/data/human/hg19/Sectioning_HindIII.txt
							FILE_enzyme_def=${HOME}/Genome/data/human/hg19/HindIII_sites.txt ;;
				MboI)	FILE_enzyme_index=${HOME}/Genome/data/human/hg19/Sectioning_MboI.txt
						FILE_enzyme_def=${HOME}/Genome/data/human/hg19/MboI_sites.txt ;;
				NA) ;;
				*)	echo "$RESTRICTION is not registered for $ORGANISM"
					exit ;;
			esac
			;;
	hg19_EBV)	BOWTIE2_INDEX=${HOME}/Genome/data/human/hg19_EBV/hg19_EBV
			CHROM_LENGTH=3157782322
			FILE_CHROME_LENGTH=${HOME}/Genome/data/human/hg19_EBV/LENGTH.txt
			case $RESTRICTION in 
				MboI)	FILE_enzyme_index=${HOME}/Genome/data/human/hg19_EBV/Sectioning_MboI.txt
						FILE_enzyme_def=${HOME}/Genome/data/human/hg19_EBV/MboI_sites.txt;;
				MboI-HinfI)	FILE_enzyme_index=${HOME}/Genome/data/human/hg19_EBV/Sectioning_MboI-HinfI.txt
						FILE_enzyme_def=${HOME}/Genome/data/human/hg19_EBV/MboI-HinfI_sites.txt ;;
				NA) ;;
				*)	echo "$RESTRICTION is not registered for $ORGANISM"
					exit ;;
			esac
			;;
	mm10)	BOWTIE2_INDEX=${HOME}/Genome/data/mouse/mm10/Bowtie2/mm10
			CHROM_LENGTH=2725537669
			FILE_CHROME_LENGTH=${HOME}/Genome/data/mouse/mm10/LENGTH.txt
			case $RESTRICTION in 
				MboI)	FILE_enzyme_index=${HOME}/Genome/data/mouse/mm10/Sectioning_MboI.txt
						FILE_enzyme_def=${HOME}/Genome/data/mouse/mm10/MboI_sites.txt ;;
				NA) ;;
				*)	echo "$RESTRICTION is not registered for $ORGANISM"
					exit ;;
			esac
			;;
	*)	echo "Please specify correct reference name"
		eixt 1 ;;
esac
