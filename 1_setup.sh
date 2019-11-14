#!/bin/bash
# Set-up environment


if hash module 2>/dev/null; then
	module load bowtie2
fi

DIR_LIB=$(dirname $0)
source ${DIR_LIB}/utils/load_setting.sh -x $REF -r $RESTRICTION


### Make bowtie2 index
# cd ${BOWTIE2_INDEXES}
# bowtie2-build --threads 12 -f ../all.fa hg19


