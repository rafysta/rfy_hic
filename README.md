# rfy_hic
This is a pipeline for automating Hi-C analysis from fastq files to make matrices.

## Run the program in the following order:
1. 2_make_map_file.sh - alignment of fastq files with bowtie2 for various filtering. map.txt.gz file is created, which can also be used to create Juicer files, too.
2. 3_make_fragment_db.sh - Count interactions for each restriction enzyme fragment. Checks for bias and creates a database for fast data access.
3. 5_matrix_generation.sh - creates three matrices: Raw, ICE, and inter-chromosome. Inter-chromosome matrices are used for ICE normalization.

## The following command summarizes the results of Hi-C's read filtering.
- 4_read_filtering_summary.sh

## The following program can be used to summarize data from technical replicate. (except for PCR duplication)
- merge_map_technical_replicate.sh


## Examples of how to use each script
All script parameters can be run with `-h` or `--help` to see the details.
For example
```
sh 2_make_map_file.sh --help
```
on the terminal. You will see a list of all parameters and a detailed explanation of how to specify them.


Specific examples are shown below.
### 2_make_map_file.sh
```
sh 2_make_map_file.sh --name ${NAME} --ref mm10 --restriction $ENZYME -d ${DIR_DATA} --f1 ${DIR_RAW}/${FILE_fastq1} --f2 ${DIR_RAW}/${FILE_fastq2}
```
If you want the threshold for self ligation to be something other than 10 kb, specify the `--threshold` option.

### 3_make_fragment_db.sh
```
sh 3_make_fragment_db.sh --directory ${DIR_DATA} --in ${FILE_IN} --name ${NAME} --ref mm10 --restriction $ENZYME
```
If you want the threshold for self ligation to be something other than 10 kb, specify the `--threshold` option.

### 4_read_filtering_summary.sh
```
sh 4_read_filtering_summary.sh -d ${DIR_DATA} -o Read_summary_sample.txt $(sqlite3 $FILE_DB "select name from sample" | xargs)
```

### 5_matrix_generation.sh
```
sh 5_matrix_generation.sh --directory ${DIR_DATA} --name ${NAME} --ref mm10 --resolution ${RESOLUTION} --intra TRUE
```
If you set the threshold for self ligation at 2 and 3 to something other than 10 kb, specify the `--threshold` option.



