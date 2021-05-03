# rfy_hic
fastqファイルから、matricesを作るためのHi-C自動解析パイプラインです。
質問等は、谷澤英樹 ( rafysta@gmail.com ) まで。

## 以下の順番で、プログラムを実行してください。
1. 2_make_map_file.sh - fastqファイルをbowtie2でalignmentして各種filteringを行います。map.txt.gzファイルができるので、このファイルを使ってJuicerファイルの作成もできます。
2. 3_make_fragment_db.sh - 制限酵素断片ごとに相互作用をカウントします。バイアスのチェック、高速データアクセスのためのデータベースの作成などを行います。
3. 5_matrix_generation.sh - Raw, ICE, inter-chromosomeの３種類のマトリックスを作ります。inter-chromosomeは、ICE normalizationで用います。

## 以下のコマンドでHi-Cのread filteringの結果をまとめることができます。
- 4_read_filtering_summary.sh

## Technical replicateのデータをまとめるには、以下のプログラムを使うことができます。(PCR duplicationを除きます
- merge_map_technical_replicate.sh


## 各スクリプトの使い方例
すべてのスクリプトのパラメータは、`-h` もしくは、`--help`をつけて実行することで、詳細を確認することができます。
例えば
```
sh 2_make_map_file.sh --help
```
とターミナル上で実行してください。すべてのパラメータリストと、指定の仕方についての詳しい説明が表示されます。


具体的な実行例を以下に示します。
### 2_make_map_file.sh
```
sh 2_make_map_file.sh --name ${NAME} --ref mm10 --restriction $ENZYME -d ${DIR_DATA} --f1 ${DIR_RAW}/${FILE_fastq1} --f2 ${DIR_RAW}/${FILE_fastq2}
```

### 3_make_fragment_db.sh
```
sh 3_make_fragment_db.sh --directory ${DIR_DATA} --in ${FILE_IN} --name ${NAME} --ref mm10 --restriction $ENZYME
```

### 4_read_filtering_summary.sh
```
sh 4_read_filtering_summary.sh -d ${DIR_DATA} -o Read_summary_sample.txt $(sqlite3 $FILE_DB "select name from sample" | xargs)
```

### 5_matrix_generation.sh
```
sh 5_matrix_generation.sh --directory ${DIR_DATA} --name ${NAME} --ref mm10 --resolution ${RESOLUTION} --intra TRUE --threshold 10000
```




