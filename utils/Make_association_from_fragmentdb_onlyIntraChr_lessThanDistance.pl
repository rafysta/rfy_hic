#!/usr/bin/perl
# 2015/07/02 filtered dbからの読み取りに変更
# 2015/03/28 dbファイルからデータを読み取って相互作用ファイルを作成する
# 2020-03-03 ある距離以下に限定したものをmatrixでなくてdata.frame形式で出力
# 2021-04-14 ある距離以上のデータの合計を1次元方式で出力。different direction readsのcut-offも導入。single fragment resolutionにも対応。

use strict;
use warnings;
use IO::File;
use Getopt::Std;
use File::Basename;
use Carp qw(croak);
$| = 0;

use DBI;

if((@ARGV != 12 and @ARGV != 14 and @ARGV != 16) or $ARGV[0] eq '--help'){
	die "Usage : $0 -i [detabase files] -o [output file name] -r [resolution. 1 for single fragment resolution] -c [target single chromosome] -b [black list of fragment] -m [maximum distance] -t [different direction reads cut-off threshold]\n";
}

my %opt;
getopts("i:o:r:c:b:m:t:", \%opt);
my $FILE_database = $opt{i};
my $FILE_out = $opt{o};
my $CHROMOSOME = $opt{c};
my $Resolution = $opt{r};
my $FILE_black = $opt{b};
my $MAX_DISTANCE = $opt{m};
my $THRESHOLD_SELF = $opt{t};
unless(defined $THRESHOLD_SELF){
	$THRESHOLD_SELF = 10000;
}

# 全データを入れる変数
my %data;

#---------------------------------------
# fragmentのblack listを読み込む
#---------------------------------------
my %Black;
if(defined $FILE_black){
	my $fh_in = IO::File->new($FILE_black ) or die "cannot open $FILE_black: $!";
	$fh_in->getline();
	while($_ = $fh_in->getline()){
		s/\r?\n//;
		my ($chr, $fragID) = split /\t/;
		$Black{"$chr\t$fragID"} = 1;
	}
	$fh_in->close();
}


my $dbh = DBI->connect("dbi:SQLite:dbname=$FILE_database");


#---------------------------------------
# collect information
#---------------------------------------
my $sth_data;

$sth_data = $dbh->prepare("select chr1, start1, end1, fragNum1, chr2, start2, end2, fragNum2, score from fragment where chr1=='$CHROMOSOME' or chr2=='$CHROMOSOME';");
$sth_data->execute();
while(my $ref = $sth_data->fetchrow_arrayref()){
	my ($chr1, $start1, $end1, $frag1, $chr2, $start2, $end2, $frag2, $score) = @$ref;

	if($chr1 eq $chr2 and $start1 > $start2){
		($start1, $end1, $frag1, $start2, $end2, $frag2) = ($start2, $end2, $frag2, $start1, $end1, $frag1);
	}

	# fragmentがblack listに含まれていたら計算しない
	if(exists $Black{"$chr1\t$frag1"}){
		next;
	}
	if(exists $Black{"$chr2\t$frag2"}){
		next;
	}

	my $middle1 = ($start1 + $end1) / 2;
	my $middle2 = ($start2 + $end2) / 2;
	my $distance = abs($middle1 - $middle2);

	# max distanceよりも離れている組み合わせ
	my $FLAG_longDistance = 0;
	if($distance > $MAX_DISTANCE or $chr1 ne $chr2){
		$FLAG_longDistance = 1;
	}

	# THRESHOLD_SELF以内の距離だった場合には、scoreを２倍にする
	# (THRESHOLD_SELF以内については、同じ向きのデータしか無いから)
	if($chr1 eq $chr2 and $distance < $THRESHOLD_SELF){
		$score = $score * 2;
	}

	my ($id1a, $id1b, $id2a, $id2b);
	if($Resolution == 1){
		$id1a = $chr1 . ":" . $start1 . ":" . $end1;
		$id1b = $chr1 . ":" . $start1 . ":" . $end1;
		$id2a = $chr2 . ":" . $start2 . ":" . $end2;
		$id2b = $chr2 . ":" . $start2 . ":" . $end2;
	}else{
		my $bin1a = int($start1/$Resolution) * $Resolution;
		my $bin1b = int($end1/$Resolution) * $Resolution;
		my $bin2a = int($start2/$Resolution) * $Resolution;
		my $bin2b = int($end2/$Resolution) * $Resolution;
		$id1a = $chr1 . ":" . $bin1a . ":" . ($bin1a + $Resolution - 1);
		$id1b = $chr1 . ":" . $bin1b . ":" . ($bin1b + $Resolution - 1);
		$id2a = $chr2 . ":" . $bin2a . ":" . ($bin2a + $Resolution - 1);
		$id2b = $chr2 . ":" . $bin2b . ":" . ($bin2b + $Resolution - 1);
	}



	# count data（既に左側が小さいということは保証されている）
	if($FLAG_longDistance == 0){
		$score = $score / 4;
		$data{"$id1a\t$id2a"} += $score;
		$data{"$id1a\t$id2b"} += $score;
		$data{"$id1b\t$id2a"} += $score;
		$data{"$id1b\t$id2b"} += $score;
	}else{
		$score = $score / 2;
		if($chr1 eq $CHROMOSOME){
			$data{"$id1a\tlong_distance"} += $score;
			$data{"$id1b\tlong_distance"} += $score;
		}
		if($chr2 eq $CHROMOSOME){
			$data{"$id2a\tlong_distance"} += $score;
			$data{"$id2b\tlong_distance"} += $score;
		}
	}

}
$sth_data->finish();
$dbh->disconnect();



#---------------------------------------
# output
#---------------------------------------
my $fh_out = IO::File->new($FILE_out, 'w') or die "cannot write $FILE_out: $!";
$fh_out->print("loc1\tloc2\tscore\n");
foreach my $key(keys %data){
	$fh_out->printf("%s\t%.2f\n", $key, $data{$key});
}
$fh_out->close();



