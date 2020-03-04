#!/usr/bin/perl
# 2015/07/02 filtered dbからの読み取りに変更
# 2015/03/28 dbファイルからデータを読み取って相互作用ファイルを作成する
# 2020-03-03 ある距離以下に限定したものをmatrixでなくてdata.frame形式で出力

use strict;
use warnings;
use IO::File;
use Getopt::Std;
use File::Basename;
use Carp qw(croak);
$| = 0;

use DBI;

if((@ARGV != 10 and @ARGV != 12 and @ARGV != 14) or $ARGV[0] eq '--help'){
	die "Usage : $0 -i [detabase files] -o [output prefix] -r [resolution] -d [distance normalize file] -b [black list of fragment] -c [chromosome list] -m [maximum distance]\n";
}

my %opt;
getopts("i:o:r:d:b:c:m:", \%opt);
my $FILE_database = $opt{i};
my $FILE_out_prefix = $opt{o};
my $Resolution = $opt{r};
my $FILE_distance = $opt{d};
my $FILE_black = $opt{b};
my @chromosomes = split /,/, $opt{c};
my $MAX_DISTANCE = $opt{m};


# 全データを入れる変数
my %data;


#---------------------------------------
# Distance curveを読み込む
#---------------------------------------
my %AVE;
if(defined $FILE_distance){
	my $fh_dis = IO::File->new($FILE_distance) or die "cannot open $FILE_distance: $!";
	while($_ = $fh_dis->getline()){
		s/\r?\n//;
		my ($chr1, $chr2, $d, $score, $probability) = split /\t/;
		$AVE{"$chr1\t$chr2\t$d"} = $probability;
	}
	$fh_dis->close();
}else{
	print "Distance curve file is required\n";
	exit;
}

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

# intra-chromosomeのデータのみを取得する
$sth_data = $dbh->prepare("select chr1, start1, end1, fragNum1, chr2, start2, end2, fragNum2, score from fragment where chr1=chr2;");
$sth_data->execute();
while(my $ref = $sth_data->fetchrow_arrayref()){
	my ($chr1, $start1, $end1, $frag1, $chr2, $start2, $end2, $frag2, $score) = @$ref;

	# 隣同士のfragmentはカウントしない
	if(abs($frag1 - $frag2) < 2){
		next;
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

	# max distanceよりも離れている組み合わせは考慮しない
	if($distance > $MAX_DISTANCE){
		next;
	}

	# 10kb以内の距離だった場合には、scoreを２倍にする
	# (10kb以内については、同じ向きのデータしか無いから)
	if($distance < 10000){
		$score = $score * 2;
	}


	my $distanceAve1 = 1;
	my $distanceAve2 = 1;
	my $distanceAve3 = 1;
	my $distanceAve4 = 1;

	if(defined $FILE_distance){
		# calculate distance average score
		my ($logDistance1, $logDistance2, $logDistance3, $logDistance4) = (-1, -1, -1, -1);

		my $distance1 = abs($start1 - $start2);
		my $distance2 = abs($start1 - $end2);
		my $distance3 = abs($end1 - $start2);
		my $distance4 = abs($end1 - $end2);

		$distance1 = int($distance1 / 100) * 100 + 50;
		$distance2 = int($distance2 / 100) * 100 + 50;
		$distance3 = int($distance3 / 100) * 100 + 50;
		$distance4 = int($distance4 / 100) * 100 + 50;

		if($distance1 < 50000){
			$logDistance1 = $distance1;
		}else{
			$logDistance1 = exp(sprintf("%.3f", log($distance1)));
			$logDistance1 = sprintf("%d", $logDistance1);
		}
		if($distance2 < 50000){
			$logDistance2 = $distance2;
		}else{
			$logDistance2 = exp(sprintf("%.3f", log($distance2)));
			$logDistance2 = sprintf("%d", $logDistance2);
		}
		if($distance3 < 50000){
			$logDistance3 = $distance3;
		}else{
			$logDistance3 = exp(sprintf("%.3f", log($distance3)));
			$logDistance3 = sprintf("%d", $logDistance3);
		}
		if($distance4 < 50000){
			$logDistance4 = $distance4;
		}else{
			$logDistance4 = exp(sprintf("%.3f", log($distance4)));
			$logDistance4 = sprintf("%d", $logDistance4);
		}

		unless(exists $AVE{"$chr1\t$chr2\t$logDistance1"}){
			die "$chr1\t$chr2\t$logDistance1\n";
		}
		unless(exists $AVE{"$chr1\t$chr2\t$logDistance2"}){
			die "$chr1\t$chr2\t$logDistance2\n";
		}
		unless(exists $AVE{"$chr1\t$chr2\t$logDistance3"}){
			die "$chr1\t$chr2\t$logDistance3\n";
		}
		unless(exists $AVE{"$chr1\t$chr2\t$logDistance4"}){
			die "$chr1\t$chr2\t$logDistance4\n";
		}

		$distanceAve1 = $AVE{"$chr1\t$chr2\t$logDistance1"};
		$distanceAve2 = $AVE{"$chr1\t$chr2\t$logDistance2"};
		$distanceAve3 = $AVE{"$chr1\t$chr2\t$logDistance3"};
		$distanceAve4 = $AVE{"$chr1\t$chr2\t$logDistance4"};
	}

	my $bin1a = int($start1/$Resolution) * $Resolution;
	my $bin1b = int($end1/$Resolution) * $Resolution;
	my $bin2a = int($start2/$Resolution) * $Resolution;
	my $bin2b = int($end2/$Resolution) * $Resolution;

	my $id1a = $chr1 . ":" . $bin1a;
	my $id1b = $chr1 . ":" . $bin1b;
	my $id2a = $chr2 . ":" . $bin2a;
	my $id2b = $chr2 . ":" . $bin2b;


	# 4つの組み合わせに均等に分配する
	$score = $score / 4;


	# count data（既に左側が小さいということは保証されている）
	$data{"$id1a\t$id2a"}{'Raw'} += $score;
	$data{"$id1a\t$id2b"}{'Raw'} += $score;
	$data{"$id1b\t$id2a"}{'Raw'} += $score;
	$data{"$id1b\t$id2b"}{'Raw'} += $score;

	$data{"$id1a\t$id2a"}{'norm'} += $score / $distanceAve1;
	$data{"$id1a\t$id2b"}{'norm'} += $score / $distanceAve2;
	$data{"$id1b\t$id2a"}{'norm'} += $score / $distanceAve3;
	$data{"$id1b\t$id2b"}{'norm'} += $score / $distanceAve4;

}
$sth_data->finish();
$dbh->disconnect();



#---------------------------------------
# output
#---------------------------------------
foreach my $chr(@chromosomes){
	my $FILE_out = $FILE_out_prefix . $chr . "_dataframe.txt";
	my $fh_out = IO::File->new($FILE_out, 'w') or die "cannot write $FILE_out: $!";

	# header
	$fh_out->print("loc1\tloc2\ttscore\tnormScore\n");

	foreach my $key(keys %data){
		my ($loc1, $loc2) = split /\t/, $key;
		my ($c1, $s1) = split /:/, $loc1;
		my ($c2, $s2) = split /:/, $loc2;
		$fh_out->printf("%s:%d:%d\t%s:%d:%d\t%f.2f\t%.3f\n", $c1, $s1, $s1 + $Resolution - 1, $c2, $s2, $s2 + $Resolution - 1, $data{$key}{'Raw'}, $data{$key}{'norm'});
	}
	$fh_out->close();
}


