#!/usr/bin/perl

use strict;
use warnings;
use IO::File;
use Getopt::Std;
use File::Basename;
use Carp qw(croak);
$| = 0;


if((@ARGV != 8 and @ARGV != 10) or $ARGV[0] eq '--help'){
	die "Usage : $0 -i [list of map file (it can be gziped) separated by cammma] -l [chromosome length] -o [file list] -m [mapQ threshold (default:30)] -e [enzyme definition file]\n";
}

my %opt;
getopts("i:l:o:m:e:", \%opt);
my @FILE_maps = split /,/, $opt{i};
my $CHROM_LENGTH = $opt{l};
my $FILE_list = $opt{o};
my ($name, $dir, $ext) = &fileparse($FILE_maps[0], '\..*');
my $TMP_PREFIX = $dir . 'tmp_database_' . $name;
my $MAPQ_threshold = $opt{m};
unless(defined $MAPQ_threshold){
	$MAPQ_threshold = 30;
}

my $FILE_ENZYME_def = $opt{e};


# Resolution to make split files
my $RESOLUTION = $CHROM_LENGTH / 100;

#---------------------------------------
# read restriction information
#---------------------------------------
my %Enzymes;
my %Chromosomes;
{
	my $fh_in = IO::File->new($FILE_ENZYME_def) or die "cannot open $FILE_ENZYME_def: $!";
	while($_ = $fh_in->getline()){
		if(m/^#/){
			next;
		}
		s/\r?\n//;
		my ($number, $chr, $pos, $before, $after) = split /\t/;

		# 番号0の断片を作っておく
		if($number == 1){
			my $id0 = $chr . ':0';
			$Enzymes{$id0} = "1\t$pos";
		}

		my $id = $chr . ':' . $number;
		my $end = $pos + $after;
		my $middle = $pos + $after/2;
		$Enzymes{$id} = "$pos\t$end\t$middle";
		$Chromosomes{$chr} = 1;
	}
	$fh_in->close();
}

my %FH;
my %NAMES;
foreach my $file (@FILE_maps){
	my $fh_in;
	if($file =~ /\.gz/){
		$fh_in = IO::File->new("gzip -dc $file |") or die "cannot open $file: $!";
	}else{
		$fh_in = IO::File->new($file) or die "cannot open $file: $!";
	}
	$fh_in->getline();
	while($_ = $fh_in->getline()){
		s/\r?\n//;
		my ($id, $chr1, $position1, $direction1, $mapQ1, $restNum1, $restLoc1, $uniq1, $chr2, $position2, $direction2, $mapQ2, $restNum2, $restLoc2, $uniq2) = split /\t/;

		if($uniq1 eq 'U' or $uniq2 eq 'U'){next;}
		if($restLoc1 eq 'NA' or $restLoc2 eq 'NA'){next;}
		if($mapQ1 < $MAPQ_threshold or $mapQ2 < $MAPQ_threshold){next;}

		if($direction1 eq '+'){
			$restNum1--;
		}
		if($direction2 eq '+'){
			$restNum2--;
		}

		my $id1 = $chr1 . ':' . $restNum1;
		my $id2 = $chr2 . ':' . $restNum2;
		my ($start1, $end1, $middle1) = split /\t/, $Enzymes{$id1};
		my ($start2, $end2, $middle2) = split /\t/, $Enzymes{$id2};

		# Remove self ligation and potential un-digested pairs (less than 10kb and different directions)
		if($chr1 eq $chr2 and abs($middle1 - $middle2) < 10000 and $direction1 ne $direction2){
			next;
		}

		my $bin = int($middle1 / $RESOLUTION);

		# もしファイルハンドルが未生成だったら作成（合計で約１００個できる可能性あり）
		unless(exists $FH{$chr1}{$bin}){
			my $file_out = $TMP_PREFIX . '_' . $chr1 . '_' . $bin;
			my $f = IO::File->new($file_out, 'w') or die "cannot write $file_out: $!";
			$FH{$chr1}{$bin} = $f;
			$NAMES{$chr1}{$bin} = $file_out;
		}
		$FH{$chr1}{$bin}->print("$chr1\t$start1\t$end1\t$restNum1\t$chr2\t$start2\t$end2\t$restNum2\n");
	}
	$fh_in->close();
}

#---------------------------------------
# make file list
#---------------------------------------
my $fh_list = IO::File->new($FILE_list, 'w') or die "cannot write $FILE_list: $!";
foreach my $chr(sort {$a cmp $b} keys %FH){
	foreach my $bin(sort {$a <=> $b} keys %{$FH{$chr}}){
		$FH{$chr}{$bin}->close();
		$fh_list->print("$NAMES{$chr}{$bin}\n");
	}
}
$fh_list->close();

