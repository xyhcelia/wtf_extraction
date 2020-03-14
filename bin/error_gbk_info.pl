#!/usr/bin/perl
use File::Basename;
use strict;



my @suffixlist = qw(.gbk);
my$id = basename($ARGV[0],@suffixlist);
my $start=0;
my $end=0;

open GBK,$ARGV[0];
while($start!=1){
    $_=<GBK>;
    chomp;
    if($_=~/mol_type="genomic DNA"/){ $start=1; }else{}
}

open OUT,">$id.outtmp";
while($end!=1){
    $_=<GBK>;
    chomp;
    $_=~s/^\s+//;
    my@arr=split;
    if($arr[0] eq "wtf"){
        $_=~/([0-9]+)..([0-9]+)/;
        print OUT "$1\t$2\twtf\n";
    }
    if($_=~/ORIGIN/){$end=1;}    
}

close GBK;


my @alpha=qw(a b c d e f g h i j k l m n o p q r s t u v w x y z);
$id=~/(\S+)-flrvs/;
my $loci=$1;
my @output=split(/\n/,`cat $id.outtmp|sort -k1,1n `);
if(@output == 1){
     my @arr=split(/\s+/,$output[0]);
     print "$id\t$arr[0]\t$arr[1]\t+\t$loci\n";
}else{
    for my $i(0..$#output){
        my @arr=split(/\s+/,$output[$i]);
        print "$id\t$arr[0]\t$arr[1]\t+\t$loci$alpha[$i]\n";
    }
}
`rm $id.outtmp`;
