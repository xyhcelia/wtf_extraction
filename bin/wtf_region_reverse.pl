#!/usr/bin/perl
use strict;
use File::Basename;

my $dir=dirname(__FILE__);

open STR,$ARGV[0];
my %strand;
while(<STR>){
     chomp;
     my @arr=split;
     $strand{$arr[3]."-fl"}=$arr[5];
}
close STR;

system("perl /$dir/fasta_length.pl $ARGV[1] >$ARGV[1].lenth");

open IN,"<$ARGV[1].lenth";
open OUT,">$ARGV[3]";
while(<IN>){
      chomp;
      my @arr=split;
      print OUT "$arr[0]\t1\t$arr[1]\t$strand{$arr[0]}\t$ARGV[2]_$arr[0]rvs\n";
}

close IN;
close OUT;

system("rm $ARGV[1].lenth");
