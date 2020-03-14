#!/usr/bion/perl
use File::Basename;
use strict;

my $dir=dirname(__FILE__);

my @input=split(/\n/,` cat  $dir/../data/wtf.*10cdsID.txt`);
my %hash;
for my$line(@input){
    my @arr=split(/\s+/,$line);
    if(exists($hash{$arr[1]})){$hash{$arr[1]}.=",".$arr[0];}else{$hash{$arr[1]}=$arr[0];}
}

open CDS,"<$dir/../data/Schizosaccharomyces_pombe.ASM294v2.45.gff3.PomBase_CDS.nowtf.noTf.link.txt";
my %posi;
while(<CDS>){
   chomp;
   my @arr=split;
   $posi{$arr[4]}=$arr[0].":".$arr[1]."-".$arr[2];
}
close CDS;

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
my $cdss;
my $cdse;
my $label;
while($end!=1){
    $_=<GBK>;
    chomp;
    $_=~s/^\s+//;
    my@arr=split;
    if($arr[0] eq "CDS" || $arr[0] eq "Tf" || $arr[0] eq "LTR"){
        $_=~/([0-9]+)..([0-9]+)/;
        $cdss=$1;
        $cdse=$2;
        ########capture the label#####
        $_=<GBK>;
        chomp;
        $_=~/label=(\S+)/;
        $label=$1;
        if(exists($hash{$label})){
           print OUT "$cdss\t$cdse\t$label\t$hash{$label}\n";
        }elsif($label eq "Tf" || $label eq "LTR"){
           print OUT "$cdss\t$cdse\t$label\n";
        }else{
           print OUT "$cdss\t$cdse\t$label\t$posi{$label}\n";
        }
    }
    if($arr[0] eq "wtf"){
        $_=~/([0-9]+)..([0-9]+)/;
        print OUT "$1\t$2\twtf\n";
    }
    if($_=~/ORIGIN/){$end=1;}    
}

close GBK;


my @wtfrow=split(/\n/,`cat $id.outtmp|sort -k1,1n|awk '{if(\$3==\"wtf\") print NR-1;}'`);
my @output=split(/\n/,`cat $id.outtmp|sort -k1,1n `);
my @wtf;
my @left;
my @right;
my @itv;
for my $i(0..$#output){
    my @arr=split(/\s+/,$output[$i]);
    if($i<$wtfrow[0]){
        if($arr[2] ne "Tf" && $arr[2] ne "LTR"){push(@left,$arr[2].":".$arr[3]);}else{ push(@left,$arr[2]);}
    }elsif($i>=$wtfrow[0] && $i<=$wtfrow[-1]){
        if($arr[2] eq "wtf"){
            push(@wtf,$arr[0]."..".$arr[1]);
        }elsif($arr[2] eq "CDS"){
            push(@itv,$arr[2].":".$arr[3]);
        }else{
            push(@itv,$arr[2]);
        }
    }else{
        if($arr[2] ne "Tf" && $arr[2] ne "LTR"){push(@right,$arr[2].":".$arr[3]);}else{ push(@right,$arr[2]);}
    }
}

print $id,"\t",join(",",@wtf),"\t",join(",",@itv),"\t",join(",",@left),"\t",join(",",@right),"\n";

`rm $id.outtmp`;
