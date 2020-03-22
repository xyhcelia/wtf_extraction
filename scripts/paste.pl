#!/usr/bin/perl
use strict;

my @loci;
my %lcontig;
my %lposi1;
my %lposi2;
my %rcontig;
my %rposi1;
my %rposi2;

open IN,$ARGV[0];
while(<IN>){
   chomp;
   my @arr=split;
   $arr[0]=~/(\S+)-[L|R]-/;
   my $iloci=$1;
   push(@loci,$iloci);
   if($arr[0]=~/-L-/){
       if($arr[1] ne "--"){
           $lcontig{$iloci}=$arr[1];
           $lposi1{$iloci}=$arr[2];
           $lposi2{$iloci}=$arr[3];
       }
   }   
   if($arr[0]=~/-R-/){
       if($arr[1] ne "--"){
           $rcontig{$iloci}=$arr[1];
           $rposi1{$iloci}=$arr[2];
           $rposi2{$iloci}=$arr[3];
       }
   }
}
close IN;

open OUT,">$ARGV[1]";
my %count;
@loci=grep{++$count{ $_ } < 2;}@loci;
my $str;
my $start;
my $end;
for my $iloci(@loci){
    my @posi=($lposi1{$iloci},$lposi2{$iloci},$rposi1{$iloci},$rposi2{$iloci});
    my @sortposi=sort{$a<=>$b}@posi;
    $start=$sortposi[0];
    $end=$sortposi[-1];
    if($lposi1{$iloci}<$lposi2{$iloci} && $rposi1{$iloci}<$rposi2{$iloci} && $lposi2{$iloci}<$rposi2{$iloci}){
    $str="+"; 
    }elsif($lposi1{$iloci}>$lposi2{$iloci} && $rposi1{$iloci}>$rposi2{$iloci} && $rposi1{$iloci}<$lposi2{$iloci}){ 
    $str="-";
    }else{
    $str="NA";
    }
    if($lcontig{$iloci} eq $rcontig{$iloci} && $end-$start<=25000 && $str ne "NA"){
         print OUT "$lcontig{$iloci}\t$start\t$end\t$str\t${iloci}-fl\n";
    }else{
         print $iloci,"\n";
    }
}
close OUT;
