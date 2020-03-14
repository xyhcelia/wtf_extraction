#!/usr/bin/perl
use strict;

open LEN,$ARGV[0];
my %len;
while(<LEN>){
     chomp;
     my @arr=split;
     $len{$arr[0]}=$arr[1];
}
close LEN;

my @alpha=qw(a b c d e f g h i j k l m n o p q r s t u v w x y z);

my $loci;
my @ends;
my %count;
open IN,$ARGV[1];
open OUT,">$ARGV[2]";
while(<IN>){
      chomp;
      my @arr=split;
      if($arr[0]=~ "#"){next;}
      if($arr[0] eq $loci){
          if($arr[8]==1){
              if(grep { $_ eq $arr[6] } @ends){
                  $count{$arr[6]}+=1;
              }else{
                  push(@ends,$arr[6]);
                  $count{$arr[6]}=1;
              } 
          }
          if($arr[9]==$len{$arr[1]}){
              if(grep { $_ eq $arr[7] } @ends){
                  $count{$arr[7]}+=1;
              }else{
                  push(@ends,$arr[7]);
                  $count{$arr[7]}=1;
              }
          }
      }elsif($arr[0] ne $loci ){
          @ends=sort{$a<=>$b}@ends;
          
          my @mergeends;
          
          my $mcend=$ends[0];
          my $mergecount=$count{$ends[0]};
          for my $i(1..@ends-1){
             if($ends[$i]-$mcend<=5){
                 if($count{$ends[$i]}<$count{$mcend}){
                       $mergecount+=$count{$ends[$i]};
                 }else{
                       $mcend=$ends[$i];
                       $mergecount+=$count{$ends[$i]};
                 }
             }else{
                 $count{$mcend}=$mergecount;
                 push(@mergeends,$mcend);
                 $mcend=$ends[$i];
                 $mergecount=$count{$ends[$i]};
             }
          }
                 $count{$mcend}=$mergecount;
                 push(@mergeends,$mcend);
          my @selectends;
          for my $iend(@mergeends){
              if($count{$iend}>=$ARGV[3]){
                   push(@selectends,$iend);
              }
          }
          @selectends=sort{ $a<=>$b }@selectends;
          $loci=~s/-flrvs//;          
          if(@selectends%2==0){
              if(@selectends==2){ 
                  print OUT "$loci-flrvs\t$selectends[0]\t$selectends[1]\t+\t$loci\n";
              }else{
                  for my $i(0..@selectends/2-1){
                      print OUT "$loci-flrvs\t$selectends[2*$i]\t$selectends[2*$i+1]\t+\t$loci$alpha[$i]\n";
                  }
              }
          }else{
               if(length($loci)>0){
               print "Warning: confusing results for $loci  in $ARGV[2]!\n";
               print "ends\t",join("\t",@selectends),"\n";
               }
          }
          
          $loci=$arr[0];
          undef @ends;
          if($arr[8]==1 ){
               push(@ends,$arr[6]);
               $count{$arr[6]}=1;
          }
          if($arr[9]==$len{$arr[1]} ){
               push(@ends,$arr[7]);
               $count{$arr[7]}=1;
          }
      }
}
          @ends=sort{$a<=>$b}@ends;
          
          my @mergeends;
          
          my $mcend=$ends[0];
          my $mergecount=$count{$ends[0]};
          for my $i(1..@ends-1){
             if($ends[$i]-$mcend<=5){
                 if($count{$ends[$i]}<$count{$mcend}){
                       $mergecount+=$count{$ends[$i]};
                 }else{
                       $mcend=$ends[$i];
                       $mergecount+=$count{$ends[$i]};
                 }
             }else{
                 $count{$mcend}=$mergecount;
                 push(@mergeends,$mcend);
                 $mcend=$ends[$i];
                 $mergecount=$count{$ends[$i]};
             }
          }
                 $count{$mcend}=$mergecount;
                 push(@mergeends,$mcend);
          my @selectends;
          for my $iend(@mergeends){
              if($count{$iend}>=$ARGV[3]){
                   push(@selectends,$iend);
              }
          }
          @selectends=sort{ $a<=>$b }@selectends;
          $loci=~s/-flrvs//;
          if(@selectends%2==0){
              if(@selectends==2){ 
                  print OUT "$loci-flrvs\t$selectends[0]\t$selectends[1]\t+\t$loci\n";
              }else{
                  for my $i(0..@selectends/2-1){
                      print OUT "$loci-flrvs\t$selectends[2*$i]\t$selectends[2*$i+1]\t+\t$loci$alpha[$i]\n";
                  }
              }
          }else{
               if(length($loci)>0){
               print "Warning: confusing results for $loci  in $ARGV[2]!\n";
               print "ends\t",join("\t",@selectends),"\n";
               }
          }


close IN;
close OUT;
