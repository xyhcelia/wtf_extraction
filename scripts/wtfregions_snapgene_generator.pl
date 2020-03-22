#!/usr/bin/perl
use strict;
use Getopt::Long;
use File::Basename;
use vars qw( $inputfa $ann $wtfdb $feature $wtfgs $dir);

GetOptions(
   "h|help"=>\&USAGE,
   "i=s"=>\$inputfa,
   "ann-wtf_LTR_Tf"=>\$ann,
   "db=s"=>\$wtfdb,
   "f=s"=>\$feature,
   "gs=s"=>\$wtfgs,
) or USAGE();


sub USAGE{
my $usage=<<EOF;
Usage:perl $0  -i input.fa(fasta) --ann-wtf_LTR_Tf  -db wtfdb
"wtfdb" means database of known wtf sequences for blast; 
or
perl $0  -i input.fa(.fasta)  -f wtf_LTR_Tf.feature.txt -gs wtfgs.feature.txt
Format example for wtf_LTR_Tf.feature.txt:
JB939_wtf9-flrvs	2058	3832	+	JB939_wtf9	wtf
JB939_wtf9-flrvs	1751	1982	-	LTR	LTR
JB939_wtf9-flrvs	4389	4711	+	LTR	LTR
Format example for wtfgs.feature.txt(gene structure file):
CDS1    JB22_wtf1       -1      287     425     CDS1
CDS2    JB22_wtf1       -1      593     838     CDS2
CDS3    JB22_wtf1       -1      886     1116    CDS3
CDS4    JB22_wtf1       -1      1158    1211    CDS4
CDS5    JB22_wtf1       -1      1266    1339    CDS5
CDS1    JB22_wtf10      -1      273     411     CDS1
CDS2    JB22_wtf10      -1      581     826     CDS2
CDS3    JB22_wtf10      -1      876     1148    CDS3
EOF
print $usage;
exit;
}

$dir=dirname(__FILE__);

######optional:blast against  wtf and LTR #######:
my $id;
my @suffixlist = qw(.fa .fasta);
$id = basename($inputfa,@suffixlist);
if(defined $ann){
    system("blastn -query $inputfa   -db $wtfdb -outfmt 7  -out query-${id}.wtf.out"
);
    `cat query-${id}.wtf.out|awk '{if(\$1!~"#") print \$0;}'|awk '{if(\$9<\$10)print \$1"\t"\$7"\t"\$8"\t"\$2"\t.\t+";else print \$1"\t"\$7"\t"\$8"\t"\$2"\t.\t-";}'|sort -k1,1 -k2,2n >query-${id}.wtf.bed`;
    if(-z "query-${id}.wtf.bed"){
        `:>${id}.wtf_LTR_Tf.feature.txt`;
    }else{
        system("bedtools merge -i query-${id}.wtf.bed  -s -c 6 -o distinct >query-${id}.wtf.merge.bed");
        `cat query-${id}.wtf.merge.bed|awk '{print \$1"\t"\$2"\t"\$3"\t"\$4"\twtf\twtf";}' >${id}.wtf_LTR_Tf.feature.txt`;
     }
    system("blastn -query $inputfa   -db $dir/../data/LTR -task blastn -outfmt 7  -out query-${id}.LTR.out");
    `cat query-${id}.LTR.out|awk '{if(\$1!~"#" && \$4>=30) print \$0;}'|awk '{if(\$9<\$10)print \$1"\t"\$7"\t"\$8"\t"\$2"\t.\t+";else print \$1"\t"\$7"\t"\$8"\t"\$2"\t.\t-";}'|sort -k1,1 -k2,2n  >query-${id}.LTR.bed`;
    if(-z "query-${id}.LTR.bed"){
    }else{
         system("bedtools merge -i query-${id}.LTR.bed  -s -c 6 -o distinct >query-${id}.LTR.merge.bed");
         `cat query-${id}.LTR.merge.bed|awk '{print \$1"\t"\$2"\t"\$3"\t"\$4"\tLTR\tLTR";}' >>${id}.wtf_LTR_Tf.feature.txt`;
    }
    system("blastn -query $inputfa   -db $dir/../data/Tf_CDS -task blastn -outfmt 7  -out query-${id}.Tf.out");
    `cat query-${id}.Tf.out|awk '{if(\$1!~"#" && \$4>=50) print \$0;}'|awk '{if(\$9<\$10)print \$1"\t"\$7"\t"\$8"\t"\$2"\t.\t+";else print \$1"\t"\$7"\t"\$8"\t"\$2"\t.\t-";}'|sort -k1,1 -k2,2n  >query-${id}.Tf.bed`;
    if(-z "query-${id}.Tf.bed"){
    }else{
         system("bedtools merge -i query-${id}.Tf.bed  -s -c 6 -o distinct >query-${id}.Tf.merge.bed");
         `cat query-${id}.Tf.merge.bed|awk '{print \$1"\t"\$2"\t"\$3"\t"\$4"\tTf\tTf";}' >>${id}.wtf_LTR_Tf.feature.txt`;
    }
    `rm -f  query-${id}.wtf.out query-${id}.LTR.out query-${id}.Tf.out query-${id}.wtf.bed query-${id}.LTR.bed query-${id}.Tf.bed query-${id}.wtf.merge.bed query-${id}.LTR.merge.bed query-${id}.Tf.merge.bed`;    

}

######blast against Pombase CDS and Pombase gene####

system("blastn -query $inputfa   -db $dir/../data/PomBase_CDS -outfmt 7  -out query-${id}.CDS.out");

######CDS feature file########
my %cdsname;
open CDSTXT,"<$dir/../data/Schizosaccharomyces_pombe.ASM294v2.45.gff3.PomBase_CDS.nowtf.noTf.link.txt";
while(<CDSTXT>){
    chomp;
    my @arr=split;
    my $posi=$arr[0].":".$arr[1]."-".$arr[2]."_"."$arr[3]";
    $cdsname{$posi}=$arr[4];
}
close CDSTXT;

open CDSOUT,"<query-${id}.CDS.out";
open CDSFEATURE,">${id}.CDS.feature.txt";
while(<CDSOUT>){
    chomp;
    my @arr=split;
    if($arr[0]=~/#/){ next; }
    if($arr[8]<$arr[9]){
       print CDSFEATURE "$arr[0]\t$arr[6]\t$arr[7]\t+\t$cdsname{$arr[1]}\tCDS\n";    
    }else{
       print CDSFEATURE "$arr[0]\t$arr[6]\t$arr[7]\t-\t$cdsname{$arr[1]}\tCDS\n";
    }
}
close CDSOUT;
close CDSFEATURE;


`rm query-${id}.CDS.out `;

#######generate snapgene  file#######
my %color = (
    "CDS"  => "#993366",
    "LTR" => "#ffe4c4",
    "Tf" => "#99cc00",
    "wtf" => "#0000ff",
);


if( defined $ann){ 
     `cat ${id}.CDS.feature.txt $id.wtf_LTR_Tf.feature.txt >$id.feature.txt`;
     `rm ${id}.CDS.feature.txt $id.wtf_LTR_Tf.feature.txt`;
}else{
     `cat ${id}.CDS.feature.txt $feature >$id.feature.txt`;
     `rm ${id}.CDS.feature.txt `;
}


     
     $/=">";
     open FA,"<$inputfa";
     <FA>;
     $_=<FA>;
     my @arr=split;
     my $locus=$arr[0];
     my $seq=join("",@arr[1..@arr-1]);
     close FA;


open OUT, ">$id.gbk";

my $time=localtime();
my @time=split(/\s+/,$time);

print OUT "LOCUS       Exported                ",length($seq)," bp ds-DNA     linear   UNA $time[2]-",uc($time[1]),"-$time[4]\
DEFINITION  .\
ACCESSION   .\
VERSION     .\
KEYWORDS    .\
SOURCE      natural DNA sequence\
  ORGANISM  Schizosaccharomyces pombe
REFERENCE   1  (bases 1 to ",length($seq),")\
  AUTHORS   .\
  TITLE     Direct Submission\
  JOURNAL   Exported $time[1] $time[2], $time[4] from SnapGene Viewer 4.1.4\
            http://www.snapgene.com\
FEATURES             Location/Qualifiers\
     source          1..",length($seq),"\
                     /organism=\"Schizosaccharomyces pombe\"\
                     /mol_type=\"genomic DNA\"\
";     

$/="\n";
my $class;
my $col;
my $direction;
my @seg;
open FEATURE,"<$id.feature.txt";
while(<FEATURE>){
   chomp;
   my @arr=split;
   $class=$arr[5];
   $col=$color{$class};
   if($arr[5] !~/wtf/ || $arr[4] eq "wtf"){
       if($arr[3] eq "+"){
            print OUT "     $class             $arr[1]..$arr[2]\n";
            $direction="RIGHT";
       }else{   
            print OUT "     $class             complement($arr[1]..$arr[2])\n";
            $direction="LEFT";
       }
       print OUT "                     /label=$arr[4]\n";
       print OUT "                     /biotype=$arr[5]\n";
       print OUT "                     /note=\"color: $col; direction: $direction\"\n";
   }
   if($arr[5] eq "wtf" && $arr[4] ne "wtf"){
       `cat $wtfgs|grep $arr[4] >$arr[4].gene_structure.txt`;          
        if($arr[3] eq "+"){
             my $start=`cat $arr[4].gene_structure.txt|sort -k2,2 -k4,4n|head -n1 |awk '{print \$4;}'`;
             print OUT "     conserved_up             $arr[1]..",$arr[1]+$start-2,"\n";
             print OUT "                     /note=\"color: #ffcc00\"\n";
             my @seg=split(/\n/,`cat  $arr[4].gene_structure.txt|sort -k2,2 -k4,4n|awk '{print \$4+$arr[1]-1".."\$5+$arr[1]-1;}'`);
             print OUT "     CDS             join(",join(",",@seg),")\n";
             print OUT "                     /label=$arr[4]\n";
             print OUT "                     /note=\"This forward directional feature has ",$#seg+1," segments:\n";
             for my $i(1..@seg-1){
             print OUT "                       $i:$seg[$i-1]\n";
             }
             print OUT "                       $#seg:$seg[$#seg-1]\"\n";
        }else{
             my $start=`cat $arr[4].gene_structure.txt|sort -k2,2 -k4,4n|head -n1 |awk '{print \$4;}'`;
             print OUT "     conserved_up             ",$arr[2]-$start+2,"..$arr[2]\n";
             print OUT "                     /note=\"color: #ffcc00\"\n";
             my @seg=split(/\n/,`cat  $arr[4].gene_structure.txt|sort -k2,2 -k5,5nr|awk '{print $arr[2]-\$5+1".."$arr[2]-\$4+1;}'`);
             print OUT "     CDS             join(",join(",",@seg),")\n";
             print OUT "                     /label=$arr[4]\n";
             print OUT "                     /note=\"This reverse directional feature has ",$#seg+1," segments:\n";
             for my $i(1..@seg-1){
             print OUT "                       $i:$seg[$i-1]\n";
             }
             print OUT "                       $#seg:$seg[$#seg-1]\"\n";
        }
   `rm $arr[4].gene_structure.txt`;  
   }
}
close FEATURE;

` rm $id.feature.txt`;

     print OUT "ORIGIN\n";
     my $n=1;
     while($n<=length($seq)){
         printf OUT ( "%9s",$n);
         print  OUT " ",substr($seq,$n-1,10)," ",substr($seq,$n+10-1,10)," ",substr($seq,$n+20-1,10)," ",substr($seq,$n+30-1,10)," ",substr($seq,$n+40-1,10)," ",substr($seq,$n+50-1,10),"\n";
         $n+=60;
     }
     print OUT "//\n";

close OUT;
