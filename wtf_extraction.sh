if [ $# -ne 3 ]
then
cat<<EOF
Usage:  sh $0 assembly.fa knownwtf.fa n
'assembly.fa':
The fasta file of a de novo assembly,from which we want to extract wtf sequences;
Notice: 1).The extension can be .fasta  or .fa and a simple file name is recommended; 2).This file should be in the current directory;
'knownwtf.fa':
A fasta file containing sequences of known wtf genes,which is from conserved_up to stop codon;
Notice: The absolute path of this file should be provided;
'n':
An integer which give a cutoff to demarcate start points and  end points of wtf sequences.
Notice: A large value is recommended,however,it should not be larger than the number of genomes that known wtf sequences are from;
EOF
     exit
fi

path=$(cd `dirname $0`; pwd)

id=`echo $1|awk '{split($1,a,".fa");print a[1];}'`
mkdir ${id}_wtf_extraction    
cd ${id}_wtf_extraction
#########use pombe protein_coding CDS to search wtf regions#######
makeblastdb -in ../$1  -out $id -dbtype nucl
blastn -query $path/data/wtf_closest_CDS.fa -db $id -outfmt 7  -out $id.wtfregions.out
perl $path/scripts/blast_1st_hits.pl $id.wtfregions.out $id.wtfregions.1st.out
perl $path/scripts/add_feature.pl -i $path/data/wtf_closest_CDS.txt -c1 1 -p $id.wtfregions.1st.out -c2 1 -c3 2,9,10 -o wtf_closest_CDS.${id}_contig.txt
perl $path/scripts/paste.pl wtf_closest_CDS.${id}_contig.txt ${id}_wtf_region.txt >${id}_unextracted.txt

rm   $id.n* $id.wtfregions.1st.out  wtf_closest_CDS.${id}_contig.txt 

#########wtf region fasta file###########
perl $path/scripts/03_extract_fasta_from_ref.pl ../$1 ${id}_wtf_region.txt ${id}_wtf_region.fa
perl $path/scripts/wtf_region_reverse.pl $path/data/wtf_gene.tandem_link.newname.bed  ${id}_wtf_region.fa  $id ${id}_wtf_region.reverse.txt
perl $path/scripts/03_extract_fasta_from_ref.pl ${id}_wtf_region.fa ${id}_wtf_region.reverse.txt ${id}_wtf_region.reverse.fa

rm  ${id}_wtf_region.reverse.txt
#########extract wtf sequences without structure variations or abnormal situations#######
makeblastdb -in $2  -out known.wtf -dbtype nucl
blastn -query ${id}_wtf_region.reverse.fa -db known.wtf -outfmt 7  -out ${id}_wtf_region.reverse.wtf.out
perl  $path/scripts/fasta_length.pl $2 >known.wtf.len.txt
perl  $path/scripts/wtfregion_wtf_ends.newname.pl  known.wtf.len.txt  ${id}_wtf_region.reverse.wtf.out ${id}_wtf_region.reverse.wtf.txt  $3 >${id}_wtf_region.reverse.wtf.error

rm ${id}_wtf_region.reverse.wtf.out known.wtf.len.txt

perl $path/scripts/03_extract_fasta_from_ref.pl  ${id}_wtf_region.reverse.fa ${id}_wtf_region.reverse.wtf.txt $id.wtf.routine.fa

##############generate gbk files for wtf loci having extraction error########
cat ${id}_wtf_region.reverse.wtf.error|grep "Warning"|awk '{print $5;}' >error.list
for loci  in ` cat error.list`;do
      cat ${id}_wtf_region.reverse.fa|awk 'BEGIN{RS=">";}{if($1~"'$loci'-flrvs")print ">"$0;}'|grep -v '^$' >$loci-flrvs.fa
      perl $path/scripts/wtfregions_snapgene_generator.pl -i $loci-flrvs.fa --ann-wtf_LTR_Tf -db known.wtf
      rm $loci-flrvs.fa
done

rm error.list


#########extract wtf sequences with structure variations or abnormal situations#######
#######search abnormal wtf regions#########
cat ${id}_wtf_region.txt |awk '{print $1"\t"$2-1"\t"$3;}' >${id}_wtf_region.bed
bedtools maskfasta -fi ../$1  -bed  ${id}_wtf_region.bed -fo  $id.normal_wtf_region_mask.fa
makeblastdb -in $id.normal_wtf_region_mask.fa  -out $id.normal_wtf_region_mask -dbtype nucl
blastn -query $2 -db $id.normal_wtf_region_mask -outfmt 7  -out $id.extra_wtf.out
extra=`cat $id.extra_wtf.out|awk '{if($1!~"#") print $0;}'|wc -l`

rm  ${id}_wtf_region.bed  $id.normal_wtf_region_mask.n*
if [ $extra -eq 0 ]
then
    rm $id.extra_wtf.out
    rm $id.normal_wtf_region_mask.fa
    rm known.wtf.*
    rm *.fai
    cd ../
    exit
fi

cat $id.extra_wtf.out|awk '{if($1!~"#") print $0;}'|awk '{if($9<$10) print $2"\t"$9-1"\t"$10"\textra_wtf\t.\t+"; else print $2"\t"$10-1"\t"$9"\t'$id'_extra_wtf\t.\t-";}' |sort -k1,1 -k2,2n |bedtools merge -d 2000 -s -c 6 -o distinct -i - |awk '{if($2-9999>=1)print $1"\t"$2-9999"\t"$3+10000"\t"$4"\t'$id'_extra_wtf-"NR; else print $1"\t1\t"$3+10000"\t"$4"\t'$id'_extra_wtf-"NR;}' >$id.extra_wtf.region.txt

perl $path/scripts/03_extract_fasta_from_ref.pl $id.normal_wtf_region_mask.fa $id.extra_wtf.region.txt $id.extra_wtf.region.fa

rm $id.normal_wtf_region_mask.fa  $id.extra_wtf.out
##########generate gbk files###########
cat  $id.extra_wtf.region.fa|grep ">"|sed 's/>//'>$id.extra_wtf.region.list
for seqname  in ` cat $id.extra_wtf.region.list`;do
      cat $id.extra_wtf.region.fa|awk 'BEGIN{RS=">";}{if($1~"'$seqname'")print ">"$0;}'|grep -v '^$' >$seqname.fa
      perl $path/scripts/wtfregions_snapgene_generator.pl -i $seqname.fa --ann-wtf_LTR_Tf -db known.wtf
      rm $seqname.fa
done
rm $id.extra_wtf.region.list

mkdir gbk_files
mv *.gbk gbk_files
rm known.wtf.*
rm *.fai
cd ../
