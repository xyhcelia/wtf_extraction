# wtf_extraction.sh
#### Project-specific pipeline for wtf sequence extraction. As wtf are repetitive sequences and simply extracting by blast can not distinguish sequences from different  loci,also can not deal with structure variations or other complex situations.

## Usage
```
sh wtf_extraction.sh assembly.fa knownwtf.fa n
'assembly.fa':
The fasta file of a de novo assembly,from which we want to extract wtf sequences;
Notice: 1).The extension can be .fasta  or .fa and a simple file name is recommended; 2).This file should be in the current directory;
'knownwtf.fa':
A fasta file containing sequences of known wtf genes,which is from conserved_up to stop codon;
Notice: The absolute path of this file should be provided;
'n':
An integer which give a cutoff to demarcate start points and  end points of wtf sequences.
Notice: A large value is recommended,however,it should not be larger than the number of genomes that known wtf sequences are from;
```

### Process schematics
![image](https://github.com/xyhcelia/Readme_images/blob/master/wtf_extraction/wtf_extraction_schematic1.png)
![image](https://github.com/xyhcelia/Readme_images/blob/master/wtf_extraction/wtf_extraction_schematic2.png)

### Outputs
![image](https://github.com/xyhcelia/Readme_images/blob/master/wtf_extraction/wtf_extraction_schematic3.png)

### Notice 1
The following scripts are also needed to successfully run the pipeline.Download them and put them into the "bin" folder.
https://github.com/xyhcelia/sequence_analysis_tookit/blob/master/03_extract_fasta_from_ref.pl
https://github.com/xyhcelia/sequence_analysis_tookit/blob/master/blast_1st_hits.pl
https://github.com/xyhcelia/sequence_analysis_tookit/blob/master/fasta_length.pl
https://github.com/xyhcelia/txt_toolkit/blob/master/add_feature.pl

### Notice 2
To deal with abnormal situations(Structure variation or confusing starts&ends positions), Genebank files have been generated to give detail information in the major pipeline.Howerver,one always expect as little manual work as possible.Another two scripts is provided in the "bin" folder.
#### 1)error_gbk_info.pl(for confusing starts&ends positions,execute the following commands in the same directory where "wtf_extraction.sh" is executed )
```Bash
find ./ -name "*wts*-flrvs.gbk" -exec ls  {} \; >error.list
for id in ` cat error.list`;do
     basename=`echo $id | gawk 'match($0, /gbk_files\/(.*).gbk/, a) {print a[1];}'`
     strain=`echo $id | gawk 'match($0, /gbk_files\/(.*)_wts/, a) {print a[1];}'`
     perl error_gbk_info.pl $id >${strain}_wtf_extraction/$basename.wtf.txt
     perl 03_extract_fasta_from_ref.pl ${strain}_wtf_extraction/${strain}_wtf_region.reverse.fa ${strain}_wtf_extraction/$basename.wtf.txt ${strain}_wtf_extraction/$basename.wtf.fa
done
```
#### 2)extrawtf_gbk_info.pl(for structure variation,also executed in the same directory as the major pipeline)
first execute the following commands:
```Bash
find ./ -name "*extra_wtf-*.gbk" -exec ls  {} \; >extra_wtf.gbk.list
echo -e "regionid\twtf_posi\tinterval\tleft\tright">extra_wtf.txt
for id in ` cat extra_wtf.gbk.list`;do
    perl extrawtf_gbk_info.pl $id >>extra_wtf.txt
done
```
* **A txt file will be produced:**

regionid  	|wtf_posi 	|interval	 |left  |right
------------|:---------:|:--------:|:-----:|:------:
JB873.flyepilon1_extra_wtf-1	|1..1545	|	|	|LTR,SPCC576.15c:wts22-L-1,SPCC576.14:wts22-L-2,SPCC576.13:wts22-L-3,SPCC576.12c:wts22-L-4,SPAC1783.08c:I:2203741-2204346,SPCC576.11:wts22-L-5
JB873.ravenpilon1_extra_wtf-2	|899..2647,3285..5043,5678..7013 |  |	LTR,LTR	 |SPCC1183.11:wts10-R-1,SPCC31H12.02c:wts10-R-2,SPCC31H12.03c:wts10-R-3
JB1206.canupilon_extra_wtf-1	|10001..11828,12472..13786	|	|SPCC1235.14:wts0403-L-2,wts5-L-8,SPCC548.07c:wts5-L-3,wts0403-R-4,SPCC548.06c:wts5-L-4,wts0403-R-3,SPCC1235.13:wts0403-L-5,SPCC1529.01:wts5-L-2,wts0403-R-5,SPCC794.01c:wts5-L-1,wts0403-R-6	|SPCC794.03:wts0403-R-7,wts5-R-1
JB1206.canupilon_extra_wtf-3	|10001..11831,12476..13791	|	 |SPCC1529.01:wts5-L-2,wts0403-R-5,SPCC794.01c:wts5-L-1,wts0403-R-6	SPCC794.03:wts0403-R-7,wts5-R-1,SPCC794.04c:wts0403-R-8,wts5-R-2,SPCC794.16:wts0403-R-9,wts5-R-3,LTR,Tf |

* **Then the locus id can be conveniently decided（extra_locus_assign.txt）：**
```
JB1206.canupilon_extra_wtf-3    wts5
JB760.canupilon1_extra_wtf-1    wts0403
JB760.canupilon2_extra_wtf-1    wts0403
JB760.ravenpilon1_extra_wtf-1   wts1920
JB760.ravenpilon2_extra_wtf-1   wts1920
JB760.wtdbgpilon1_extra_wtf-1   wts1920
JB760.wtdbgpilon2_extra_wtf-1   wts1920
JB873.canupilon1_extra_wtf-1    wts22
JB873.canupilon2_extra_wtf-1    wts22
JB873.flyepilon1_extra_wtf-1    wts22
```
* **Then also take usage of the "error_gbk_info.pl"**
```Bash
while read extra wts
do
    strain=` echo $extra|awk '{split($1,a,"_"); print a[1];}'`
    ln -s  $PWD/${strain}_wtf_extraction/gbk_files/$extra.gbk   ${strain}_wtf_extraction/gbk_files/${strain}_$wts-flrvs.gbk
    echo ">${strain}_$wts-flrvs" >./${strain}_wtf_extraction/${strain}_$wts-flrvs.fa
    cat  ${strain}_wtf_extraction/${strain}.extra_wtf.region.fa|awk 'BEGIN{RS=">";}{if($1~"'$extra'")print ">"$0;}'| grep -v '^$'|grep -v ">" >>./${strain}_wtf_extraction/${strain}_$wts-flrvs.fa
    perl error_gbk_info.pl  ./${strain}_wtf_extraction/gbk_files/${strain}_$wts-flrvs.gbk >./${strain}_wtf_extraction/${strain}_$wts-flrvs.wtf.txt
    perl 03_extract_fasta_from_ref.pl  ./${strain}_wtf_extraction/${strain}_$wts-flrvs.fa ./${strain}_wtf_extraction/${strain}_$wts-flrvs.wtf.txt ./${strain}_wtf_extraction/${strain}_$wts-flrvs.wtf.fa
    rm ./${strain}_wtf_extraction/${strain}_$wts-flrvs.fa ./${strain}_wtf_extraction/${strain}_$wts-flrvs.fa.fai ./${strain}_wtf_extraction/${strain}_$wts-flrvs.wtf.txt
done<extra_locus_assign.txt
```
