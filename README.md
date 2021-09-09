
# Core for "De novo circular RNA represents a novel class of HCMV transcripts"

Human cytomegalovirus (HCMV) infects 40-100% of the human population globally and may be a life-threatening pathogen in immunocompromised individuals. A number of HCMV non-coding RNAs are involved in regulating viral gene expression and the virus life cycle. Circular RNAs (circRNAs) are found in most eukaryotes and play critical roles in several human diseases including cancer, diabetes, neurological disorders and viral pathogenesis. Here for the first time, we report circRNAs as a new class of HCMV transcripts. We bioinformatically predicted 704 unique potential circRNAs for HCMV TB40E strain and 230 circRNAs for HCMV HAN strain, and experimentally identified 324 circRNAs for HCMV Towne strain from 426 clones. We systematically compared the features of HCMV circRNAs with those of host, other DNA virus and RNA virus circRNAs, including consensus sequence in breakpoint, strand-preference, length distribution, and isoform numbers and demonstrated that the unique characteristics of viral circRNAs were related to their genome type. Competitive endogenous RNA network analysis suggested that HCMV circRNAs play important roles in viral DNA synthesis and immune response via circRNA-miRNA-mRNA networks. We further annotated the potential function of all bioinformatically and experimentally identified circRNAs in HCMV growth properties and virus mediated manipulation of the host cell. Our findings highlight circRNAs as an important new class transcriptome component of HCMV that may contribute to viral replication and pathogenesis.  

![image](https://github.com/ShaominYang/Circular-RNAs-represent-a-novel-class-of-HCMV-transcripts/blob/main/HCMV-2021-5-3-2.jpg)
### <p align="center"> The structure diagram of HCMV (CINEMA 4D R20 for Windows 64 bit) </p>


# Requirements

### Computer

Two Intel W-3175X CPUs with 128 GB memory running Ubuntu system (version 18.04)

### Tools
- 1.Bash (Ubuntu, version 18.04)
- 2.Perl [https://www.perl.org](https://www.perl.org/)
- 3.Java [https://javadl.oracle.com](https://javadl.oracle.com/)
- 3.BWA [ttp://bio-bwa.sourceforge.net](http://bio-bwa.sourceforge.net/)
- 4.CIRI2 [https://sourceforge.net/projects/ciri/files/CIRI2/](https://sourceforge.net/projects/ciri/files/CIRI2/) 
- 5.SAMtools [http://www.htslib.org/](http://www.htslib.org/)
- 6.sratoolkit [https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software/)
- 7.aspera [https://www.ibm.com/products/aspera](https://www.ibm.com/products/aspera/)
- 8.ViReMa [https://sourceforge.net/projects/virema/](https://sourceforge.net/projects/virema/)

# Data preparation

###  Data collection
- 1.HCMV HAN strain infected HELF cells at 72 hours post-infection (PRJNA644588) [https://www.ebi.ac.uk/ena/browser/view/PRJNA577553](https://www.ebi.ac.uk/ena/browser/view/PRJNA577553)
- 2.HCMV TB40/E strain infected primary fibroblasts (HFF), endothelial cells (EC), and neural progenitors (NPCs) derived from embryonic stem cells at 48 hpi and 96 hpi (PRJNA299678)
- 3.KSHV BCBL1 strain infected B-cell lymphoma cells with RNase R-treatment (PRJNA483204)
- 4.EBV Akata strain infected B-cell lymphoma cells with RNase R-treatment (PRJNA479852)
###  Download RNA_seq data from ebi 
```Shell
for i in vol1/fastq/SRR102/087/SRR10277187/SRR10277187_1.fastq.gz vol1/fastq/SRR102/087/SRR10277187/SRR10277187_2.fastq.gz
do
echo $i
ascp -QT -l 300m -P33001 \
-i ~/miniconda2/envs/rna/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/${i} .
done
```
### multi-threading to speed up the extraction of fastq from SRA-accessions

```Shell
for i in your_SRR*_id
do
fasterq-dump ${i} -3 -e 16 -p
echo ${i}
done
rm SRR*
```

### Generating genome contained human and HCMV(HAN)

```Shell
cat hg19.fa KJ426589.1.fasta > hg19_HAN.fa
cat hg19.gtf KJ426589.1.gtf > hg19_HAN.gtf
```

### Build  BWA index

```Shell
bwa index hg19_HAN.fa
```

# running CIRI2 and circ-full pipeline

![image](https://github.com/ShaominYang/Circular-RNAs-represent-a-novel-class-of-HCMV-transcripts/blob/main/Illustration.jpg)
### <p align="center"> Illustration of CIRI2-based identification of circRNAs </p>

  
```Shell



for i in SRR10277187
do
echo $i
mkdir ${i}_output
bwa mem -t 52 hg19_HAN.fa ${i}_1.fastq.gz ${i}_2.fastq.gz >${i}_output/${i}.sam
perl ~/CIRI2/CIRI_v2.0.6/CIRI2.pl -I ${i}_output/${i}.sam -O ${i}_output/${i}.ciri -F hg19_HAN.fa -A hg19_HAN.gtf -T 24
## Reconstructed HCMV circRNAs circ-full
perl ~/CIRI2/CIRI_AS/CIRI_AS_v1.2.pl -S ${i}_output/${i}.sam -C ${i}_output/${i}.ciri -F hg19_HAN.fa -A hg19_HAN.gtf -O ${i}_output/${i} -D yes
java -jar ~/CIRI2/CIRI-full_v2.0/CIRI-full.jar RO1 -1 ${i}_1.fastq.gz -2 ${i}_2.fastq.gz -o ${i}_output/${i}
bwa mem -t 52 hg19_HAN.fa ${i}_output/${i}_ro1.fq > ${i}_output/${i}_ro1.sam
java -jar ~/CIRI2/CIRI-full_v2.0/CIRI-full.jar RO2 -r hg19_HAN.fa -s ${i}_output/${i}_ro1.sam -l 300 -o ${i}_output/${i}RO2
java -jar ~/CIRI2/CIRI-full_v2.0/CIRI-full.jar Merge -c ${i}_output/${i}.ciri -as ${i}_output/${i}_jav.list -ro ${i}_output/${i}RO2_ro2_info.list -a hg19_HAN.gtf -r hg19_HAN.fa -o ${i}_output/${i}
unset DISPLAY
java -jar ~/CIRI2/CIRI_vis/CIRI-vis_v1.4.jar -i ${i}_output/${i}_merge_circRNA_detail.anno -l ${i}_output/${i}_library_length.list -r hg19_HAN.fa -min 1
echo ${i}_CircRNA_full_finished
done
  
```

# running ViReMa pipeline

```Shell

bowtie-build KJ426589.1.fasta KJ426589.1.fasta

for i in SRR10277187
do
trim_galore --phred33 --gzip --trim-n -o trim_galore_out_dir --paired ${i}_1.fastq.gz ${i}_2.fastq.gz
echo ${i}_trim_galore
cat trim_galore_out_dir/${i}_1_val_1.fq.gz trim_galore_out_dir/${i}_2_val_2.fq.gz > ${i}.fq.gz
python2 ViReMa.py VZV_HSV1.fasta trim_galore_out_dir/${i}.fq.gz ${i}.sam --Output_Dir ${i}_virema/ --Output_Tag ${i} -BED --p 50 --MicroInDel_Length 5 --Defuzz 0 -FuzzEntry
done

```

# running statistics pipeline

### mapping statistics

```Shell
for i in SRR10277187
do
echo $i
samtools view -bS ${i}_output/${i}.sam > ${i}_output/${i}.bam
rm ${i}_output/${i}.sam
samtools sort ${i}_output/${i}.bam -o ${i}_output/${i}_sorted.bam -@ 42
qualimap bamqc -bam ${i}_output/${i}_sorted.bam -oc count.matrix -outdir ${i}_output/${i}_bamqc -outformat PDF:HTML --java-mem-size=50G
rm ${i}_output/${i}.bam
done
```
### genome coverage statistics

```Shell
for i in SRR10277187
do
samtools index -b ${i}_output/${i}_sort.bam ${i}_output/${i}_sort.bam.bai
samtools idxstats ${i}_output/${i}_sort.bam
samtools depth -a -m 0 ${i}_output/${i}_sort.bam > ${i}_output/${i}_coverage.txt
done
```


# Citations


>1.  Yang S, Zhou H, Cruz-Cosme R, Liu M, Xu J, Niu X, Li Y, Xiao L, Wang Q, Zhu H, Tang Q. Circular RNA profiling reveals abundant and diverse circRNAs of SARS-CoV-2, SARS-CoV and MERS-CoV origin. bioRxiv [Preprint]. 2020 Dec 8:2020.12.07.415422. doi: 10.1101/2020.12.07.415422. PMID: 33330860; PMCID: PMC7743059.
>2.  Yuan Gao†, Jinfeng Wang† and Fangqing Zhao*. CIRI: an efficient and unbiased algorithm for de novo circular RNA identification. Genome Biology (2015) 16:4.
>3.  Yuan Gao, Jinyang Zhang and Fangqing Zhao*. Circular RNA identification based on multiple seed matching. Briefings in Bioinformatics (2017) DOI: 10.1093/bib/bbx014.
>4.  Routh A, Johnson JE. Discovery of functional genomic motifs in viruses with ViReMa-a Virus Recombination Mapper-for analysis of next-generation sequencing data. Nucleic Acids Res. 2014 Jan;42(2):e11. doi: 10.1093/nar/gkt916. Epub 2013 Oct 16. PMID: 24137010; PMCID: PMC3902915.

