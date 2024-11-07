![image](https://github.com/user-attachments/assets/76787cea-94a2-41c4-9106-07206aec3926)
# A structural variant filtering and prioritization tool for long-read sequencing data

## INTRODUCTION

🚧 ** **_needLR  is currently a beta version and actively under construction_** ** 🚧

*NOTE: needLR_3.2 has been updated to needLR_3.3 as of November 7th, 2024. The major change is that needLR_3.3 runs Truvari instead of Jasmine for SV merging. This streamlines the input (just 1 backend vcf) and reduces the amount of space needLR needs*

*needLR_trio_3.3 is also now available! Click [here](https://github.com/jgust1/needLR/blob/main/needLR_trio) for detailed instructions*

*Please contact jgust1@uw.edu with issues or suggestions.*

needLR is a command line tool that uses Truvari merging to compare a query structural variant (SV) vcf to our collection of 1000 Genomes Project (1KGP) samples sequenced by Oxford Nanopore Technologies long-read sequencing (ONT LRS). The output is a .txt file with detailed annotations about the genomic context, OMIM phenotype association, and ancestry-specific allele frequencies of each of the SVs in the query vcf. 

There are 3 key concepts that drive this project:  
* More than half of suspected Mendelian conditions remain molecularly unsolved after current clinical testing methods.  
* Basepair-for-basepair, there is more genetic diversity associated with structural variants (SVs) than SNVs and indels combined.  
* Traditional short read sequencing technologies are underpowered to resolve up to half of the SVs per genome.  

Thus, we need a comprehensive catalog of long-read sequencing (LRS)-based SV calls from healthy individuals in order to filter for rare potentially disease-causing SVs in clinical cases. 

The Miller Lab is actively using Oxford Nanopore Technologies (ONT) LRS to sequence samples from the 1000 Genomes Project (1KGP). Our recent [preprint](https://pubmed.ncbi.nlm.nih.gov/38496498/) describes the first 100 samples in the cohort. This version of needLR incorporates SV calls made by Sniffles2 v2.2 for 150 1KGP samples. 

#### Please cite our 2024 preprint:  
<sup>*Gustafson JA, Gibson SB, Damaraju N, et al. Nanopore sequencing of 1000 Genomes Project samples to build a comprehensive catalog of human genetic variation. Preprint. medRxiv. 2024;2024.03.05.24303792. Published 2024 Mar 7. doi:10.1101/2024.03.05.24303792*</sup>

## WORKFLOW

#### needLR performs the following steps on an input query vcf:

1. Preprocess query sample vcf to match 1KGP vcf format (SVs >=50bp, FILTER=PASS, full chromosomes)
2. Run Truvari on the query vcf and a pre-merged vcf of 1KGP samples
3. Isolate merged SVs seen in query sample (common and unique)
4. Assign ancestry aware allele frequencies to each query SV based on 1KGP sample input
5. Annotate query sample SVs with genomic context, OMIM phenotype association, and Hardy-Weinberg equilibrium check

## INSTALLATION AND SET UP

#### needLR requires an environment with the following dependencies:

[truvari v4.2.2](https://github.com/acenglish/truvari/wiki) <sup>1</sup>  
[bedtools v2.31.1](https://github.com/arq5x/bedtools2/) <sup>2</sup>  
[bcftools v1.19](https://github.com/samtools/bcftools/) <sup>3</sup>  

<sup>1</sup> <sub>*English, Adam C et al. “Truvari: refined structural variant comparison preserves allelic diversity.” Genome biology vol. 23,1 271. 27 Dec. 2022, doi:10.1186/s13059-022-02840-6*</sub>    
<sup>2</sup> <sub>*Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 2010;26(6):841-842. doi:10.1093/bioinformatics/btq033*</sub>  
<sup>3</sup> <sub>*Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools and BCFtools. Gigascience. 2021;10(2):giab008. doi:10.1093/gigascience/giab008*</sub>  


#### Download the needLR_local directory from AWS into a parent directory you want to run needLR from
```
wget https://s3.amazonaws.com/1000g-ont/needLR/needLR_local.tar.gz
```

#### Extract the directory 
```
tar -xvzf needLR_local.tar.gz
```
This will create a working directory called needLR_local. Everything needed to run needLR is inside.


#### Navigate into the working directory
```
cd needLR_local
```
This is what you should see:
```
backend_files/
needLR_3.3.sh 
needLR_output/
EXAMPLE/
```

## RUN  NEEDLR (standard) - go [here](https://github.com/jgust1/needLR/blob/main/needLR_trio) for needLR_trio

#### The needLR_3.3.sh command takes 2 required arguments:

| Flag | Description |
| :------------ |:-------------|
|-f| A .txt file that lists the full file path(s) to the query vcf(s) (The vcfs must be gzipped (*.vcf.gz) and have an index in the same directory as the vcf)|
|-g| A fasta file for a reference genome (we use the hg38 reference recommended [here](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use))|

#### To run needLR:
```
needLR_3.3.sh -f /file/path/to/list.txt -g /file/path/to/reference/genome.fa -t 20
```
## OUTPUT

Each query vcf will spawn an output subdirectory within the `/needLR_local/needLR_output` directory. 
The output files will be:

| Output | Description |
|:------------|:-------------|
|{SAMPLE_ID}_RESULTS.txt| Annotated query SVs (the main output)|
|preprocessed_{SAMPLE_ID}.vcf| Preprocessed query vcf that is used in the Truvari merge|
|{SAMPLE_ID}.truvari.vcf| Raw Truvari output from the query sample + 1KGP sample merge|

>[!NOTE]
>needLR generates many temporary files when running, this can add up to ~2G

## ANALYZING OUTPUT

The `{SAMPLE_ID}_RESULTS.txt` file can easily be opened in Excel. 
_Be sure that Excel is set up to delimit columns by tab (and only tab)_

These are the columns which can be sorted at will:

| Column Name           | Column Description                                                                 |
|:-----------------------|:-----------------------------------------------------------------------------------|
| Chr                   | Query SV chromosome                                                               |
| Start_Pos             | Query SV start coordinate (hg38)                                                   |
| End_Pos               | Query SV end coordinate (hg38)                                                     |
| SV_Length             | Query SV length                                                                    |
| SV_Type               | Query SV type                                                                      |
| 1KGP_support          | Samples in the 1000 Genomes Project samples that share the SV with the query sample (n=150) |
| Genotype              | Query SV genotype                                                                  |
| Alt_reads             | Number of reads in the query sample supporting the SV                              |
| Ref_reads             | Number of reference reads in the query sample at the SV locus                      |
| Total_reads           | Total number of reads at the SV locus (in the query sample)                        |
| Allele_Freq_ALL       | Allele frequency of SV in 1KGP samples (n=150)                                     |
| Genes                 | Genes that the SV intersects with (canonical hg38 coordinates, per gencode)        |
| OMIM                  | OMIM phenotypes associated with any gene the SV intersects (OMIM 8/2023)           |
| Exonic                | If the SV intersects with a canonical exon (hg38 coordinates, per gencode)         |
| Centromeric           | If the SV intersects with a centromere (UCSC hg38)                                 |
| Pericentromeric       | If the SV intersects with a pericentromeric region (+/-5Mb on either side of UCSC-defined centromere) |
| Telomeric             | If the SV intersects with a telomere (5Mb of either end of a chromosome)           |
| STR                   | If the SV intersects with a Short Tandem Repeat region (vamos original motifs, n=148) |
| VNTR                  | If the SV intersects with a Variable Number Tandem Repeat region (vamos original motifs, individuals=148) |
| Segdup                | If the SV intersects with a segmental duplication (Genome in a Bottle v3.3)        |
| Repeat                | If the SV intersects with a repeat region (UCSC hg38 repeat masker)                |
| Gap                   | If the SV intersects with an hg38 gap region (UCSC hg38 mapping and sequencing: gap) |
| HiConf                | If the SV is fully contained within a high confidence region (Genome in a Bottle T2TQ100-V1.0_stvar) |
| Pop_Count_AFR         | How many 1KGP AFR ancestry samples have the SV (n=51)                              |
| Pop_Freq_AFR          | Frequency (%) of 1KGP AFR ancestry samples with SV (n=51)                          |
| Pop_Count_AMR         | How many 1KGP AMR ancestry samples have the SV (n=18)                              |
| Pop_Freq_AMR          | Frequency (%) of 1KGP AMR ancestry samples with SV (n=18)                          |
| Pop_Count_EAS         | How many 1KGP EAS ancestry samples have the SV (n=25)                              |
| Pop_Freq_EAS          | Frequency (%) of 1KGP EAS ancestry samples with SV (n=25)                          |
| Pop_Count_EUR         | How many 1KGP EUR ancestry samples have the SV (n=24)                              |
| Pop_Freq_EUR          | Frequency (%) of 1KGP EUR ancestry samples with SV (n=24)                          |
| Pop_Count_SAS         | How many 1KGP SAS ancestry samples have the SV (n=32)                              |
| Pop_Freq_SAS          | Frequency (%) of 1KGP SAS ancestry samples with SV (n=32)                          |
| Pop_Count_ALL         | How many 1KGP SAS ancestry samples have the SV (n=150)                             |
| Pop_Freq_ALL          | Frequency (%) of 1KGP SAS ancestry samples with SV (n=150)                         |
| Allele_Count_AFR      | How many 1KGP AFR ancestry alleles have the SV (n=102)                             |
| Allele_Freq_AFR       | Frequency (%) of 1KGP AFR ancestry alleles with SV (n=102)                         |
| Allele_Count_AMR      | How many 1KGP AMR ancestry alleles have the SV (n=36)                              |
| Allele_Freq_AMR       | Frequency (%) of 1KGP AMR ancestry alleles with SV (n=36)                          |
| Allele_Count_EAS      | How many 1KGP EAS ancestry alleles have the SV (n=50)                              |
| Allele_Freq_EAS       | Frequency (%) of 1KGP EAS ancestry alleles with SV (n=50)                          |
| Allele_Count_EUR      | How many 1KGP EUR ancestry alleles have the SV (n=48)                              |
| Allele_Freq_EUR       | Frequency (%) of 1KGP EUR ancestry alleles with SV (n=48)                          |
| Allele_Count_SAS      | How many 1KGP SAS ancestry alleles have the SV (n=64)                              |
| Allele_Freq_SAS       | Frequency (%) of 1KGP SAS ancestry alleles with SV (n=64)                          |
| Allele_Count_ALL      | How many 1KGP SAS ancestry alleles have the SV (n=150)                             |
| Allele_Freq_ALL       | Frequency (%) of 1KGP SAS ancestry alleles with SV (n=150)                         |
| GT_homWT              | Number of 1KGP samples that are homozygous for the reference allele at the locus (n=150) |
| GT_het                | Number of 1KGP samples that are heterozygous for the SV allele at the locus (n=150) |
| GT_homVAR             | Number of 1KGP samples that are homozygous for the SV allele at the locus (n=150)  |
| HWE-p                 | Reference allele frequency in 1KGP samples (n=150)                                 |
| HWE-q                 | SV allele frequency in 1KGP samples (n=150)                                        |
| HWE                   | Is the SV in Hardy-Weinberg equilibrium in the 1KGP samples (n=150)                |

Column headers are also defined in `/needLR_local/needLR_output/column_key.txt`  

## EXAMPLE

In the downloaded zip file there is a directory `EXAMPLE/` which contains an example vcf and an example `input_list.txt`.  
If eveything is sut up correctly, you should be able to run the following command from inside the `/needLR_local` directory
```
./needLR_3.3.sh -f EXAMPLE/input_list.txt -g /path/to/reference.fasta
```
and find the output here: `needLR_local/needLR_output/EXAMPLE_HG01511_sniffles2_needLR_3.3`  
The results should match the example output available [here](https://docs.google.com/spreadsheets/d/1-GKLSyWqRhMJeR7RFZOmtTysYUgt1MUh/edit?gid=311322986#gid=311322986) in Excel format.

**Considerations/shortcomings:**

All of the allele frequencies are currently based on the number of autosomes in the 1KGP sample set, rendering allele frequencies (other than 0) for SVs on chrX and chrY inaccurate.

## COMING SOON

#### More 1KGP samples
Our goal is to have 500 1KGP samples sequenced and analyzed by the begining of 2025

#### chm13 compatibility 

#### Additional phenotype annotations (ie.HPO terms)
