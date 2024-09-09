![image](https://github.com/user-attachments/assets/76787cea-94a2-41c4-9106-07206aec3926)
# A structural variant filtering and prioritization tool for long-read sequencing data

## INTRODUCTION

🚧 ** **_needLR_3.2  is currently a beta version and actively under construction_** ** 🚧

*Please contact jgust1@uw.edu with issues or suggestions.*

needLR is a command line tool that uses Jasmine merging to compare a query structural variant (SV) vcf to our collection of 1000 Genomes Project (1KGP) samples sequenced by Oxford Nanopore Technologies long-read sequencing (ONT LRS). The output is a .txt file with detailed annotations about the genomic context, OMIM phenotype association, and ancestry-specific allele frequencies of each of the SVs in the query vcf. 

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
2. Run Jasmine merge on 150 1KGP samples and 1 query sample
3. Isolate merged SVs seen in query sample (common and unique)
4. Assign ancestry aware allele frequencies to each query SV based on 1KGP sample input
5. Annotate query sample SVs with genomic context, OMIM phenotype association, and Hardy-Weinberg equilibrium check

## INSTALLATION AND SET UP

#### needLR requires an environment with the following dependencies:

[jasmine v1.1.5](https://github.com/mkirsche/Jasmine/) <sup>1</sup>  
[bedtools v2.31.1](https://github.com/arq5x/bedtools2/) <sup>2</sup>  
[bcftools v1.19](https://github.com/samtools/bcftools/) <sup>3</sup>  

<sup>1</sup> <sub>*Kirsche M, Prabhu G, Sherman R, et al. Jasmine and Iris: population-scale structural variant comparison and analysis. Nat Methods. 2023;20(3):408-417. doi:10.1038/s41592-022-01753-3*</sub>    
<sup>2</sup> <sub>*Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 2010;26(6):841-842. doi:10.1093/bioinformatics/btq033*</sub>  
<sup>3</sup> <sub>*Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools and BCFtools. Gigascience. 2021;10(2):giab008. doi:10.1093/gigascience/giab008*</sub>  


#### Download the needLR_local directory from AWS into a parent directory you want to run needLR from
```
wget https://s3.amazonaws.com/1000g-ont/needLR/needLR_local.tar.gz
```
>[!NOTE]
>The zip file is ~3.2G, the unzipped directory is ~9.5G


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
needLR_3.2.sh 
needLR_output/
EXAMPLE/
```

## RUN  NEEDLR

#### The needLR_3.2.sh command takes 3 required arguments:

| Flag | Description |
| :------------ |:-------------|
|-f| A .txt file that lists the full file path(s) to the query vcf(s) (The vcfs must be gzipped (*.vcf.gz) and have an index in the same directory as the vcf)|
|-g| A fasta file for a reference genome (we use the hg38 reference recommended [here](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use)) - This fasta is included in the backend_files download|
|-t| The number of threads to use in the Jasmine merging step|

#### To run needLR:
```
needLR_3.2.sh -f /file/path/to/list.txt -g /file/path/to/reference/genome.fa -t 20
```
As needLR is running, you will see the Jasmine output building - it will look like this:
<img width="878" alt="image" src="https://github.com/user-attachments/assets/7f2972ff-670c-4042-861d-172f704fc744">

Then it will look like this:
<img width="878" alt="image" src="https://github.com/user-attachments/assets/c6b03af3-bc68-4174-952d-c11243c46d02">

And eventually it will look like this:
<img width="878" alt="image" src="https://github.com/user-attachments/assets/40ae6ecb-087c-4d76-a1d5-42695b4ff063">  
(Sometimes this step can take quite a while)

## OUTPUT

Each query vcf will spawn an output subdirectory within the `/needLR_local/needLR_output` directory. 
The output files will be:

| Output | Description |
|:------------|:-------------|
|{SAMPLE_ID}_RESULTS.txt| Annotated query SVs (the main output)|
|preprocessed_{SAMPLE_ID}.vcf| Preprocessed query vcf that is used in the Jasmine merge|
|{SAMPLE_ID}_jasmine.vcf| Raw Jasmine output from the query sample + 1KGP sample merge|

>[!NOTE]
>needLR generates many temporary files when running, this can add up to ~2G

>[!NOTE]
>Jasmine will automatically make a seperate  `/needLR_local/output` directory (this can be ignored)

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
./needLR_3.2.sh -f EXAMPLE/input_list.txt -g /backend_files/hg38_fasta/hg38.no_alt.fa -t 20
```
and find the output here: `needLR_local/needLR_output/EXAMPLE_HG01511_sniffles2_needLR_3.2`  
The results should match the example output available [here](https://docs.google.com/spreadsheets/d/1N_OOIQIQJJLoDX_BcjiF02pvJQ7u4QU3/edit?gid=1980402699#gid=1980402699) in Excel format.

**Considerations/shortcomings:**

All of the allele frequencies are currently based on the number of autosomes in the 1KGP sample set, rendering allele frequencies (other than 0) for SVs on chrX and chrY inaccurate.

Some SVs are called and pass all filters but have a genotype of “./.” This can lead to a SV with an allele frequency of 0, but a population frequency of >0.


## PROOF OF CONCEPT

We ran needLR on sniffles2 vcfs from 8 individuals with pathogenic SVs that were unsolved by short-read sequencing.

needLR reduced the number of SVs per individual from an average of 22,583 total SVs to an average of 565 unique SVs:
<img width="768" alt="image" src="https://github.com/user-attachments/assets/ce6d86ec-0e90-4770-8499-5d84a11de208">

Of the unique SVs, an average of 144 intersect with hg38 annotated genes:
<img width="768" alt="image" src="https://github.com/user-attachments/assets/2c203333-e461-4a22-ac6a-249d2d08b9a2">

Of the unique, genic SVs, even fewer intersect with canonical exons and/or are in genes associated with OMIM phenotypes:
<img width="776" alt="image" src="https://github.com/user-attachments/assets/1649d44d-8e37-4355-8206-63d532099f56">

## ADITIONAL RESOURCES

[Here](https://s3.amazonaws.com/1000g-ont/needLR/UWONT150_needLR_3.2_RESULTS.txt.zip) is a "quick reference" annotated needLR output for just the first 150 1KGP samples (no query sample) - This has allele frequencies and genomic context info for each of the SVs seen in our first 150 samples.

## COMING SOON

#### More 1KGP samples
Our goal is to have 500 1KGP samples sequenced and analyzed by the begining of 2025

#### needLR_Trio  
For annotating inheritance and potential *de novo* SVs when trio data is available

#### chm13 compatibility 

#### Additional phenotype annotations (ie.HPO terms)
