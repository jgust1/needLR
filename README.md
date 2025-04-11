![image](https://github.com/user-attachments/assets/76787cea-94a2-41c4-9106-07206aec3926)
# A structural variant filtering and prioritization tool for long-read sequencing data

## INTRODUCTION

🚧 ** **_needLR  is currently a beta version and actively under construction_** ** 🚧

needLR_3.4 has replaced needLR_3.3 as of April 10th, 2024. Major changes include:
* More control samples! (450 samples from the 1KGP ONT Consortium)
* ALT and REF allele sequences in output
* More functionality (needLR_basic, needLR_duo, needLR_trio, needLR_custom_controls, needLR_annotate_multisample)
* VCF and txt file output formats

Access the depricated needLR_3.3 README [here](https://github.com/jgust1/needLR/blob/main/needLR_3.3_README_depricated_20250409.md)

*Please contact jgust1@uw.edu with issues or suggestions.*

needLR is a command line tool that uses Truvari merging to compare query structural variant (SV) vcfs to our collection of 1000 Genomes Project (1KGP) samples sequenced by Oxford Nanopore Technologies long-read sequencing (ONT LRS). The output includes .vcf and .txt files with detailed annotations about the genomic context, OMIM phenotype association, and ancestry-specific allele frequencies of each of the SVs in the query vcf. 

There are 3 key concepts that drive this project:  
* More than half of suspected Mendelian conditions remain molecularly unsolved after current clinical testing methods.  
* Basepair-for-basepair, there is more genetic diversity associated with structural variants (SVs) than SNVs and indels combined.  
* Traditional short read sequencing technologies are underpowered to resolve up to half of the SVs per genome.  

Thus, we need a comprehensive catalog of long-read sequencing (LRS)-based SV calls from healthy individuals in order to filter for rare potentially disease-causing SVs in clinical cases. 

The Miller Lab is actively using Oxford Nanopore Technologies (ONT) LRS to sequence samples from the 1000 Genomes Project (1KGP). Our recent [preprint](https://pubmed.ncbi.nlm.nih.gov/38496498/) describes the first 100 samples in the cohort. This version of needLR incorporates SV calls made by Sniffles2 v2.5.2 for 450 1KGP samples. 

#### Please cite our 2024 manuscript:  
<sup>*Gustafson, J. A., Gibson, S. B., Damaraju, N., Zalusky, M. P. G., Hoekzema, K., Twesigomwe, D., Yang, L., Snead, A. A., Richmond, P. A., De Coster, W., Olson, N. D., Guarracino, A., Li, Q., Miller, A. L., Goffena, J., Anderson, Z. B., Storz, S. H. R., Ward, S. A., Sinha, M., Gonzaga-Jauregui, C., … Miller, D. E. (2024). High-coverage nanopore sequencing of samples from the 1000 Genomes Project to build a comprehensive catalog of human genetic variation. Genome research, 34(11), 2061–2073. https://doi.org/10.1101/gr.279273.124*</sup>

## WORKFLOW AND FUNCTIONALITY

#### needLR performs the following steps on an input query vcf:

1. Preprocess query sample vcf to match 1KGP vcf format (SVs >=50bp, FILTER=PASS, full chromosomes)
2. Run Truvari on the query vcf and a pre-merged vcf of 1KGP samples
3. Isolate merged SVs seen in query sample (common and unique)
4. Assign ancestry aware allele frequencies to each query SV based on 1KGP sample input
5. Annotate query sample SVs with genomic context, OMIM phenotype association, and Hardy-Weinberg equilibrium check

needLR_3.4 offers multiple functionalities:
* [needLR_basic](#run-needlr_34_basic): Compares one query vcf to a pre-merged, multisample vcf of 450 1KGP samples and annotates the SVs in the query individual.
* [needLR_duo](#run-needlr_34_duo): Compares a query sample and parental sample to a pre-merged, multisample vcf of 450 1KGP samples and annotates the SVs in the query individual. This function uniquely annotates the SVs from the query vcf as being "inherited" or "uncertain" based on SVs from the parental vcf.
* [needLR_trio](#run-needlr_34_trio): Compares a query sample and two parental samples (maternal and paternal) to a pre-merged, multisample vcf of 450 1KGP samples and annotates the SVs in the query individual. This function uniquely annotates the SVs from the query vcf as being "inherited" or "de novo" based on SVs from the parental vcfs.
* [needLR_custom_controls](#run-needlr_34_custom_controls): Compares one query vcf to a pre-merged, multisample vcf of a cohort of samples defined by the user and annotates the SVs in the query individual.
* [needLR_annotate_multisample](#run-needlr_34_annotate_multisample): Annotates any multisample vcf (from bcftools merge or Trivari). Does not require a specific query sample

>[!NOTE]
>needLR is currently optimized for Sniffles2 SV calling and all backend annotation data is based on the GRCh38 reference genome

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
wget https://s3.amazonaws.com/1000g-ont/needLR/needLR_3.4_local.tar.gz
```

#### Extract the directory 
```
tar -xvzf needLR_3.4_local.tar.gz
```
This will create a working directory called needLR_3.4_local. Everything needed to run needLR is inside.


#### Navigate into the working directory
```
cd needLR_3.4_local
```
This is what you should see:
```
backend_files/
needLR_3.4_annotate_multisample.sh
needLR_3.4_basic.sh
needLR_3.4_custom_controls.sh
needLR_3.4_duo.sh
needLR_3.4_trio.sh
needLR_output/
EXAMPLE/
```

## RUN needLR_3.4_basic

#### needLR_3.4_basic takes 2 required arguments:

| Flag | Description |
| :------------ |:-------------|
|-g| A fasta file for a reference genome (we use the hg38 reference recommended [here](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use))|
|-l| A .txt file that lists the full file path(s) to the query vcf(s) (The vcfs must be gzipped (*.vcf.gz) and have an index in the same directory as the vcf)|

#### Example:
```
./needLR_3.4_basic -g /file/path/to/reference/genome.fa -l /file/path/to/list.txt
```

## RUN needLR_3.4_duo

#### needLR_3.4_duo takes 4 required arguments:

| Flag | Description |
| :------------ |:-------------|
|-g| A fasta file for a reference genome (we use the hg38 reference recommended [here](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use))|
|-b| The full file path to the query/proband vcf (The vcf must be gzipped (*.vcf.gz) and have an index in the same directory as the vcf) |
|-r| The full file path to the parental vcf (The vcf must be gzipped (*.vcf.gz) and have an index in the same directory as the vcf) |

#### Example:
```
./needLR_3.4_duo -g /file/path/to/reference/genome.fa -b /file/path/to/proband.vcf.gz -r /file/path/to/parent.vcf.gz
```

>[!NOTE]
>needLR_duo is imperfect in predicting inherited vs. _de novo_ SVs. The annotation is fully dependent on how well the SVs were merged. We recommend using needLR as a starting point and then manually inspecting inheritance in IGV.

## RUN needLR_3.4_trio

#### needLR_3.4_trio takes 5 required arguments:

| Flag | Description |
| :------------ |:-------------|
|-g| A fasta file for a reference genome (we use the hg38 reference recommended [here](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use))|
|-b| The full file path to the query/proband vcf (The vcf must be gzipped (*.vcf.gz) and have an index in the same directory as the vcf) |
|-m| The full file path to the maternal vcf (The vcf must be gzipped (*.vcf.gz) and have an index in the same directory as the vcf) |
|-p| The full file path to the paternal vcf (The vcf must be gzipped (*.vcf.gz) and have an index in the same directory as the vcf) |

#### Example:
```
./needLR_3.4_trio -g /file/path/to/reference/genome.fa -b /file/path/to/proband.vcf.gz -m /file/path/to/maternal.vcf.gz -p/file/path/to/paternal.vcf.gz
```

>[!NOTE]
>needLR_trio is imperfect in predicting inherited vs. _de novo_ SVs. The annotation is fully dependent on how well the SVs were merged. We recommend using needLR as a starting point and then manually inspecting inheritance in IGV.

## RUN needLR_3.4_custom_controls

#### needLR_3.4_custom_controls takes 5 required arguments:

| Flag | Description |
| :------------ |:-------------|
|-g| A fasta file for a reference genome (we use the hg38 reference recommended [here](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use))|
|-l| A .txt file that lists the full file path(s) to the query vcf(s) (The vcfs must be gzipped (*.vcf.gz) and have an index in the same directory as the vcf)|
|-t| A Truvari-merged, sorted, and tabix'd vcf of the the user's control samples (see below for instructions on generating this)  |
|-n| A .txt file with the names of the samples in the exact order they are in in the Truvari-merged vcf |
|-s| The total number of samples in the Truvari-merged vcf|

#### To prepare a custom cohort to be used in needLR:
#### 1) Prepare all of the individual vcfs to include SVs that are >=50bp, "PASS" the filtering caller's criteria, and are on full-length chromosomes 1-22, X, Y, and M

```
bcftools view -i '(INFO/SVTYPE="BND") || (INFO/SVTYPE="INS" || INFO/SVTYPE="DEL" || INFO/SVTYPE="DUP" || INFO/SVTYPE="INV") && (INFO/SVLEN > 49 || INFO/SVLEN < -49)' -f PASS -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM -o preprocessed.vcf original.vcf
```
#### 2) Merge the preprocessed vcfs using bcftools merge and tabix the output (where sample_path_list.txt is a list of all of the sample vcfs to merge in the order they should be in). This merge only merges exact variant matches and outputs a multisample vcf
```
bcftools merge -m none --force-samples -l sample_path_list.txt -Oz -o bcftools_merged.vcf.gz

tabix bcftools_merged.vcf.gz
```

#### 3) Use Truvari to further merge the bcftools-merged vcf. Below are the parameters used in needLR itself (these can be customized). Sort, gzip, and tabix the Truvari merged output

```
truvari collapse -i bcftools_merged.vcf.gz -o truvari_merged.vcf -c truvari_collapsed.vcf -f reference_gemome.fa -k common -r 2000 -p 0 -B 50 -P 0.2 -O 0.2 --chain -s 50 -S 10000000 --passonly

bcftools sort truvari_merged.vcf > truvari_merged_sorted.vcf

bgzip truvari_merged_sorted.vcf

tabix truvari_merged_sorted.vcf.gz
```

#### Example:
```
./needLR_3.4_custom_controls -g /file/path/to/reference/genome.fa -l /file/path/to/list.txt -t /file/path/to/truvari_merged_sorted.vcf.gz -n /file/path/to/sample_names.txt -s 100
```

## RUN needLR_3.4_annotate_multisample

#### needLR_3.4_annotate_multisample takes 4 required arguments:

| Flag | Description |
| :------------ |:-------------|
|-g| A fasta file for a reference genome (we use the hg38 reference recommended [here](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use))|
|-u| A multipsample vcf - must be gzipped (*.vcf.gz) and have an index in the same directory as the vcf |
|-n| A .txt file with the names of the samples in the exact order they are in in the multisample vcf |
|-s| The total number of samples in the multisample vcf|


#### Example:
```
./needLR_3.4_annotate_multisample -g /file/path/to/reference/genome.fa -u /file/path/to/multisample.vcf.gz -n /file/path/to/sample_names.txt -s 100
```

## OUTPUT

Each needLR instance will spawn an output subdirectory within the `/needLR_local/needLR_output` directory. 
The output files will be:

| Output | Description |
|:------------|:-------------|
|{SAMPLE_ID}_RESULTS.txt| Annotated query or cohort SVs (the main output)|
|{SAMPLE_ID}_RESULTS_unique.txt| Annotated SVs that are unique to the query sample (not seen in the control cohort) |
|{SAMPLE_ID}_RESULTS_0.01.txt| Annotated SVs with an AF <=0.01  (in realation to the control cohort)|
|{SAMPLE_ID}_RESULTS.vcf.gz| Annotated query SVs in vcf format |
|{SAMPLE_ID}_RESULTS.vcf.gz.tbi| Tabix'd vcf|


>[!NOTE]
>needLR generates many temporary files when running, this can add up to ~2G

## ANALYZING OUTPUT

The `{SAMPLE_ID}_RESULTS*.txt` files can easily be opened in Excel. 
_Be sure that Excel is set up to delimit columns by tab (and only tab)_

Below are the output columns. Some are specific to the needLR function used. "Query SV" is used to refer to the SV identified in the proband/query sample _or_, in the case of needLR_3.4_annotate_multisample, it is any SV. 
>[!NOTE]
>The SV start, end, length, REF, and ALT are the data from the most common SV of the merged SVs at that locus. Thus, this data in "query only" SVs (SVs seen in the query sample but not the control samples) will exactly match the query original vcf. If the query sample has an SV that is shared with control samples, it will be annotated with the data from the most common SV in that merge. 

| Column Name           | Column Description                                                                 |
|:-----------------------|:-----------------------------------------------------------------------------------|
| Chr                   | Query SV chromosome                                                               |
| Start_Pos             | Query SV start coordinate (hg38)                                                   |
| End_Pos               | Query SV end coordinate (hg38)                                                     |
| REF             | Reference allele associated with the query SV                                                                    |
| ALT             | Alternate allele associated with the query SV                                                                    |
| SV_Length             | Query SV length                                                                    |
| SV_Type               | Query SV type                                                                      |
| Query ID               | Unique ID of the query sample                                                                     |
| Sample_support          | Samples in the control dataset that share the SV with the query sample |
| Genotype              | Query SV genotype                                                                  |
| Alt_reads             | Number of reads in the query sample supporting the SV                              |
| Ref_reads             | Number of reference reads in the query sample at the SV locus                      |
| Total_reads           | Total number of reads at the SV locus (in the query sample)                        |
| Allele_Freq_ALL       | Allele frequency of SV in control dataset                                   |
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
| Pop_Count_AFR         | How many 1KGP AFR ancestry samples have the SV                             |
| Pop_Freq_AFR          | Frequency (%) of 1KGP AFR ancestry samples with SV                          |
| Pop_Count_AMR         | How many 1KGP AMR ancestry samples have the SV                              |
| Pop_Freq_AMR          | Frequency (%) of 1KGP AMR ancestry samples with SV                          |
| Pop_Count_EAS         | How many 1KGP EAS ancestry samples have the SV                              |
| Pop_Freq_EAS          | Frequency (%) of 1KGP EAS ancestry samples with SV                          |
| Pop_Count_EUR         | How many 1KGP EUR ancestry samples have the SV                              |
| Pop_Freq_EUR          | Frequency (%) of 1KGP EUR ancestry samples with SV                          |
| Pop_Count_SAS         | How many 1KGP SAS ancestry samples have the SV                              |
| Pop_Freq_SAS          | Frequency (%) of 1KGP SAS ancestry samples with SV                          |
| Pop_Count_ALL         | How many control samples have the SV                             |
| Pop_Freq_ALL          | Frequency (%) of control samples with SV                         |
| Allele_Count_AFR      | How many 1KGP AFR ancestry alleles have the SV                             |
| Allele_Freq_AFR       | Frequency (%) of 1KGP AFR ancestry alleles with SV                      |
| Allele_Count_AMR      | How many 1KGP AMR ancestry alleles have the SV                            |
| Allele_Freq_AMR       | Frequency (%) of 1KGP AMR ancestry alleles with SV                           |
| Allele_Count_EAS      | How many 1KGP EAS ancestry alleles have the SV                             |
| Allele_Freq_EAS       | Frequency (%) of 1KGP EAS ancestry alleles with SV                         |
| Allele_Count_EUR      | How many 1KGP EUR ancestry alleles have the SV                               |
| Allele_Freq_EUR       | Frequency (%) of 1KGP EUR ancestry alleles with SV                          |
| Allele_Count_SAS      | How many 1KGP SAS ancestry alleles have the SV                               |
| Allele_Freq_SAS       | Frequency (%) of 1KGP SAS ancestry alleles with SV                          |
| Allele_Count_ALL      | How many control samples have the SV                              |
| Allele_Freq_ALL       | Frequency (%) of control samples alleles with SV                          |
| GT_homWT              | Number of control samples that are homozygous for the reference allele at the locus  |
| GT_het                | Number of control that are heterozygous for the SV allele at the locus |
| GT_homVAR             | Number of control samples that are homozygous for the SV allele at the locus  |
| HWE                   | Is the SV in Hardy-Weinberg equilibrium in the 1KGP samples                |


## EXAMPLE

In the downloaded zip file there is a directory `EXAMPLE/` which contains an example query vcf, input list, and results.  
If eveything is sut up correctly, you should be able to run the following command from inside the `/needLR_local` directory
```
./needLR_3.4_basic.sh -g /file/path/to/reference/genome.fa -l EXAMPLE/EXAMPLE_list.txt
```
and find the output here: `needLR_local/needLR_output/EXAMPLE_HG02555_needLR_3.4_basic`  
The results should match the example output in `EXAMPLE/`

## NOTES ON TRUVARI PARAMETERS

By design, needLR uses extremely relaxed merging parameters, tending more toward over-merging than under-merging. We have shown that this is suitable for all of our test cases. However, if the user wants to adjust the Truvari parameters, they are defined at the very top of the script for easy hacking. 

## CONSIDERATIONS AND LIMITATIONS

All of the allele frequencies are currently based on the number of autosomes in the 1KGP sample set, rendering allele frequencies (other than 0) for SVs on chrX and chrY inaccurate.

BNDs and SVs >=10Mb are filtered out (huge variants are not efficient to annotate this way).

Currently, needLR is not effective in identifying candidate pathogenic tandem repeats due to the nature of the relaxed SV merging parameters, but a future version will incorporate a tandem repeat annotation tool.

## COMING SOON

* More control samples! We will continue to increase the number of 1KGP controls as they are sequenced/processed.

* Optimization with PacBio data and a diversity of SV callers (both alignment- and assembly-based)

* Additional annotations 

* Compatibility with the chm13 reference genome
