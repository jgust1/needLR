### needLR_trio is deigned to annotate a query proband SV with inheritance information when parental vcfs are available. 

>[!NOTE]
>Unlike needLR (standard), needLR_trio can only run one trio at a time.

>[!NOTE]
>SV caller limitations, filtering of low-quality SVs, and merging errors  lead to a vast overestimation of de novo SVs. We reccommend using needLR_trio results as a starting point and always visualizing confirming candidate de novo SVs in a data viewer like IGV.

### RUN  NEEDLR_TRIO

#### The needLR_trio_3.3.sh command takes 4 required arguments:

| Flag | Description |
| :------------ |:-------------|
|-p| The full file path to the proband vcf (The vcfs must be gzipped (*.vcf.gz) and have an index in the same directory as the vcf)|
|-m| The full file path to the maternal vcf (The vcfs must be gzipped (*.vcf.gz) and have an index in the same directory as the vcf)|
|-f| The full file path to the paternal vcf (The vcfs must be gzipped (*.vcf.gz) and have an index in the same directory as the vcf)|
|-g| A fasta file for a reference genome (we use the hg38 reference recommended [here](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use))|

#### To run needLR_trio:
```
needLR_3.3.sh -f /path/to/proband.vcf -m /path/to/maternal.vcf -f /path/to/paternal.vcf -g /file/path/to/reference/genome.fa
```
## OUTPUT

An output subdirectory named for the proband will be initiated within the `/needLR_local/needLR_output` directory. 
The output files will be:

| Output | Description |
|:------------|:-------------|
|{SAMPLE_ID}_RESULTS.txt| Annotated proband SVs (the main output)|
|preprocessed_{SAMPLE_ID}.vcf| Preprocessed proband vcf that is used in the Truvari merge|
|{SAMPLE_ID}_truvari.vcf| Raw Truvari output from the proband sample + 1KGP sample merge|

>[!NOTE]
>needLR generates many temporary files when running, this can add up to ~2G

## ANALYZING OUTPUT

The `{SAMPLE_ID}_RESULTS.txt` file can easily be opened in Excel. 
_Be sure that Excel is set up to delimit columns by tab (and only tab)_

These are the columns which can be sorted at will:

| Column Name           | Column Description                                                                 |
|:-----------------------|:-----------------------------------------------------------------------------------|
| Chr                   | Proband SV chromosome                                                               |
| Start_Pos             | Proband SV start coordinate (hg38)                                                   |
| End_Pos               | Proband SV end coordinate (hg38)                                                     |
| SV_Length             | Proband SV length                                                                    |
| SV_Type               | Proband SV type                                                                      |
| 1KGP_support          | Samples in the 1000 Genomes Project samples that share the SV with the proband sample (n=150) |
| Genotype              | Proband SV genotype                                                                  |
| Alt_reads             | Number of reads in the proband sample supporting the SV                              |
| Ref_reads             | Number of reference reads in the proband sample at the SV locus                      |
| Maternal_genotye	|	Maternal genotype for proband SV|
| Paternal_genotye	|	Paternal genotype for proband SV|
| Inheritance	|	If the SV appears to be inherited from a parent (inherited) or de novo|
| Total_reads           | Total number of reads at the SV locus (in the proband sample)                        |
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
