![image](https://github.com/user-attachments/assets/76787cea-94a2-41c4-9106-07206aec3926)
# A structural variant filtering and prioritization tool for long-read sequencing data

## INTRODUCTION

🚧 ** **_needLR  is actively under construction_** ** 🚧

needLR_v4.0 has replaced needLR_v3.5 as of April 3rd, 2026. Major changes include:
* needLR modes are now subcommands -- please see updated usage
* Custom control sets may be used with any subcommand
* Input can be provided as a single VCF, a list of files, or a merged multi sample VCF
* Analysis can be restricted to a smaller chromosomal region
* Annotations can be selected a la carte
* OMIM and GenCC annotations are now multicolumn for easier sorting
* needLR trio and duo preserve parental read support levels
* A reference fasta is no longer required -- currently needLR only works for hg38
* vcfs do not need to be bgzipped and indexed

Changes from 3.4 -> 3.5
* Increase to 500 control samples from the 1KGP-LRSC
* New annotations: pLI (Probability of Loss-of-function Intolerance) scores and [ORegAnno](https://www.bcgsc.ca/resources/software/oreganno) annotation
* Gene annotation buffer is increased from 1kbp to  5kbp
* UTRs and coding exons (CDS) are annotated seperately
* Vamos STR and VNTR are combined into one annotation

Access the deprecated needLR_v3.4 README [here](https://github.com/millerlaboratory/needLR/blob/main/docs/needLR_v3.4_README_depricated_20251213.md)

Access the deprecated needLR_v3.5 README [here](https://github.com/millerlaboratory/needLR/blob/main/docs/needLR_v3.4_README_deprecated_20260403.md)

#### Please cite our 2025 needLR preprint:  
<sup>*Gustafson JA, Lin J, Zalusky MPG, Eichler EE, Miller DE. needLR: Long-read structural variant annotation with population-scale frequency estimation. arXiv preprint arXiv:2512.08175. 2025 Dec 9.*</sup>

*Please contact jgust1@uw.edu with issues or suggestions.*

needLR is a command line tool that uses Truvari merging to compare query structural variant (SV) vcfs to our collection of 1000 Genomes Project (1KGP) samples sequenced by Oxford Nanopore Technologies long-read sequencing (ONT LRS). The output includes .vcf and .tsv files with detailed annotations about the genomic context, OMIM phenotype association, and ancestry-specific allele frequencies of each of the SVs in the query vcf. 

There are 3 key concepts that drive this project:  
* More than half of suspected Mendelian conditions remain molecularly unsolved after current clinical testing methods.  
* Basepair-for-basepair, there is more genetic diversity associated with structural variants (SVs) than SNVs and indels combined.  
* Traditional short read sequencing technologies are underpowered to resolve up to half of the SVs per genome.  

Thus, we need a comprehensive catalog of long-read sequencing (LRS)-based SV calls from healthy individuals in order to filter for rare potentially disease-causing SVs in clinical cases. 

The Miller Lab is actively using Oxford Nanopore Technologies (ONT) LRS to sequence samples from the 1000 Genomes Project (1KGP). Our 2024 [publication](https://pubmed.ncbi.nlm.nih.gov/39358015/) describes the first 100 samples in the cohort. 


This version of needLR incorporates SV calls made by Sniffles_v2.6.2 for 500 1KGP samples using the following parameters. (The tandem repeat file was sourced from PacBio PBSV annotation [here](https://github.com/PacificBiosciences/pbsv/tree/master/annotations)) - we recommend generating SVs for user query/input files using the exact same parameters for optimal merging efficiency.  

```
        sniffles_v2.6.2 \
            --input sample.bam \
            --reference hg38_reference.fa \
            --output-rnames \
            --vcf sample.vcf \
            --allow-overwrite \
            --tandem-repeats human_GRCh38_no_alt_analysis_set.trf.bed
```

## WORKFLOW AND FUNCTIONALITY

#### needLR performs the following steps on an input query vcf:

1. Preprocess query sample vcf to match 1KGP vcf format (SVs >=50bp, full chromosomes)
2. Run Truvari on the query vcf and a pre-merged vcf of 1KGP samples
3. Isolate merged SVs seen in query sample (common and unique)
4. Assign ancestry aware allele frequencies to each query SV based on 1KGP sample input
5. Annotate query sample SVs with genomic context, OMIM phenotype association, and Hardy-Weinberg equilibrium check

needLR_v4.0 has three subcommands:
* [annotate](#subcommand-annotate): Compares one or more query vcfs to a pre-merged, multisample vcf of 500 1KGP samples and annotates the SVs in the query individual. A custom control set may optionally be provided. If a multisample vcf is provided, SVs that are present in one more affected individuals in the cohort are annotated.
* [comparator](#subcommand-comparator): Compares a single query sample and one or two parental samples to a pre-merged, multisample vcf of 500 1KGP samples and annotates the SVs in the query individual. This function uniquely annotates the SVs from the query vcf as being "inherited", "maternal", "paternal", "de_novo", or "not_inherited" based on SVs from the parental vcf(s).
* [bed](#subcommand-bed): Annotates any sorted bed file with needLR annotations

[Follow these steps to make a custom cohort for use either as a query or a control.](https://github.com/millerlaboratory/needLR/blob/main/docs/custom_cohort.md)

>[!NOTE]
>needLR is currently optimized for sniffles_v2.6.2 SV calling, truvari_v4.2.2 merging, and all backend annotation data is based on the GRCh38 reference genome

## INSTALLATION AND SET UP

Please install needLR using conda or Docker/podman/apptainer. 

>[!NOTE]
>needLR relies heavily on bash scripts and was developed in bash version 4.4-release. It was also tested using bash version 5.0-relase. There is a known incompatibility with bash versions earlier than 4.4-release. If using an earlier version of bash, please try the docker/podman/apptainer container installation method.

### Conda

Build an environment for needLR to run in like so:

```
conda create -n needLR-4.0 -c bioconda -c conda-forge needlr=4.0
conda activate needLR-4.0
```


Alternatively, you can make a custom conda installation following these steps:


1. Build a conda environment using the `.yaml` file: `envs/needLR-4.0.yaml`
2. Clone this repository.
3. Copy or make a sm link of `needLR` and `src/` in this repository to `${CONDA_PREFIX}/bin`
4. Download the backend files required to run needLR from AWS: 

```
wget https://s3.amazonaws.com/1000g-ont/needLR/needLR-v4.0-backend-files.tar.gz
tar -xvzf needLR_v4.0_backend_files.tar.gz
```

If using a custom installation, you **must** use flag `-B` with your path to the `backend_files` folder downloaded in step 4 above.


### Docker

Pull needLR from dockerhub:

```
docker pull miragale/needlr:latest
```

When running with docker, your command must be modified to mount input and output locations, and needLR **must** be given an output location. Please follow this pattern:

```
INPUTDIR=/full/path/to/input/loc
OUTPUTDIR=/full/path/to/output/loc

docker run -v ${INPUTDIR}:/mnt/inputs \
  -v ${OUTPUTDIR}:/mnt/outputs \
  miragale/needlr:latest needLR {subcommand} -O /mnt/outputs/needLR_output <more options> /mnt/inputs/{your.vcf.gz}
```

### Dependencies included in package:

[truvari v4.2.2](https://github.com/acenglish/truvari/wiki) <sup>1</sup>  
[bedtools v2.31.1](https://github.com/arq5x/bedtools2/) <sup>2</sup>  
[bcftools v1.23.1](https://github.com/samtools/bcftools/) <sup>3</sup>  

<sup>1</sup> <sub>*English, Adam C et al. “Truvari: refined structural variant comparison preserves allelic diversity.” Genome biology vol. 23,1 271. 27 Dec. 2022, doi:10.1186/s13059-022-02840-6*</sub>    
<sup>2</sup> <sub>*Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 2010;26(6):841-842. doi:10.1093/bioinformatics/btq033*</sub>  
<sup>3</sup> <sub>*Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools and BCFtools. Gigascience. 2021;10(2):giab008. doi:10.1093/gigascience/giab008*</sub>  


## USAGE

>[!NOTE] 
>Options and commands have been updated since needLR v3.5

`needLR {subcommand} <options> {input.vcf.gz} {input.2.vcf.gz} {more.input.vcf.gz}...`

  Global options:
  ``` -B                     : [ path to folder containing backend files for non-conda custom installation ]
   -T                     : [ additional CPU threads to pass to bcftools ]
   -R                     : [ restrict analysis to region (e.g. chr1:23456-34567) ]
   -B                     : [ the full path to needLR's back end file folder, only used in custom installations ]
   -H                     : [ help ]

   ```

Alternatively, you can run needLR through our nextflow wrapper, which can speed up running many samples at once when a server or cluster with many cpus is available. [See the documentation here.](https://github.com/millerlaboratory/needLR/blob/nextflow/nextflow/nextflow_readme.md)




### Subcommand: annotate

`annotate` requires VCF input, either as positional arguments, with option `-Q` to provide a multisample VCF, or with flag `-L` to provide a list of input files. 

Additional options:

| Option | Description |
| :------------ |:-------------|
|-T| additional CPU threads to pass to bcftools |
|-C| full path to a control cohort VCF, merged with Truvari |
|-Q| full path to a query cohort VCF, merged with Truvari |
|-O| output directory name (default is needLR_output relative to current directory) |
|-R| restrict analysis to a region (e.g. chr22:12345-23456 **or** chr22)|
|-L| A .txt file that lists the full file path(s) to the query vcf(s) |
See above for recommended sniffles2 version/parameters |

General Annotation options
| | |
| :------------ |:-------------|
|--all| annotate VCF with all available options listed below (default TRUE) |

A la Carte Annotation options
| | |
| :------------ |:-------------|
|--omim| annotate VCF with OMIM phenotypes and modes of inheritance (default FALSE) |
|--hpo| annotate VCF with OMIM phenotype associated HPO terms (default FALSE) |
|--gencc| annotate VCF with GENCC phenotypes, level of support, and modes of inheritance (default FALSE) |
|--pli| annotate VCF with probability of loss-of-function intolerance (pLI) scores from gnomADg v4.1 (default FALSE) |
|--utr| annotate SV UTR overlap (gencodev45) (default FALSE) |
|--cds| annotate SV coding exon overlap (gencodev45) (default FALSE) |
|--oreganno| annotate SV regulatory element overlap (ORegAnno)(default FALSE) |
|--tre| annotate SV tandem repeat overlap (VAMOS) (default FALSE) |
|--mapflags| annotate SV with overlap of difficult to map regions (centromeres, telomeres, repeats,segdups, gaps, homopolymers) (default FALSE) |
|--hiconf| annotate SVs fully contained in high-confidence regions (default FALSE) |

#### Examples

Compare a single query VCF to the 500 1KGP database and apply all available annotations

```
needLR annotate examples/inputs/single_genome_example_chr22.vcf.gz
```

>[!NOTE]
> This example willl run much more quickly (and equivalently) if option `-R chr22` is included, since the input\
VCF is limited to SVs on chr22.

Output for this example: `examples/outputs/single_genome_example_chr22_needLR_1kg_v4.0/`

This example is included in the docker image and can be run like so:

```
OUTPUTDIR=/full/path/to/output/loc

docker run -v ${OUTPUTDIR}:/mnt/outputs \
  miragale/needlr:latest needLR annotate -O /mnt/outputs/needLR_output /mnt/needLR_examples/inputs/single_genome_example_chr22.vcf.gz
```


Compare a list of query VCFs to a different merged VCF and annotate with only OMIM and hiconfidence regions
```
needLR annotate -L examples/inputs/list_of_samples.txt -C examples/inputs/merged_cohort_chr22.vcf.gz --omim --hiconf
```

Outputs for this example:
```
examples/outputs/single_genome_example_chr22_needLR_customControl_v4.0/
examples/outputs/HG005_Pb_hantrio_sniffles_chr22_needLR_customControl_v4.0/
```

Compare the SVs in a merged cohort VCF to the 500 1KGP database and limit analysis to a smaller region
```
needLR annotate -Q examples/inputs/merged_cohort_chr22.vcf.gz -R chr22:10731900-11588324
```

Output for this example `examples/outputs/merged_cohort_chr22_needLR_1kg_v4.0/`


### Subcommand: comparator

`comparator` requires a single query VCF and at least one parental VCF provided with option `-P`

All options

| Option | Description |
| :------------ |:-------------|
|-P| full path to one or two parental VCFs[gz], separated by commas. Maternal vcf is assumed to be first vcf |
|-T| additional CPU threads to pass to bcftools |
|-C| full path to a control cohort VCF, merged with Truvari |
|-R| restrict analysis to a region (e.g. chr22:12345-23456 **or** chr22)|

All annotation options available in `annotate` are also available in `comparator`.

#### Examples

Compare a proband VCF to two parental VCFs along with the 500 1KGP database and apply all available annotations

```
needLR comparator -P examples/inputs/trio/HG007_Mo_hantrio_sniffles_chr22.vcf.gz,examples/inputs/trio/HG006_Fa_hantrio_sniffles_chr22.vcf.gz examples/inputs/trio/HG005_Pb_hantrio_sniffles_chr22.vcf.gz
```

Output for this example `examples/outputs/HG005_Pb_hantrio_sniffles_chr22_needLR_TRIO_1kg_v4.0/`


Compare a proband VCF to a single parent VCF along with a custom control VCF and apply only gencc annotations
```
needLR comparator -P examples/inputs/trio/HG007_Mo_hantrio_sniffles_chr22.vcf.gz -C examples/inputs/merged_cohort_chr22.vcf.gz --gencc examples/inputs/trio/HG005_Pb_hantrio_sniffles_chr22.vcf.gz
```

Output for this example `examples/outputs/HG005_Pb_hantrio_sniffles_chr22_needLR_DUO_customControl_v4.0/`



>[!NOTE]
>`needLR comparator` is imperfect in predicting inherited vs. _de novo_ SVs. The annotation is fully dependent on how well the SVs were merged. We recommend using needLR as a starting point and then manually inspecting inheritance in IGV. Observe extra precaution in VNTR regions.

### Subcommand: bed

`bed` requires only the path to a sorted bed file. Annotation columns are appended to the end of any existing columns

All options

| Option | Description |
| :------------ |:-------------|
|-L| .txt file list of full file paths to query beds |
|-O| output directory name (default is needLR_output relative to current directory) |
|-R| restrict analysis to a region (e.g. chr22:12345-23456 **or** chr22)|

All annotation options available in `annotate` are also available in `bed`.

#### Examples

Annotate kinnex data, limit to chr22, apply all available annotations

```
needLR bed -R chr22 examples/inputs/kinnex_example_chr22.bed
```
Output for this example `examples/outputs/kinnex_example_chr22_needLR_bed_v4.0/`


## OUTPUT

Output for each VCF is contained in a subdirectory of `needLR_output` by default. Modify this parent output directory with flag `-O`

Output:

| Output | Description |
|:------------|:-------------|
|{SAMPLE_ID}_RESULTS.tsv| Annotated query or cohort SVs (the main output)|
|{SAMPLE_ID}_RESULTS_unique.tsv| Annotated SVs that are unique to the query sample or cohort (not seen in the control cohort) |
|{SAMPLE_ID}_RESULTS_0.01.tsv| Annotated SVs with an AF <=0.01  (in relation to the control cohort)|
|{SAMPLE_ID}_RESULTS.vcf.gz| Annotated query SVs in vcf format |
|{SAMPLE_ID}_RESULTS_denovo.tsv| Annotated SVs that are unique to the proband (`comparator` trios only) |
|{SAMPLE_ID}_RESULTS_not_inherited.tsv| Annotated SVs that are not present in the comparator (`comparator` duos only) |


>[!NOTE]
>needLR generates many temporary files when running, this can add up to ~20Gb for ONT WGS datasetsless 

## ANALYZING OUTPUT

The `{SAMPLE_ID}_RESULTS*.tsv` files can easily be opened in Excel. 
_Be sure that Excel is set up to delimit columns by tab (and only tab)_

Below are the output columns. Some are specific to the needLR subcommand used. "Query SV" is used to refer to the SV identified in the proband/query sample _or_, in the case of `bed`, it is any SV. 
>[!NOTE]
>The SV start, end, length, REF, and ALT are the data from the most common SV of the merged SVs at that locus. Thus, this data in "query only" SVs (SVs seen in the query sample but not the control samples) will exactly match the query original vcf. If the query sample has an SV that is shared with control samples, it will be annotated with the data from the most common SV in that merge. 

### Columns always output when using VCF input

| Column Name           | Column Description                                                                 |
|:-----------------------|:-----------------------------------------------------------------------------------|
| Chr                   | Query SV chromosome                                                               |
| Start_Pos             | Query SV start coordinate (hg38)                                                   |
| End_Pos               | Query SV end coordinate (hg38)                                                     |
| REF             | Reference allele associated with the query SV  |
| ALT             | Alternate allele associated with the query SV   |
| SV_Length             | Query SV length                                                                    |
| SV_Type               | Query SV type                                                                      |
| Query ID               | Unique ID of the query sample                                                                     |
| Ctrl_support          | Samples in the control dataset that share the SV with the query sample |
| Allele_Freq_ALL       | Allele frequency of SV in control dataset                                   |
| Pop_Count_ALL         | How many control samples have the SV                             |
| Pop_Freq_ALL          | Frequency (%) of control samples with SV                         |
| Allele_Count_ALL      | How many control samples have the SV                              |
| Allele_Freq_ALL       | Frequency (%) of control samples alleles with SV                          |
| GT_homWT              | Number of control samples that are homozygous for the reference allele at the locus  |
| GT_het                | Number of control that are heterozygous for the SV allele at the locus |
| GT_homVAR             | Number of control samples that are homozygous for the SV allele at the locus  |
| HWE                   | Is the SV in Hardy-Weinberg equilibrium in the 1KGP samples                |

### Annotation columns
| Column Name           | Column Description                                                                 |
|:-----------------------|:-----------------------------------------------------------------------------------|
| Genes                 | Genes that the SV intersects with (gencode v45, 5kbp slop)        |
| OMIM                  | OMIM phenotypes associated with any gene the SV intersects (OMIM 01/2026)           |
| OMIM_MOI               | Mode of inheritance associated with phenotypes and genes the SV overlaps (OMIM 01/2026)           |
| GenCC                  | GenCC phenotypes associated with any gene the SV intersects (GenCC 20250510)      |
| GenCC_Support          | Support level(s) for GenCC phenotypes for overlapping genes (GenCC 20250510)  |
| GenCC_MOI               | Mode of inheritance for GenCC phenotypes for overlapping genes (GenCC 20250510) |
| HPO                  | HPO phenotypes associated with any gene the SV intersects           |
| pLI                  | pLI score for each gene the SV intersects           |
| UTR                | If the SV intersects with a UTR (gencode v45)         |
| CDS                | If the SV intersects with a coding exon (gencode v45)         |
| ORegAnno                | If the SV intersects with an ORegAnno region          |
| Centromeric           | If the SV intersects with a centromere (UCSC hg38)                                 |
| Pericentromeric       | If the SV intersects with a pericentromeric region (+/-5Mb on either side of UCSC-defined centromere) |
| Telomeric             | If the SV intersects with a telomere (5Mb of either end of a chromosome)           |
| Vamos                   | If the SV intersects with an STR/VNTR (vamos original motifs, n=148) |
| Segdup                | If the SV intersects with a segmental duplication (Genome in a Bottle v3.3)        |
| Repeat                | If the SV intersects with a repeat region (UCSC hg38 repeat masker)                |
| Gap                   | If the SV intersects with an hg38 gap region (UCSC hg38 mapping and sequencing: gap) |
| Homopolymer                   | If the SV intersects with an homopolymer >50bp |
| HiConf                | If the SV is fully contained within a high confidence region (Genome in a Bottle T2TQ100-V1.0_stvar) |

### Columns output when using 1KGP dataset as a control reference
| Column Name           | Column Description                                                                 |
|:-----------------------|:-----------------------------------------------------------------------------------|
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

### Columns output when using a proband or single vcf
| Column Name           | Column Description                                                                 |
|:-----------------------|:-----------------------------------------------------------------------------------|
| Genotype              | Query SV genotype                                                                  |
| Alt_reads             | Number of reads in the query sample supporting the SV                              |
| Ref_reads             | Number of reference reads in the query sample at the SV locus                      |
| Total_reads           | Total number of reads at the SV locus (in the query sample)                        |

### Columns output when using a query cohort vcf
| Column Name           | Column Description                                                                 |
|:-----------------------|:-----------------------------------------------------------------------------------|
| Cohort_support        | Sample ID of the samples in the input query cohort with the SV |
| Cohort_Pop_Count  | How many query cohort samples have the SV |
| Cohort_Pop_Freq | Frequency (%) of query cohort samples with SV |
| Cohort_Allele_Count    | How many query cohort alleles have the SV |
| Cohort_Allele_Freq    | Frequency (%) of query cohort alleles with SV |

### Columns output when using comparator
| Column Name           | Column Description                                                                 |
|:-----------------------|:-----------------------------------------------------------------------------------|
| Maternal_genotype  |   Genotype of SV in first parental VCF  (Parental_genotype when one VCF supplied) |
| Maternal_Alt_Reads  |   Count of reads supporting SV in first parental VCF  (Parental_genotype when one VCF supplied) |
| Maternal_Ref_Reads  |   Count of reads supporting reference in first parental VCF  (Parental_genotype when one VCF supplied) |
| Maternal_Total_Reads  |   Total count of reads at SV locus in first parental VCF  (Parental_genotype when one VCF supplied) |
| Paternal_genotype  |   Genotype of SV in second parental VCF  (not included when one VCF supplied) |
| Maternal_Alt_Reads  |   Count of reads supporting SV  in second parental VCF  (Parental_genotype when one VCF supplied) |
| Maternal_Ref_Reads  |   Count of reads supporting reference in second parental VCF  (Parental_genotype when one VCF supplied) |
| Maternal_Total_Reads  |   Total count of reads at SV locus in second parental VCF  (Parental_genotype when one VCF supplied) |
| Inheritance   |   Predicted inheritence pattern |


## RECOMMENDATIONS FOR ANALYSIS

In our experience using sniffles_v2.6.2 vcfs as queries (ONT sequenced, ~20-30X depth of coverage), we tend to see ~200-500 unique SVs per sample (ie. 1KGP allele frequency = 0). To further prioritize candidate pathogenic SVs, we recommend the following considerations:

Highest priority SVs:
* In CDS or UTR of a gene associated with a phenotype (OMIM, HPO, GenCC)
* In an annotated ORegAnno region
* >10 supporting reads
* pLI > 0.9

Regions of lower confidence:
* Inside tandem repeat or segmental duplication regions (especially INS)
* Within centromeres/telomeres


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
