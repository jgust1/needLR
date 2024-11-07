needLR_trio is deigned to annotate a query SV with inheritance information when parental vcfs are available. Unlike needLR (standard), needLR_trio can only run one trio at a time.

## RUN  NEEDLR_TRIO

#### The needLR_trio_3.3.sh command takes 4 required arguments:

| Flag | Description |
| :------------ |:-------------|
|-p| The full file path to the proband vcf (The vcfs must be gzipped (*.vcf.gz) and have an index in the same directory as the vcf)|
|-m| The full file path to the maternal vcf (The vcfs must be gzipped (*.vcf.gz) and have an index in the same directory as the vcf)|
|-f| The full file path to the paternal vcf (The vcfs must be gzipped (*.vcf.gz) and have an index in the same directory as the vcf)|
|-g| A fasta file for a reference genome (we use the hg38 reference recommended [here](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use))|

#### To run needLR:
```
needLR_3.3.sh -f /path/to/proband.vcf -m /path/to/maternal.vcf -f /path/to/paternal.vcf -g /file/path/to/reference/genome.fa
```
