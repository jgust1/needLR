# Nextflow wrapper for NeedLR 4.0

A wrapper around needLR to make running multiple VCFs concurrently more convenient.

Currently implemented for:

    * single VCFs
    * merged cohort VCFs
    * custom controls and 1kg controls
    * bed annotations

Currently **not** available for:

    * trios
    * duos

## Usage

```
nextflow run main.nf --subcommand annotate --fileList /path/to/your/file/list.txt --cpus 50 --maxSimultaneous 5
```

## Installation

Clone this repository and run `main.nf` from the `nextflow` subdirectory of `needLR`. You will need to run from an environment with both `needlR=4.0` and `nextflow=25.10.4` available.

Make one with conda:

```conda create env -n needLRNextflow
conda install -n needLRNextflow -c bioconda needlr=4.0 nextflow=25.10.4
conda activate needLRNextflow
```

## All options

*Options that are necessary for the nextflow to run are bolded*

|Option|Description|
|------|-----------|
| **--subcommand** | needLR subcommand to run. Valid values are `annotate` and `bed`|
| **--fileList** | path to your input file list|
| --control_vcf| path to a Truvari merged vcf[.gz] to use as population controls |
| --subpops| path to a .csv or .tsv file listing sample names and subpopulation membership in provided `control_vcf` |
| --keep1kg | run 1kg controls along with provided `control_vcf` |
| --melt1kg| merge provided subpopulation counts with 1kg subpopulation counts where labels overlap |
| --region | a genomic region to restrict analysis to (e.g. chrX:12345-23456)|
| --globalcpus | total number of CPUs for nextflow to use (default: 1) |
| --maxSimultaneous | number of needLR instances to run simultaneously (default: 1)|
| --merged | input files are merged VCFs (true or false)|
| --publish_dir | output directory for needLR files (default: results)|
| --no_annotations | do not add annotations, population stats only (default false)|
| --omim |include omim annotations (default true)|
| --gencc |include gencc annotations (default true)|
| --hpo |include hpo term annotations (default true)|
| --pli |include pLI annotations (default true)|
| --utr |include UTR annotations (default true)|
| --cds |include CDS annotations (default true)|
| --oreganno |include ORegAnno annotations (default true)|
| --mapflags |include difficult mapping region flags (default true)|
| --hiconf |annotate overlap with high confidence regions (default true)|
| --example | only include this if you are running the example. Allows for paths relative to launchdir in `fileList`|

See the [Nextflow Documentation](https://www.nextflow.io/docs/stable/reference/cli.html) for more general command line options. This is especially relevant if running the nextflow on a cluster.


## Examples

One example of an input file list is included in `needLR/examples/input/nextflow_input/list_ofsamples.txt`

```
nextflow run main.nf --subcommand annotate --fileList ../examples/nextflow_input/list_of_samples.txt --example=true --cpus 20 --maxSimultaneous 2
```



