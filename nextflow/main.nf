nextflow.enable.dsl=2

log.info    """
NeedLR v 4.0 Nextflow Wrapper
-----------------------------

This wrapper assists in running needLR v4.0 on more than one sample concurrently

// Input parameters:

Region of Interest: ${params.region}
CPUS: ${params.cpus}
Merged VCFS:    ${params.merged}


// Run type: ${params.subcommand}

// Annotations:
Genes:                      true
OMIM:                       ${params.omim}
GenCC:                      ${params.gencc}
HPO Terms:                  ${params.hpo}
pLI:                        ${params.pli}
UTR:                        ${params.utr}
CDS:                        ${params.cds}
ORegAnno:                   ${params.oreganno}
Mapping Flags:              ${params.mapflags}
High Confidence Regions:    ${params.hiconf}

// Reference parameters
Control Reference Dataset: ${params.control_vcf}
Subpopulation Counts: ${params.subpops}
Include 1KG: ${params.keep1kg}
Merge 1KG with provide reference dataset: ${params.melt1kg}
"""

include { run_needLR_bed } from "${launchDir}/workflows/annotate_bed.nf"
include { run_needLR_annotate; run_needLR_annotate_custom_controls } from "${launchDir}/workflows/annotate_svs.nf"
//include { run_needLR_duo; run_needLR_duo_custom_controls; run_needLR_trio; run_needLR_trio_custom_controls } from "${launchDir}/workflows/comparators.nf"

workflow {
    def id_list
    if(params.example){
        def inputPath=file("${launchDir}")
        def filenames=file(params.fileList).readLines().collect { x -> inputPath.resolve(x) }
        ch_inputs=Channel.fromList(filenames )
        id_list = file(params.fileList).readLines().collect { x -> x.split('/')[-1]}
        ch_ids=Channel.fromList(id_list)

    }else{
        ch_inputs=Channel.fromPath(params.fileList)
                        .splitText()
        id_list = file(params.fileList).readLines().collect { x -> x.split('/')[-1]}
        ch_ids=Channel.fromList(id_list)
    }

    //adjust task.cpus based on how many instances of needlr are requested and params.maxChunk (default 5)
    
    annotation_string=params.annotation_string
    if(params.subcommand=="annotate"){ 
        if(params.control_vcf==null)
        {
            run_needLR_annotate(
                ch_inputs,
                ch_ids,
                params.region,
                params.annotation_string,
                params.merged )
        } else {
            ch_controls = Channel.fromPath(params.control_vcf)
            if(params.subpops!=null){
                ch_subpops = Channel.fromPath(params.subpops)
            }else{
                ch_subpops = Channel.fromPath("${projectDir}/assets/NO_FILE")
            }
            run_needLR_annotate_custom_controls(
                ${ch_inputs},
                ${ch_controls},
                ${ch_subpops},
                ${ch_ids},
                ${params.region},
                ${params.annotation_string},
                ${params.merged}
            )
        }
    }
    else if(params.subcommand=="bed"){
        run_needLR_bed(
                ${ch_inputs},
                ${ch_ids},
                ${params.region},
                ${annotation_string},
            )
    } else {
        throw new IllegalArgumentException("Invalid subcommand: ${params.subcommand}")
    }

}