process run_needLR_annotate {
    publishDir "${params.publish_dir}/${sample_id}", mode: 'copy'
    
    input:
        path( vcf )
        val( sample_id )
        val( region )
        val( annotations )
        val( isMerged )

    output:
        path "${sample_id}_needLR_1kg_v4.0", emit: results_folder

    script:
        """
        expectedOutput="${sample_id}"
        if [[ \${expectedOutput} =~ "gz" ]]
        then
            actualOutput=\${expectedOutput%*.vcf.gz}
        else
            actualOutput=\${expectedOutput%*.vcf}
        fi
        argstopass=()
        annotationString=${annotations}
        annotationList=( \$( echo \${annotationString}| tr ',' ' ') )
        if [[ ! \${annotationString} =~ "all" ]]
        then
            for anno in \${annotationList[@]}
            do
                argstopass+=('--'\$anno)
            done
        fi
        if [[ ${region} != "none" ]]
        then
            argstopass+=( -R ${region} )
        fi
        if [[ ${isMerged} == TRUE ]]
        then
            needLR annotate -T ${task.cpus} -Q ${vcf} \${argstopass[@]}
        else
            if [[ \${#argstopass[@]} -gt 0 ]]
            then
                needLR annotate -T ${task.cpus} \${argstopass[@]} ${vcf}
            else
                needLR annotate -T ${task.cpus} ${vcf}
            fi
        fi
        mv needLR_output/\${actualOutput}_needLR_1kg_v4.0 ${sample_id}_needLR_1kg_v4.0
        """
}

process run_needLR_annotate_custom_controls {
    publishDir "${params.publish_dir}/${sample_id}", mode: 'copy'
    
    input:
        path( vcf )
        path( control_vcf )
        path( subpopulations )
        val( sample_id )
        val( region )
        val( annotations )
        val( isMerged )

    output:
        path "${sample_id}_needLR_customControl_v4.0", emit: results_folder

    script:
        """
        expectedOutput="${sample_id}"
        if [[ \${expectedOutput} =~ "gz" ]]
        then
            actualOutput=\${expectedOutput%*.vcf.gz}
        else
            actualOutput=\${expectedOutput%*.vcf}
        fi
        argstopass=()
        annotationString=${annotations}
        annotationList=( \$( echo \${annotationString}| tr ',' ' ') )
        if [[ ! \${annotationString} =~ "all" ]]
        then
            for anno in \${annotationList[@]}
            do
                argstopass+=('--'\$anno)
            done
        fi
        if [[ ${region} != "none" ]]
        then
            argstopass+=( -R ${region} )
        fi
        if [[ ${control_vcf.name} != "NO_FILE" ]]
        then
            argstopass+=( -S ${subpopulations} )
        fi
        if [[ ${isMerged} == TRUE ]]
        then
            needLR annotate -T ${task.cpus} -Q ${vcf} -C ${control_vcf} \${argstopass[@]}
        else
            if [[ \${#argstopass[@]} -gt 0 ]]
            then
                needLR annotate -T ${task.cpus} -C ${control_vcf} \${argstopass[@]} ${vcf}
            else
                needLR annotate -T ${task.cpus} -C ${control_vcf} ${vcf}
            fi
        fi
        mv needLR_output/\${actualOutput}_needLR_customControl_v4.0 ${sample_id}_needLR_customControl_v4.0
        """
}





