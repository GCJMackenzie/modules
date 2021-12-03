process GATK4_VARIANTRECALIBRATOR {
    tag "$meta.id"
    label 'process_low'

        conda (params.enable_conda ? "bioconda::gatk4=4.2.3.0" : null)
        container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/gatk4:4.2.3.0--hdfd78af_0' :
            'quay.io/biocontainers/gatk4:4.2.3.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf) , path(tbi)
    path fasta
    path fai
    path dict
    val allelespecific
    tuple path(resvcfs), path(restbis), val(reslabels)
    val annotation
    val mode
    val create_rscript

    output:
    tuple val(meta), path("*.recal")   , emit: recal
    tuple val(meta), path("*.idx")     , emit: idx
    tuple val(meta), path("*.tranches"), emit: tranches
    tuple val(meta), path("*plots.R")  , optional:true, emit: plots
    path "versions.yml"                , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    refCommand = fasta ? "-R ${fasta} " : ''
    alleleSpecificCommand = allelespecific ? '-AS' : ''
    resourceCommand = '--resource:' + reslabels.join( ' --resource:')
    annotationCommand = '-an ' + annotation.join( ' -an ')
    modeCommand = mode ? "--mode ${mode} " : 'SNP'
    rscriptCommand = create_rscript ? "--rscript-file ${prefix}.plots.R" : ''

    """
    gatk VariantRecalibrator \\
        ${refCommand} \\
        -V ${vcf} \\
        ${alleleSpecificCommand} \\
        ${resourceCommand} \\
        ${annotationCommand} \\
        ${modeCommand} \\
        -O ${prefix}.recal \\
        --tranches-file ${prefix}.tranches \\
        ${rscriptCommand}\\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}