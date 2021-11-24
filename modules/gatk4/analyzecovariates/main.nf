// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_ANALYZECOVARIATES {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::gatk4=4.2.3.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.2.3.0--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/gatk4:4.2.3.0--hdfd78af_0"
    }

    input:
    tuple val(meta), path(firstpass), path(secondpass), path(other)
    val csvout
    val ignoretimewarning

    output:
    tuple val(meta), path("*.pdf"), emit: plots
    tuple val(meta), path("*.csv"), optional:true, emit: csv
    path "versions.yml"           , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def inputsCommand = ''
    if(other) {
    inputsCommand = "-before ${firstpass} -after ${secondpass} -bqsr ${other}"
    } else if(secondpass) {
    inputsCommand =  "-before ${firstpass} -after ${secondpass}"
    } else {
    inputsCommand = "-bqsr ${firstpass}"
    }
    ignoreTimeCommand = ignoretimewarning ? "--ignore-last-modification-times" : ''
    csvCommand = csvout ? "-csv ${prefix}.csv" : ''

    """
    echo ${PATH}
    gatk AnalyzeCovariates  \\
        ${inputsCommand} \\
        ${ignoreTimeCommand} \\
        -plots ${prefix}.pdf \\
        ${csvCommand} \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
