// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process OPTITYPE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::optitype=1.3.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/optitype:1.3.5--0"
    } else {
        container "quay.io/biocontainers/optitype:1.3.5--0"
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path('*.tsv'), emit: log
    tuple val(meta), path('*.pdf'), emit: pdf
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    # Create a config for OptiType on a per sample basis with options.args2
    configbuilder --max-cpus $task.cpus $options.args2 > config.ini

    # Run the actual OptiType typing with options.args
    OptiTypePipeline.py -i ${bam} -c config.ini --${meta.seq_type} $options.args --outdir .

    cat \$(which OptiTypePipeline.py) 2>&1 ${software}.version.txt
    """
}
