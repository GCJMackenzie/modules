//
// Run GATK mutect2 in tumor only mode, getepileupsummaries, calculatecontamination and filtermutectcalls
//

params.mutect2_options     = [:]
params.getpileup_options   = [:]
params.calccontam_options  = [:]
params.filtercalls_options = [suffix: '_filtered']

include { GATK4_MUTECT2                as MUTECT2 }                  from '../../../modules/gatk4/mutect2/main'                addParams( options: params.mutect2_options )
include { GATK4_GETPILEUPSUMMARIES     as GETPILEUPSUMMARIES }       from '../../../modules/gatk4/getpileupsummaries/main'     addParams( options: params.getpileup_options )
include { GATK4_CALCULATECONTAMINATION as CALCULATECONTAMINATION }   from '../../../modules/gatk4/calculatecontamination/main' addParams( options: params.calccontam_options )
include { GATK4_FILTERMUTECTCALLS      as FILTERMUTECTCALLS }        from '../../../modules/gatk4/filtermutectcalls/main'      addParams( options: params.filtercalls_options )

workflow GATK_TUMOR_ONLY_SOMATIC_VARIANT_CALLING {
    take:
    ch_mutect2_in             // channel: [ val(meta), [ input ], [ input_index ], [] ]
    fasta                     // channel: /path/to/reference/fasta
    fastaidx                  // channel: /path/to/reference/fasta/index
    dict                      // channel: /path/to/reference/fasta/dictionary
    germline_resource         // channel: /path/to/germline/resource
    germline_resource_idx     // channel: /path/to/germline/index
    panel_of_normals          // channel: /path/to/panel/of/normals
    panel_of_normals_idx      // channel: /path/to/panel/of/normals/index
    interval_file             // channel: /path/to/interval/file


    main:
    ch_versions = Channel.empty()
    input = channel.from(ch_mutect2_in)

    //
    //Perform variant calling using mutect2 module in tumor single mode.
    //
    MUTECT2 ( input , true , false , false , [] , fasta , fastaidx , dict , germline_resource , germline_resource_idx , panel_of_normals , panel_of_normals_idx )
    ch_versions = ch_versions.mix(MUTECT2.out.versions)

    //
    //Generate pileup summary table using getepileupsummaries.
    //
    pileup_input = channel.from(ch_mutect2_in)
    pileup_input = pileup_input.map {
        meta, input_file, input_index, which_norm ->
        [meta, input_file[0], input_index[0]]
    }
    GETPILEUPSUMMARIES ( pileup_input , germline_resource , germline_resource_idx , interval_file )
    ch_versions = ch_versions.mix(GETPILEUPSUMMARIES.out.versions)

    //
    //Contamination and segmentation tables created using calculatecontamination on the pileup summary table.
    //
    ch_pileup = GETPILEUPSUMMARIES.out.table.collect()
    ch_pileup.add([])
    CALCULATECONTAMINATION ( ch_pileup, true )
    ch_versions = ch_versions.mix(CALCULATECONTAMINATION.out.versions)

    //
    //Mutect2 calls filtered by filtermutectcalls using the contamination and segmentation tables.
    //
    ch_vcf =           MUTECT2.out.vcf.collect()
    ch_tbi =           MUTECT2.out.tbi.collect()
    ch_stats =         MUTECT2.out.stats.collect()
    ch_stats.add([])
    ch_segment =       CALCULATECONTAMINATION.out.segmentation.collect()
    ch_contamination = CALCULATECONTAMINATION.out.contamination.collect()
    ch_contamination.add([])
    ch_filtermutect_in = ch_vcf.combine(ch_tbi, by: 0).combine(ch_stats, by: 0).combine(ch_segment, by: 0).combine(ch_contamination, by: 0)
    FILTERMUTECTCALLS ( ch_filtermutect_in, fasta, fastaidx, dict )
    ch_versions = ch_versions.mix(FILTERMUTECTCALLS.out.versions)

    emit:
    mutect2_vcf            = MUTECT2.out.vcf.collect()                             // channel: [ val(meta), [ vcf ] ]
    mutect2_index          = MUTECT2.out.tbi.collect()                             // channel: [ val(meta), [ tbi ] ]
    mutect2_stats          = MUTECT2.out.stats.collect()                           // channel: [ val(meta), [ stats ] ]

    pileup_table           = GETPILEUPSUMMARIES.out.table.collect()                // channel: [ val(meta), [ table ] ]

    contamination_table    = CALCULATECONTAMINATION.out.contamination.collect()    // channel: [ val(meta), [ contamination ] ]
    segmentation_table     = CALCULATECONTAMINATION.out.segmentation.collect()     // channel: [ val(meta), [ segmentation ] ]

    filtered_vcf           = FILTERMUTECTCALLS.out.vcf.collect()                   // channel: [ val(meta), [ vcf ] ]
    filtered_index         = FILTERMUTECTCALLS.out.tbi.collect()                   // channel: [ val(meta), [ tbi ] ]
    filtered_stats         = FILTERMUTECTCALLS.out.stats.collect()                 // channel: [ val(meta), [ stats ] ]

    versions               = ch_versions                                           // channel: [ versions.yml ]
}
