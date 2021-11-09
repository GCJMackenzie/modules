//
// Run GATK mutect2 in tumour only mode, getepileupsummaries, calculatecontamination and filtermutectcalls
//

params.mutect2_options     = [:]
params.getpileup_options   = [:]
params.calccontam_options  = [:]
params.filtercalls_options = [:]

include { GATK4_MUTECT2                     } from '../../../modules/gatk4/mutect2/main'                addParams( options: params.mutect2_options )
include { GATK4_GETPILEUPSUMMARIES          } from '../../../modules/gatk4/getepileupsummaries/main'    addParams( options: params.getpileup_options )
include { GATK4_CALCULATECONTAMINATION      } from '../../../modules/gatk4/calculatecontamination/main' addParams( options: params.calccontam_options )
include { GATK4_FILTERMUTECTCALLS           } from '../../../modules/gatk4/filtermutectcalls/main'      addParams( options: params.filtercalls_options )

workflow GATK_TUMOUR_ONLY_SOMATIC_VARIANT_CALLING {
    take:
    ch_mutect2_in       // channel: [ val(meta), [ input ], [ input_index ], [] ]
    fasta               // channel: /path/to/reference/fasta
    fastaidx            // channel: /path/to/reference/fasta/index
    dict                // channel: /path/to/reference/fasta/dictionary
    pon_name            // channel: name for panel of normals
    interval_file       // channel: /path/to/interval/file

    // tuple val(meta), path(vcf), path(tbi), path(intervalfile), [], []


    // tuple val(meta), path(genomicsdb)

    main:
    ch_versions      = Channel.empty()
    input = channel.from(ch_mutect2_in)
    //
    //Perform variant calling for each sample using mutect2 module in panel of normals mode.
    //
    GATK4_MUTECT2 ( input , false , true, false , [] , fasta , fastaidx , dict , [], [] , [] , [] )
    ch_versions = ch_versions.mix(GATK4_MUTECT2.out.versions.first())

    //
    //Convert all sample vcfs into a genomicsdb workspace using genomicsdbimport.
    //
    ch_vcf = GATK4_MUTECT2.out.vcf.collect{it[1]}.toList()
    ch_index = GATK4_MUTECT2.out.tbi.collect{it[1]}.toList()
    // ch_mutect2_out = ch_vcf.combine([ch_index])
    gendb_input = Channel.of([[ id:pon_name ]]).combine(ch_vcf).combine(ch_index).combine([interval_file]).combine(['']).combine([dict])
    GATK4_GENOMICSDBIMPORT ( gendb_input, false, false, false )
    ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions.first())

    //
    //Panel of normals made from genomicsdb workspace using createsomaticpanelofnormals.
    //
    GATK4_GENOMICSDBIMPORT.out.genomicsdb.view()
    GATK4_CREATESOMATICPANELOFNORMALS ( GATK4_GENOMICSDBIMPORT.out.genomicsdb, fasta, fastaidx, dict )
    ch_versions = ch_versions.mix(GATK4_CREATESOMATICPANELOFNORMALS.out.versions.first())

    emit:
    mutect2_vcf      = GATK4_MUTECT2.out.vcf.collect()                     // channel: [ val(meta), [ vcf ] ]
    mutect2_index    = GATK4_MUTECT2.out.tbi.collect()                     // channel: [ val(meta), [ tbi ] ]
    mutect2_stats    = GATK4_MUTECT2.out.stats.collect()                   // channel: [ val(meta), [ stats ] ]

    genomicsdb       = GATK4_GENOMICSDBIMPORT.out.genomicsdb     // channel: [ val(meta), [ genomicsdb ] ]

    pon_vcf          = GATK4_CREATESOMATICPANELOFNORMALS.out.vcf // channel: [ val(meta), [ vcf.gz ] ]
    pon_index        = GATK4_CREATESOMATICPANELOFNORMALS.out.tbi // channel: [ val(meta), [ tbi ] ]

    versions         = ch_versions                               // channel: [ versions.yml ]
}
