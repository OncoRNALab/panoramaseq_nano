/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { SEQTK_SAMPLE           } from '../modules/nf-core/seqtk/sample/main'
include { CHOPPER                } from '../modules/nf-core/chopper/main'
include { RESTRANDER             } from '../modules/local/restrander/main'
include { EXTRACT_BARCODE        } from '../modules/local/extract_barcode/main'
include { SPLIT_READS            } from '../modules/local/split_reads/main'
include { QUIK_STARSOLO          } from '../modules/local/quik_starsolo/main'
include { MINIMAP2_INDEX         } from '../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN         } from '../modules/nf-core/minimap2/align/main'
include { TAG_BAM                } from '../modules/local/tag_bam/main'
include { QUANTIFICATION_EPI2ME  } from '../subworkflows/local/quantification_epi2me'
include { QUANTIFICATION_ISOQUANT } from '../subworkflows/local/quantification_isoquant'
include { QUANTIFICATION_OARFISH  } from '../subworkflows/local/quantification_oarfish'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_starlight_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow STARLIGHT {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    //
    // MODULE: Run FastQC on the full (unsubsampled) input
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // Optional subsampling: keep only the first N reads when --n_reads is set
    //
    if (params.n_reads) {
        SEQTK_SAMPLE (
            ch_samplesheet.map { meta, reads -> [ meta, reads, params.n_reads ] }
        )
        ch_reads = SEQTK_SAMPLE.out.reads
        // versions_seqtk uses topic: versions — collected automatically by Channel.topic("versions")
    } else {
        ch_reads = ch_samplesheet
    }

    //
    // MODULE: Optional Chopper quality trim/filter (before Restrander)
    // Pass --chopper_enabled true to activate; when false, reads go directly to Restrander.
    //
    if (params.chopper_enabled) {
        ch_chopper_fasta = params.chopper_contaminants
            ? channel.fromPath(params.chopper_contaminants, checkIfExists: true)
            : channel.value([])

        CHOPPER (
            ch_reads.map { meta, reads ->
                [ meta, reads instanceof List ? reads[0] : reads ]
            },
            ch_chopper_fasta.first()
        )
        ch_reads_for_restrand = CHOPPER.out.fastq
    }
    else {
        ch_reads_for_restrand = ch_reads
    }

    //
    // MODULE: Orient reads with Restrander (mandatory)
    //
    ch_restrander_config = channel.fromPath(params.restrander_config, checkIfExists: true)

    RESTRANDER (
        ch_reads_for_restrand,
        ch_restrander_config.first()
    )
    ch_oriented = RESTRANDER.out.reads
    ch_versions = ch_versions.mix(RESTRANDER.out.versions.first())

    //
    // MODULE: Extract Spatial Barcodes and UMIs using parasail SW alignment
    //
    EXTRACT_BARCODE (
        ch_oriented
    )
    ch_versions = ch_versions.mix(EXTRACT_BARCODE.out.versions.first())

    //
    // MODULE: Split oriented reads into R1 (barcode) and R2 (cDNA) using per-read positions
    //
    ch_split_input = ch_oriented.join(EXTRACT_BARCODE.out.tags)

    SPLIT_READS (
        ch_split_input.map { meta, reads, tags -> [ meta, reads ] },
        ch_split_input.map { meta, reads, tags -> [ meta, tags  ] }
    )
    ch_versions = ch_versions.mix(SPLIT_READS.out.versions.first())

    //
    // MODULE: QUIK barcode correction — R1 is the 36 bp barcode at position 0, R2 is clean cDNA
    //
    ch_barcode_file = params.barcode_whitelist
        ? channel.fromPath(params.barcode_whitelist, checkIfExists: true)
        : channel.empty()

    ch_quik_input = SPLIT_READS.out.r1.join(SPLIT_READS.out.r2)

    QUIK_STARSOLO (
        ch_quik_input,
        ch_barcode_file.first(),
        params.barcode_length
    )
    ch_versions = ch_versions.mix(QUIK_STARSOLO.out.versions.first())

    //
    // Oarfish: scnanoseq-style path — R2 FASTQ → transcriptome align (no genome alignment)
    //
    if (params.gene_quant_mode == 'oarfish' && params.genome_gtf && params.genome_fasta) {
        ch_gtf          = channel.fromPath(params.genome_gtf, checkIfExists: true)
        ch_genome_fasta = file(params.genome_fasta, checkIfExists: true)

        QUANTIFICATION_OARFISH (
            QUIK_STARSOLO.out.r2,
            EXTRACT_BARCODE.out.tags,
            ch_genome_fasta,
            ch_gtf.first()
        )
        ch_versions = ch_versions.mix(QUANTIFICATION_OARFISH.out.versions)
    }

    //
    // MINIMAP2: build index (once) then align filtered R2 (cDNA) reads — epi2me / isoquant only
    //
    if (params.gene_quant_mode in ['epi2me', 'isoquant']) {
        if (params.genome_mmi) {
            ch_genome_index = Channel.value([ [id: 'genome'], file(params.genome_mmi, checkIfExists: true) ])
        } else if (params.genome_fasta) {
            MINIMAP2_INDEX (
                Channel.value([ [id: 'genome'], file(params.genome_fasta, checkIfExists: true) ])
            )
            ch_genome_index = MINIMAP2_INDEX.out.index
            // versions_minimap2 uses topic: versions — collected automatically by Channel.topic("versions")
        } else {
            ch_genome_index = Channel.empty()
        }

        if (params.genome_fasta || params.genome_mmi) {
            MINIMAP2_ALIGN (
                QUIK_STARSOLO.out.r2,
                ch_genome_index.first(),
                true,   // bam_format
                "bai",  // bam_index_extension
                false,  // cigar_paf_format
                false   // cigar_bam
            )
            // versions_minimap2 uses topic: versions — collected automatically by Channel.topic("versions")

            ch_tag_bam_input = MINIMAP2_ALIGN.out.bam.join(EXTRACT_BARCODE.out.tags)

            TAG_BAM (
                ch_tag_bam_input.map { meta, bam, bc_tags -> [ meta, bam ] },
                ch_tag_bam_input.map { meta, bam, bc_tags -> [ meta, bc_tags ] }
            )
            ch_versions = ch_versions.mix(TAG_BAM.out.versions.first())

            if (params.genome_gtf && params.genome_fasta) {
                ch_gtf          = channel.fromPath(params.genome_gtf, checkIfExists: true)
                ch_genome_fasta = file(params.genome_fasta, checkIfExists: true)

                if (params.gene_quant_mode == 'epi2me' && params.stringtie_enabled) {
                    QUANTIFICATION_EPI2ME (
                        TAG_BAM.out.bam,
                        QUIK_STARSOLO.out.r2,
                        ch_gtf.first(),
                        ch_genome_fasta
                    )
                    ch_versions = ch_versions.mix(QUANTIFICATION_EPI2ME.out.versions)
                }
                else if (params.gene_quant_mode == 'isoquant') {
                    QUANTIFICATION_ISOQUANT (
                        TAG_BAM.out.bam,
                        TAG_BAM.out.index,
                        ch_genome_fasta,
                        ch_gtf.first()
                    )
                    ch_versions = ch_versions.mix(QUANTIFICATION_ISOQUANT.out.versions)
                }
            }
        }
    }

    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'starlight_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
