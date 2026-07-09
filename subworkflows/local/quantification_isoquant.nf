/*
 * Classical quantification path:
 * UMI-tools dedup on genome BAM → IsoQuant per-barcode quantification
 */

include { FILTER_PRIMARY_BAM } from '../../modules/local/filter_primary_bam/main'
include { UMITOOLS_DEDUP     } from '../../modules/nf-core/umitools/dedup/main'
include { SAMTOOLS_INDEX     } from '../../modules/local/samtools_index/main'
include { SAMTOOLS_FAIDX     } from '../../modules/local/samtools_faidx/main'
include { ISOQUANT           } from '../../modules/local/isoquant/main'

workflow QUANTIFICATION_ISOQUANT {

    take:
    ch_tagged_bam    // channel: [ meta, bam ]
    ch_tagged_index  // channel: [ meta, bai ]
    ch_genome_fasta  // path(fasta)
    ch_genome_gtf    // channel: path(gtf)

    main:
    ch_versions = channel.empty()

    ch_tagged = ch_tagged_bam
        .map { meta, bam -> [ meta.id, meta, bam ] }
        .join(ch_tagged_index.map { meta, bai -> [ meta.id, meta, bai ] })
        .map { sample_id, meta, bam, m2, bai -> [ meta, bam, bai ] }

    if (params.isoquant_use_primary) {
        FILTER_PRIMARY_BAM (
            ch_tagged.map { meta, bam, bai -> [ meta, bam ] }
        )
        ch_versions = ch_versions.mix(FILTER_PRIMARY_BAM.out.versions.first())

        ch_dedup_in = FILTER_PRIMARY_BAM.out.bam
            .map { meta, bam -> [ meta.id, meta, bam ] }
            .join(FILTER_PRIMARY_BAM.out.index.map { meta, bai -> [ meta.id, meta, bai ] })
            .map { sample_id, meta, bam, m2, bai -> [ meta, bam, bai ] }
    }
    else {
        ch_dedup_in = ch_tagged
    }

    UMITOOLS_DEDUP (
        ch_dedup_in,
        false
    )
    // umitools version is published via topic: versions (handled in starlight.nf)

    SAMTOOLS_INDEX (
        UMITOOLS_DEDUP.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    SAMTOOLS_FAIDX (
        Channel.value([ [ id: 'genome' ], ch_genome_fasta ])
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())

    ch_reference = SAMTOOLS_FAIDX.out.fasta.first()
        .combine(ch_genome_gtf)
        .map { meta, fasta, fai, gtf -> [ fasta, fai, gtf ] }

    ch_isoquant_in = SAMTOOLS_INDEX.out.bam_index
        .map { meta, bam, bai -> [ meta.id, meta, bam, bai ] }
        .combine(ch_reference)
        .map { sample_id, meta, bam, bai, fasta, fai, gtf ->
            [ meta, bam, bai, fasta, fai, gtf ]
        }

    ISOQUANT (
        ch_isoquant_in,
        params.isoquant_read_group
    )
    ch_versions = ch_versions.mix(ISOQUANT.out.versions.first())

    emit:
    grouped_gene_mtx            = ISOQUANT.out.grouped_gene_mtx
    grouped_gene_mtx_barcodes   = ISOQUANT.out.grouped_gene_mtx_barcodes
    grouped_gene_mtx_features   = ISOQUANT.out.grouped_gene_mtx_features
    grouped_transcript_mtx      = ISOQUANT.out.grouped_transcript_mtx
    grouped_transcript_mtx_barcodes = ISOQUANT.out.grouped_transcript_mtx_barcodes
    grouped_transcript_mtx_features = ISOQUANT.out.grouped_transcript_mtx_features
    gene_counts                 = ISOQUANT.out.gene_counts
    transcript_counts           = ISOQUANT.out.transcript_counts
    read_assignments            = ISOQUANT.out.read_assignments
    versions                    = ch_versions
}
