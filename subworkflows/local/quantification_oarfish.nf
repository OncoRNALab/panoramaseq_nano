/*
 * Oarfish quantification path (scnanoseq-style):
 * R2 FASTQ → reference transcriptome align → tag → UMI dedup → CB sort → Oarfish → gene aggregation
 *
 * Uses Oarfish 0.9.4 (not the older scnanoseq pin).
 */

include { GFFREAD                                } from '../../modules/nf-core/gffread/main'
include { MINIMAP2_INDEX as MINIMAP2_INDEX_TXOME } from '../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_TXOME } from '../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_FILTER_MAPPED                 } from '../../modules/local/samtools_filter_mapped/main'
include { TAG_BAM                                } from '../../modules/local/tag_bam/main'
include { UMITOOLS_DEDUP                         } from '../../modules/nf-core/umitools/dedup/main'
include { SORT_BAM_CB                            } from '../../modules/local/sort_bam_cb/main'
include { OARFISH                                } from '../../modules/local/oarfish/main'
include { OARFISH_AGGREGATE                      } from '../../modules/local/oarfish_aggregate/main'

workflow QUANTIFICATION_OARFISH {

    take:
    ch_cdna_fastq   // channel: [ meta, fastq ]
    ch_bc_tags      // channel: [ meta, tags ]
    ch_genome_fasta // path(fasta)
    ch_genome_gtf   // channel: path(gtf)

    main:
    ch_versions = channel.empty()

    ch_ref_gtf = ch_genome_gtf.map { gtf ->
        [ [ id: 'reference', strandedness: params.strandedness ], gtf ]
    }

    GFFREAD (
        ch_ref_gtf,
        ch_genome_fasta
    )

    MINIMAP2_INDEX_TXOME (
        GFFREAD.out.gffread_fasta
    )

    ch_txome_index = MINIMAP2_INDEX_TXOME.out.index.first()

    ch_txome_align = ch_cdna_fastq
        .combine(ch_txome_index)
        .map { meta, reads, idx_meta, index ->
            [ meta, reads, idx_meta, index ]
        }

    MINIMAP2_ALIGN_TXOME (
        ch_txome_align.map { meta, reads, idx_meta, index -> [ meta, reads ] },
        ch_txome_align.map { meta, reads, idx_meta, index -> [ idx_meta, index ] },
        true,   // bam_format
        "bai",  // bam_index_extension
        false,  // cigar_paf_format
        false   // cigar_bam
    )

    SAMTOOLS_FILTER_MAPPED (
        MINIMAP2_ALIGN_TXOME.out.bam
            .join(MINIMAP2_ALIGN_TXOME.out.index, by: 0)
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FILTER_MAPPED.out.versions.first())

    ch_tag_input = SAMTOOLS_FILTER_MAPPED.out.bam
        .map { meta, bam -> [ meta.id, meta, bam ] }
        .join(ch_bc_tags.map { meta, tags -> [ meta.id, meta, tags ] })
        .map { sample_id, meta, bam, m2, tags -> [ meta, bam, tags ] }

    TAG_BAM (
        ch_tag_input.map { meta, bam, tags -> [ meta, bam ] },
        ch_tag_input.map { meta, bam, tags -> [ meta, tags ] }
    )
    ch_versions = ch_versions.mix(TAG_BAM.out.versions.first())

    ch_dedup_in = TAG_BAM.out.bam
        .map { meta, bam -> [ meta.id, meta, bam ] }
        .join(TAG_BAM.out.index.map { meta, bai -> [ meta.id, meta, bai ] })
        .map { sample_id, meta, bam, m2, bai -> [ meta, bam, bai ] }

    UMITOOLS_DEDUP (
        ch_dedup_in,
        false
    )

    SORT_BAM_CB (
        UMITOOLS_DEDUP.out.bam
    )
    ch_versions = ch_versions.mix(SORT_BAM_CB.out.versions.first())

    OARFISH (
        SORT_BAM_CB.out.bam
    )
    ch_versions = ch_versions.mix(OARFISH.out.versions.first())

    ch_aggregate_in = OARFISH.out.features
        .join(OARFISH.out.barcodes, by: 0)
        .join(OARFISH.out.mtx, by: 0)
        .combine(ch_genome_gtf)
        .map { meta, features, barcodes, mtx, gtf ->
            [ meta, features, barcodes, mtx, gtf ]
        }

    OARFISH_AGGREGATE (
        ch_aggregate_in
    )
    ch_versions = ch_versions.mix(OARFISH_AGGREGATE.out.versions.first())

    emit:
    transcript_matrix = OARFISH_AGGREGATE.out.transcript_matrix
    gene_matrix       = OARFISH_AGGREGATE.out.gene_matrix
    stats             = OARFISH_AGGREGATE.out.stats
    meta_info         = OARFISH.out.meta_info
    versions          = ch_versions
}
