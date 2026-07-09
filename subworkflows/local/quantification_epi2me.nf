/*
 * epi2me / wf-single-cell quantification path:
 * StringTie → txome align → gffcompare → assign_features → create_matrix
 */

include { FILTER_PRIMARY_BAM         } from '../../modules/local/filter_primary_bam/main'
include { STRINGTIE_STRINGTIE        } from '../../modules/nf-core/stringtie/stringtie/main'
include { GFFREAD                    } from '../../modules/nf-core/gffread/main'
include { MINIMAP2_INDEX as MINIMAP2_INDEX_TRANSCRIPTOME } from '../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_TRANSCRIPTOME } from '../../modules/nf-core/minimap2/align/main'
include { GFFCOMPARE                 } from '../../modules/nf-core/gffcompare/main'
include { SAMTOOLS_FAIDX             } from '../../modules/local/samtools_faidx/main'
include { SORT_BAM_NAMES             } from '../../modules/local/sort_bam_names/main'
include { EXTRACT_BAM_TAGS           } from '../../modules/local/extract_bam_tags/main'
include { ASSIGN_FEATURES            } from '../../modules/local/assign_features/main'
include { CREATE_MATRIX              } from '../../modules/local/create_matrix/main'

workflow QUANTIFICATION_EPI2ME {

    take:
    ch_tagged_bam   // channel: [ meta, bam ]
    ch_quik_r2      // channel: [ meta, reads ]
    ch_gtf          // channel: path(gtf)
    ch_genome_fasta // path(fasta)

    main:
    ch_versions = channel.empty()

    FILTER_PRIMARY_BAM (
        ch_tagged_bam
    )
    ch_versions = ch_versions.mix(FILTER_PRIMARY_BAM.out.versions.first())

    ch_stringtie_bam = FILTER_PRIMARY_BAM.out.bam.map { meta, bam ->
        [ meta + [ strandedness: params.strandedness ], bam ]
    }

    STRINGTIE_STRINGTIE (
        ch_stringtie_bam,
        ch_gtf
    )

    GFFREAD (
        STRINGTIE_STRINGTIE.out.transcript_gtf,
        ch_genome_fasta
    )

    MINIMAP2_INDEX_TRANSCRIPTOME (
        GFFREAD.out.gffread_fasta
    )

    ch_txome_align = ch_quik_r2
        .map { meta, reads -> [ meta.id, meta, reads ] }
        .join(
            MINIMAP2_INDEX_TRANSCRIPTOME.out.index
                .map { meta, index -> [ meta.id, meta, index ] }
        )
        .map { sample_id, meta, reads, idx_meta, index ->
            [ meta, reads, idx_meta, index ]
        }

    MINIMAP2_ALIGN_TRANSCRIPTOME (
        ch_txome_align.map { meta, reads, idx_meta, index -> [ meta, reads ] },
        ch_txome_align.map { meta, reads, idx_meta, index -> [ idx_meta, index ] },
        true,   // bam_format
        "bai",  // bam_index_extension
        false,  // cigar_paf_format
        false   // cigar_bam
    )

    SAMTOOLS_FAIDX (
        Channel.value([ [ id: 'genome' ], ch_genome_fasta ])
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())

    ch_ref_gtf = ch_gtf.map { gtf -> [ [ id: 'genome' ], gtf ] }

    GFFCOMPARE (
        STRINGTIE_STRINGTIE.out.transcript_gtf
            .map { meta, gtf -> [ meta, [ gtf ] ] },
        SAMTOOLS_FAIDX.out.fasta.first(),
        ch_ref_gtf.first()
    )

    SORT_BAM_NAMES (
        MINIMAP2_ALIGN_TRANSCRIPTOME.out.bam
    )
    ch_versions = ch_versions.mix(SORT_BAM_NAMES.out.versions.first())

    EXTRACT_BAM_TAGS (
        ch_tagged_bam
    )
    ch_versions = ch_versions.mix(EXTRACT_BAM_TAGS.out.versions.first())

    ch_assign_input = SORT_BAM_NAMES.out.bam
        .map { meta, bam -> [ meta.id, meta, bam ] }
        .join(GFFCOMPARE.out.tmap.map { meta, tmap -> [ meta.id, meta, tmap ] })
        .join(EXTRACT_BAM_TAGS.out.mapq_tags.map { meta, tags -> [ meta.id, meta, tags ] })
        .map { sample_id, meta, bam, t2, tmap, t3, mapq_tags ->
            [ meta, bam, tmap, mapq_tags ]
        }

    ASSIGN_FEATURES (
        ch_assign_input.map { meta, bam, tmap, mapq_tags -> [ meta, bam ] },
        ch_assign_input.map { meta, bam, tmap, mapq_tags -> [ meta, tmap ] },
        ch_gtf,
        ch_assign_input.map { meta, bam, tmap, mapq_tags -> [ meta, mapq_tags ] }
    )
    ch_versions = ch_versions.mix(ASSIGN_FEATURES.out.versions.first())

    ch_matrix_input = EXTRACT_BAM_TAGS.out.barcode_tags
        .map { meta, tags -> [ meta.id, meta, tags ] }
        .join(ASSIGN_FEATURES.out.features.map { meta, features -> [ meta.id, meta, features ] })
        .map { sample_id, meta, barcode_tags, m2, features ->
            [ meta, barcode_tags, features ]
        }

    CREATE_MATRIX (
        ch_matrix_input.map { meta, barcode_tags, features -> [ meta, barcode_tags ] },
        ch_matrix_input.map { meta, barcode_tags, features -> [ meta, features ] }
    )
    ch_versions = ch_versions.mix(CREATE_MATRIX.out.versions.first())

    emit:
    gene_matrix        = CREATE_MATRIX.out.gene_matrix
    transcript_matrix  = CREATE_MATRIX.out.transcript_matrix
    stats              = CREATE_MATRIX.out.stats
    versions           = ch_versions
}
