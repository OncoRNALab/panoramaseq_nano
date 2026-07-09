process CREATE_MATRIX {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pysam_umi_tools_editdistance_h5py_pruned:105616b830f5346a' :
        'community.wave.seqera.io/library/pysam_umi_tools_editdistance_h5py_pruned:105616b830f5346a' }"

    input:
    tuple val(meta), path(barcode_tags)
    tuple val(meta2), path(features)

    output:
    tuple val(meta), path("${prefix}_gene_bc_matrix"),        emit: gene_matrix, optional: true
    tuple val(meta), path("${prefix}_transcript_bc_matrix"), emit: transcript_matrix, optional: true
    tuple val(meta), path("${prefix}.matrix_stats.json"),    emit: stats, optional: true
    path  "versions.yml",                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p hdf_chunks

    workflow-glue create_matrix \\
        ${params.matrix_chrom} \\
        ${barcode_tags} \\
        ${features} \\
        --tsv_out ${prefix}.tags.tsv.zst \\
        --sa_tags_out ${prefix}.sa_tags.tsv.zst \\
        --hdf_out hdf_chunks \\
        --stats ${prefix}.matrix_stats.json \\
        --chunk_size ${params.matrix_chunk_size} \\
        --ref_interval ${params.matrix_ref_interval} \\
        ${params.matrix_umi_length ? "--umi_length ${params.matrix_umi_length}" : ''} \\
        ${params.matrix_skip_umi_clustering ? '--skip_umi_clustering' : ''} \\
        ${args}

    aggregate_matrix.py \\
        hdf_chunks \\
        . \\
        --prefix ${prefix} \\
        --features gene transcript

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        workflow_glue: 0.0.1
        umi_tools: \$( python -c "import umi_tools; print(umi_tools.__version__)" )
        scipy: \$( python -c "import scipy; print(scipy.__version__)" )
        h5py: \$( python -c "import h5py; print(h5py.__version__)" )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_gene_bc_matrix ${prefix}_transcript_bc_matrix
    touch ${prefix}_gene_bc_matrix/barcodes.tsv.gz
    touch ${prefix}_gene_bc_matrix/features.tsv.gz
    touch ${prefix}_gene_bc_matrix/matrix.mtx.gz
    touch ${prefix}_transcript_bc_matrix/barcodes.tsv.gz
    touch ${prefix}_transcript_bc_matrix/features.tsv.gz
    touch ${prefix}_transcript_bc_matrix/matrix.mtx.gz
    touch ${prefix}.matrix_stats.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        workflow_glue: 0.0.1
        umi_tools: 1.1.6
        scipy: 1.18.0
        h5py: 3.16.0
    END_VERSIONS
    """
}
