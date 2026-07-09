process OARFISH_AGGREGATE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pysam_umi_tools_editdistance_h5py_pruned:105616b830f5346a' :
        'community.wave.seqera.io/library/pysam_umi_tools_editdistance_h5py_pruned:105616b830f5346a' }"

    input:
    tuple val(meta), path(features), path(barcodes), path(mtx), path(gtf)

    output:
    tuple val(meta), path("${prefix}_transcript_bc_matrix"), emit: transcript_matrix
    tuple val(meta), path("${prefix}_gene_bc_matrix"),        emit: gene_matrix
    tuple val(meta), path("${prefix}.aggregation_stats.json"), emit: stats
    path  "versions.yml",                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    oarfish_aggregate_genes.py \\
        --features ${features} \\
        --barcodes ${barcodes} \\
        --matrix ${mtx} \\
        --gtf ${gtf} \\
        --output-dir . \\
        --prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed 's/Python //' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_transcript_bc_matrix ${prefix}_gene_bc_matrix
    touch ${prefix}_transcript_bc_matrix/features.tsv.gz
    touch ${prefix}_transcript_bc_matrix/barcodes.tsv.gz
    touch ${prefix}_transcript_bc_matrix/matrix.mtx.gz
    touch ${prefix}_gene_bc_matrix/features.tsv.gz
    touch ${prefix}_gene_bc_matrix/barcodes.tsv.gz
    touch ${prefix}_gene_bc_matrix/matrix.mtx.gz
    echo '{}' > ${prefix}.aggregation_stats.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}
