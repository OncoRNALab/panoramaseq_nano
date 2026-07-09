process EXTRACT_BAM_TAGS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_polars_pysam_zstandard_pruned:233385ddf666eb9c' :
        'community.wave.seqera.io/library/python_polars_pysam_zstandard_pruned:233385ddf666eb9c' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.mapq_tags.tsv.zst"),    emit: mapq_tags
    tuple val(meta), path("${prefix}.barcode_tags.tsv.zst"), emit: barcode_tags
    path  "versions.yml",                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    extract_bam_tags.py \\
        ${bam} \\
        --prefix ${prefix} \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: \$( python -c "import pysam; print(pysam.__version__)" )
        pandas: \$( python -c "import pandas; print(pandas.__version__)" )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mapq_tags.tsv.zst
    touch ${prefix}.barcode_tags.tsv.zst

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: 0.24.0
        pandas: 3.0.3
    END_VERSIONS
    """
}
