process ASSIGN_FEATURES {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_polars_pysam_zstandard_pruned:233385ddf666eb9c' :
        'community.wave.seqera.io/library/python_polars_pysam_zstandard_pruned:233385ddf666eb9c' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(tmap)
    path gtf
    tuple val(meta3), path(mapq_tags)

    output:
    tuple val(meta), path("${prefix}.feature_assigns.tsv.zst"), emit: features
    path  "versions.yml",                                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    workflow-glue assign_features \\
        ${bam} \\
        ${tmap} \\
        ${gtf} \\
        ${mapq_tags} \\
        ${prefix}.feature_assigns.tsv.zst \\
        --min_mapq ${params.gene_assign_min_mapq} \\
        --min_tr_coverage ${params.gene_assign_min_tr_coverage} \\
        --min_read_coverage ${params.gene_assign_min_read_coverage} \\
        --chunksize ${params.gene_assign_chunksize} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        workflow_glue: 0.0.1
        pandas: \$( python -c "import pandas; print(pandas.__version__)" )
        polars: \$( python -c "import polars; print(polars.__version__)" )
        pyarrow: \$( python -c "import pyarrow; print(pyarrow.__version__)" )
        pysam: \$( python -c "import pysam; print(pysam.__version__)" )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.feature_assigns.tsv.zst

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        workflow_glue: 0.0.1
        pandas: 3.0.3
        polars: 1.42.0
        pyarrow: 18.0.0
        pysam: 0.24.0
    END_VERSIONS
    """
}
