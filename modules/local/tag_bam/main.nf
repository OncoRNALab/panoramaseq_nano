process TAG_BAM {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pysam_samtools_python:d57d14fab94eb674' :
        'community.wave.seqera.io/library/pysam_samtools_python:d57d14fab94eb674' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(bc_tags)

    output:
    tuple val(meta), path("${prefix}.tagged.bam"),     emit: bam
    tuple val(meta), path("${prefix}.tagged.bam.bai"), emit: index
    path  "versions.yml",                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    """
    tag_bam.py \\
        ${bam} \\
        ${bc_tags} \\
        --prefix ${prefix} \\
        --barcode_length ${params.barcode_length} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: \$( python -c "import pysam; print(pysam.__version__)" )
        samtools: \$( samtools --version 2>&1 | head -n1 | sed 's/^samtools //' )
        python: \$( python --version 2>&1 | sed 's/Python //' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tagged.bam
    touch ${prefix}.tagged.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: 0.23.0
        samtools: 1.21
        python: 3.11.0
    END_VERSIONS
    """
}
