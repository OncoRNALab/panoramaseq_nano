process SORT_BAM_CB {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pysam_samtools_python:d57d14fab94eb674' :
        'community.wave.seqera.io/library/pysam_samtools_python:d57d14fab94eb674' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.cb_sorted.bam"), emit: bam
    path  "versions.yml",                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    samtools sort \\
        -@ ${task.cpus} \\
        -t CB \\
        -o ${prefix}.cb_sorted.bam \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$( samtools --version 2>&1 | head -n1 | sed 's/^samtools //' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cb_sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.21
    END_VERSIONS
    """
}
