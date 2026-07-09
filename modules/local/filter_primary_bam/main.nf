process FILTER_PRIMARY_BAM {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pysam_samtools_python:d57d14fab94eb674' :
        'community.wave.seqera.io/library/pysam_samtools_python:d57d14fab94eb674' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.primary.bam"),     emit: bam
    tuple val(meta), path("${prefix}.primary.bam.bai"), emit: index
    path  "versions.yml",                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    def flags  = task.ext.view_flags ?: '0x904'
    """
    samtools view \\
        -@ ${task.cpus} \\
        -F ${flags} \\
        -b ${bam} \\
        | samtools sort -@ ${task.cpus} -o ${prefix}.primary.bam -

    samtools index ${prefix}.primary.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$( samtools --version 2>&1 | head -n1 | sed 's/^samtools //' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.primary.bam
    touch ${prefix}.primary.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.21
    END_VERSIONS
    """
}
