process SAMTOOLS_INDEX {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pysam_samtools_python:d57d14fab94eb674' :
        'community.wave.seqera.io/library/pysam_samtools_python:d57d14fab94eb674' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), path("${bam.name}.bai"), emit: bam_index
    path  "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools index -@ ${task.cpus} ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$( samtools --version 2>&1 | head -n1 | sed 's/^samtools //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${bam.name}.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.21
    END_VERSIONS
    """
}
