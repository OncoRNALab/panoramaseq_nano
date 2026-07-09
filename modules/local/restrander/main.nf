process RESTRANDER {
    tag "${meta.id}"
    label 'process_medium'

    conda null
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://registry-1.docker.io/francops1722/restrander:v1.1.3-sif' :
        'francops1722/restrander:v1.1.3' }"

    input:
    tuple val(meta), path(reads)
    path config

    output:
    tuple val(meta), path("${prefix}.restranded.fastq.gz"), emit: reads
    tuple val(meta), path("${prefix}.stats.json"),         emit: stats
    path  "versions.yml",                                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    restrander \\
        ${reads} \\
        ${prefix}.restranded.fastq.gz \\
        ${config} \\
        > ${prefix}.stats.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        restrander: "1.1.3"
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.restranded.fastq.gz
    echo '{"stats":{"totalReads":0}}' > ${prefix}.stats.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        restrander: "1.1.3"
    END_VERSIONS
    """
}
