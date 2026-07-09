process MINIMAP2_ALIGN_BAM {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/37671219cfd244eb9b33db9345d3543ffd83037419a1c57f4648aace493ec2c2/data' :
        'community.wave.seqera.io/library/minimap2_samtools:b09096fc890429ce' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("${prefix}.txome.bam"),     emit: bam
    tuple val(meta), path("${prefix}.txome.bam.bai"), emit: index
    path  "versions.yml",                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def tags   = task.ext.copy_tags ?: 'CB,UB,UR,CR,CY,UY'
    prefix     = task.ext.prefix ?: "${meta.id}"
    """
    samtools fastq \\
        -T ${tags} \\
        ${bam} \\
    | minimap2 \\
        ${args} \\
        -t ${task.cpus} \\
        ${index} \\
        - \\
    | samtools sort -@ ${task.cpus} -o ${prefix}.txome.bam -

    samtools index ${prefix}.txome.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$( minimap2 --version )
        samtools: \$( samtools --version 2>&1 | head -n1 | sed 's/^samtools //' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txome.bam
    touch ${prefix}.txome.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: 2.28
        samtools: 1.21
    END_VERSIONS
    """
}
