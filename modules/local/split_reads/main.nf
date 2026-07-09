process SPLIT_READS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/parasail-python_seqkit_vsearch_editdistance_python:7fd2e78b16bceadc' :
        'community.wave.seqera.io/library/parasail-python_seqkit_vsearch_editdistance_python:7fd2e78b16bceadc' }"

    input:
    tuple val(meta), path(oriented_reads)
    tuple val(meta2), path(bc_tags)

    output:
    tuple val(meta), path("${prefix}_R1.fastq.gz"), emit: r1
    tuple val(meta), path("${prefix}_R2.fastq.gz"), emit: r2
    path  "versions.yml",                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    split_reads.py \\
        ${oriented_reads} \\
        ${bc_tags} \\
        --barcode_length ${params.barcode_length} \\
        --min_cdna_len ${params.min_cdna_len} \\
        --prefix ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version 2>&1 | sed 's/Python //' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_R1.fastq.gz
    touch ${prefix}_R2.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}
