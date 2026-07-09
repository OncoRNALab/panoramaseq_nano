process OARFISH {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::oarfish=0.9.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/oarfish:0.9.4--h7f5d12c_0' :
        'biocontainers/oarfish:0.9.4--h7f5d12c_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("features.tsv.gz") , emit: features
    tuple val(meta), path("barcodes.tsv.gz") , emit: barcodes
    tuple val(meta), path("matrix.mtx.gz")   , emit: mtx
    tuple val(meta), path("${meta.id}.meta_info.json"), emit: meta_info
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    oarfish \\
        --output ${prefix} \\
        --alignments ${bam} \\
        -j ${task.cpus} \\
        ${args}

    if [[ ! -f "${prefix}.features.txt" ]]; then
        echo "ERROR: Oarfish did not produce ${prefix}.features.txt (is --single-cell set?)" >&2
        exit 1
    fi

    cp ${prefix}.features.txt features.tsv
    cp ${prefix}.barcodes.txt barcodes.tsv

    grep '^%' ${prefix}.count.mtx > matrix.mtx
    grep -v '^%' ${prefix}.count.mtx | awk '{print \$2" "\$1" "\$3}' >> matrix.mtx

    gzip -f features.tsv barcodes.tsv matrix.mtx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oarfish: \$(oarfish --version | sed 's/oarfish //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch features.tsv.gz
    touch barcodes.tsv.gz
    touch matrix.mtx.gz
    touch ${prefix}.meta_info.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oarfish: 0.9.4
    END_VERSIONS
    """
}
