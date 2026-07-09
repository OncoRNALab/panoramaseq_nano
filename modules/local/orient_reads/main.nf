// DEPRECATED: replaced by modules/local/restrander. Retained for reference only;
// not included in the STARLIGHT workflow.
process ORIENT_READS {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/parasail-python_seqkit_vsearch_editdistance_python:7fd2e78b16bceadc' :
        'community.wave.seqera.io/library/parasail-python_seqkit_vsearch_editdistance_python:7fd2e78b16bceadc' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.oriented.fastq.gz"), emit: reads
    tuple val(meta), path("${prefix}.orient_summary.tsv"), emit: summary, optional: true
    path  "versions.yml",                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    def adapter_args = params.adapter_fasta
        ? "--adapters ${params.adapter_fasta}"
        : "--tso_sequence \"${params.tso_sequence}\" --rt_adapter \"${params.rt_adapter}\""
    """
    orient_reads.py \\
        ${reads} \\
        ${adapter_args} \\
        --min_adapter_id ${params.orient_min_adapter_id} \\
        --threads ${task.cpus} \\
        --output ${prefix}.oriented.fastq.gz \\
        --summary ${prefix}.orient_summary.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$( vsearch --version 2>&1 | head -1 | sed 's/vsearch v\\([^ ]*\\).*/\\1/' )
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
        python: \$( python --version 2>&1 | sed 's/Python //' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.oriented.fastq.gz
    touch ${prefix}.orient_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: 2.28.1
        seqkit: 2.8.2
        python: 3.11.0
    END_VERSIONS
    """
}
