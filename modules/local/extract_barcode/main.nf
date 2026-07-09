process EXTRACT_BARCODE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/parasail-python_seqkit_vsearch_editdistance_python:7fd2e78b16bceadc' :
        'community.wave.seqera.io/library/parasail-python_seqkit_vsearch_editdistance_python:7fd2e78b16bceadc' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.bc_tags.tsv.gz"),    emit: tags
    tuple val(meta), path("${prefix}.barcode_start.txt"), emit: barcode_start
    path  "versions.yml",                                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    """
    extract_barcode.py \\
        ${reads} \\
        --rt_adapter "${params.rt_adapter}" \\
        --adapter_anchor "${params.adapter_anchor}" \\
        --barcode_length ${params.barcode_length} \\
        --polya_length ${params.barcode_polya_length} \\
        --window ${params.barcode_window} \\
        --umi_window ${params.umi_window} \\
        --gap_open ${params.barcode_gap_open} \\
        --gap_extend ${params.barcode_gap_extend} \\
        --max_anchor_ed ${params.barcode_max_anchor_ed} \\
        --max_adapter_ed ${params.max_adapter_ed} \\
        --min_barcode_qv ${params.barcode_min_qv} \\
        ${params.min_barcode_len != null ? "--min_barcode_len ${params.min_barcode_len}" : ''} \\
        --output ${prefix}.bc_tags.tsv.gz \\
        --consensus_output ${prefix}.barcode_start.txt \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parasail: \$( python -c "import parasail; print(parasail.__version__)" )
        editdistance: \$( python -c "import importlib.metadata; print(importlib.metadata.version('editdistance'))" )
        python: \$( python --version 2>&1 | sed 's/Python //' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bc_tags.tsv.gz
    echo "0" > ${prefix}.barcode_start.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parasail: 1.3.4
        editdistance: 0.8.1
        python: 3.11.0
    END_VERSIONS
    """
}
