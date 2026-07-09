process QUIK_STARSOLO {
    tag "${meta.id}"
    label 'use_gpu'
    // QUIK is GPU-only and requires a Singularity/Apptainer container.
    // Run with -profile singularity (or apptainer) and the gpu profile to
    // pass --nv (CUDA) automatically via containerOptions in modules.config.
    conda null
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/francoaps/quik-runtime-compile:runtime-compile-v1' :
        'quay.io/francoaps/quik-runtime-compile:runtime-compile-v1' }"

    input:
    tuple val(meta), path(r1), path(r2)
    path barcode_file
    val  barcode_length

    output:
    tuple val(meta), path("*_R1_filtered.fastq.gz"),                          emit: r1
    tuple val(meta), path("*_R2_filtered.fastq.gz"),                          emit: r2
    tuple val(meta), path("*_whitelist.txt"),                                  emit: whitelist
    tuple val(meta), path("*_barcode_calling_stats.txt"),                      emit: stats
    tuple val(meta), path("*_R1_rejected.fastq.gz"), optional: true,           emit: r1_rejected
    tuple val(meta), path("*_R2_rejected.fastq.gz"), optional: true,           emit: r2_rejected
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def strategy          = params.strategy
    def distance_measure  = params.distance_measure
    def rejection_threshold = params.rejection_threshold ?: (barcode_length * 0.25).toInteger()

    """
    echo "=== QUIK Runtime Compilation for ONT Workflow ==="
    echo "Parameters: barcode_length=${barcode_length}, rejection_threshold=${rejection_threshold}"
    echo "Strategy: ${strategy}, Distance: ${distance_measure}, Barcode start: 0 (fixed after split)"
    echo "====================================="

    # Step 1: Compile QUIK from source inside the container.
    # Nextflow mounts the task work directory so build artefacts persist on the host.
    mkdir -p quik_build
    cmake /opt/quik \
        -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=\${PWD}/quik_build \
        -B \${PWD}/quik_build
    make -j${task.cpus} -C \${PWD}/quik_build
    QUIK_EXEC="\${PWD}/quik_build/single_strategy_benchmark_fastq_paired"
    echo "QUIK compiled successfully: \${QUIK_EXEC}"

    # Step 2: Decompress input FASTQ files
    echo "Decompressing input FASTQ files..."
    gunzip -c ${r1} > input_R1.fastq
    gunzip -c ${r2} > input_R2.fastq

    # Step 3: Extract barcode sequences from the whitelist CSV
    echo "Extracting barcode sequences..."
    tail -n +2 ${barcode_file} | cut -d',' -f1 > barcodes_only.txt

    # Step 4: Run QUIK barcode calling
    \${QUIK_EXEC} \\
        barcodes_only.txt \\
        input_R1.fastq \\
        input_R2.fastq \\
        0 \\
        ${barcode_length} \\
        ${strategy} \\
        ${distance_measure} \\
        ${rejection_threshold} \\
        ${prefix}_R1_filtered.fastq \\
        ${prefix}_R2_filtered.fastq \\
        > ${prefix}_barcode_calling_stats.txt 2>&1

    # Step 5: Whitelist from reference barcodes
    echo "Creating whitelist from reference barcodes..."
    cp barcodes_only.txt ${prefix}_whitelist.txt

    # Step 6: Validate whitelist
    whitelist_count=\$(wc -l < ${prefix}_whitelist.txt)
    echo "Whitelist contains \${whitelist_count} reference barcodes"

    # Step 7: Compress output FASTQ files
    gzip ${prefix}_R1_filtered.fastq
    gzip ${prefix}_R2_filtered.fastq

    # Step 8: Extract rejected reads (optional, for Columba barcode rescue)
    if [ "${params.enable_columba_rescue}" = "true" ]; then
        echo "Extracting rejected reads for Columba barcode rescue..."

        echo "Extracting filtered read IDs..."
        zcat ${prefix}_R1_filtered.fastq.gz | awk 'NR % 4 == 1 {
            read_id = substr(\$1, 2)
            split(read_id, parts, "_calledidx_")
            print parts[1]
        }' > filtered_ids.txt

        filtered_count=\$(wc -l < filtered_ids.txt)
        echo "Filtered read count: \${filtered_count}"

        echo "Extracting rejected R1 reads..."
        awk 'NR==FNR{ids[\$1]=1; next}
             FNR % 4 == 1 {
                 hdr = \$0
                 read_id = substr(\$1, 2)
                 split(read_id, parts, " ")
                 id = parts[1]
             }
             FNR % 4 == 2 {seq = \$0}
             FNR % 4 == 3 {plus = \$0}
             FNR % 4 == 0 {
                 qual = \$0
                 if (!(id in ids)) {
                     print hdr"\\n"seq"\\n"plus"\\n"qual
                 }
             }' filtered_ids.txt input_R1.fastq | gzip > ${prefix}_R1_rejected.fastq.gz

        echo "Extracting rejected R2 reads..."
        awk 'NR==FNR{ids[\$1]=1; next}
             FNR % 4 == 1 {
                 hdr = \$0
                 read_id = substr(\$1, 2)
                 split(read_id, parts, " ")
                 id = parts[1]
             }
             FNR % 4 == 2 {seq = \$0}
             FNR % 4 == 3 {plus = \$0}
             FNR % 4 == 0 {
                 qual = \$0
                 if (!(id in ids)) {
                     print hdr"\\n"seq"\\n"plus"\\n"qual
                 }
             }' filtered_ids.txt input_R2.fastq | gzip > ${prefix}_R2_rejected.fastq.gz

        rejected_count=\$(zcat ${prefix}_R1_rejected.fastq.gz | wc -l)
        rejected_reads=\$((rejected_count / 4))
        echo "Rejected R1 reads: \${rejected_reads}"
        rm filtered_ids.txt
    fi

    # Step 9: Clean up
    rm -f input_R1.fastq input_R2.fastq barcodes_only.txt
    rm -rf quik_build

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quik: \$(echo "2.0-runtime-flex")
        cuda: \$(nvcc --version 2>/dev/null | grep release | cut -d' ' -f6 | cut -d',' -f1 || echo "unknown")
        barcode_length: \$(echo "${barcode_length}")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_R1_filtered.fastq.gz
    touch ${prefix}_R2_filtered.fastq.gz
    touch ${prefix}_whitelist.txt
    touch ${prefix}_barcode_calling_stats.txt
    touch ${prefix}_R1_rejected.fastq.gz
    touch ${prefix}_R2_rejected.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quik: \$(echo "2.0-runtime-flex")
        barcode_length: \$(echo "${barcode_length}")
    END_VERSIONS
    """
}
