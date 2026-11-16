#!/usr/bin/env nextflow


workflow {
    Channel
        .fromList(params.sample_ids)
        .map { sample_id ->
            sample_lanes = params.samples[sample_id] // multiple lanes per sample
            tuple(sample_id, sample_lanes) // [[sample1, [lane1, lane2, lane3, ...]], [...]]
        }
        .set { sample_lanes_ch }

    CAT_FASTQS(sample_lanes_ch)
    SEQKIT_STATS(
        CAT_FASTQS.out.sample_fastq
            .map { it[1] }
            .collect()
    )
    KRAKEN(
        params.kraken_db,
        CAT_FASTQS.out.sample_fastq
    )
    BRACKEN(
        params.bracken_db,
        KRAKEN.out.sample_report
    )
    KREPORT2KRONA_KREPORT(KRAKEN.out.sample_report)
    KREPORT2KRONA_BREPORT(BRACKEN.out.sample_report)
    KTIMPORTTEXT_KKRONA(KREPORT2KRONA_KREPORT.out.sample_krona)
    KTIMPORTTEXT_BKRONA(KREPORT2KRONA_BREPORT.out.sample_krona)
    MULTIQC(
        KRAKEN.out.sample_report
            .map { it[1] }
            .collect(),
        BRACKEN.out.sample_report
            .map { it[1] }
            .collect()
    )
}

// merge FASTQ files by sample
process CAT_FASTQS {
    tag "${id}"
    // publishDir "${params.output_dir}/${params.project_id}/${id}", mode: 'copy'
    debug true
    cpus 1
    time { 10.min * task.attempt }
    memory { 1.GB * task.attempt }

    input:
    tuple val(id), path(fastqs)

    output:
    tuple val(id), path("${id}_merged.fastq.gz"), emit: sample_fastq

    script:
    """
    cat ${fastqs} > ${id}_merged.fastq.gz
    """
}

// get FASTQ statistics
process SEQKIT_STATS {
    tag "${params.project_id}"
    publishDir "${params.output_dir}/${params.project_id}", mode: 'copy'
    cpus 8
    time { 10.min * task.attempt }
    memory { 6.GB * task.attempt }
    container "${params.seqkit_docker}"
    scratch true // I/O intensive, run on SSD (/tmp)
    stageInMode 'copy'

    input:
    path fastqs

    output:
    path "${params.project_id}_seqkit_stats.tsv"

    script:
    """
    seqkit stats \
        --all \
        --tabular \
        --out-file ${params.project_id}_seqkit_stats.tsv \
        --threads ${task.cpus * 2} \
        ${fastqs}
    """
}

// classify reads by taxons
process KRAKEN {
    tag "${id}"
    publishDir "${params.output_dir}/${params.project_id}/${id}", mode: 'copy'
    cpus 8
    time { 180.min * task.attempt }
    memory '48.GB'
    container "${params.kraken_docker}"
    scratch true // I/O intensive, run on SSD (/tmp)
    stageInMode 'copy'
    maxForks 2 // /tmp has limited capacity, run a few instances at a time

    input:
    path 'database/*'
    tuple val(id), path(fastq)

    output:
    tuple val(id), path("${id}.kraken"), emit: sample_kraken
    tuple val(id), path("${id}.kreport"), emit: sample_report

    script:
    """
    kraken2 \
        --db database \
        --threads ${task.cpus * 2} \
        --report ${id}.kreport \
        --output ${id}.kraken \
        --minimum-hit-groups 3 \
        --confidence 0.1 \
        --report-minimizer-data \
        --memory-mapping \
        ${fastq}
    """
}

process BRACKEN {
    tag "${id}"
    publishDir "${params.output_dir}/${params.project_id}/${id}", mode: 'copy'
    cpus 1
    time { 1.min * task.attempt }
    memory { 1.GB * task.attempt }
    container "${params.bracken_docker}"

    input:
    path 'database/*'
    tuple val(id), path(report)

    output:
    tuple val(id), path("${id}.bracken"), emit: sample_bracken
    tuple val(id), path("${id}.breport"), emit: sample_report

    script:
    """
    bracken \
        -d database \
        -i ${report} \
        -o ${id}.bracken \
        -w ${id}.breport \
        -r 75 \
        -l S \
        -t 20
    """
}

process KREPORT2KRONA_KREPORT {
    tag "${id}"
    publishDir "${params.output_dir}/${params.project_id}/${id}", mode: 'copy'
    cpus 1
    time { 30.min * task.attempt }
    memory { 1.GB * task.attempt }

    input:
    tuple val(id), path(report)

    output:
    tuple val(id), path("${id}.kkrona"), emit: sample_krona

    script:
    """
    kreport2krona.py \
        --report ${report} \
        --output ${id}.kkrona
    """
}

process KREPORT2KRONA_BREPORT {
    tag "${id}"
    publishDir "${params.output_dir}/${params.project_id}/${id}", mode: 'copy'
    cpus 1
    time { 1.min * task.attempt }
    memory { 1.GB * task.attempt }

    input:
    tuple val(id), path(report)

    output:
    tuple val(id), path("${id}.bkrona"), emit: sample_krona

    script:
    """
    kreport2krona.py \
        --report ${report} \
        --output ${id}.bkrona
    """
}

process KTIMPORTTEXT_KKRONA {
    tag "${id}"
    publishDir "${params.output_dir}/${params.project_id}/${id}", mode: 'copy'
    cpus 1
    time { 1.min * task.attempt }
    memory { 1.GB * task.attempt }
    container "${params.krona_docker}"

    input:
    tuple val(id), path(krona)

    output:
    tuple val(id), path("${id}_kkrona.html")

    script:
    """
    ktImportText \
        -o ${id}_kkrona.html \
        ${krona}
    """
}

process KTIMPORTTEXT_BKRONA {
    tag "${id}"
    publishDir "${params.output_dir}/${params.project_id}/${id}", mode: 'copy'
    cpus 1
    time { 1.min * task.attempt }
    memory { 1.GB * task.attempt }
    container "${params.krona_docker}"

    input:
    tuple val(id), path(krona)

    output:
    tuple val(id), path("${id}_bkrona.html")

    script:
    """
    ktImportText \
        -o ${id}_bkrona.html \
        ${krona}
    """
}

process MULTIQC {
    tag "${params.project_id}"
    publishDir "${params.output_dir}/${params.project_id}", mode: 'copy'
    cpus 1
    time { 10.min * task.attempt }
    memory { 6.GB * task.attempt }
    container "${params.multiqc_docker}"

    input:
    path kreport
    path breport

    output:
    path "${params.project_id}.html"
    path "${params.project_id}_data/*"

    script:
    """
    multiqc \
        --filename ${params.project_id} \
        .
    """
}
