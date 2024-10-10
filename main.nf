#! /usr/bin/env nextflow

// Description
// MAG & genome quality checks, dereplication

nexflow.enable.dsl = 2

params.threads = 8
params.outdir = './results'
params.dbdir = "databases"

log.info """\

CHECK GENOME ASSEMBLY QUALITY AND DEREPLICATE SETS OF GENOMES
=============================================================
input_mags: ${params.input_mags}
mag_metadata: ${params.mag_metadata}
input_isolates: ${params.input_isolates}
dbdir: ${params.dbdir}
outdir: ${params.outdir}
threads: ${params.threads}
"""

// Define channels
mag_fastas = Channel.fromPath("${params.input_mags}/*.fa")
    .map { file ->
        def baseName = file.getBaseName()
        return [file, baseName]
    }

isolate_fastas = Channel.fromPath("${params.input_isolates}/*.fa")
    .map { file ->
        def baseName = file.getBaseName()
        return [file, baseName]
    }

metadata_ch = channel.fromPath(params.metadata)
database_directory = channel.fromPath(params.dbdir)

workflow {
    // get assembly stats with quast
    mag_quast_stats = quast_stats(mag_fastas)
    isolate_quast_stats = quast_stats(isolate_fastas)

    // run checkM on bacDive isolates
    download_checkm_db(database_directory)
    checkm_db = download_checkm_db.out.checkm2_db
    run_checkm(isolate_fastas.map{ it[0] }.collect())

    // combine quality stats for all genomes

    // run dRep on all genomes at 95% identity
}

process quast_stats {
    tag "${genome_name}_quast_stats"
    publishDir "${params.outdir}/quast", mode: 'copy', pattern:"*.tsv"

    conda "envs/quast.yml"

    input: 
    tuple path(fasta), val(genome_name)

    output:
    path("*.tsv"), emit: quast_tsv

    script:
    """
    quast *.fa --output-dir ${genome_name} -t 1
    ln -s ${genome_name}/${genome_name}.report.tsv
    ln -s ${genome_name}/${genome_name}.transposed_report.tsv
    """
}

process download_checkm_db {
    tag "download_checkm_db"
    
    conda "envs/checkm.yml"

    input:
    path(database_directory)

    output: 
    path("*"), emit: checkm2_db

    script:
    """
    checkm2 database --download --path ${database_directory}
    """
}

process run_checkm {
    tag "run_checkm"

    conda "envs/checkm.yml"

    input: 
    path(fasta_directory)
    path(database_directory)

    output:
    path("*"), emit: checkm_results

    script:
    """
    checkm2 predict -i ${fasta_directory}*.fa -o checkm_results --database_path ${database_directory}
    """
}
