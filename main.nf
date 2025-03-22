#! /usr/bin/env nextflow

// Description
// MAG & genome quality checks, dereplication

nextflow.enable.dsl = 2

params.threads = 5
params.outdir = './genomequal_results'

log.info """\

CHECK GENOME ASSEMBLY QUALITY AND DEREPLICATE SETS OF GENOMES
=============================================================
input_genomes: ${params.input_genomes}
outdir: ${params.outdir}
threads: ${params.threads}
"""

// Define channels
genome_fasta_dir = Channel.fromPath(params.input_genomes)
genome_fastas = Channel.fromPath("${params.input_genomes}/*.fa")
    .map { file ->
        def baseName = file.getBaseName()
        return [file, baseName]
    }

workflow {
    // get assembly stats with quast
    genome_quast_stats = quast_stats(genome_fastas)

}

process quast_stats {
    tag "${genome_name}_quast_stats"
    publishDir "${params.outdir}/quast", mode: 'move', pattern:"*.tsv"

    conda "envs/quast.yml"

    input: 
    tuple path(fasta), val(genome_name)

    output:
    path("*stats.tsv"), emit: quast_tsv

    script:
    """
    quast *.fa --output-dir ${genome_name} -t 1
    mv ${genome_name}/transposed_report.tsv ${genome_name}.transposed_report.tsv
    tail -n +2 ${genome_name}.transposed_report.tsv > ${genome_name}.stats.tsv
    """
}