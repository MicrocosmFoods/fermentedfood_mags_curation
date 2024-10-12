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
genome_metadata: ${params.genome_metadata}
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

metadata_ch = channel.fromPath(params.genome_metadata)

workflow {
    // get assembly stats with quast
    genome_quast_stats = quast_stats(genome_fastas)

    // run dRep on all genomes at 95% identity
    drep_results = dereplicate_genomes(genome_fasta_dir, metadata_ch)
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

process dereplicate_genomes {
    tag "dereplicate_genomes"
    publishDir "${params.outdir}/drep", mode: 'copy', pattern:"*.csv"

    conda "envs/drep.yml"

    input:
    path(genome_fasta_dir)
    path(genome_metadata)

    output:
    path("*"), emit: drep_results

    script:
    """
    dRep dereplicate drep_results -g ${genome_fasta_dir}/*.fa -sa 0.95 -p ${params.threads} --genomeInfo ${genome_metadata}
    """
}