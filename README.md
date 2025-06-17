# Genomes from Microbes of Fermented Foods
This repository documents curating genomes and metadata from publicly available studies of microbes from different fermented foods. 

The full database of ~13,500 microbial genomes and associated curated metadata can be accessed on [Zenodo](https://zenodo.org/records/15588452). We have also made a subset of these genomes available as a Narrative on KBase. We clustered the full set of 13,500 genomes at 99% average nucleotide identity (ANI) to obtain ~4,300 "strain"-representative genomes. You can access the static KBase Narrative [here](https://kbase.us/n/218406/47/). To access the KBase platform to explore the database and run your own analyses, you will need to create a KBase account. 

<img src="https://github.com/MicrocosmFoods/fermentedfood_metadata_curation/blob/main/figures/kbase_narrative/metadata_summary.png?raw=true">

The most up-to-date corresponding metadata is available [here](https://raw.githubusercontent.com/MicrocosmFoods/fermentedfood_metadata_curation/refs/heads/main/data/2025-05-21-genome-metadata-food-taxonomy.tsv). 

## Accessed Datasets and Repositories
For curating the set of microbial genomes from diverse fermented foods, we accessed metagenome-assembled genomes (MAGs) and isolates from publicly available sources. 

- Carlino et al. 2024 [Unexplored microbial diversity from 2,500 food metagenomes and links with the human microbiome](https://www.cell.com/cell/fulltext/S0092-8674(24)00833-X), data downloaded from [https://zenodo.org/doi/10.5281/zenodo.10891046](https://zenodo.org/doi/10.5281/zenodo.10891046)
- Carlino et al. 2024 cFMD v1.2.1 [Unexplored microbial diversity from 2,500 food metagenomes and links with the human microbiome](https://www.cell.com/cell/fulltext/S0092-8674(24)00833-X) second release of food MAGs, data downloaded from [https://zenodo.org/records/14906604](https://zenodo.org/records/14906604), excluding MAGs from Saak et al. 2023 since we already incorporate those separately
- Caffrey et al. 2024 [MiFoDB, a workflow for microbial food metagenomic characterization, enables high-resolution analysis of fermented food microbial dynamics](https://www.biorxiv.org/content/10.1101/2024.03.29.587370v1.full), data downloaded from [https://zenodo.org/records/13830159](https://zenodo.org/records/13830159) for MiFoDB_beta_v3
- Rappaport et al. 2024 [Genomics and synthetic community experiments uncover the key metabolic roles of acetic acid bacteria in sourdough starter microbiomes](https://journals.asm.org/doi/10.1128/msystems.00537-24), data available at [PRJNA589612](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA589612/) and metadata accessed from authors
- Du et al. 2023 [Metagenomics reveals the habitat specificity of biosynthetic potential of secondary metabolites in global food fermentations](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-023-01536-8), data downloaded from [Github](https://github.com/durubing-jn/food-fermentation-mategenome)
- Saak et al. 2023 [Longitudinal, Multi-Platform Metagenomics Yields a High-Quality Genomic Catalog and Guides an In Vitro Model for Cheese Communities](https://journals.asm.org/doi/10.1128/msystems.00701-22), data downloaded from [Dryad](https://datadryad.org/stash/dataset/doi:10.5061/dryad.bg79cnpd8)
- BacDive Genbank accessions from fermented foods. [BacDive](https://bacdive.dsmz.de/) accessions from fermented foods were collected in `metadata/raw_metadata/isolate_genomes/BacDive_fermented_food_filtered_list.tsv` and matched to available Genbank assembly records in `metadata/raw_metadata/isolate_genomes/bacdive/2024-10-08-parse-bacdive-accessions.tsv`. Genomes were downloaded with `ncbi-genome-download` with: 
    ```
    ncbi-genome-download --section genbank \\
    --assembly-accessions metadata/raw_metadata/isolate_genomes/bacdive/bacdive-accessions-download.txt \\
    -m bacdive-metadata.txt \\
    --format "fasta" \\
    -p 3 \\
    bacteria
    ```

## Environment Setup
After [installing conda for your OS](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html), you can create a conda environment with all the dependencies required for running the scripts with: 
```
conda env create -n fermented_foods envs/dev.yml
```

## Metadata Curation 
Metadata associated with each genome including sample accession, food information and taxonomy, and further curation of the genome set including dereplication and GTDB-tk taxonomic assignment is documented in the [fermentedfood_metadata_curation](https://github.com/MicrocosmFoods/fermentedfood_metadata_curation) repository. 

## Repository Structure & Files
The repository is split both for scripts and directories for handling genomes from MAG datasets or collections of isolates. The cleaned, curated metadata for the MAG datasets and bacdive isolates is in the main `metadata` directory and copied in the subdirectories. The subdirectories contain the raw files for curating metadata from different sources together.

``` 
- metadata/ - Most of this are now intermediate files used to create the final metadata files in the [fermentedfood_metadata_curation](https://github.com/MicrocosmFoods/fermentedfood_metadata_curation) repository.
- scripts/
    - batch_fasta_files.py - Helper script to create batches of fasta files and associated samplesheet templates for uploading to KBase.
- envs/
 - dev.yml
 - quast.yml
```