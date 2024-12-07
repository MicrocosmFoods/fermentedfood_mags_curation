library(tidyverse)

## Organizing metadata for BacDive isolates from fermented foods that have Genbank assemblies
# bacdive accessions
bacdive_accessions <- read_tsv("metadata/raw_metadata/bacdive_genomes/2024-10-08-parse-bacdive-accessions.tsv")
# bacdive fermented food records
bacdive_ff_records <- read_tsv("metadata/raw_metadata/bacdive_genomes/BacDive_fermented_food_filtered_list.tsv")

# join together
bacdive_ff_accessions <- left_join(bacdive_ff_records, bacdive_accessions) %>% 
  drop_na(GCA) %>% 
  select(GCA, Species, 'Isolation source', 'Category 3')

colnames(bacdive_ff_accessions) <- c("genbank_id", "isolate", "isolation_source", "broad_category")

bacdive_ff_metadata <- bacdive_ff_accessions %>% 
  mutate(genbank_accession = paste(genbank_id, ".1", sep=""))

write_tsv(bacdive_ff_metadata, "metadata/cleaned_metadata/bacdive_genomes/2024-10-08-bacdive-ff-isolates-metadata.tsv")

## join with CheckM stats downloaded from NCBI at https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/CheckM_report_prokaryotes.txt
genbank_checkm_metadata <- read_tsv("metadata/raw_metadata/bacdive_genomes/CheckM_report_prokaryotes.txt") %>% 
  select("genbank-accession", "checkm-completeness", "checkm-contamination", "species-name")

colnames(genbank_checkm_metadata) <- c("genbank_accession", "completeness", "contamination", "ncbi_species_name")

bacdive_checkm_metadata <- left_join(bacdive_ff_metadata, genbank_checkm_metadata)

missing_bacdive_records <- bacdive_checkm_metadata[is.na(bacdive_checkm_metadata$completeness), ] %>% 
  select(-genbank_accession) %>% 
  mutate(genbank_accession = paste(genbank_id, ".2", sep="")) %>% 
  select(-completeness, -contamination, -ncbi_species_name)

missing_bacdive_records_genbankids <- missing_bacdive_records %>% 
  pull(genbank_accession)

missing_bacdive_genbank_checkm_metadata <- left_join(missing_bacdive_records, genbank_checkm_metadata)

write.table(missing_bacdive_records_genbankids, "metadata/raw_metadata/bacdive_genomes/missing-bacdive-genbank-records.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

complete_bacdive_genbank_checkm_metadata <- bind_rows(bacdive_checkm_metadata, missing_bacdive_genbank_checkm_metadata) %>% 
  drop_na(completeness)

## quast stats
bacdive_combined_quast_stats <- read_tsv("metadata/raw_metadata/all_quast_stats.tsv", col_names = c("mag_id", "contigs_0", "contigs_1000", "contigs_5000", "contigs_10000", "contigs_25000", "contigs50000", "total_length_0", "total_length_1000", "total_length_5000", "total_length_10000", "total_length_25000", "total_length_50000", "contigs", "largest_contig", "total_length", "gc", "n50", "n90", "auN", "L50", "L90", "n_per_100kbp")) %>% 
  select(mag_id, contigs, total_length, gc, n50) %>% 
  filter(grepl("^GCA", mag_id)) %>%
  mutate(genbank_id = sub("\\..*$", "", mag_id))

complete_bacdive_genbank_metadata <- left_join(complete_bacdive_genbank_checkm_metadata, bacdive_combined_quast_stats) %>% 
  drop_na(contigs) %>% 
  drop_na(completeness) %>% 
  select(genbank_id, isolate, isolation_source, broad_category, completeness, contamination, contigs, total_length, gc, n50)

complete_bacdive_genbank_metadata %>% 
  filter(completeness > 90) %>% 
  filter(contigs < 50) %>% 
  count() # ~ 100 isolates of good quality

hq_bacdive_genbank_metadata <- complete_bacdive_genbank_metadata %>% 
  filter(completeness > 90) %>% 
  filter(contigs < 50)

# write intermediate metadata to manually curate and match the mag metadata categories
write_tsv(complete_bacdive_genbank_metadata, "metadata/raw_metadata/manually_curated_metadata/2024-11-04-complete-bacdive-ff-genbank-metadata.tsv")

# manually curated metadata for fermented food and hierarchical categories
curated_bacdive_metadata <- read_tsv("metadata/raw_metadata/manually_curated_metadata/2024-12-04-manually-curated-bacdive-metadata.tsv") %>% 
  select(genbank_id, sample_description,fermented_food,	specific_substrate, substrate_category,	category)
bacdive_genbank_accession_ids <- read_tsv("metadata/raw_metadata/bacdive_genomes/bacdive-genbank-accessions.txt", col_names = c("genbank_accession"))
bacdive_genbank_accession_list <- bacdive_genbank_accession_ids %>% 
  mutate(genbank_accession = gsub(".fa", "", genbank_accession)) %>% 
  mutate(genbank_id = sub("\\..*$", "", genbank_accession))

# join with manually curated metadata and the correct genbank accession number with the ID
bacdive_full_metadata <- complete_bacdive_genbank_metadata %>% 
  left_join(curated_bacdive_metadata, by = "genbank_id") %>% 
  left_join(bacdive_genbank_accession_list) %>% 
  select(-genbank_id) %>% 
  select(genbank_accession, everything()) %>% 
  drop_na(fermented_food) 

bacdive_workflow_metadata <- bacdive_full_metadata %>% 
  mutate(mag_id = genbank_accession) %>% 
  select(mag_id, fermented_food) %>% 
  mutate(fermented_food = gsub("[^a-zA-Z0-9]", "_", fermented_food))

# write final curated sets of metadata
write_tsv(bacdive_full_metadata, "metadata/cleaned_metadata/2024-12-06-bacdive-accessions-curated-metadata.tsv")
write_tsv(bacdive_workflow_metadata, "metadata/cleaned_metadata/2024-12-06-bacdive-workflow-metadata.tsv")
