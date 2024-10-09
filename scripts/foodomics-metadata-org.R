library(tidyverse)

# metadata sheets
metabolomics_metadata <- read.csv("metabolomics/GlobalFoodOmics/foodomics_metabolomics_sample_metadata.csv")
amplicon_metadata <- read_tsv("metabolomics/GlobalFoodOmics/foodomics_16S_amplicon_runs_metadata.tsv")

# select metabolomics data for fermented foods
fermented_foods_metabolomics <- metabolomics_metadata %>% 
  select(filename, sample_name, description, dna_extracted, fermented, fermentation_day) %>% 
  filter(fermented == "yes")

fermented_foods_dna <- fermented_foods_metabolomics %>% 
  filter(dna_extracted == TRUE)

# join with ENA amplicon metadata
fermented_foods_amplicon_metabolomic_metadata <- amplicon_metadata %>% 
  mutate(sample_name = sample_title) %>% 
  select(sample_accession, experiment_accession, read_count, fastq_ftp, sample_name) %>% 
  left_join(fermented_foods_dna) %>% 
  filter(read_count > 5000) %>% 
  drop_na()

write_csv(fermented_foods_amplicon_metabolomic_metadata, "cleaned_metadata/metabolomics/2024-09-17-fermented-foods-amplicon-metabolomics-metadata.csv")
