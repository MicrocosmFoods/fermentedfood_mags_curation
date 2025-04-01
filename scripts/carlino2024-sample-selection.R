library(tidyverse)

# all curated MAG metadata
curated_mag_metadata <- read_tsv("metadata/2025-03-24-all-ff-mag-metadata-cleaned-curated.tsv")

# all Carlino 2024 cFMD v1.0.0 sample metadata
carlino_sample_metadata <- read_tsv("metadata/raw_metadata/mag_datasets/CellFoodMetagenomics/2024-11-04-Carlino-sample-metadata.tsv") %>% 
  mutate(source = paste0(dataset_name, "__", sample_id))

# curated Carlino 2024 MAGs joined with sample metadata
carlino_mag_sample_info <- curated_mag_metadata %>% 
  filter(study_catalog == "Carlino2024") %>% 
  left_join(carlino_sample_metadata)

filtered_carlino_samples <- carlino_mag_sample_info %>% 
  filter(library_layout == "PAIRED") %>% 
  filter(sequencing_platform == "Illumina_NovaSeq_6000" | sequencing_platform == "Illumina_HiSeq_4000" | sequencing_platform == "Illumina_HiSeq_2500" | sequencing_platform == "NextSeq_500" | sequencing_platform == "Illumina_HiSeq_2000" | sequencing_platform == "Illumina_HiSeq_1500") %>% 
  filter(database_origin == "NCBI" | database_origin == "ENA") %>% 
  filter(n_of_reads > 10000000)

filtered_fermented_foods <- filtered_carlino_samples %>% 
  group_by(fermented_food) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  pull(fermented_food)

filtered_food_samples <- filtered_carlino_samples %>% 
  filter(fermented_food %in% filtered_fermented_foods) %>% 
  filter(completeness > 90 & contamination < 10)
