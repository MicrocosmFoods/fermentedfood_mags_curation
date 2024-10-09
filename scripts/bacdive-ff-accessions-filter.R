library(tidyverse)

# bacdive accessions
bacdive_accessions <- read_tsv("raw_metadata/isolate_genomes/bacdive/2024-10-08-parse-bacdive-accessions.tsv")
# bacdive fermented food records
bacdive_ff_records <- read_tsv("raw_metadata/isolate_genomes/bacdive/BacDive_fermented_food_filtered_list.tsv")

# join together
bacdive_ff_accessions <- left_join(bacdive_ff_records, bacdive_accessions) %>% 
  drop_na(GCA) %>% 
  select(GCA, Species, 'Isolation source', 'Category 3')

colnames(bacdive_ff_accessions) <- c("genbank_id", "isolate", "isolation_source", "broad_category")

bacdive_ff_metadata <- bacdive_ff_accessions %>% 
  mutate(genbank_accession = paste(genbank_id, ".1", sep=""))

write_tsv(bacdive_ff_metadata, "cleaned_metadata/isolate_genomes/2024-10-08-bacdive-ff-isolates-metadata.tsv")
