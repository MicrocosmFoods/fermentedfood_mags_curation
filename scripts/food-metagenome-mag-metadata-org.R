library(tidyverse)

# Organizing metadata for food MAGs from multiple studies

# Du 2023 BGC MAGs
du2023_metadata <- read.csv("raw_metadata/mag_datasets/FermentedFoodBGCs/Du2023_MAG_metadata.csv") %>% 
  select(MAG_id, Food.fermentation, Genome.size..bp., Completeness...., Contamination...., Domain, Phylum, Class, Order, Family, Genus, Species) %>% 
  mutate(study_catalog = c("Du2023"))

colnames(du2023_metadata) <- c("mag_id", "substrate", "genome_size", "completeness", "contamination", "domain", "phylum", "class", "order", "family", "genus", "species", "study_catalog")

du2023_metadata_cleaned <- du2023_metadata %>% 
  mutate(domain = gsub("d__", "", domain)) %>% 
  mutate(phylum = gsub("p__", "", phylum)) %>% 
  mutate(class = gsub("c__", "", class)) %>% 
  mutate(order = gsub("o__", "", order)) %>% 
  mutate(family = gsub("f__", "", family)) %>% 
  mutate(genus = gsub("g__", "", genus)) %>% 
  mutate(species = gsub("s__", "", species)) %>% 
  select(-genome_size)

# Carlino 2024 food metagenomes
carlino2024_mag_metadata <- read.csv("raw_metadata/mag_datasets/CellFoodMetagenomics/2024-cell-food-metagenomics-mag-metadata.csv") %>% 
  mutate(dataset_sample_id = paste(dataset_id, sample_id, sep="__")) %>% 
  select(MAG_id, dataset_sample_id, completeness, contamination, superkingdom, phylum, class, order, family, genus, species) %>% 
  mutate(study_catalog = c("Carlino2024"))
colnames(carlino2024_mag_metadata) <- c("mag_id", "dataset_sample_id", "completeness", "contamination", "domain", "phylum", "class", "order", "family", "genus", "species", "study_catalog")
carlino2024_sample_metadata <- read.csv("raw_metadata/mag_datasets/CellFoodMetagenomics/2024-cell-food-metagenomics-sample-metadata.csv") %>% 
  select(sample_id, category, type)
colnames(carlino2024_sample_metadata) <- c("dataset_sample_id", "substrate", "specific_substrate")

carlino2024_complete_metadata <- left_join(carlino2024_mag_metadata, carlino2024_sample_metadata) %>% 
  mutate(source = dataset_sample_id) %>% 
  select(mag_id, substrate, source, completeness, contamination, domain, phylum, class, order, family, genus, species, study_catalog)

# Caffrey 2024 MiFoDB
caffrey2024_metadata <- read.csv("raw_metadata/mag_datasets/MiFoDB/MiFoDB_beta_v3_AllReferences.csv") %>% 
  select(genome, substrate, source, completeness, contamination, domain, phylum, class, order, family, genus, species) %>% 
  drop_na() %>% 
  mutate(mag_id = genome) %>% 
  select(-genome) %>% 
  mutate(study_catalog = c("Caffrey2024"))

# Saak 2023 cheese MAGs
saak2023_metadata <- read.csv("raw_metadata/mag_datasets/Saak2023CheeseMAGs/saak2023_cheese_mags_isolates_metadata.csv") %>% 
  select(Bin.Name, Completeness, Contamination, gtdbtk.Classifaction) %>% 
  drop_na()

saak2023_metadata_cleaned <- saak2023_metadata %>% 
  separate(gtdbtk.Classifaction,
           into = c("domain", "phylum", "class", "order", "family", "genus", "species"),
           sep = ";") %>%
  mutate(across(c(domain, phylum, class, order, family, genus, species), ~ sub("^[a-z]__", "", .))) %>% 
  mutate(mag_id = Bin.Name) %>% 
  mutate(completeness = Completeness) %>% 
  mutate(contamination = Contamination) %>% 
  mutate(substrate = c("dairy")) %>% 
  mutate(study_catalog = c("Saak2023_CheeseMAGs")) %>% 
  select(mag_id, substrate, completeness, contamination, domain, phylum, class, order, family, genus, species, study_catalog)

# Msystems Sourdough AABs
sourdough_aabs_metadata <- read.csv("raw_metadata/mag_datasets/msystemsAAB/SourdoughAAB-MAGs-metadata-modified.csv") %>% 
  filter(Source_environment == "sourdough") %>% 
  select(assembly_name, Species, Source_category, Completeness, Contam)

colnames(sourdough_aabs_metadata) <- c("mag_id", "species", "source_category", "completeness", "contamination")

# taxonomy information from GTDB for Acetobacter spp.
sourdough_aabs_metadata_cleaned <- sourdough_aabs_metadata %>% 
  mutate(substrate = c("grains")) %>% 
  mutate(domain = c("Bacteria")) %>% 
  mutate(phylum = c("Pseudomonadota")) %>% 
  mutate(class = c("Alphaproteobacteria")) %>% 
  mutate(order = c("Acetobacterales")) %>% 
  mutate(family = c("Acetobacteraceae")) %>% 
  select(-source_category) %>% 
  mutate(genus = sub(" .*", "", species)) %>% 
  mutate(study_catalog = c("Rappaport2024_sourdough")) %>% 
  select(mag_id, substrate, completeness, contamination, domain, phylum, class, order, family, genus, species, study_catalog)
  

# combine all dataset metadata
caffrey_carlino_dbs <- bind_rows(carlino2024_complete_metadata, caffrey2024_metadata) %>% 
  filter(contamination < 10)

write.csv(caffrey_carlino_dbs, "cleaned_metadata/mag_datasets/caffrey_carlino_datasets.csv", quote = FALSE, row.names = FALSE)

all_food_mags_metadata <- bind_rows(du2023_metadata_cleaned, carlino2024_complete_metadata, caffrey2024_metadata, saak2023_metadata_cleaned, sourdough_aabs_metadata_cleaned) %>% 
  select(-source) %>% 
  filter(contamination < 10) %>% 
  mutate(substrate = tolower(substrate))

all_food_mags_metadata %>% 
  ggplot(aes(x=contamination, y=completeness)) +
  geom_point(aes(color=phylum))

all_food_mags_metadata %>% 
  group_by(substrate) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  print(n=40)

all_food_mags_metadata %>% 
  group_by(substrate) %>% 
  filter(study_catalog == "Caffrey2024") %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  select(-n)

# substrate categories reorganized to match broad categories in Caffrey 2024
# below the first word listed is the broad category, after the arrow are entries that were changed to the broad category. 
# dairy -> cheese, milk kefir, nunu, koumiss
# grains  -> fermented_grains, sourdough
# salt     
# soy -> soy sauce
# legumes  -> fermented_legumes
# sugar    
# seafood  -> fermented_fish, fish
# meat  -> fermented_meat  
# seed -> fermented_seeds
# However fruit and vegetable separate cateogries were combined as the Carlino2024 study has these grouped together. Below are the new custom categories added on top of Caffrey2024, or merged into categories from other studies
# fruits_and_vegetables -> fruit, vegetables, fermented_fruits_and_vegetables, kimchi
# fermented_beverages -> alcohol, chinese liquor, wine, kombucha
# \\<pattern\\>

all_food_mags_metadata_cleaned <- all_food_mags_metadata %>% 
  mutate(substrate = gsub("cheese", "dairy", substrate)) %>% 
  mutate(substrate = gsub("dairy ", "dairy", substrate)) %>% 
  mutate(substrate = gsub("fermented_fruits_and_vegetables", "fruits_and_vegetables", substrate)) %>% 
  mutate(substrate = gsub("fermented_meat", "meat", substrate)) %>% 
  mutate(substrate = gsub("fermented_grains", "grains", substrate)) %>% 
  mutate(substrate = gsub("fermented_legumes", "legumes", substrate)) %>% 
  mutate(substrate = gsub("fermented_seeds", "seed", substrate)) %>% 
  mutate(substrate = gsub("\\<fruit\\>", "fruits_and_vegetables", substrate)) %>% 
  mutate(substrate = gsub("\\<vegetable\\>", "fruits_and_vegetables", substrate)) %>% 
  mutate(substrate = gsub("fermented_fish", "seafood", substrate)) %>% 
  mutate(substrate = gsub("milk kefir", "dairy", substrate)) %>% 
  mutate(substrate = gsub("chinese liquor", "fermented_beverages", substrate)) %>% 
  mutate(substrate = gsub("wine", "fermented_beverages", substrate)) %>% 
  mutate(substrate = gsub("fish", "seafood", substrate)) %>% 
  mutate(substrate = gsub("sourdough", "grains", substrate)) %>% 
  mutate(substrate = gsub("soy sauce", "soy", substrate)) %>% 
  mutate(substrate = gsub("nunu", "dairy", substrate)) %>% 
  mutate(substrate = gsub("kimchi", "fruits_and_vegetables", substrate)) %>% 
  mutate(substrate = gsub("kombucha", "fermented_beverages", substrate)) %>% 
  mutate(substrate = gsub("alcohol", "fermented_beverages", substrate)) %>% 
  mutate(substrate = gsub("koumiss", "dairy", substrate))

all_food_mags_metadata_cleaned <- all_food_mags_metadata_cleaned %>% 
  mutate(group = case_when(
    phylum == "Pseudomonadota" ~ "Proteobacteria",
    str_detect(phylum, "Firmicutes") ~ "Firmicutes",
    phylum == "Actinobacteriota" ~ "Actinobacteria",
    phylum == "Actinomycetota" ~ "Actinobacteria",
    phylum == "Bacteroidota" ~ "Bacteroidetes",
    phylum == "Bacillota" ~ "Firmicutes",
    phylum == "Halobacteriota" ~ "Euryarchaea",
    phylum == "Acidobacteriota" ~ "Acidobacteria", 
    phylum == "Armatimonadota" ~ "Armatimonadetes",
    phylum == "Deinococcota" ~ "Deinococcus_Thermus",
    phylum == "Thermoproteota" ~ "Euryarchaea",
    phylum == "Euryarchaeota" ~ "Euryarchaea",
    phylum == "Balneolaeota" ~ "Bacteroidetes",
    TRUE ~ phylum  # Keeps the original phylum name if no rule matches
  ))

all_food_mags_metadata_cleaned %>% 
  group_by(group) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  print(n=30)
  

all_food_mags_metadata_cleaned %>% 
  group_by(substrate) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  print(n=20)

all_food_mags_metadata_cleaned %>% 
  filter(completeness > 90) %>% 
  filter(domain == "Bacteria") %>% 
  count() # ~5900 HQ bacterial MAGs

hq_bac_food_mags <- all_food_mags_metadata_cleaned %>% 
  filter(completeness > 90) %>% 
  filter(domain == "Bacteria" )

all_bac_food_mags <- all_food_mags_metadata_cleaned %>% 
  filter(domain == "Bacteria")

all_euk_food_mags <- all_food_mags_metadata_cleaned %>% 
  filter(domain == "Eukaryota")

write.csv(all_food_mags_metadata_cleaned, "cleaned_metadata/mag_datasets/2024-10-09-all-food-mag-metadata-cleaned.csv", row.names = FALSE, quote = FALSE)

write.csv(hq_bac_food_mags, "cleaned_metadata/mag_datasets/2024-10-08-HQ-bacterial-food-mags-metadata.csv", row.names = FALSE, quote = FALSE)

write.csv(all_bac_food_mags, "cleaned_metadata/mag_datasets/2024-10-11-all-bac-food-mags-metadata.csv", row.names = FALSE, quote = FALSE)

write.csv(all_euk_food_mags, "cleaned_metadata/mag_datasets/2024-10-11-all-euk-food-mags-metadata.csv", row.names = FALSE, quote = FALSE)
