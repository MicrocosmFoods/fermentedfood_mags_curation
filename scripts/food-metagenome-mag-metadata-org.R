library(tidyverse)

# Organizing metadata for food MAGs from multiple studies

# Du 2023 BGC MAGs
du2023_metadata <- read.csv("metadata/raw_metadata/mag_datasets/FermentedFoodBGCs/Du2023_MAG_metadata.csv") %>% 
  select(MAG_id, Food.fermentation, Genome.size..bp., Completeness...., Contamination...., Domain, Phylum, Class, Order, Family, Genus, Species) %>% 
  mutate(study_catalog = c("Du2023")) %>% 
  mutate(substrate_category = Food.fermentation) %>% 
  mutate(source = NA)

colnames(du2023_metadata) <- c("mag_id", "specific_substrate", "genome_size", "completeness", "contamination", "domain", "phylum", "class", "order", "family", "genus", "species", "study_catalog", "substrate_category", "source")

du2023_metadata_cleaned <- du2023_metadata %>% 
  mutate(domain = gsub("d__", "", domain)) %>% 
  mutate(phylum = gsub("p__", "", phylum)) %>% 
  mutate(class = gsub("c__", "", class)) %>% 
  mutate(order = gsub("o__", "", order)) %>% 
  mutate(family = gsub("f__", "", family)) %>% 
  mutate(genus = gsub("g__", "", genus)) %>% 
  mutate(species = gsub("s__", "", species)) %>% 
  select(-genome_size) %>% 
  select(mag_id, substrate_category, specific_substrate, source, everything())

# Carlino 2024 food metagenomes
carlino2024_mag_metadata <- read.csv("metadata/raw_metadata/mag_datasets/CellFoodMetagenomics/2024-cell-food-metagenomics-mag-metadata.csv") %>% 
  mutate(dataset_sample_id = paste(dataset_id, sample_id, sep="__")) %>% 
  select(MAG_id, dataset_sample_id, completeness, contamination, superkingdom, phylum, class, order, family, genus, species) %>% 
  mutate(study_catalog = c("Carlino2024"))

colnames(carlino2024_mag_metadata) <- c("mag_id", "dataset_sample_id", "completeness", "contamination", "domain", "phylum", "class", "order", "family", "genus", "species", "study_catalog")

carlino2024_sample_metadata <- read.csv("metadata/raw_metadata/mag_datasets/CellFoodMetagenomics/2024-cell-food-metagenomics-sample-metadata.csv") %>% 
  select(sample_id, category, type)
colnames(carlino2024_sample_metadata) <- c("dataset_sample_id", "substrate", "specific_substrate")

carlino2024_complete_metadata <- left_join(carlino2024_mag_metadata, carlino2024_sample_metadata) %>% 
  mutate(source = dataset_sample_id) %>% 
  mutate(substrate_category = substrate) %>% 
  select(mag_id, substrate_category, specific_substrate, source, completeness, contamination, domain, phylum, class, order, family, genus, species, study_catalog)

# Updated v1.2.1 Carlino 2024 samples
updated_carlino2024_datasets <- read_tsv("metadata/raw_metadata/mag_datasets/CellFoodMetagenomics_v1.2.1/cFMD_datasets.tsv") %>% 
  mutate(dataset_id = Dataset) %>% 
  select(-Dataset)

updated_carlino2024_sample_metadata <- read_tsv("metadata/raw_metadata/mag_datasets/CellFoodMetagenomics_v1.2.1/cFMD_metadata.tsv") %>% 
  mutate(dataset_id = dataset_name) %>% 
  mutate(dataset_sample_id = paste(dataset_name, sample_id, sep="__")) %>% 
  select(-dataset_name)

updated_carlino2024_mag_metadata <- read_tsv("metadata/raw_metadata/mag_datasets/CellFoodMetagenomics_v1.2.1/cFMD_mags_list.tsv") %>% 
  mutate(dataset_sample_id = paste(dataset_id, sample_id, sep="__"))

updated_carlino2024_sample_dataset_info <- left_join(updated_carlino2024_sample_metadata, updated_carlino2024_datasets) %>% 
  filter(Version == "cFMDv1.2.1")

fermented_updated_carlino2024_sample_dataset_info <- updated_carlino2024_sample_dataset_info %>% 
  filter(`fermented/non-fermented` == "F") %>% 
  select(dataset_sample_id, macrocategory, category, type, subtype, sample_accession, project_accession, Version, Reference)

fermented_updated_carlino2024_mag_metadata <- left_join(updated_carlino2024_mag_metadata, fermented_updated_carlino2024_sample_dataset_info, by="dataset_sample_id", relationship = "many-to-many") %>% 
  filter(Version == "cFMDv1.2.1") %>% 
  filter(dataset_id != "SaakC_2023") # filter out Saak2023 because we already include those in our DB

# Caffrey 2024 MiFoDB
caffrey2024_metadata <- read.csv("metadata/raw_metadata/mag_datasets/MiFoDB/MiFoDB_beta_v3_AllReferences.csv") %>% 
  select(genome, description, substrate, source, completeness, contamination, domain, phylum, class, order, family, genus, species) %>% 
  drop_na() %>% 
  mutate(mag_id = genome) %>% 
  mutate(specific_substrate = description) %>% 
  mutate(substrate_category = substrate) %>% 
  select(-genome, -description, -substrate) %>% 
  select(mag_id, substrate_category, specific_substrate, source, everything()) %>% 
  mutate(study_catalog = c("Caffrey2024"))

# Saak 2023 cheese MAGs
saak2023_metadata <- read.csv("metadata/raw_metadata/mag_datasets/Saak2023CheeseMAGs/saak2023_cheese_mags_isolates_metadata.csv") %>% 
  select(Cheese, Bin.Name, Completeness, Contamination, gtdbtk.Classifaction) %>% 
  drop_na()

saak2023_metadata_cleaned <- saak2023_metadata %>% 
  separate(gtdbtk.Classifaction,
           into = c("domain", "phylum", "class", "order", "family", "genus", "species"),
           sep = ";") %>%
  mutate(across(c(domain, phylum, class, order, family, genus, species), ~ sub("^[a-z]__", "", .))) %>% 
  mutate(mag_id = Bin.Name) %>% 
  mutate(completeness = Completeness) %>% 
  mutate(contamination = Contamination) %>% 
  mutate(specific_substrate = c("cheese")) %>% 
  mutate(substrate_category = c("cheese")) %>% 
  mutate(source = paste("cheese_", Cheese)) %>% 
  mutate(study_catalog = c("Saak2023_CheeseMAGs")) %>% 
  select(mag_id, substrate_category, specific_substrate, source, completeness, contamination, domain, phylum, class, order, family, genus, species, study_catalog)

# Msystems Sourdough AABs
sourdough_aabs_metadata <- read.csv("metadata/raw_metadata/mag_datasets/msystemsAAB/SourdoughAAB-MAGs-metadata-modified.csv") %>% 
  filter(Source_environment == "sourdough") %>% 
  select(assembly_name, Species, Source_category, Completeness, Contam, Source)

colnames(sourdough_aabs_metadata) <- c("mag_id", "species", "source_category", "completeness", "contamination", "source")

# taxonomy information from GTDB for Acetobacter spp.
sourdough_aabs_metadata_cleaned <- sourdough_aabs_metadata %>%
  mutate(substrate_category = c("grains")) %>% 
  mutate(specific_substrate = c("sourdough")) %>% 
  mutate(domain = c("Bacteria")) %>% 
  mutate(phylum = c("Pseudomonadota")) %>% 
  mutate(class = c("Alphaproteobacteria")) %>% 
  mutate(order = c("Acetobacterales")) %>% 
  mutate(family = c("Acetobacteraceae")) %>% 
  select(-source_category) %>% 
  mutate(genus = sub(" .*", "", species)) %>% 
  mutate(study_catalog = c("Rappaport2024_sourdough")) %>% 
  mutate(source = c("Rappaport2024_sourdough")) %>% 
  select(mag_id, substrate_category, specific_substrate, source, completeness, contamination, domain, phylum, class, order, family, genus, species, study_catalog)
  

# combine all dataset metadata
caffrey_carlino_dbs <- bind_rows(carlino2024_complete_metadata, caffrey2024_metadata) %>%
  filter(contamination < 10)

write.csv(caffrey_carlino_dbs, "metadata/cleaned_metadata/mag_datasets/caffrey_carlino_datasets.csv", quote = FALSE, row.names = FALSE)

all_food_mags_metadata <- bind_rows(du2023_metadata_cleaned, carlino2024_complete_metadata, caffrey2024_metadata, saak2023_metadata_cleaned, sourdough_aabs_metadata_cleaned) %>% 
  filter(contamination < 10) %>% 
  mutate(substrate_category = tolower(substrate_category)) %>% 
  mutate(specific_substrate = tolower(specific_substrate))

# substrate categories cleaned up as best as possible before manual curation
# below the first word listed is the broad category, after the arrow are entries that were changed to the broad category. 
# phyla are simplified into "groups" based on previous NCBI major groupings
# dairy -> cheese, nunu, koumiss
# grains  -> fermented_grains, sourdough
# salt     
# soy -> soy sauce
# legumes  -> fermented_legumes
# sugar    
# seafood  -> fermented_fish, fish
# meat  -> fermented_meat  
# seed -> fermented_seeds
# However fruit and vegetable separate categories were combined as the Carlino2024 study has these grouped together. Below are the new custom categories added on top of Caffrey2024, or merged into categories from other studies
# fruits_and_vegetables -> fruit, vegetables, fermented_fruits_and_vegetables, kimchi
# fermented_beverages -> alcohol, chinese liquor, wine, kombucha
# \\<pattern\\>

all_food_mags_metadata_cleaned <- all_food_mags_metadata %>%
  mutate(
    # Trim whitespace
    original_substrate = specific_substrate,
    specific_substrate = trimws(specific_substrate, which = "both"),
    
    # Assign specific_substrate based on keywords in original_substrate
    specific_substrate = case_when(
      grepl("cheese", specific_substrate, ignore.case = TRUE) ~ "cheese",
      grepl("kefir", specific_substrate, ignore.case = TRUE) ~ "kefir",
      grepl("yoghurt|yogurt", specific_substrate, ignore.case = TRUE) ~ "yogurt",
      grepl("kombucha", specific_substrate, ignore.case = TRUE) ~ "kombucha",
      grepl("rice", specific_substrate, ignore.case = TRUE) ~ "rice",
      grepl("beer", specific_substrate, ignore.case = TRUE) ~ "beer",
      grepl("sauerkraut", specific_substrate, ignore.case = TRUE) ~ "sauerkraut",
      grepl("wine", specific_substrate, ignore.case = TRUE) ~ "wine",
      grepl("sourdough", specific_substrate, ignore.case = TRUE) ~ "sourdough",
      grepl("sausage", specific_substrate, ignore.case = TRUE) ~ "sausage",
      grepl("kimchi", specific_substrate, ignore.case = TRUE) ~ "kimchi",
      grepl("cocoa", specific_substrate, ignore.case = TRUE) ~ "cocoa",
      grepl("maize", specific_substrate, ignore.case = TRUE) ~ "maize",
      TRUE ~ specific_substrate  # Default to the original specific_substrate
    ),
    
    # Apply replacements to categorize substrate
    substrate_category = substrate_category %>%
      gsub("cheese", "dairy", .) %>%
      gsub("dairy ", "dairy", .) %>%
      gsub("fermented_fruits_and_vegetables", "fruits_and_vegetables", .) %>%
      gsub("fermented_meat", "meat", .) %>%
      gsub("fermented_grains", "grains", .) %>%
      gsub("fermented_legumes", "legumes", .) %>%
      gsub("fermented_seeds", "seed", .) %>%
      gsub("\\<fruit\\>", "fruits_and_vegetables", .) %>%
      gsub("\\<vegetable\\>", "fruits_and_vegetables", .) %>%
      gsub("fermented_fish", "seafood", .) %>%
      gsub("chinese liquor", "fermented_beverages", .) %>%
      gsub("wine", "fermented_beverages", .) %>%
      gsub("fish", "seafood", .) %>%
      gsub("sourdough", "grains", .) %>%
      gsub("soy sauce", "soy", .) %>%
      gsub("nunu", "dairy", .) %>%
      gsub("kimchi", "fruits_and_vegetables", .) %>%
      gsub("kombucha", "fermented_beverages", .) %>%
      gsub("alcohol", "fermented_beverages", .) %>%
      gsub("koumiss", "dairy", .),
    
    # Assign group based on phylum
    group = case_when(
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
      TRUE ~ phylum  # Default to original phylum if no match
    ) 
  ) %>% 
  select(mag_id, original_substrate, substrate_category, everything())

# write out this intermediate table for manually curating categories
write_tsv(all_food_mags_metadata_cleaned, "metadata/raw_metadata/2024-11-01-raw-metadata-for-curation.tsv")

####################
# Manually curated metadata 
# Sample description - original description from the study of origin
# Fermented Food - the food type
# specific_substrate - the food substrate, such as cheese is the fermented food but the specific substrate is milk
# substrate_category - larger category of the substrate, such as cheese is dairy
# General category - broadest category, such as fermented_beverages could have several types of drinks in them
####################

# read in the curated metadata table to join with assembly quality information
manually_curated_metadata <- read_tsv("metadata/raw_metadata/manually_curated_metadata/2024-11-04-manually-curated-mag-metadata.tsv")

colnames(manually_curated_metadata) <- c("mag_id", "sample_description", "fermented_food", "specific_substrate", "substrate_category", "general_category", "source", "completeness", "contamination", "domain", "phylum", "class", "order", "family", "genus", "species", "study_catalog", "phylo_group")

carlino_sample_info <- read_tsv("metadata/raw_metadata/mag_datasets/CellFoodMetagenomics/2024-11-04-Carlino-sample-metadata.tsv") %>% 
  mutate(source = paste(dataset_name, sample_id, sep="__")) %>% 
  select(source, "fermented/non-fermented")

colnames(carlino_sample_info) <- c("source", "fermented_status")

manually_curated_metadata_cleaned <- manually_curated_metadata %>% 
  left_join(carlino_sample_info) %>% 
  filter(!is.na(substrate_category)) %>% 
  filter(is.na(fermented_status) | fermented_status == "F") %>% 
  select(-fermented_status) %>%
  mutate(across(
    c("fermented_food", "specific_substrate", "substrate_category", "general_category"),
    ~ str_trim(.) %>%  
      str_replace_all("[,/]", " ") %>%
      str_replace_all("\\s+", "_"))) %>% 
  filter(!is.na(fermented_food))

## QUAST stats on genome assembly quality
## combine the manually curated metadata with assembly quality
## Then split out bacteria, archaea, and eukaryote assemblies since mostly focusing on processing bacterial assemblies of varying quality
combined_quast_stats <- read_tsv("metadata/raw_metadata/all_quast_stats.tsv", col_names = c("mag_id", "contigs_0", "contigs_1000", "contigs_5000", "contigs_10000", "contigs_25000", "contigs50000", "total_length_0", "total_length_1000", "total_length_5000", "total_length_10000", "total_length_25000", "total_length_50000", "contigs", "largest_contig", "total_length", "gc", "n50", "n90", "auN", "L50", "L90", "n_per_100kbp")) %>% 
  select(mag_id, contigs, total_length, gc, n50)

all_food_mags_metadata_cleaned_curated <- left_join(manually_curated_metadata_cleaned, combined_quast_stats) %>% 
  drop_na(contigs) %>% 
  select(mag_id, sample_description, fermented_food, specific_substrate, substrate_category, general_category, completeness, contamination, contigs, total_length, gc, n50, domain, phylum, class, order, family, genus, species, phylo_group, study_catalog, source)

hq_bac_food_mags <- all_food_mags_metadata_cleaned_curated %>% 
  filter(completeness > 90) %>% 
  filter(contigs < 100)

all_food_mags_metadata_cleaned_curated %>% 
  filter(completeness > 50) %>% 
  filter(contamination < 10) %>% 
  filter(contigs < 200) %>% 
  count() # 5550 MQ genomes

mq_bac_food_mags <- all_food_mags_metadata_cleaned_curated %>% 
  filter(domain == 'Bacteria') %>% 
  filter(completeness > 50) %>% 
  filter(contamination < 10) %>% 
  filter(contigs < 200)

# cleaned metadata
write_tsv(all_food_mags_metadata_cleaned_curated, "metadata/cleaned_metadata/2025-03-05-all-bac-mag-metadata-cleaned.tsv")
write_tsv(all_food_mags_metadata_cleaned, "metadata/cleaned_metadata/2025-03-05-all-mags-metadata-cleaned.tsv")


# mag stats
all_food_mags_metadata_cleaned_curated %>% 
  filter(completeness > 90) %>% 
  filter(contigs < 100) %>% 
  count() # ~2000 "HQ" bacterial MAGs - 90% completeness, less than 100 contigs

all_food_mags_metadata_cleaned_curated %>% 
  filter(completeness > 90) %>% 
  filter(contigs < 100) %>% 
  group_by(phylo_group, substrate_category) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  print(n=41)

all_food_mags_metadata_cleaned_curated %>% 
  filter(completeness > 90) %>% 
  filter(contigs < 100) %>% 
  group_by(study_catalog) %>% 
  count()

hqmq_mags_list <- all_food_mags_metadata_cleaned_curated %>% 
  filter(completeness > 50) %>% 
  filter(contamination < 10) %>% 
  filter(n50 > 5000) %>% 
  pull(mag_id)

eukaryotic_mags <- all_food_mags_metadata_cleaned %>% 
  filter(domain=="Eukaryota") # 770 eukaryotic genomes


all_food_mags_metadata_cleaned %>% 
  filter(domain == "Archaea") # 17 archaeal genomes

# plot stats of substrate & group counts
substrate_categories_plot <- all_food_mags_metadata_cleaned_curated %>%
  filter(completeness > 90) %>% 
  filter(contigs < 100) %>%
  group_by(substrate_category) %>%
  mutate(substrate_count = n()) %>% 
  ungroup() %>%
  mutate(substrate_category = fct_lump_min(substrate_category, min = 20, other_level = "other")) %>%  
  ggplot(aes(y=fct_infreq(substrate_category))) +   
  geom_bar(aes(fill=substrate_category)) + 
  scale_fill_brewer(palette = "Set3") +   
  theme_classic() +                   
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12), 
    axis.text.y = element_text(size = 12), 
    text = element_text(size = 12), 
    legend.position = "none", 
    legend.text = element_text(size = 12),    
    legend.title = element_text(size = 12),
    title = element_text(size = 12, face = "bold")
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  xlab("Number of Bacterial Genomes") +
  ylab("Food Substrate Category") + 
  ggtitle("High-Quality Genomes across Food Substrate Categories")

substrate_categories_plot

ggsave("figs/hq-genomes-substrate-categories-plot.png", substrate_categories_plot, width=11, height=8, units=c("in"))

# mag stats of different metadata categories
length(unique(all_food_mags_metadata_cleaned$fermented_food))
length(unique(all_food_mags_metadata_cleaned$specific_substrate))
length(unique(all_food_mags_metadata_cleaned$substrate_category))
length(unique(all_food_mags_metadata_cleaned$general_category))

# write list of modified HQMQ MAGs to process
write_lines(hqmq_mags_list, "metadata/2025-03-05-modified-hqmq-bac-mags-list.txt")

# write out eukaryotic MAGs to metadata since there are a sizeable amount of them
write_tsv(eukaryotic_mags, "metadata/cleaned_metadata/2025-03-05-all-euk-mags-metadata.tsv")
