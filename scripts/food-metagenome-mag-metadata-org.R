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

# substrate categories reorganized to match broad categories in Caffrey 2024
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

## QUAST stats on genome assembly quality
combined_quast_stats <- read_tsv("metadata/raw_metadata/all_quast_stats.tsv", col_names = c("mag_id", "contigs_0", "contigs_1000", "contigs_5000", "contigs_10000", "contigs_25000", "contigs50000", "total_length_0", "total_length_1000", "total_length_5000", "total_length_10000", "total_length_25000", "total_length_50000", "contigs", "largest_contig", "total_length", "gc", "n50", "n90", "auN", "L50", "L90", "n_per_100kbp")) %>% 
  select(mag_id, contigs, total_length, gc, n50)

all_food_mags_metadata_cleaned <- left_join(all_food_mags_metadata_cleaned, combined_quast_stats) %>% 
  drop_na(contigs) %>% 
  select(mag_id, substrate, completeness, contamination, contigs, total_length, gc, n50, domain, phylum, class, order, family, genus, species, group, study_catalog)

hq_bac_food_mags <- all_food_mags_metadata_cleaned %>% 
  filter(completeness > 90) %>% 
  filter(contigs < 100) %>% 
  filter(domain == "Bacteria" )

all_bac_food_mags <- all_food_mags_metadata_cleaned %>% 
  filter(domain == "Bacteria")

all_euk_food_mags <- all_food_mags_metadata_cleaned %>% 
  filter(domain == "Eukaryota")

# cleaned metadata copied in both subdirectory and main metadata directory
write.csv(all_food_mags_metadata_cleaned, "metadata/cleaned_metadata/mag_datasets/2024-10-09-all-food-mag-metadata-cleaned.csv", row.names = FALSE, quote = FALSE)

write_tsv(all_food_mags_metadata_cleaned, "metadata/all-food-mags-metadata.tsv")

# HQ bac MAGS metadata copied in both subdirectory and main metadata directory
write.csv(hq_bac_food_mags, "metadata/cleaned_metadata/mag_datasets/2024-10-08-HQ-bacterial-food-mags-metadata.csv", row.names = FALSE, quote = FALSE)

write_tsv(hq_bac_food_mags, 'metadata/all-hq-bac-food-mag-metadata.tsv')

write.csv(all_bac_food_mags, "metadata/cleaned_metadata/mag_datasets/2024-10-11-all-bac-food-mags-metadata.csv", row.names = FALSE, quote = FALSE)

write.csv(all_euk_food_mags, "metadata/cleaned_metadata/mag_datasets/2024-10-11-all-euk-food-mags-metadata.csv", row.names = FALSE, quote = FALSE)

# plot mag stats with main group
all_food_mags_metadata_cleaned %>% 
  filter(domain == "Bacteria") %>% 
  filter(completeness > 90) %>% 
  filter(contigs < 100) %>% 
  count() # ~2000 "HQ" bacterial MAGs - 90% completeness, less than 100 contigs

all_food_mags_metadata_cleaned %>% 
  filter(domain == "Bacteria") %>% 
  filter(completeness > 90) %>% 
  filter(contigs < 100) %>% 
  group_by(group, substrate) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  print(n=41)

all_food_mags_metadata_cleaned %>% 
  filter(domain == "Bacteria") %>% 
  filter(completeness > 90) %>% 
  filter(contigs < 100) %>% 
  group_by(study_catalog) %>% 
  count()

# plot stats of substrate & group counts
substrate_categories_plot <- all_food_mags_metadata_cleaned %>%
  filter(domain == "Bacteria") %>% 
  filter(completeness > 90) %>% 
  filter(contigs < 100) %>%
  group_by(substrate) %>%
  mutate(substrate_count = n()) %>% 
  ungroup() %>%
  mutate(substrate = fct_lump_min(substrate, min = 20, other_level = "other")) %>%  
  ggplot(aes(y=fct_infreq(substrate))) +   
  geom_bar(aes(fill=substrate)) + 
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

phylo_groups_plot <- all_food_mags_metadata_cleaned %>%
  filter(domain == "Bacteria") %>% 
  filter(completeness > 90) %>% 
  filter(contigs < 100) %>%
  ggplot(aes(y=fct_infreq(group))) +   
  geom_bar(aes(fill=group)) + 
  scale_fill_brewer(palette = "Set2") +   
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
  ylab("Phylogenetic Group") + 
  ggtitle("Number of HQ Bacterial Genomes across Phylogenetic Groups")

ggsave("figs/hq-genomes-substrate-categories-plot.png", substrate_categories_plot, width=11, height=8, units=c("in"))
ggsave("figs/hq-genomes-phylo-groups-plot.png", phylo_groups_plot, width=11, height=8, units=c("in"))
