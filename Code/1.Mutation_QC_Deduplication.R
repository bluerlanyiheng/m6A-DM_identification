#data:correlation of mutation and tumor samples
#Quality Control，Data Deduplication, DNA point substitution mutations screening

library(readr)
library(maftools)
library(tidyverse)
library(dplyr)
library(stringr)

setwd("/data/sdd/workplace")

#merged mutation datasets
#cosmic_data
cosmic_data <- read_tsv("CosmicMutantExport.tsv") 
cosmic_subset <- cosmic_data %>% select("Gene.name", "Gene.CDS.length", "Sample.name", "ID_sample", "Primary.site", 
                                        "Genome.wide.screen", "Mutation.CDS", "Mutation.AA",
                                        "Mutation.Description", "Mutation.genome.position","Sample.Type", "HGVSG")
cosmic_subset2 <- cosmic_subset %>% filter(!grepl("TCGA", Sample.name))


# Filter rows where "Sample.name" contains "TCGA"
tcga_data <-read.maf(maf = "TCGA-mc3.maf")
tcga_subset <- tcga_data %>% select("Hugo_Symbol", "Chromosome","Start_position","End_Position", 
                                    "Variant_Classification","Variant_Type", "Tumor_Sample_Barcode") 
tcga_subset2 <- tcga_subset %>%mutate(Mutation.genome.position = paste(Chromosome, Start_position, sep = ":") 
                                      %>%paste(End_Position, sep = "-", collapse = ""))
tcga_subset3 <- cosmic_subset %>% filter(Mutation.genome.position %in% tcga_subset2$Mutation.genome.position)
tcga_subset4 <- tcga_subset3 %>%filter(str_detect(Sample.name, "TCGA"))

# Merge the two datasets and remove duplicates
mutation_set <- bind_rows(cosmic_subset2, tcga_subset4) %>%
mutation_set <- distinct(mutation_set)

#QC & remove duplicates
# 1. Remove duplicates based on "Gene.name", "HGVSG", and "Sample.name" while keeping all information
mutation_set <- mutation_set %>%
  distinct(Gene.name, HGVSG, Sample.name, .keep_all = TRUE)

# 2. Remove rows where "Sample.Type" contains "cell-line", "xenograft", or "organoid" strings
mutation_set <- mutation_set %>%
  filter(!grepl("cell-line|xenograft|organoid", Sample.Type))

# 3. Remove rows where "Genome.wide.screen" contains "n" string
mutation_set <- mutation_set %>%
  filter(!grepl("n", Genome.wide.screen))

# 4. Keep rows where "Mutation.Description" contains "Substitution - coding silent" or "Substitution - Missense" strings
mutation_set <- mutation_set %>%
  filter(grepl("Substitution - coding silent|Substitution - Missense", Mutation.Description))

# Create synonymous_mutation_set, containing rows with "Substitution - coding silent" in Mutation.Description
synonymous_mutation_set <- mutation_set %>%
  filter(Mutation.Description == "Substitution - coding silent")

# Create missense_mutation_set, containing rows with "Substitution - Missense" in Mutation.Description
missense_mutation_set <- mutation_set %>%
  filter(Mutation.Description == "Substitution - Missense")


