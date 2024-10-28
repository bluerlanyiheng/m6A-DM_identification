#mDM identification and definition
library(readr)
library(dplyr)
library(GenomicRanges)
library(maftools)
library(tidyverse)
library(stringr)
library(tidyr)
setwd("/data/sdd/workplace")

#################################################################################################################################
# merged mutation-sample association data
mutation_set <- mutation_set %>%
  distinct(Gene.name, HGVSG, Sample.name, .keep_all = TRUE)
mutation_set <- mutation_set %>%
  filter(grepl("Substitution", Mutation.Description))

## Recurrence screening 
mutation_set <- mutation_set %>%
  count(HGVSG) %>%
  mutate(mutation_count = n())
mDMs <- mutation_set %>%
  filter(mutation_count > 1)

## m6A screening
mDMs <- mDMs %>%
  filter(grepl("A>", Mutation.CDS))

mDM <- mDM %>%
  separate(Mutation.genome.position, into = c("Seqname", "Start", "End"), sep = "[:-]")

# Convert Seqname, Start, and End columns to GenomicRanges1 objects
GenomicRanges1 <- with(mDM, GRanges(Chromosome = Seqname, IRanges(Start, End)))

# Split m6A_range into three columns
m6A_range_split <- separate(m6A_range, into = c("m6A_Seqname", "m6A_Start", "m6A_End"), sep = "[:-]")

# Convert m6A_Seqname, m6A_Start, and m6A_End columns to GenomicRanges2 objects
GenomicRanges2 <- with(m6A_range_split, GRanges(Chromosome = m6A_Seqname, IRanges(m6A_Start, m6A_End)))

# Find intersection of intervals between GenomicRanges1 and GenomicRanges2
GenomicRanges3 <- GenomicRanges1[GenomicRanges1 %over% GenomicRanges2]

# Convert GenomicRanges3 to "Seqname:Start-End" format
overlap_pos <- with(GenomicRanges3, paste0(seqnames(GenomicRanges3), ":", start(GenomicRanges3), "-", end(GenomicRanges3)))

# Filter rows in mDM where Mutation.genome.position values are present in overlap_pos
mDM_filtered <- mDM[with(mDM, paste0(Seqname, ":", Start, "-", End) %in% overlap_pos), ]

## Reliability screening
# Generate the HGVSG2 column
mDM_filtered$HGVSG2 <- gsub(":g\\.", "-", mDM_filtered$HGVSG)
mDM_filtered$HGVSG2 <- gsub(">", "-", mDM_filtered$HGVSG2)
mDM_filtered$HGVSG2 <- gsub("([[:alpha:]])([[:digit:]])", "\\1-\\2", mDM_filtered$HGVSG2)

# Initialize an empty mDM_final data frame
mDM_final <- data.frame()

# Iterate through each row of mDM_filtered
for (i in 1:nrow(mDM_filtered)) {
  # Check if the HGVSG2 value exists in the mutation_iden column of gnomAD_data
  if (mDM_filtered$HGVSG2[i] %in% gnomAD_data$mutation_iden) {
    # Retrieve the corresponding INFO value from gnomAD_data
    info_value <- gnomAD_data$INFO[gnomAD_data$mutation_iden == mDM_filtered$HGVSG2[i]]
    
    # Check if the INFO value satisfies the condition
    if (info_value < 0.05) {
      # Add the row to mDM_final if the condition is met
      mDM_final <- rbind(mDM_final, mDM_filtered[i, ])
    }
  }
  
  # Check if the Gene.name value exists in WER_association
  if (mDM_filtered$Gene.name[i] %in% WER_association[, 2]) {
  
    mDM_final <- rbind(mDM_final, mDM_filtered[i, ])
  }
  
  # Check if the Mutation.genome.position value exists in the pos column of atoi_edit_pos
  if (!(mDM_filtered$Mutation.genome.position[i] %in% atoi_edit_pos$pos)) {
    
    mDM_final <- rbind(mDM_final, mDM_filtered[i, ])
  }
}

# Output the results to mDM_final
mDM_final

# Create empty data frames smDM and mmDM
smDM <- data.frame()
mmDM <- data.frame()

# Classify based on the values in Mutation.Description column
for (i in 1:nrow(mDM_final)) {
  if (mDM_final[i, "Mutation.Description"] == "Substitution - coding silent") {
    # Append the rows with "Substitution - coding silent" to smDM
    smDM <- rbind(smDM, mDM_final[i, ])
  } else if (mDM_final[i, "Mutation.Description"] == "Substitution - Missense") {
    # Append the rows with "Substitution - Missense" to mmDM
    mmDM <- rbind(mmDM, mDM_final[i, ])
  }
}
# Write mDM_final to CSV
write.csv(mDM_final, "mDM_final.csv", row.names = FALSE)

# Write smDM to CSV
write.csv(smDM, "smDM.csv", row.names = FALSE)

# Write mmDM to CSV
write.csv(mmDM, "mmDM.csv", row.names = FALSE)


