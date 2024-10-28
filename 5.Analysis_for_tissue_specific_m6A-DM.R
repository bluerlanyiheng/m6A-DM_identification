library(readr)
library(dplyr)
library(GenomicRanges)
library(maftools)
library(tidyverse)
library(stringr)
library(tidyr)
library(ggplot2)
library(readxl)
#################################################################################
# Get tissue-specific peak regions
bed_files <- list.files(pattern = "merged\\.bed$", full.names = TRUE)
if (length(bed_files) != 11) {
  stop("Expected 11 BED files, but found ", length(bed_files))
}

# Read BED files and create GRanges objects
gr_list <- lapply(bed_files, function(file) {
  bed_data <- read.table(file, header = FALSE, stringsAsFactors = FALSE)
  gr <- GRanges(
    seqnames = Rle(bed_data$V1),
    ranges = IRanges(start = bed_data$V2, end = bed_data$V3),
    strand = Rle("*")  # Assume no strand information, use "*" for unknown
    
  )
  return(gr)
})

# Name the GRanges objects based on tissue names
names(gr_list) <- sapply(bed_files, function(file) {
  basename(file) %>% sub("_.*", "", .)
})
names(gr_list) <- gsub("Central Nervous System", "Central_Nervous_System", names(gr_list))
names(gr_list) <- gsub("Large Intestine", "Large_Intestine", names(gr_list))

# Add GRanges objects to the global environment
list2env(gr_list, envir = .GlobalEnv)
print(names(gr_list))

# Read mutation data

smDM <- read.table("input.txt", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
# Determine if each mutation overlaps with the corresponding tissue's peak regions

smDM$InPeak <- sapply(1:nrow(smDM), function(i) {
  chrom <- smDM$Chrom[i]
  start <- smDM$Start[i]
  end <- smDM$End[i] + 1
  tissue <- smDM$Tissue[i]
  
  mutation_gr <- GRanges(
    seqnames = Rle(chrom),
    ranges = IRanges(start = start, end = end),
    strand = Rle("*")
  )
  
  tissue_gr <- gr_list[[tissue]]
  
  if (is.null(tissue_gr)) {
    return(FALSE)
  }
  
  overlaps <- findOverlaps(mutation_gr, tissue_gr)
  has_overlap <- length(overlaps) > 0
  
  return(has_overlap)
})

# Write the results to a new file
write.table(smDM, file = "smDM_in_tissue-specific_peaks.txt", row.names = FALSE, quote = FALSE, sep = "\t")

#################################################################################
# Get smDM data in pancancer m6A peaks

bed_file <- "all_merged.bed"
bed_data <- read.table(bed_file, header = FALSE, stringsAsFactors = FALSE)
all_merged_gr <- GRanges(
  seqnames = Rle(bed_data$V1),
  ranges = IRanges(start = bed_data$V2, end = bed_data$V3),
  strand = Rle("*")  # Assume no strand information, use "*" for unknown
  
)

# Read mutation data
smDM2 <- read.table("input.txt", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)

# Determine if each mutation overlaps with all peak regions
smDM2$InallPeak <- sapply(1:nrow(smDM2), function(i) {
  chrom <- smDM$Chrom[i]
  start <- smDM$Start[i]
  end <- smDM$End[i] + 1
  mutation_gr <- GRanges(
    seqnames = Rle(chrom),
    ranges = IRanges(start = start, end = end),
    strand = Rle("*")
  )
  
  overlaps <- findOverlaps(mutation_gr, all_merged_gr)
  has_overlap <- length(overlaps) > 0
  
  return(has_overlap)
})

# Write the results to a new file
write.table(smDM2, file = "smDM_in_all_peaks.txt", row.names = FALSE, quote = FALSE, sep = "\t")

#################################################################################
# Create a ranked scatter plot
# Set working directory and read data

data <- read_excel("Tissue-specific-results-figures.xlsx", sheet = 1, col_names = TRUE)
# Extract highlighted x values
highlighted_x <- data %>% filter(highlight) %>% pull(x)
# Create data for vertical lines
lines_data <- data.frame(
  x = highlighted_x,
  y = -0.0001,  # Start y-coordinate for lines
  
  xend = highlighted_x,
  yend = -0.001  # End y-coordinate for lines
  
)

# Create and save the plot
pdf("Lymphoid_rank_plot.pdf", width = 5, height = 4)
ggplot(data, aes(x = x, y = y)) +
  geom_point(aes(color = highlight), size = 2) +
  scale_color_manual(values = c("FALSE" = "#3A70B7", "TRUE" = "#E83A37")) +
  theme_minimal() +
  labs(title = "Ranked Scatter Plot", x = "X Value", y = "Y Value") +
  geom_segment(data = lines_data, aes(x = x, y = y, xend = xend, yend = yend), color = "red", size = 0.000000001) +
  coord_cartesian(ylim = c(-0.002, 0.015)) +  # Adjust y-axis range to show short red lines
  
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    
    axis.line = element_line(color = "black"),  # Set axis lines to black
    
    axis.ticks = element_line(color = "black"),  # Set axis ticks to black
    
    axis.text = element_text(color = "black"),  # Set axis text to black
    
    axis.title = element_text(color = "black")  # Set axis titles to black
    
  )
dev.off()