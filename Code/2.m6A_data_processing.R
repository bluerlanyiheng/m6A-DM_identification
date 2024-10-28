library(readr)
library(dplyr)
library(GenomicRanges)
library(vcfR)
#################################################################################################################################
# Data collection and processing for high-resolution m6A-data & meRIP-peak data

cancer_gene_census <- read_tsv("Cosmic_CancerGeneCensus_v98_GRCh38.tsv")

# Filter the rows where "TIER" column is equal to "1" and keep all the information
CGC_set <- cancer_gene_census %>% filter(TIER == "1")

#high-resolution m6A-data
# Read the file
m6A_data <- read_tsv("m6A-Atlas2_HighRes_human.txt")

# Remove the "chr" string prefix from each value in the "Seqname" column
m6A_data$Seqname <- gsub("chr", "", m6A_data$Seqname)

# Create a new column "m6A_position"
m6A_data <- mutate(m6A_data, m6A_position = paste(Seqname, Position, sep = ":"))

# Append "-Position" to each value in the "m6A_position" column
m6A_data$m6A_position <- paste(m6A_data$m6A_position, "-", m6A_data$Position, sep = "")

# Assign the column data to the variable "HighRes_m6A"
HighRes_m6A <- m6A_data$m6A_position

#meRIP-peak data
# Read the file
m6A_data <- read_tsv("m6A=sites=species=human=hg38_REPIC.txt")

# Remove the "chr" string prefix from each value in the "pos" column
m6A_data$pos <- gsub("chr", "", m6A_data$pos)
m6A_data$pos <- gsub("[-]|[+]", "", m6A_data$pos)

# Output the "pos" column data to the variable "meRIP_peak"
meRIP_peak <- m6A_data$pos

# Merge HighRes_m6A and meRIP_peak into m6A_range
m6A_range <- union(HighRes_m6A, meRIP_peak)

#################################################################################################################################
#Merge m6A data into non-redundant peaks by tissues

bed_files <- list.files(path, pattern = "*.bed", full.names = TRUE)

for (file in bed_files) {
  bed <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  # Extract the file name (excluding the path and extension)
  file_name <- gsub("\\.bed$", "", basename(file))
  # Store the BED data in a separate object, named after the file name
  assign(file_name, bed)
}
# Get all object names
object_names <- gsub("\\.bed$", "", basename(bed_files))
# For each BED data, merge peaks, calculate peak lengths, and output the merged BED dat
for (file_name in object_names) {
  # Get the BED data
  bed <- get(file_name)
  # Convert the data frame to a GRanges object, specifying the columns for seqnames, start, and end
  gr <- makeGRangesFromDataFrame(bed, 
                                 keep.extra.columns = TRUE, 
                                 seqnames.field = "V1",  # Assume the seqnames column is the first column
                                 
                                 start.field = "V2",     # Assume the start column is the second column
                                 
                                 end.field = "V3")       # Assume the end column is the third column
  
  
  # Merge peaks
  merged_bed <- reduce(gr)
  # Calculate peak lengths
  peak_lengths <- width(merged_bed)
  # Add peak lengths to a new column
  merged_bed$length <- peak_lengths
  # Convert the GRanges object back to a data frame
  merged_df <- as.data.frame(merged_bed)
  # Construct the output file name
  output_file <- file.path(path, paste0(file_name, "_merged.bed"))
  # Output the merged BED data
  write.table(merged_df, file = output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  # Count the number of lines in the merged BED file
  num_lines <- length(readLines(output_file))
  # Print the number of lines
  cat(paste0("Number of lines in ", output_file, ": ", num_lines, "\n"))
}

###############################################################################################################################################
##merge peaks by using all tissues

merged_bed_files <- list.files(path, pattern = "_merged.bed$", full.names = TRUE)

all_merged_peaks <- GRanges()

for (file in merged_bed_files) {

  bed <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  # Convert the data frame to a GRanges object
  gr <- makeGRangesFromDataFrame(bed, 
                                 keep.extra.columns = TRUE, 
                                 seqnames.field = "V1",  # Assume the seqnames column is the first column
                                 start.field = "V2",     # Assume the start column is the second column
                                 end.field = "V3")       # Assume the end column is the third column
  # Combine the current GRanges object with the all_merged_peaks GRanges object
  all_merged_peaks <- c(all_merged_peaks, gr)
}
class(all_merged_peaks)
# Merge all peaks
all_merged_peaks <- reduce(all_merged_peaks)
# Calculate peak lengths
peak_lengths <- width(all_merged_peaks)
# Add peak lengths to a new column
all_merged_peaks$length <- peak_lengths
# Convert the GRanges object back to a data frame
all_merged_df <- as.data.frame(all_merged_peaks)
# Construct the output file name
output_file <- file.path(path, "all_merged.bed")
#Output the merged BED data
write.table(all_merged_df, file = output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
# Count the number of lines in the all_merged.bed file
num_lines <- length(readLines(output_file))
# Print the number of lines
cat(paste0("Number of lines in ", output_file, ": ", num_lines, "\n"))


#############################################################################################################################
##counting overlapping peaks by tissues

overlap_counts <- data.frame(peak_id = seq_along(all_merged_peaks), stringsAsFactors = FALSE)
for (file in merged_bed_files) {
  tissue_name <- gsub("\\.bed$", "", basename(file))
  tissue_name <- gsub("_merged", "", tissue_name)
  overlap_counts[[tissue_name]] <- 0
}

for (file in merged_bed_files) {
  # Read the _merged.bed file
  bed <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  # Convert the data frame to a GRanges object
  gr_tissue <- makeGRangesFromDataFrame(bed, 
                                        keep.extra.columns = TRUE, 
                                        seqnames.field = "V1",  # Assume the seqnames column is the first column
                                        
                                        start.field = "V2",     # Assume the start column is the second column
                                        
                                        end.field = "V3")       # Assume the end column is the third column
 
  overlaps <- findOverlaps(all_merged_peaks, gr_tissue)
 
# Count overlaps
  tissue_name <- gsub("\\.bed$", "", basename(file))
  tissue_name <- gsub("_merged", "", tissue_name)
  overlap_counts[[tissue_name]][queryHits(overlaps)] <- 1
  
}

overlap_counts$sum <- rowSums(overlap_counts[, -1])

output_file_overlap_counts <- file.path(path, "all_merged_overlap_counts.txt")

write.table(overlap_counts, file = "output_file_overlap_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

num_lines_overlap_counts <- length(readLines(output_file_overlap_counts))
cat(paste0("Number of lines in ", output_file_overlap_counts, ": ", num_lines_overlap_counts, "\n"))



















