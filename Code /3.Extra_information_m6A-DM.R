#data collection for CGCgenes,meRIP-peaks and high-resolution m6A data
library(readr)
library(dplyr)
library(GenomicRanges)
library(vcfR)
#################################################################################################################################
#meRIP-peak length for genes
# Create new columns "Seqname", "Start", and "End"
meRIP_peak <- data.frame(meRIP_peak, stringsAsFactors = FALSE)
meRIP_peak <- separate(meRIP_peak, col = meRIP_peak, into = c("Seqname", "Start", "End"),sep = "[:-]")
# Read the hg38 genome annotation file (GTF format)
gtf <- read.delim("hg38.gtf", header = FALSE, comment.char = "#", sep = "\t", stringsAsFactors = FALSE)
# Extract the CDS positions
cds <- gtf[gtf$V3 == "CDS", c("V1", "V4", "V5")]
# Create GenomicRanges objects
peak_gr <- with(meRIP_peak, GRanges(Chromosome = Seqname, IRanges(Start, End)))
cds_gr <- with(cds, GRanges(Chromosome = V1, IRanges(V4, V5)))
# Filter peaks on CDS
cds_peaks <- subsetByOverlaps(peak_gr, cds_gr)
# Calculate the length of each peak on CDS
peak_lengths <- width(cds_peaks)
# Create a new data frame with peak lengths for each gene CDS
cds_peak_lengths <- data.frame(Gene = unique(cds$V9), PeakLength = tapply(peak_lengths, factor(as.character(queryHits(cds_peaks))), sum))
cds_peak_lengths <- unique(cds_peak_lengths)
cds_peak_lengths$Gene <- sub('.*?"gene_name":"(.*?)".*', '\\1', cds_peak_lengths$Gene)
# Calculate the total length of peaks on the entire CDS for each gene
gene_peak_lengths <- aggregate(PeakLength ~ Gene, data = cds_peak_lengths, FUN = sum)

##############################################################################################################################################
# Extra calculation of the mutational features
mut <- read.table("CosmicMutantExport.tsv", header=T, sep="\t")
mut <- data.frame(mut)
mut <- subset(mut,grepl("^.*y.*$",Genome.wide.screen))
synonymous_mutation <- subset(mut,mut$Mutation.Description=="Substitution - coding silent")
missense_mutation<- subset(mut,mut$Mutation.Description=="Substitution - Missense")
synonymous_mutation <- distinct(synonymous_mutation, HGVSG, ID_sample, .keep_all= TRUE)
missense_mutation <- distinct(missense_mutation, HGVSG, ID_sample, .keep_all= TRUE)

# Count for #synonymous mutations per gene and #missense mutations per gene
count_result_synonymous_mutation <- count(synonymous_mutation, HGNC.ID)
count_result_missense_mutation <- count(missense_mutation, HGNC.ID)

#Count the number of synonymous and missense mutations
df_syn <- distinct(synonymous_mutation, HGVSG, .keep_all= TRUE)
nrow(df_syn)
df_mis <- distinct(missense_mutation, HGVSG, .keep_all= TRUE)
nrow(df_mis)

#Count for the genome-wide screening samples
df_sample <- distinct(a, ID_sample, .keep_all= TRUE)
nrow(df_sample)

#Count for mutation type
df_muttype <- distinct(a, Mutation.Description, .keep_all= TRUE)
df_muttype

#Count for CNV data in TCGA
cnv <- read.table("CosmicCompleteCNA.tsv", header = TRUE, sep = "\t", na.strings = TRUE, fill = TRUE) 
cnv <- data.frame(cnv)
cnv <- subset(cnv,grepl("^.*TCGA.*$",SAMPLE_NAME))
cnv2 <- distinct(cnv, CNV_ID, ID_SAMPLE, .keep_all= TRUE)
cnv
cnv_countpersample <- count(cnv2, ID_SAMPLE)

#Collection of CNVs within CGC genes
CGC <- CGC_set$GENE_SYMBOL
CGC<-as.vector(as.matrix(CGC))
cgc_cnv <- cnv %>% filter(cnv$gene_name %in% CGC)
############################################################################################################################
#collection of gnomAD data, A-to-I data and WER-data
#AF information from gnomAD4.0.0
gnomad.exomes.r4.0.0.sites.vcf <-"bcftools annotate -x ^INFO/AF gnomad.exomes.r4.0.0.sites.vcf.bgz"
# Read the VCF file
vcf_data <- readVcf("gnomad.exomes.r4.0.0.sites.vcf")
vcf_data$CHROM <- gsub("^chr", "", vcf_data$CHROM)
vcf_data$mutation_iden <- paste(vcf_data$CHROM, vcf_data$POS, vcf_data$REF, vcf_data$ALT, sep = "-")#chr-pos-ref-alt
gnomAD_data <- vcf_data[, c("mutation_iden", "INFO")]#chr-pos-ref-alt AF

#WER-data
data <- read.table("hg38_Human_m6A_ClipSeqRipSeqResult.txt", header = F)
WER_association <- unique(data[, c(3, 4)])

#A-to-I data
data <- read.table("TABLE1_hg38.txt", header = F)
atoi_edit_pos <- unique(data[, 1:2])
atoi_edit_pos$V1 <- gsub("chr", "", atoi_edit_pos$V1)
atoi_edit_pos$pos <- paste0(atoi_edit_pos$V1, ":", atoi_edit_pos$V2, "-", atoi_edit_pos$V2)


