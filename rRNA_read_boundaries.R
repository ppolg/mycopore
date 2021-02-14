####################################################################

# Most (or whole) of script shamelessly stolen from Felix Grunberger
# https://github.com/felixgrunberger/Native_RNAseq_Microbes
# (check out the guy, he is pretty cool)

####################################################################

####################################################################
#     Somehow plot the 5' and 3' end of reads compared to map?     #
####################################################################

#------#
# Init #
#______#

# Get libraries (probs more than needed?)

library(here)
source(here("seq/R/load.libraries.R"))

# Load functions
# Load bam file
get_bam_input <- function(org){
  input_bam_file <- paste(here::here("seq/R/data/"),org, ".bam", sep = "")
  allReads <- readGAlignments(input_bam_file, use.names = T, param = ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE), tag=c("NM"), what=c("mapq", "flag")))
  allReads_table <- GenomicAlignments::as.data.frame(allReads) %>%
    mutate(minion_read_name = names(allReads)) 
  
  # > get read sequences
  param <- ScanBamParam(
    flag=scanBamFlag(isUnmappedQuery=FALSE),
    what="seq")
  res0 <- scanBam(input_bam_file,param = param)[[1]] # always list-of-lists
  allReads_sequence <- res0[["seq"]]                 # query widths
  allReads_sequence_table <- as.list(as.character(allReads_sequence))
  allReads_table$sequence <- unlist(allReads_sequence_table)
  allReads_table$n_char   <- nchar(allReads_table$sequence[1:length(allReads_table$sequence)])
  left  <- paste(str_split_fixed(string = allReads_table$cigar, pattern = "M", n = 2)[,1],"M", sep = "")
  right <- paste(str_split_fixed(string = allReads_table$cigar, pattern = "M", n = 2)[,2],"1M", sep = "")
  
  #................................calculate cigar tables / SOFT AND HARD CLIPPING!!!
  allReads_table$soft_l <- as_tibble(cigarOpTable(left))$S
  allReads_table$hard_l <- as_tibble(cigarOpTable(left))$H
  allReads_table$soft_r <- as_tibble(cigarOpTable(right))$S
  allReads_table$hard_r <- as_tibble(cigarOpTable(right))$H
  return(allReads_table)
  
}

# Set plotting colours
two_color_npg <- rev(c("#A3280A", "#E3812B"))

# Modify the bam table with relative position to rRNA start and end sites
# One of the toe rrn operons are currently picked manually
# Since Msmeg rrn is on the negative strand, end and start need to be switched
# Circularisation is not a thing, but due to Nanopore's shortcomings need to move by 12

mod_bam_all <- function(input_table, set_name){
  return(input_table %>%
           mutate(UTR5_16 = ifelse(hard_l == 0 & hard_r == 0 & soft_l < 500, gff$end[3] - (end + 12),gff$end[3] - (end)),
                  UTR3_16 = (gff$start[3] - start)) %>%
           mutate(UTR5_23 = ifelse(hard_l == 0 & hard_r == 0 & soft_l < 500, gff$end[2] - (end + 12),gff$end[2] - (end)),
                  UTR3_23 = (gff$start[2] - start)) %>%
           mutate(UTR5_5 = ifelse(hard_l == 0 & hard_r == 0 & soft_l < 500, gff$end[1] - (end + 12),gff$end[1] - (end)),
                  UTR3_5 = (gff$start[1] - start)) %>%
           mutate(sample = set_name)
  )
}

#------#
# MAIN #
#______#

# Load data
# Genome (GFF, FASTA)
gff <- read.gff(here("seq/R/data/msmeg.gff")) %>%
  dplyr::filter(type == "rRNA")

fasta <- readDNAStringSet(here("seq/R/data/msmeg.fasta"))

# Get data
RUN1_data <- get_bam_input("RUN1")
RUN2_data <- get_bam_input("RUN2")
# Modify data to include relative position to rrn
RUN1_data_mod  <- mod_bam_all(RUN1_data, "RUN1")
RUN2_data_mod  <- mod_bam_all(RUN2_data, "RUN2")

# Filter to remove long DNA contamination
RUN1_data_f <- RUN1_data_mod %>%
  dplyr::filter(width<9000)

RUN2_data_f <- RUN2_data_mod %>%
  dplyr::filter(width<9000)

# Filter to remove reads not mapping rrn
RUN1_data_f2 <- RUN1_data_f %>%
  dplyr::filter(end < (gff$end[3] + 400),     # Between -400 from 5' start site of precursor...
                end > (gff$start[1] - 100),   # ...and 100+ sites past the TSS of 5S
                start < (gff$end[3] + 400),
                start > (gff$start[1] - 100))

RUN2_data_f2 <- RUN2_data_f %>%
  dplyr::filter(end < (gff$end[3] + 400),    
                end > (gff$start[1] - 100),   
                start < (gff$end[3] + 400),
                start > (gff$start[1] - 100))

# Merge runs
all_data <- rbind(RUN1_data_f2, RUN2_data_f2)

# Order data by run
all_data$sample <-  factor(all_data$sample, levels = rev(c("RUN1", "RUN2")))

#------#
# Plot #
#______#

# Plot 16S TSS
pdf(here("seq/R/Out/16S_TSS.pdf"), 
    width = 8, height = 5, paper = "special",onefile=FALSE)
ggplot(data = all_data, aes(x = UTR5_16 , fill = sample, y = sample)) +
  geom_density_ridges2(aes(height = ..ndensity..),alpha = 1, size = 0.1, scale = 0.95,
                       stat = "binline", draw_baseline = FALSE, bins = 150, color = NA) +
  scale_x_continuous(limits = c(-310,20), 
                     breaks = c(-302,-152,0),
                     labels = c("precursor TSS", "RNase III", 0),
                     expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_Publication_white() +
  ylab("") +
  xlab("5' end [nt]") +
  ggtitle("") +
  scale_fill_manual(values = two_color_npg) +
  scale_color_manual(values = two_color_npg) +
  guides(fill = F, color = F) +
  theme(axis.text.y = element_text(face = "italic")) 
dev.off() 

# Plot 16S TTS
pdf(here("seq/R/Out/16S_TTS.pdf"), 
    width = 8, height = 5, paper = "special",onefile=FALSE)
ggplot(data = all_data, aes(x = UTR3_16 , fill = sample, y = sample)) +
  geom_density_ridges2(aes(height = ..ndensity..),alpha = 1, size = 0.1, scale = 0.95, stat = "binline", draw_baseline = FALSE, bins = 200, color = NA) +
  scale_x_continuous(limits = c(-100,100), 
                     breaks = c(-50,0,50), 
                     expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_Publication_white() +
  ylab("") +
  xlab("3' end [nt]") +
  ggtitle("") +
  scale_fill_manual(values = two_color_npg) +
  scale_color_manual(values = two_color_npg) +
  guides(fill = F, color = F) +
  theme(axis.text.y = element_text(face = "italic")) 
dev.off()

# Plot 23S TSS
pdf(here("seq/R/Out/23S_TSS.pdf"), 
    width = 8, height = 5, paper = "special",onefile=FALSE)
ggplot(data = all_data, aes(x = UTR5_23 , fill = sample, y = sample)) +
  geom_density_ridges2(aes(height = ..ndensity..),alpha = 1, size = 0.1, scale = 0.95, stat = "binline", draw_baseline = FALSE, bins = 150, color = NA) +
  scale_x_continuous(limits = c(-250,20), 
                     breaks = c(-217,-156,-128,0), 
                     labels = c("RNase III",-156,-128,0),
                     expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_Publication_white() +
  ylab("") +
  xlab("5' end [nt]") +
  ggtitle("") +
  scale_fill_manual(values = two_color_npg) +
  scale_color_manual(values = two_color_npg) +
  guides(fill = F, color = F) +
  theme(axis.text.y = element_text(face = "italic")) 
dev.off()

# Plot 23S TTS
pdf(here("seq/R/Out/23S_TTS.pdf"), 
    width = 8, height = 5, paper = "special",onefile=FALSE)
ggplot(data = all_data, aes(x = UTR3_23 , fill = sample, y = sample)) +
  geom_density_ridges2(aes(height = ..ndensity..),alpha = 1, size = 0.1, scale = 0.95, stat = "binline", draw_baseline = FALSE, bins = 150, color = NA) +
  scale_x_continuous(limits = c(-100,100), 
                     breaks = c(-50, 0, 50), 
                     expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_Publication_white() +
  ylab("") +
  xlab("3' end [nt]") +
  ggtitle("") +
  scale_fill_manual(values = two_color_npg) +
  scale_color_manual(values = two_color_npg) +
  guides(fill = F, color = F) +
  theme(axis.text.y = element_text(face = "italic")) 
dev.off()

# Plot 5S TSS
pdf(here("seq/R/Out/5S_TSS.pdf"), 
    width = 8, height = 5, paper = "special",onefile=FALSE)
ggplot(data = all_data, aes(x = UTR5_5 , fill = sample, y = sample)) +
  geom_density_ridges2(aes(height = ..ndensity..),alpha = 1, size = 0.1, scale = 0.95, stat = "binline", draw_baseline = FALSE, bins = 150, color = NA) +
  scale_x_continuous(limits = c(-100,50), 
                     breaks = c(-50,0), 
                     expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_Publication_white() +
  ylab("") +
  xlab("5' end [nt]") +
  ggtitle("") +
  scale_fill_manual(values = two_color_npg) +
  scale_color_manual(values = two_color_npg) +
  guides(fill = F, color = F) +
  theme(axis.text.y = element_text(face = "italic")) 
dev.off()

# Plot 5S TTS
pdf(here("seq/R/Out/5S_TTS.pdf"), 
    width = 8, height = 5, paper = "special",onefile=FALSE)
ggplot(data = all_data, aes(x = UTR3_5 , fill = sample, y = sample)) +
  geom_density_ridges2(aes(height = ..ndensity..),alpha = 1, size = 0.1, scale = 0.95, stat = "binline", draw_baseline = FALSE, bins = 150, color = NA) +
  scale_x_continuous(limits = c(-100,100), 
                     breaks = c(-50, 0, 50), 
                     expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_Publication_white() +
  ylab("") +
  xlab("3' end [nt]") +
  ggtitle("") +
  scale_fill_manual(values = two_color_npg) +
  scale_color_manual(values = two_color_npg) +
  guides(fill = F, color = F) +
  theme(axis.text.y = element_text(face = "italic")) 
dev.off()

# Plot the region cause where did all my reads go? :(
pdf(here("seq/R/Out/Allreads_TSS.pdf"), 
    width = 8, height = 5, paper = "special",onefile=FALSE)
ggplot(data = all_data, aes(x =  end, fill = sample, y = sample)) +
  geom_density_ridges2(aes(height = ..ndensity..),alpha = 1, size = 0.1, scale = 0.95, stat = "binline", draw_baseline = FALSE, bins = 150, color = NA) +
  scale_y_discrete(expand = c(0,0)) +
  theme_Publication_white() +
  scale_x_reverse(breaks = NULL) +
  geom_vline(xintercept = c(gff$start[1],gff$end[1],
                            gff$start[2],gff$end[2],
                            gff$start[3],gff$end[3]), linetype = "dashed", alpha = 0.5) +
  annotate(geom = "text",
           label = c("5S end","5S start","23S end","23S start","16S end","16S start"),
           x = c((gff$start[1]-20),(gff$end[1]-20),(gff$start[2]-20),(gff$end[2]-20),(gff$start[3]-20),(gff$end[3]-20)),
           y = 2.7,
           size = 2.5,
           angle = 90,
           vjust = 1) +
  ylab("") +
  xlab("5' end [nt]") +
  ggtitle("") +
  scale_fill_manual(values = two_color_npg) +
  scale_color_manual(values = two_color_npg) +
  guides(fill = F, color = F) +
  theme(axis.text.y = element_text(face = "italic")) 
dev.off()

pdf(here("seq/R/Out/Allreads_TTS.pdf"), 
    width = 8, height = 5, paper = "special",onefile=FALSE)
ggplot(data = all_data, aes(x =  start - 12, fill = sample, y = sample)) +
  geom_density_ridges2(aes(height = ..ndensity..),alpha = 1, size = 0.1, scale = 0.95, stat = "binline", draw_baseline = FALSE, bins = 150, color = NA) +
  scale_y_discrete(expand = c(0,0)) +
  theme_Publication_white() +
  scale_x_reverse(breaks = NULL) +
  geom_vline(xintercept = c(gff$start[1],gff$end[1],
                            gff$start[2],gff$end[2],
                            gff$start[3],gff$end[3]), linetype = "dashed", alpha = 0.5) +
  annotate(geom = "text",
           label = c("5S end","5S start","23S end","23S start","16S end","16S start"),
           x = c((gff$start[1]-20),(gff$end[1]-20),(gff$start[2]-20),(gff$end[2]-20),(gff$start[3]-20),(gff$end[3]-20)),
           y = 2.7,
           size = 2.5,
           angle = 90,
           vjust = 1) +
  ylab("") +
  xlab("3' end [nt]") +
  ggtitle("") +
  scale_fill_manual(values = two_color_npg) +
  scale_color_manual(values = two_color_npg) +
  guides(fill = F, color = F) +
  theme(axis.text.y = element_text(face = "italic")) 
dev.off()

# Plot histogram for read length

pdf(here("seq/R/Out/RUN1_lengths.pdf"), 
    width = 8, height = 5, paper = "special",onefile=FALSE)
ggplot(data = RUN1_data_f2) +
  geom_bar(aes(x =  width),alpha = 1, size = 0.1, color = two_color_npg[2]) +
  scale_x_continuous(limits = c(0,5000), 
                     breaks = c(0,1000,2000,3000,4000), 
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(0,3000),
                     breaks = c(0,3000),
                     expand = c(0,0)) +
  theme_Publication_white() +
  ylab("RUN1") +
  xlab("") +
  ggtitle("") +
  theme(axis.text.y = element_text(face = "italic"))
dev.off()

pdf(here("seq/R/Out/RUN2_lengths.pdf"), 
    width = 8, height = 5, paper = "special",onefile=FALSE)
ggplot(data = RUN2_data_f2) +
  geom_bar(aes(x =  width),alpha = 1, size = 0.1, color = two_color_npg[1]) +
  scale_x_continuous(limits = c(0,5000), 
                     breaks = c(0,1000,2000,3000,4000), 
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(0,3000),
                     breaks = c(0,3000),
                     expand = c(0,0)) +
  theme_Publication_white() +
  ylab("RUN2") +
  xlab("read length [nt]") +
  ggtitle("") +
  theme(axis.text.y = element_text(face = "italic"))
dev.off()
