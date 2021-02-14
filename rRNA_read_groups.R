####################################################################

# Parts of script shamelessly stolen from Felix Grunberger
# https://github.com/felixgrunberger/Native_RNAseq_Microbes
# (check out the guy, he is pretty cool)

####################################################################



####################################################################
#    Somehow make this pick the start and end of read together?    #
####################################################################

#------#
# Init #
#______#

#Get libraries (probs more than needed?)

library(here)
source(here("seq/R/load.libraries.R"))

# Load functions
# Bam file loader
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

# Get colour scheme
colours <- rev(c("#A3280A", "#8A5122", "#E0D253"))

#------#
# MAIN #
#______#

# Load data
# Genome (GFF, FASTA)
gff <- read.gff(here("seq/R/data/msmeg.gff")) %>%
  dplyr::filter(type == "rRNA")

fasta <- readDNAStringSet(here("seq/R/data/msmeg.fasta"))

# Get data and filter out everything over 9000 (looking at you, contaminating DNA)
RUN1_data <- get_bam_input("RUN1") %>%
  dplyr::filter(width<9000)

RUN2_data <- get_bam_input("RUN2") %>%
  dplyr::filter(width<9000)

all_data <- rbind(RUN1_data, RUN2_data)

# Filter data to matching rrn operon
all_data_f <- all_data %>%
  dplyr::filter(end < (gff$end[3] + 500),     # Between -400 from 5' start site of precursor...
                end > (gff$start[1] - 200),   # ...and 100+ sites past the TSS of 5S
                start < (gff$end[3] + 500),
                start > (gff$start[1] - 200))

# Define window size for processing sites (standard)
window_s = 3

# Define processing sites
range_1 <- ((gff$end[3] + 302 - window_s):(gff$end[3] + 302 + window_s)) # pre-START
range_2 <- ((gff$end[3] + 152 - window_s):(gff$end[3] + 152 + window_s)) # RNase III @ -152 (16S)
range_3 <- ((gff$end[3] - window_s):(gff$end[3] + window_s))             # 16S START
range_4 <- ((gff$start[3] - window_s):(gff$start[3] + 40 + window_s))    # 16S END
range_5 <- ((gff$end[2] + 217 - window_s):(gff$end[2] + 217 + window_s)) # RNase III @ -217 (23S)
range_6 <- ((gff$end[2] + 128 - window_s):(gff$end[2] + 156 + window_s)) # Any peaks in the region there?
range_7 <- ((gff$end[2] - window_s -3):(gff$end[2] + window_s +3))       # 23S START
range_8 <- ((gff$start[2] - window_s):(gff$start[2] + 70 + window_s))    # 23S END
range_9 <- ((gff$end[1] - window_s):(gff$end[1] + window_s))             # 5S START
range_0 <- ((gff$start[1] - window_s):(gff$start[1] + window_s))         # 5S END
range_x <- ((gff$end[3] + 310):(gff$end[3] + 500))                       # Anything before pre-START
range_y <- ((gff$start[1] - 200):(gff$start[1] - 10))                    # Anything past 5S TSS

# Add column about type of start site
all_data_c1 <- all_data_f %>%
  dplyr::mutate(fivegroup = ifelse((end+12) %in% range_1, "1",
                             ifelse((end+12) %in% range_2, "2",
                              ifelse((end+12) %in% range_3, "3",
                               ifelse((end+12) %in% range_4, "4",
                                ifelse((end+12) %in% range_5, "5",
                                 ifelse((end+12) %in% range_6, "6",
                                  ifelse((end+12) %in% range_7, "7",
                                   ifelse((end+12) %in% range_8, "8",
                                    ifelse((end+12) %in% range_9, "9",
                                     ifelse((end+12) %in% range_0, "0",
                                      ifelse((end+12) %in% range_x, "x", "rest")))))))))))) %>%
  dplyr::filter(fivegroup != "rest")


# Add column about type of end site
all_data_c2 <- all_data_c1 %>%
  dplyr::mutate(threegroup = ifelse(start %in% range_1, "1",
                              ifelse(start %in% range_2, "2",
                               ifelse(start %in% range_3, "3",
                                ifelse(start %in% range_4, "4",
                                 ifelse(start %in% range_5, "5",
                                  ifelse(start %in% range_6, "6",
                                   ifelse(start %in% range_7, "7",
                                    ifelse(start %in% range_8, "8",
                                     ifelse(start %in% range_9, "9",
                                      ifelse(start %in% range_0, "0",
                                       ifelse(start %in% range_y, "y", "rest")))))))))))) %>%
  dplyr::filter(threegroup != "rest")

# Calculate transcript group from start and end type
all_data_c3 <- all_data_c2 %>%
  dplyr::mutate(grouptype = paste(fivegroup,threegroup,sep=":")) %>%
  group_by(grouptype)

# Get numbers of elements in groups
all_data_groups <- all_data_c3 %>%
  group_by(grouptype) %>%
  summarise(n = n())

# Filter for reasonable amount of reads
count_limit=20
all_data_groups_f <- all_data_groups %>%
  dplyr::filter(n > count_limit)

# Export CSV to put in Excel for presentations
write_csv(all_data_groups, here("seq/R/Out/read_groups.csv"), col_names = TRUE )
write_csv(all_data_groups_f, here("seq/R/Out/read_groups_filtered.csv"), col_names = TRUE )

#------#
# PLOT #
#______#

# Plot 40 reads (like Felix's stuff)
#--------------------------------------

# Set 1: Pre-processed 16S (2:5)
set1 <- all_data_c3 %>%
  dplyr::filter(grouptype == "2:5") %>%
  mutate(set = "set1", set_type = "intermediate") %>%
  arrange(desc(width)) %>%
  head(40)

# Set 2: Mature 16S (3:4)
set2 <- all_data_c3 %>%
  dplyr::filter(grouptype == "3:4") %>%
  mutate(set = "set2", set_type = "mature") %>%
  arrange(desc(width)) %>%
  head(40)

# Set 3: Weird thing between 16S 3' and RNaseII site @ 1977  (4:5)
set3 <- all_data_c3 %>%
  dplyr::filter(grouptype == "4:5") %>%
  mutate(set = "set3", set_type = "by-product") %>%
  arrange(desc(width)) %>%
  head(40)

# Set 4: Mature 23S  (7:8)
set4 <- all_data_c3 %>%
  dplyr::filter(grouptype == "7:8") %>%
  mutate(set = "set4", set_type = "mature") %>%
  arrange(desc(width)) %>%
  head(40)

# Set 5: Mature 5S  (9:0)
set5 <- all_data_c3 %>%
  dplyr::filter(grouptype == "9:0") %>%
  mutate(set = "set5", set_type = "mature") %>%
  arrange(desc(width)) %>%
  head(40)

# Set 6: Extended 5S  (9:y)
set6 <- all_data_c3 %>%
  dplyr::filter(grouptype == "9:y") %>%
  mutate(set = "set6", set_type = "intermediate") %>%
  arrange(desc(width)) %>%
  head(40)

# SETS ONLY AFTER THE SECOND RUN
# Only do 20

set1_t <- set1 %>%
  arrange(desc(width)) %>%
  head(20)

set2_t <- set2 %>%
  arrange(desc(width)) %>%
  head(20)

set3_t <- set3 %>%
  arrange(desc(width)) %>%
  head(20)

set4_t <- set4 %>%
  arrange(desc(width)) %>%
  head(20)

set5_t <- set5 %>%
  arrange(desc(width)) %>%
  head(20)

set6_t <- set6 %>%
  arrange(desc(width)) %>%
  head(20)

# Set 7: 16S 5' to 5S 3' (3:0)
set7 <- all_data_c3 %>%
  dplyr::filter(grouptype == "3:0") %>%
  mutate(set = "set7", set_type = "intermediate") %>%
  arrange(desc(width)) %>%
  head(20)

# Set 8: 16S 5' to 23S 3' (3:8)
set8 <- all_data_c3 %>%
  dplyr::filter(grouptype == "3:8") %>%
  mutate(set = "set8", set_type = "intermediate") %>%
  arrange(desc(width)) %>%
  head(20)

# Set 9: 16S 3' to 23S 3' (4:8)
set9 <- all_data_c3 %>%
  dplyr::filter(grouptype == "4:8") %>%
  mutate(set = "set9", set_type = "intermediate") %>%
  arrange(desc(width)) %>%
  head(20)

# Set 10: RNase III (2040) to 23S 3' (6:8)
set10 <- all_data_c3 %>%
  dplyr::filter(grouptype == "6:8") %>%
  mutate(set = "set10", set_type = "intermediate") %>%
  arrange(desc(width)) %>%
  head(20)

# Combine
all_sets <- rbindlist(list(set1, set2, set3, set4, set5, set6), fill = T) %>%
  arrange(end) %>%
  group_by(set) %>%
  arrange(desc(width)) %>%
  ungroup() %>%
  mutate(lfd_n = 1:n()) 

all_sets_20 <- rbindlist(list(set1_t, set2_t, set3_t, set4_t, set5_t, set6_t, set7, set8, set9, set10), fill = T) %>%
  arrange(end) %>%
  group_by(set) %>%
  arrange(desc(width)) %>%
  ungroup() %>%
  mutate(lfd_n = 1:n()) 

# Plot
pdf(here("seq/R/Out/coocurring_reads.pdf"), 
    width = 16, height = 10, paper = "special",onefile=FALSE)
ggplot(data = all_sets) +
  geom_segment(aes(x = start - 12, xend = end, yend = lfd_n, y = lfd_n, color = set_type), size = 0.7) +
  scale_x_reverse() +
  theme_Publication_white() +
  ylab("") +
  xlab("") +
  ggtitle("") +
  scale_color_manual(values = colours) +
  geom_vline(xintercept = c(gff$start[1],gff$end[1],
                            gff$start[2],gff$end[2],
                            gff$start[3],gff$end[3]), linetype = "dashed", alpha = 0.5) +
  annotate(geom = "text",
           label = c("5S end","5S start","23S end","23S start","16S end","16S start"),
           x = c((gff$start[1]-20),(gff$end[1]-20),(gff$start[2]-20),(gff$end[2]-20),(gff$start[3]-20),(gff$end[3]-20)),
           y = 260,
           angle = 90,
           vjust = 1) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_blank()) 

dev.off()

pdf(here("seq/R/Out/coocurring_reads_2runs.pdf"), 
    width = 16, height = 10, paper = "special",onefile=FALSE)
ggplot(data = all_sets_20) +
  geom_segment(aes(x = start - 12, xend = end, yend = lfd_n, y = lfd_n, color = set_type), size = 0.7) +
  scale_x_reverse() +
  theme_Publication_white() +
  ylab("") +
  xlab("") +
  ggtitle("") +
  scale_color_manual(values = colours) +
  geom_vline(xintercept = c(gff$start[1],gff$end[1],
                            gff$start[2],gff$end[2],
                            gff$start[3],gff$end[3]), linetype = "dashed", alpha = 0.5) +
  annotate(geom = "text",
           label = c("5S end","5S start","23S end","23S start","16S end","16S start"),
           x = c((gff$start[1]-20),(gff$end[1]-20),(gff$start[2]-20),(gff$end[2]-20),(gff$start[3]-20),(gff$end[3]-20)),
           y = 220,
           angle = 90,
           vjust = 1) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_blank()) 

dev.off()

# Plot histogram for ratio between species??
  
pdf(here("seq/R/Out/ratio_histogram.pdf"), 
    width = 16, height = 10, paper = "special",onefile=FALSE)
ggplot(data = all_data_groups_f,aes(x=grouptype,y=n)) +
  geom_bar(stat="identity", width = 0.5, color = "black", fill = colours[2]) +
  theme_Publication_white() +
  scale_x_discrete(labels = c("a","b","c","d","e","f","g","h","i","j")) +
  scale_y_log10(expand = c(0,0),
                breaks = c(0,10,100,1000,10000)) +
  ylab("") +
  xlab("") +
  ggtitle("") +
  theme(axis.text.y = element_text(face = "italic"), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())

dev.off()


















