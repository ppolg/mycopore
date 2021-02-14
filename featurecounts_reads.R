####################################################################

# Most of script shamelessly stolen from Felix Grunberger
# https://github.com/felixgrunberger/Native_RNAseq_Microbes
# (check out the guy, he is pretty cool)

####################################################################



####################################################################
#   Somehow generate nice QC plots from the featurecounts output?  #
####################################################################

#------#
# Init #
#______#

# Get libraries (probs more than needed?)

library(here)
source(here("seq/R/load.libraries.R"))

# Load functions

category_calculator <- function(id_table, identifier){
  dataset <- id_table %>%
    mutate(group = ifelse(mapped_type == "CDS", "CDS", 
                          ifelse(mapped_type == "rRNA", locus_name, NA))) %>% 
    dplyr::filter(!is.na(group)) %>%
    distinct(minion_read_name, .keep_all = T) %>%
    mutate(total = n()) %>%
    group_by(group) %>%
    summarise(percentage = n()/max(total)*100,
              total_count = n()) %>% 
    mutate(sequencing_set = identifier)
  
  return(dataset)
}


calc_stats <- function(id_table, identifier){
  dataset <- id_table %>%
    mutate(group = ifelse(mapped_type == "CDS", "CDS", 
                          ifelse(mapped_type == "rRNA", locus_name, NA))) %>% 
    dplyr::filter(!is.na(group)) %>%
    group_by(minion_read_name) %>%
    mutate(max_identity = max(identity)) %>%
    dplyr::filter(identity == max(identity)) %>%
    ungroup() %>%
    distinct(minion_read_name, .keep_all = T)  %>% 
    mutate(sequencing_set = identifier)
  
  return(dataset)
}

# Make init data.frames
gene_table_counts <- data.frame()
stats_counts <- data.frame()

# Colours
colour4 <-rev(c("#E3812B","#A3280A","#8A5122","#E0D253"))
colour4.2 <-rev(c("#E3812B","#A3280A","#4A749E","#612882"))
colour4.3 <-rev(c("#195928","#E3812B","#4A749E","#612882"))

colour2 <-rev(c("#858482", "#A3280A"))

heat_color_npg <- c(pal_npg()(10)[4],
                    pal_npg()(10)[6],
                    pal_npg()(10)[7],
                    pal_npg()(10)[1])

two_c_npg <- c("grey70",pal_npg()(10)[1])

#------#
# Main #
#______#

# Get gene files
gene_files <- paste(here("seq/R/data/tidy_data/"), list.files(here("seq/R/data/tidy_data/"),pattern = "_id_table"), sep = "")

# Get sample names
sample_names <- unlist(lapply(gene_files, FUN=function(x){str_split_fixed(str_split_fixed(x, "_id", 2)[1],"tidy_data/",2)[2]}))

# Read files and modify them

for (i in seq_along(sample_names)){
  
  # > set sample fasta, gff and input gene table file with counts information
  working_directory   <- paste(here(), "/seq/R/data", sep = "")
  sample_name         <- sample_names[i]
  input_fasta         <- paste(working_directory, "/msmeg.fasta", sep = "")
  input_gff           <- paste(working_directory, "/msmeg3.gff", sep = "")
  input_gene_table    <- gene_files[i]
  
  # > load saved R file
  load(input_gene_table)
  
  # > add gff information to count tables
  gene_table_counts <- rbind(gene_table_counts, category_calculator(full_id_table, identifier = sample_name))
  stats_counts <- rbind(stats_counts, calc_stats(full_id_table, identifier = sample_name))
}

# Reorder
gene_table_counts$sequencing_set <-  factor(gene_table_counts$sequencing_set, 
                                            levels = rev(c("RUN1","RUN2")))

#Write tsv

fwrite(gene_table_counts, file = here("seq/R/Out/gene_table_counts.tsv"))

#------#
# Plot #
#______#


# Mapped percentages
pdf(here("seq/R/Out/mapped_percent.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = gene_table_counts, aes(x = sequencing_set, y = percentage, fill = group)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = colour4.3) +
  theme_Publication_white() +
  xlab("") +
  ylab("Ratio of reads mapped to features (%)") +
  ggtitle("") +
  labs(fill = "") +
  theme(axis.text.y = element_text(face = "italic")) +
  theme(panel.grid.major.y = element_blank()) +
  scale_y_continuous(expand = c(0, 0))
dev.off()

# Mapping boxplot lengths
pdf(here("seq/R/Out/boxplot_length.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = stats_counts, aes(x = sequencing_set, y = aligned_reads, fill = group)) +
  geom_boxplot(outlier.color = NA, alpha = 1, notch = F) +
  theme_Publication_white() +
  xlab("") +
  ylab("Read length [nt]") +
  ggtitle("") +
  scale_y_continuous(limits = c(0,3000)) +
  scale_fill_manual(values = colour4.3) +
  theme(axis.text.y = element_text(face = "italic")) +
  theme(panel.grid.major.x = element_blank()) 
dev.off()

# Boxplot identity
pdf(here("seq/R/Out/boxplot_identity.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = stats_counts, aes(x = sequencing_set, y = identity, fill = group == "CDS")) +
  geom_half_boxplot(data = subset(stats_counts, group == "CDS"),position = position_nudge(x = .05, y = 0), 
                    outlier.color = NA,width = 0.4, side = "r", errorbar.draw = F) +
  geom_half_boxplot(data = subset(stats_counts, group != "CDS"),position = position_nudge(x = -.05, y = 0), 
                    outlier.color = NA,width = 0.4, side = "l", errorbar.draw = F) +
  scale_y_continuous(limits = c(70,100), expand = c(0,0)) +
  theme_Publication_white() +
  xlab("") +
  ylab("Mapping identitiy (%)") +
  ggtitle("") +
  scale_fill_manual(values = colour2) +
  theme(axis.text.y = element_text(face = "italic")) +
  theme(panel.grid.major.x = element_blank()) 
dev.off()

# Total mapped
pdf(here("seq/R/Out/total_mapped.pdf"),  
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = gene_table_counts, aes(x = sequencing_set, y = total_count, fill = group, color = group)) +
  geom_bar(stat="identity", 
           position=position_dodge(width=0.5), width=0.1) +
  geom_point(position=position_dodge(width = 0.5), 
             mapping = aes(group = group), size = 3, fill = "white", shape = 21, stroke = 2) +
  scale_fill_manual(values = colour4.3) +
  scale_color_manual(values = colour4.3) +
  theme_Publication_white() +
  xlab("") +
  ylab("Mapped reads") +
  ggtitle("") +
  labs(fill = "") +
  coord_flip() +
  guides(color = FALSE, fill = guide_legend("")) +
  theme(axis.text.y = element_text(face = "italic")) +
  theme(panel.grid.major.y = element_blank()) +
  scale_y_log10(expand = c(0,0, 0.1, 0), breaks = c(0, 100, 1000, 10000, 100000, 1000000), labels = c("0","100","1000","1e4","1e5",""))
dev.off()
