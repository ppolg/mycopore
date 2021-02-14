####################################################################

# Parts of script shamelessly stolen from Felix Grunberger
# https://github.com/felixgrunberger/Native_RNAseq_Microbes
# (check out the guy, he is pretty cool)

####################################################################



####################################################################
#      Somehow generate nice QC plots from a variety of files?     #
####################################################################

#------#
# Init #
#______#

# Get libraries (probs more than needed?)

library(here)
source(here("seq/R/load.libraries.R"))

# Load functions

#Read and tag guppy table
get_guppy_table <- function(table_file,name){
  input_table<- paste(here::here("seq/R/data/"),table_file, ".txt", sep = "")
  table <- fread(input_table) %>%
    dplyr::mutate(run = name)
  
  return(table)
}

# Set plotting colours
colours <- rev(c("#E3812B","#A3280A"))
colours2 <- rev(c("#8A5122","#E0D253")) 

#------#
# Main # 
#______#

# Read guppy tables

gtable_RUN1 <- get_guppy_table("RUN1_summary","RUN1")
gtable_RUN2 <- get_guppy_table("RUN2_summary","RUN2")

# Merge guppy tables

gtable_all <- rbind(gtable_RUN1,gtable_RUN2) %>%
  group_by(run) %>%
  dplyr::mutate(
    run = as.factor(run),
    start_time = as.numeric(start_time), 
    sequence_length_template = as.numeric(sequence_length_template),
    mean_qscore_template = as.numeric(mean_qscore_template),
    type = ifelse(calibration_strand_identity > 0.7, "enolase", "genome")) # No calibration score???

# Summarise gtable

gtable_summary <- gtable_all %>%
  group_by(run) %>%
  summarise(number_of_reads = n(),
            number_of_bases = sum(sequence_length_template, na.rm = TRUE),
            median_rawread_length = median(sequence_length_template, na.rm = TRUE),
            median_rawread_quality = median(mean_qscore_template, na.rm = TRUE))

# Regroup(?)
gtable_all$group <-  factor(gtable_all$run, 
                               levels = rev(c("RUN1","RUN2")))

# Filter for passed only
gtable_passed <- gtable_all %>%
  dplyr::filter(passes_filtering == "TRUE")
  

# Generate yield data
gtable_all_yield <- gtable_all %>%
  mutate(timeindex = round((start_time + 3600)/1800)) %>%
  group_by(run, timeindex) %>% 
  summarise(yield = sum(sequence_length_template)/1000000000) %>%
  mutate(cumulative_yield = cumsum(yield))



#------#
# Plot # 
#______#

# Yield over time

pdf(here("seq/R/Out/yield_plot.pdf"), 
    width = 7, height = 7, paper = "special",onefile=FALSE)
ggplot(gtable_all_yield, aes(x=timeindex/2, y = cumulative_yield, color = run, fill = run)) +
  geom_line(alpha = 0.8, size = 3 ) +
  xlab("Run time (h)") +
  ylab("Cumulative yield [Gb]") +
  theme_Publication_white() +
  scale_color_manual(values = colours) +
  scale_fill_manual(values = colours) +
  scale_linetype_manual(values = c(4,1)) +
  scale_y_continuous(expand = c(0,0, 0.1, 0), breaks = c(0.0,0.2,0.4,0.6,0.8)) +
  scale_x_continuous(limits = c(0,72), expand = c(0,0), breaks = c(24,48,72)) +
  guides(color = guide_legend(title = ""), fill = guide_legend(title = ""), linetype = guide_legend(title = ""))
dev.off()

# Raw read length vs control

pdf(here("seq/R/Out/length_vs_control.pdf"), 
    width = 7, height = 7, paper = "special",onefile=FALSE)
ggplot(data = gtable_all, aes(x = sequence_length_template, y = run, fill = type, color = type)) +
  geom_density_ridges2(aes(height =..ndensity..), scale = 0.9, alpha = 0.8, size = 1) +
  theme_Publication_white() +
  scale_x_continuous(trans = "log10", limits = c(50,6000), breaks = c(100,500,1000,5000),expand = c(0, 0)) +
  geom_vline(xintercept = 1314, linetype = "dashed", alpha = 0.5) +
  ylab("") +
  scale_y_discrete(expand = c(0,0, 0.2, 0)) +
  xlab("Read length [nt]") +
  scale_fill_manual(values = colours2) +
  scale_color_manual(values = colours2) +
  guides(group = F, fill= F, color = F) +
  theme(axis.text.y = element_text(face = "italic"))  
dev.off()

pdf(here("seq/R/Out/length_vs_control_pass.pdf"), 
    width = 7, height = 7, paper = "special",onefile=FALSE)
ggplot(data = gtable_passed, aes(x = sequence_length_template, y = run, fill = type, color = type)) +
  geom_density_ridges2(aes(height =..ndensity..), scale = 0.9, alpha = 0.8, size = 1) +
  theme_Publication_white() +
  scale_x_continuous(trans = "log10", limits = c(50,6000), breaks = c(100,500,1000,5000),expand = c(0, 0)) +
  geom_vline(xintercept = 1314, linetype = "dashed", alpha = 0.5) +
  ylab("") +
  scale_y_discrete(expand = c(0,0, 0.2, 0)) +
  xlab("Read length [nt]") +
  scale_fill_manual(values = colours2) +
  scale_color_manual(values = colours2) +
  guides(group = F, fill= F, color = F) +
  theme(axis.text.y = element_text(face = "italic"))  
dev.off()

# Quality vs control

pdf(here("seq/R/Out/quality_vs_control.pdf"), 
    width = 7, height = 7, paper = "special",onefile=FALSE)
ggplot(data = gtable_all, aes(x = mean_qscore_template, y = run, fill = type, color = type)) +
  geom_density_ridges2(aes(height =..ndensity..), scale = 0.9, alpha = 0.8, size = 1) +
  theme_Publication_white() +
  geom_vline(xintercept = 7, linetype = "dashed", alpha = 0.5) +
  annotate(geom = "text",
           label = "Threshold                       ",
           x = 7.1,
           y = 3,
           size = 2.5,
           angle = 90,
           vjust = 1) +
  scale_x_continuous(limits = c(0,16), expand = c(0, 0)) +
  ylab("") +
  scale_y_discrete(expand = c(0,0, 0.2, 0)) +
  xlab("Read quality") +
  scale_fill_manual(values = colours2) +
  scale_color_manual(values = colours2) +
  guides(group = F, fill= F, color = F) +
  theme(axis.text.y = element_text(face = "italic"))  
dev.off()

pdf(here("seq/R/Out/quality_vs_control_pass.pdf"), 
    width = 7, height = 7, paper = "special",onefile=FALSE)
ggplot(data = gtable_passed, aes(x = mean_qscore_template, y = run, fill = type, color = type)) +
  geom_density_ridges2(aes(height =..ndensity..), scale = 0.9, alpha = 0.8, size = 1) +
  theme_Publication_white() +
  geom_vline(xintercept = 7, linetype = "dashed", alpha = 0.5) +
  annotate(geom = "text",
           label = "Threshold                       ",
           x = 7.1,
           y = 3,
           size = 2.5,
           angle = 90,
           vjust = 1) +
  scale_x_continuous(limits = c(0,16), expand = c(0, 0)) +
  ylab("") +
  scale_y_discrete(expand = c(0,0, 0.2, 0)) +
  xlab("Read quality") +
  scale_fill_manual(values = colours2) +
  scale_color_manual(values = colours2) +
  guides(group = F, fill= F, color = F) +
  theme(axis.text.y = element_text(face = "italic"))  
dev.off()


