####################################################################

# This script was shamelessly stolen from Felix Grunberger
# https://github.com/felixgrunberger/Native_RNAseq_Microbes
# (check out the guy, he is pretty cool)

####################################################################


# Get libraries (probs more than needed?)

library(here)
source(here("seq/R/load.libraries.R"))

# Load functions

# Felix's unlist bam files function

.unlist <- function (x){
  x1 <- x[[1L]]
  if (is.factor(x1)){
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}

# Felix's wrapper function

counts_wrapper <- function (input_bam_file, input_fasta_file, input_gff_file){
  
  # > read in fasta file
  fasta <- readDNAStringSet(filepath = input_fasta_file)
  
  # > use featurecounts to calculate counts of mapped reads to features (CDS, tRNA, rRNA) 
  datalist  <-  list()
  datalist2 <- list()
  datafile <- str_split(input_bam_file, "/")
  interesting_list <- c("CDS", "rRNA", "tRNA")
  
  for (i in seq_along(interesting_list)){
    name <- interesting_list[i]
    dir.create(paste(here("seq/R/data/featurecounts_data_"),name, sep = ""),showWarnings = FALSE)
    datalist[[i]] <- featureCounts(allowMultiOverlap = T, files = input_bam_file, annot.ext = input_gff_file, isGTFAnnotationFile = T, GTF.featureType = name, GTF.attrType = "ID", isLongRead = T,nthreads = 8, reportReads = "CORE", reportReadsPath = paste(here("seq/R/data/featurecounts_data_"),name, sep = ""))
    datalist2[[i]] <- fread(paste(here("seq/R/data/featurecounts_data_"),name, "/",datafile[[1]][length(datafile[[1]])], ".featureCounts", sep = "")) %>%
      dplyr::rename(id = V1, gene = V4) %>%
      dplyr::select(id, gene) %>%
      dplyr::filter(!is.na(gene)) %>%
      mutate(mapped_type = name)
  }
  return(list(datalist, datalist2))
}

# Felix's code on combining files
wrapper_bam_to_table <- function (input_bam_file, input_gff_file, input_fasta_file, datalist_input, output = c("read_ids", "gene_ids")){
  
  # > read in fasta file
  fasta <- readDNAStringSet(filepath = input_fasta_file)
  
  # > read in gff file and grep for feature numbers
  interesting_list <- c("CDS", "rRNA", "tRNA")
  gff_table <- read.gff(input_gff_file) %>%
    as_tibble() %>%
    mutate(start_gene = start, end_gene = end,strand_gene = strand) %>%
    dplyr::filter(type %in% interesting_list) %>%
    mutate(id_name = str_split_fixed(str_split_fixed(attributes, ";Parent=", 2)[,1], "ID=", 2)[,2],
           locus_name = ifelse(type == "CDS", str_split_fixed(str_split_fixed(attributes, ";product=", 2)[,2], ";", 2)[,1],
                               ifelse(type == "rRNA", str_split_fixed(str_split_fixed(attributes, ";product=", 2)[,2], " ", 2)[,1], 
                                      ifelse(type == "tRNA", str_split_fixed(attributes, ";product=", 2)[,2], NA )))) %>%
    dplyr::select(id_name, locus_name, start_gene, end_gene, strand_gene)
  
  # > for featurecounts calculated transcript abundacies for each gene
  if(output == "gene_ids"){
    feature_list <- c(rep(interesting_list[1],length(datalist_input[[1]]$annotation$Chr)),
                      rep(interesting_list[2],length(datalist_input[[2]]$annotation$Chr)),
                      rep(interesting_list[3],length(datalist_input[[3]]$annotation$Chr)))
    
    big_data_all <- cbind(datalist_input[[1]]$annotation, datalist_input[[1]]$counts) %>%
      rbind(cbind(datalist_input[[2]]$annotation, datalist_input[[2]]$counts)) %>%
      rbind(cbind(datalist_input[[3]]$annotation, datalist_input[[3]]$counts)) %>%
      as_tibble() %>%
      mutate(type = feature_list) %>%
      dplyr::rename(counts = 7) %>%
      rowwise() %>%
      dplyr::filter(Strand == "+" | Strand == "-") 
    
    listofdfs <- list()
    # > enable calculation for different chromosomes
    for(i in 1:length(names(fasta))){
      names(fasta) <-  str_split_fixed(names(fasta), " ", 2)[,1]
      used_chr <- str_split_fixed(names(fasta[i]), " ", 2)[,1]
      
      
      df <- big_data_all %>%
        dplyr::filter(Chr == used_chr) %>%
        rowwise() %>%
        mutate(seq = ifelse(Strand == "+" & End < length(fasta[names(fasta) == used_chr][[1]]), as.character(fasta[names(fasta) == used_chr][[1]][Start:End]),
                            ifelse(Strand == "-" & End < length(fasta[names(fasta) == used_chr][[1]]),as.character(reverseComplement(fasta[names(fasta) == used_chr][[1]][Start:End])), NA)))
      listofdfs[[i]] <- df
    }
    
    full_table_big_data <- data.frame(Reduce(rbind, listofdfs))
    big_data_all_seq_names <- left_join(full_table_big_data, gff_table, by = c("GeneID" = "id_name"))
    return(big_data_all_seq_names)
  }
  # > single read table output
  if(output == "read_ids"){
    
    assigned_features <- do.call(rbind,datalist_input)
    
    # >calculate correct start end end positions of mapped reads | read in BAM file with NM tag
    allReads <- readGAlignments(input_bam_file, use.names = T, param = ScanBamParam(tag=c("NM"), what="mapq"))
    
    allReads_table <- GenomicAlignments::as.data.frame(allReads) %>%
      mutate(minion_read_name = names(allReads)) %>%
      left_join(summary_table, by = c("minion_read_name" = "read_id")) 
    
    # > calculate number of aligned reads based on CIGAR operations (M,I)
    allReads_table$aligned_reads <- NA
    allReads_table$aligned_reads <- unlist(lapply(explodeCigarOpLengths(allReads_table$cigar, ops = c("M", "I")), function(x) sum(x)))
    
    # > join featurecounts table | calculate mapping identity | factorize strands
    allReads_table_filtered <- allReads_table %>%
      left_join(assigned_features, by = c("minion_read_name" = "id")) %>%
      mutate(identity = (1 - NM/aligned_reads)*100,
             length_read = qwidth)  %>%
      separate_rows(gene, sep = ",") %>%
      left_join(gff_table, by = c("gene" = "id_name")) %>%
      mutate(shortest_distance_to_gene = ifelse(abs(start-start_gene) <= max(sequence_length_template) & abs(end-end_gene)<=max(sequence_length_template), T, F)) %>%
      dplyr::filter(shortest_distance_to_gene == T) %>%
      mutate(strand = factor(strand, levels = c("+","-"))) %>%
      dplyr::filter(strand == strand_gene) %>%
      dplyr::select(seqnames, strand, qwidth, start, end, width, mapq, NM, minion_read_name, 
                    sequence_length_template, mean_qscore_template, aligned_reads, identity, gene, mapped_type, length_read,
                    start_gene, end_gene, locus_name, strand_gene) 
    
    return(allReads_table_filtered)
  }
}

#Read and tag guppy table
get_guppy_table <- function(table_file,name){
  input_table<- paste(here::here("seq/R/data/"),table_file, ".txt", sep = "")
  table <- fread(input_table) %>%
    dplyr::mutate(run = name)
}

#------#
# MAIN #
#______#

# Load data
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# FILES
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Get summary files
summary_files <- paste(here("seq/R/data/summary_data/"), list.files(here("seq/R/data/summary_data/")), sep = "")

# Felix's code kept the file extensions, so overriding manually
sample_names <- c("RUN1","RUN2") 
  
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD AND WRITE TO R OBJECT
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

for (i in seq_along(sample_names)){
  
  #...................................set names
  working_directory   <- paste(here(), "/seq/R/data", sep = "")
  sample_name         <- sample_names[i]
  input_summary       <- summary_files[i]
  where_to_store      <- paste(working_directory,"/tidy_data/", sample_name, sep = "")
  input_fasta         <- paste(working_directory, "/msmeg.fasta", sep = "")
  input_gff           <- paste(working_directory, "/msmeg3.gff", sep = "")
  input_bam           <- paste(working_directory, "/", sample_name,".bam", sep = "")
  type_features       <- c("All", "rRNA", "CDS", "tRNA")
  
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # COMPUTE SUMMARY TABLE
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  summary_table <- fread(input_summary)
  
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # COMPUTE BAM TABLE
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  
  # > for BAM plotting (output from featurecounts)
  full_counts_table <- counts_wrapper(input_bam_file = input_bam,
                                      input_fasta_file = input_fasta,
                                      input_gff_file = input_gff)
  
  # > for BAM plotting (identity calculation, ...)
  full_id_table <- wrapper_bam_to_table(input_bam_file = input_bam, 
                                        input_gff_file = input_gff,
                                        input_fasta_file = input_fasta,
                                        datalist_input = full_counts_table[[2]],
                                        output = "read_ids")
  
  # > for BAM plotting (names of genes, ...)
  full_gene_table <- wrapper_bam_to_table(input_bam_file = input_bam, 
                                          input_gff_file = input_gff,
                                          input_fasta_file = input_fasta,
                                          datalist_input = full_counts_table[[1]],
                                          output = "gene_ids")
  #.........save as R files
  
  # 1) single read matrix with metadata information
  save(full_id_table, file = paste(where_to_store, "_id_table", sep = ""))
  # 2) count information for each gene
  save(full_gene_table, file = paste(where_to_store, "_gene_table", sep = ""))
  
  #.........save as .tsv
  
  # 1) single read matrix with metadata information
  fwrite(full_id_table, file = paste(where_to_store, "_id_table.tsv", sep = ""))
  # 2) count information for each gene
  fwrite(full_gene_table, file = paste(where_to_store, "_gene_table.tsv", sep = ""))
  }

