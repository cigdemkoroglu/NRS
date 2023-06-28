lift_RNA <- function(RNA_file, insertions_file){
  
  # sample call
  # lift_RNA('hg38_refGeneMrna_cutout.txt','insertions_v5.xlsx')
  
  library(stringr)
  
  fasta_lines <- readLines(RNA_file)
  insertions_df <- as.data.frame(readxl::read_xlsx(insertions_file))
  
  insertions_df$deleted_length <- insertions_df$END - insertions_df$START + 1
  insertions_df$insertion_size <- insertions_df$insert_length - insertions_df$deleted_length
  
  modified_lines_df <- data.frame('chromosome'=character(),
                                  'line_number'=integer(),
                                  'original_position'=character(),
                                  'new_position'=character(),
                                  'move_amount'=integer()
  )
  k = 1
  
  for(i in 1:length(fasta_lines)){
    
    chr_string <- str_extract(fasta_lines[i], 'chr([0-9]+|X|Y):[0-9]+')
    
    if (!is.na(chr_string)){
      
      split_chr_string <- str_split(chr_string, ':')
      chromosome <- split_chr_string[[1]][1]
      chromosome_number <- str_extract(chromosome, '[0-9]+|X|Y')
      position <- as.integer(split_chr_string[[1]][2])
      
      chr_filter <- insertions_df$CHROM == chromosome_number[1]
      position_filter <- insertions_df$START <= position
      
      move_amount <- sum(insertions_df[chr_filter&position_filter, 'insertion_size'])
      
      new_position <- as.character(position + move_amount)
      
      new_chr_string <- paste0(chromosome, ':', new_position)
      
      fasta_lines[i] <- str_replace(fasta_lines[i], chr_string, new_chr_string)
      if (move_amount > 0){
        modified_lines_df[k,] <-c(chromosome, i, as.character(position), new_position, move_amount)
        k = k + 1
      }
    }
  }
  
  
  writeLines(fasta_lines, 'fasta_modified.txt')
  write.csv(modified_lines_df, 'modified_rnas.csv')
  
}