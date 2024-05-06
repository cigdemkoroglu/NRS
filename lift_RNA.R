lift_RNA <- function(RNA_file, insertions_file){
  
  # sample call
  # lift_RNA('hg38_refGeneMrna_cutout.txt','insertions_v5.xlsx')
  
  library(stringr)
  
  fasta_lines <- readLines(RNA_file)
  insertions_df <- as.data.frame(readxl::read_xlsx(insertions_file))

  # add deletion length and insertion size
  insertions_df$deleted_length <- insertions_df$END - insertions_df$START + 1
  insertions_df$insertion_size <- insertions_df$insert_length - insertions_df$deleted_length

  # initialize empty output file
  modified_lines_df <- data.frame('chromosome'=character(),
                                  'line_number'=integer(),
                                  'original_position'=character(),
                                  'new_position'=character(),
                                  'move_amount'=integer()
  )
  k = 1
  
  for(i in 1:length(fasta_lines)){

    # read the chromosome string
    chr_string <- str_extract(fasta_lines[i], 'chr([0-9]+|X|Y):[0-9]+')
    
    if (!is.na(chr_string)){

      # get the chromosome number and position from the chromosome string
      split_chr_string <- str_split(chr_string, ':')
      chromosome <- split_chr_string[[1]][1]
      chromosome_number <- str_extract(chromosome, '[0-9]+|X|Y')
      position <- as.integer(split_chr_string[[1]][2])

      # filter the insertions that happen before the position
      chr_filter <- insertions_df$CHROM == chromosome_number[1]
      position_filter <- insertions_df$START <= position

      # Find the total move amount by addding the insertion sizes
      move_amount <- sum(insertions_df[chr_filter&position_filter, 'insertion_size'])

      # add the original position to find the new position
      new_position <- as.character(position + move_amount)

      # create new chr string with the new position in the original format
      new_chr_string <- paste0(chromosome, ':', new_position)

      # write the new position to the corresponding line
      fasta_lines[i] <- str_replace(fasta_lines[i], chr_string, new_chr_string)
      
      # add to the modified lines output if position is modified
      if (move_amount > 0){
        modified_lines_df[k,] <-c(chromosome, i, as.character(position), new_position, move_amount)
        k = k + 1
      }
    }
  }
  
  # write to the output files
  writeLines(fasta_lines, 'fasta_modified.txt')
  write.csv(modified_lines_df, 'modified_rnas.csv')
  
}
