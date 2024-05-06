lift_database <- function(database, insertions_df){
  
  library(data.table)
  library(magrittr)
  
  # read the input files
  database_dt <- fread(database)
  names(database_dt)[1]  <- 'CHROM'
  insertions_df <- as.data.frame(readxl::read_xlsx(insertions_df))

  # add deletion sizes
  insertions_df$deleted_length <- insertions_df$END - insertions_df$START + 1

  # remove the rows which intersect with an insertion
  database_dt[, remove_row := FALSE]
  for (i in 1:nrow(insertions_df)){
    chromosome <- paste0('chr', insertions_df$CHROM[i])
    database_dt[CHROM == chromosome & POS >= insertions_df$START[i] & POS <= insertions_df$END[i], remove_row := TRUE]
  }
  
  database_dt <- database_dt[remove_row == FALSE,]
  database_dt[,remove_row := NULL]

  # record the original positions
  database_dt[, old_POS := POS]

  # reflect the changes after insertions and deletions to the positions
  for (i in 1:nrow(insertions_df)){
    chromosome <- paste0('chr', insertions_df$CHROM[i])
    database_dt[CHROM == chromosome & old_POS >= insertions_df$START[i],
                POS := POS + insertions_df$insert_length[i] - insertions_df$deleted_length[i]]
  }
  
  database_dt[ ,old_POS := NULL]

  # write the results
  names(database_dt)[1] <- '#CHROM'
  fwrite(database_dt, './database_modified.vcf.txt', sep = '\t')
  return(database_dt)
  
}
