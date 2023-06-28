lift_database <- function(database, insertions_df){
  
  library(data.table)
  library(magrittr)
  
  database_dt <- fread(database)
 
  names(database_dt)[1:3] <- c('chromosome', 'START', 'END')
  
  insertions_df <- as.data.frame(readxl::read_xlsx(insertions_df))
  
  insertions_df$deleted_length <- insertions_df$END - insertions_df$START + 1
  
  
  
  database_dt[, remove_row := FALSE]
  for (i in 1:nrow(insertions_df)){
    database_dt[chromosome == insertions_df$CHROM[i] & START >= insertions_df$START[i] & START <= insertions_df$END[i], remove_row := TRUE]
    database_dt[chromosome == insertions_df$CHROM[i] & END >= insertions_df$START[i] & END <= insertions_df$END[i], remove_row := TRUE]
    database_dt[chromosome == insertions_df$CHROM[i] & START <= insertions_df$START[i] & END >= insertions_df$END[i], remove_row := TRUE]
  }
  
  database_dt <- database_dt[remove_row == FALSE,]
  
  database_dt[,remove_row := NULL]
  
  database_dt[, old_START := START]
  database_dt[, old_END := END]
  for (i in 1:nrow(insertions_df)){
    database_dt[chromosome == insertions_df$CHROM[i] & old_START >= insertions_df$START[i],
                START := START + insertions_df$insert_length[i] - insertions_df$deleted_length[i]]
    database_dt[chromosome == insertions_df$CHROM[i] & old_START >= insertions_df$START[i],
                END := END + insertions_df$insert_length[i] - insertions_df$deleted_length[i]]
  }
  
  # database_dt[ ,old_START := NULL]
  # database_dt[ ,old_END := NULL]
  
  fwrite(database_dt, './database_modified.txt', sep = '\t')
  return(database_dt)
  
}