lift_RefSeq <- function(database_file, insertions_file, names_file){
  
  library(data.table)
  
  # sample call
  # lift_RefSeq('./lift_RefSeq/hg38_refGene.txt', './lift_RefSeq/insertions_v5.xlsx', './lift_RefSeq/RefSeq_changes.xlsx')
  
  database_dt <- fread(database_file)
  df_for_names <- as.data.frame(readxl::read_xlsx(names_file))
  insertions_df <- as.data.frame(readxl::read_xlsx(insertions_file))
  
  original_names <- names(database_dt) # save original column names to convert back at the end
  names(database_dt)  <- names(df_for_names) # use actual column names from the other file
  database_dt[,move_amount := 0] # move_amount column =
  database_dt[,change_flag := character()]
  
  insertions_df$deleted_length <- insertions_df$END - insertions_df$START + 1
  
  valid_chrs <- c(paste0('chr', 1:23), 'chrX', 'chrY')
  database_dt <- database_dt[chrom %in% valid_chrs,]
  
  columns_to_modify <- c('txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonStarts', 'exonEnds')
  
  database_dt[, (paste0('old_', columns_to_modify)) := .SD, .SDcols = columns_to_modify]
  modified_gene_names <- c()
  
  for (i in 1:nrow(insertions_df)){
    
    chromosome <- paste0('chr', insertions_df$CHROM[i])
    start <- insertions_df$START
    end <- insertions_df$END
    insertion_size = insertions_df$insert_length[i] - insertions_df$deleted_length[i]
    
    database_dt[chrom == chromosome & txStart >= insertions_df$START[i],
                move_amount := move_amount + insertion_size]
    
    database_dt[chrom == chromosome & txStart <= insertions_df$START[i] & txEnd >= insertions_df$START[i], change_flag := i]
    
  }
  
  move_exons <- function(exons, total_move_amount, insertion_size=0, start=0){
    exons <- as.numeric(unlist(strsplit(exons, ',')))
    for (i in 1:length(exons)){
      exons[i] = exons[i] + total_move_amount
      if (exons[i] > start){
        exons[i] = exons[i] + insertion_size
      }
    }
    return(paste0(paste0(as.character(exons), collapse=','),','))
  }
  
  percent_cut = nrow(database_dt)%/%100
  
  for (i in 1:nrow(database_dt)){
    
    if (i %% percent_cut == 0){
      percent = i / percent_cut
      {cat('\r', percent, '%')}
    }
    database_dt[i, txStart := txStart + move_amount]
    database_dt[i, txEnd := txEnd + move_amount]
    database_dt[i, cdsStart := cdsStart + move_amount]
    database_dt[i, cdsEnd := cdsEnd + move_amount]
    
    exon_starts <- database_dt$exonStarts[i]
    exon_ends <- database_dt$exonEnds[i]
    move_amount <- database_dt$move_amount[i]
    
    database_dt[i, exonStarts := move_exons(exon_starts, move_amount)]
    database_dt[i, exonEnds := move_exons(exon_ends, move_amount)]
    
    if(!is.na(database_dt$change_flag[i])){
      flag <- as.integer(database_dt$change_flag[i])
      start <- insertions_df$START[flag]
      insertion_size <- insertions_df$insert_length[flag] - insertions_df$deleted_length[flag]
      
      database_dt[i, exonStarts := move_exons(exon_starts, move_amount, insertion_size, start)]
      database_dt[i, exonEnds := move_exons(exon_ends, move_amount, insertion_size, start)]
      database_dt[i, txEnd:= txEnd + insertion_size]
      if(database_dt$old_cdsStart[i] >= insertions_df$START[flag]){
        database_dt[i, cdsStart:=cdsStart + insertion_size]
      }
      if(database_dt$old_cdsEnd[i] >= insertions_df$START[flag]){
        database_dt[i, cdsEnd:=cdsEnd + insertion_size]
      }
    }
  }
  
  database_dt[, (paste0('old_', columns_to_modify)) := NULL]
  
  names(database_dt)[1] <- '#CHROM'
  fwrite(database_dt, './database_modified.txt', sep = '\t')
  
}