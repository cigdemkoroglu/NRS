find_gained_lost <- function(database1, database_insertion_modified, insertions_file, chr_number){
  
  library(data.table)

  # deneme
  
  # sample call
  # find_gained_lost ('chr10_hg38_2023_multianno.txt', 'chr10_hg38_NRS_multianno.txt', 'insertions_v5.xlsx', '10')
  
  database1_dt <- fread(database1)
  database1_dt[, original_start := Start]
  database1_dt[, original_end := End]
  
  database_insertion_modified_dt <- fread(database_insertion_modified)
  insertions_df <- as.data.frame(readxl::read_xlsx(insertions_file))
  insertions_df <- insertions_df[insertions_df$CHROM == chr_number,]
  insertions_df$deleted_length <- insertions_df$END - insertions_df$START + 1
  
  for (i in 1:length(insertions_df)){
    start <- insertions_df$START[i]
    insertion_size <- insertions_df$insert_length[i] - insertions_df$deleted_length[i]
    database1_dt[Start>=start, ':='(Start = Start + insertion_size, End = End + insertion_size)]
  }
  
  common_rs <- intersect(database1_dt$avsnp150, database_insertion_modified_dt$avsnp150)
  common_rs <- common_rs[common_rs!='.']
  common_starts <- intersect(database1_dt[avsnp150=='.', Start], database_insertion_modified_dt[avsnp150=='.', Start])
  
  common_dt <- merge.data.table(database1_dt, database_insertion_modified_dt, by=c('Chr', 'Start', 'End', 'Ref', 'Alt', 'avsnp150'))
  
  lost_dt <- database1_dt[!(avsnp150 %in% common_rs | (avsnp150 == '.' & Start %in% common_starts)),]
  lost_dt[,Start := original_start]
  lost_dt[,End := original_end]
  lost_dt[,original_start := NULL]
  lost_dt[,original_end := NULL]
  gained_dt <- database_insertion_modified_dt[!(avsnp150 %in% common_rs | (avsnp150 == '.' & Start %in% common_starts)),]
  
  
  print(paste0('For ','chr',chr_number,': ','first database - lost - common = ', nrow(database1_dt) - nrow(lost_dt) - nrow(common_dt)))
  
  print(paste0('For ','chr',chr_number,': ','second database - gained - common = ', nrow(database_insertion_modified_dt) - nrow(gained_dt) - nrow(common_dt)))
  
  fwrite(lost_dt, 'lost.txt', sep = '\t')
  fwrite(gained_dt, 'gained.txt', sep = '\t')
}
