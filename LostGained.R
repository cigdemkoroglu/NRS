find_gained_lost <- function(annotated_variant, annotated_variant_NRS, insertions_file, chr_number){
  
  library(data.table)

  # sample call
  # find_gained_lost ('chr10_hg38_2023_multianno.txt', 'chr10_hg38_NRS_multianno.txt', 'insertions_v5.xlsx', '10')

  # read chromosome based annotation files for reference
  
  annotated_variant_dt <- fread(annotated_variant) 
  annotated_variant_dt[, original_start := Start]
  annotated_variant_dt[, original_end := End]

  # read chromosome based annotation files for reference-NRS
  
  annotated_variant_NRS_dt <- fread(annotated_variant_NRS)

  # read insertions file, filter for chromosome, calculate length of deletions
  
  insertions_df <- as.data.frame(readxl::read_xlsx(insertions_file))
  insertions_df <- insertions_df[insertions_df$CHROM == chr_number,]
  insertions_df$deleted_length <- insertions_df$END - insertions_df$START + 1

  # For each insertion in insertions file, increase start and end coordinates 
  # of the effected x as much as the insertion size
  
  for (i in 1:length(insertions_df)){
    start <- insertions_df$START[i]
    insertion_size <- insertions_df$insert_length[i] - insertions_df$deleted_length[i]
    annotated_variant_dt[Start>=start, ':='(Start = Start + insertion_size, End = End + insertion_size)]
  }

  # 
  
  common_rs <- intersect(annotated_variant_dt$avsnp150, annotated_variant_NRS_dt$avsnp150)
  common_rs <- common_rs[common_rs!='.']
  common_starts <- intersect(annotated_variant_dt[avsnp150=='.', Start], annotated_variant_NRS_dt[avsnp150=='.', Start])
  
  common_dt <- merge.data.table(annotated_variant_dt, 
                                annotated_variant_NRS_dt, 
                                by=c('Chr', 'Start', 'End', 'Ref', 'Alt', 'avsnp150')
                               )
  
  lost_dt <- annotated_variant_dt[!(avsnp150 %in% common_rs | (avsnp150 == '.' & Start %in% common_starts)),]
  lost_dt[,Start := original_start]
  lost_dt[,End := original_end]
  lost_dt[,original_start := NULL]
  lost_dt[,original_end := NULL]
  
  gained_dt <- annotated_variant_NRS_dt[!(avsnp150 %in% common_rs | (avsnp150 == '.' & Start %in% common_starts)),]
  
  
  print(paste0('For ','chr',chr_number,': ','first database - lost - common = ', nrow(annotated_variant_dt) - nrow(lost_dt) - nrow(common_dt)))
  
  print(paste0('For ','chr',chr_number,': ','second database - gained - common = ', nrow(annotated_variant_NRS_dt) - nrow(gained_dt) - nrow(common_dt)))
  
  fwrite(lost_dt, 'lost.txt', sep = '\t')
  fwrite(gained_dt, 'gained.txt', sep = '\t')
}
