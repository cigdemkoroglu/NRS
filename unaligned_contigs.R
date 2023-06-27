f <- function(file_location) {
  
  library(data.table)
  library(magrittr)
  
  alignment <- fread(file_location)
  
  names(alignment) <- c('st_ref', 'end_ref', 'st_query', 'end_query', 'length_aln_ref', 'length_aln_query', 'percent_identity', 'length_ref',  'length_query', 'tag_ref', 'tag_query')
  
  # delete unused columns
  cols_to_delete <- c('length_aln_ref', 'length_aln_query', 'percent_identity')
  alignment[, (cols_to_delete) := NULL]
  
  # get start and end coordinates of aligned parts in numerical order such that start < end
  alignment[, st_actual := amin(st_query, end_query)]
  alignment[, end_actual := amax(st_query, end_query)]
  
  # extract all tags
  tags <- unique(alignment[,tag_query])
  
  # initialize output data.table
  output <- data.table(Tag=character(), Start=numeric(), End=numeric())
  
  # loop over tags to find unaligned parts in each tag
  for (tag in tags){
    
    # filter alignment data.table to get rows with tag
    alignment_tag <- alignment[tag_query==tag,]
    tag_length <- alignment_tag$length_query[1]
    
    # pull all aligned coordinates in query
    aligned_coords <- unlist(mapply(function(x,y) {x:y}, alignment_tag$st_actual, alignment_tag$end_actual)) %>%
      sort() %>%
      unique()
    
    # write unaligned part at the beginning of query
    if (aligned_coords[1] > 100){
      output <- rbind(output, list(tag, 1, aligned_coords[1]-1))
    }
    
    # write unaligned parts in the middle of query
    for (i in 1:(length(aligned_coords)-1)){
      if (aligned_coords[i+1] - aligned_coords[i] > 100){
        output <- rbind(output, list(tag, aligned_coords[i]+1, aligned_coords[i+1]-1))
      }
    }
    
    # write unaligned part at the end of query
    if (alignment_tag$length_query[1] - aligned_coords[length(aligned_coords)] > 100){
      output <- rbind(output, list(tag, aligned_coords[length(aligned_coords)]+1, tag_length))
    }
    
  }
  
  # calculate total unaligned coordinates for each tag
  output_summary <- output[, .(total_unaligned_coords=sum(End-Start)), by = c('Tag')]
  
  # write outputs to file
  fwrite(output_summary, './output_alignment_summary.csv')
  fwrite(output, './output_alignment.csv')
  
}