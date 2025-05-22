#' Make group contrasts
#'
#' @param groupX The wanted group
#' @param groups All the groups
#'
#' @return
#' A contrast
#' @export
#'
#' @examples
#' groupX <- 'A'
#' groups <- c('A', 'B', 'C')
#' make_group_contrasts_for_one(groupX, groups)
make_group_contrasts_for_one <- function(groupX,
                                         groups){

  unique.groups<-unique(groups)
  results <- rep("", length(unique.groups)-1)
  counter <- 1
  for(groupY in unique.groups)
    if(groupX!=groupY){
      results[counter] <- paste(groupX,"-", groupY, sep="")
      counter <- counter+1
    }
  return(results)

}
