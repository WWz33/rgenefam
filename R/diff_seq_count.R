diff_seq_count <- function(a = NULL, b = NULL){
  a_set <- names(a)
  b_set <- names(b)
  a_uniq <- a_set[!a_set %in% b_set]
  b_uniq <- b_set[!b_set %in% a_set]
  diff_summary <- list(a_uniq = a_uniq, b_uniq = b_uniq)
  return(diff_summary)
}

sumy <- diff_seq_count(test, new_test)
