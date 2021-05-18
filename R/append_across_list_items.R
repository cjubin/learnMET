append_across_list_items <- function(x, item) {
  names(item) <- deparse(substitute(item))
  x <- append(x, item)
  return(x)
}