#' Apply PCA on a split: list object containing a training and test datasets. 
#' 
#' @description 
#' Fit PCA on the training set and apply the same transformation to the test set.
#' The goal is to use principal components in prediction models as a data 
#' dimensionality reduction method.
#'
#' @param
#'
#'
#'
#'


apply_pca <- function(split,geno) {
  
  geno$geno_ID = row.names(geno)
  
  col_to_keep = colnames(geno)
  
  geno_training = merge(split[[1]], geno, by = 'geno_ID', all.x = T)[, col_to_keep]
  geno_training = unique(geno_training)
  
  geno_test =  merge(split[[2]], geno, by = 'geno_ID', all.x = T)[, col_to_keep]
  geno_test = unique(geno_test)
  
  
  rec <- recipe(geno_ID ~ . ,
                data = geno_training) %>%
    update_role(geno_ID, new_role = 'outcome') %>%
    step_nzv(all_predictors()) %>%
    step_pca(
      all_predictors(),
      num_comp = num_pcs,
      options = list(center = T, scale. = T)
    )
  
  norm_obj <- prep(rec, training = geno_training)
  
  training_pca <- bake(norm_obj, geno_training)
  test_pca <- bake(norm_obj, geno_test)
  
  test_pca$geno_ID = geno_test$geno_ID
  
  pc_values <- rbind(training_pca, test_pca)
  pc_values <- unique(pc_values)
  pc_values <- as.data.frame(pc_values)
  
  
  return(pc_values)
  
}