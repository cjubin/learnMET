#' Compute variable importance according to the machine learning algorithm used
#'
#' @description
#' Variable importance calculated with permutation approach
#'
#'
#' @name permutation_based_vip
#'
#' @param object an object of class `res_fitted_split`
#'
#' @param type `model_specific` or `model_agnostic`
#'
#' @param permutations = 20
#'
#' @param unseen_data
#'
#' @author Cathy C. Westhues \email{cathy.jubin@@uni-goettingen.de}
#' @references
#'

permutation_based_vip <- function(model,
                                  x,
                                  y,
                                  permutations,
                                  predictors,
                                  path_plot) {
  pred_wrapper <- function(model, newdata)  {
    results <- model %>% predict(new_data = newdata) %>%
      as.vector()
  }
  
  pred_x_train <-
    as.data.frame(pred_wrapper(model = model, newdata = x))[, '.pred']
  
  loss_original <-
    sqrt(mean((y - pred_x_train) ^ 2))
  
  
  permutation_loss_one_variable <- function(variable, x) {
    df_shuffled <- x
    df_shuffled[, variable] <- sample(df_shuffled[, variable])
    
    pred_x_train_permuted <-
      as.data.frame(pred_wrapper(model = model, newdata = df_shuffled))[, '.pred']
    
    #calculate rmse with predictions from permuted dataset
    loss_permuted <-
      sqrt(mean((y - pred_x_train_permuted) ^ 2))
    
    # return difference loss
    diff <- loss_permuted - loss_original
    names(diff) <- variable
    return(diff)
  }
  
  # VIP for one permutation for all features
  results_permutation <- vector(mode = 'list', length = permutations)
  
  for (j in 1:permutations) {
    results_permutation[[j]] <-
      as.data.frame(t(unlist(lapply(predictors, function(v) {
        permutation_loss_one_variable(variable = v, x = x)
      }))))
  }
  
  table_all_permutations <- as.data.frame(data.table::rbindlist(results_permutation))
  long <- reshape2::melt(table_all_permutations)
  long2 <- long %>% group_by(variable) %>% dplyr::summarize(Mean = mean(value, na.rm=TRUE))
  long2 <- long2 %>% slice_max(Mean, n = 40)
  
  ggplot(long[long$variable%in%long2$variable,], aes(x=variable, y=value)) + 
    geom_boxplot(color="black") +
    ylab(expression(vip[diff](~e[perm] - e[original]) )) + 
    xlab('Predictor variable') +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(size = 1),
      legend.position = 'none',
      axis.text.x = element_text(size = 13),
      axis.text.y = element_text(size = 13),
      axis.title.y = element_text(size = 16),
      axis.title.x = element_text(size = 16)
    )+coord_flip()
  
  
  ggsave(
    file = paste0(path_plot,
                  '/vip_permutation_plot.pdf'),
    width = 25,
    height = 17,
    units = "cm"
  )
  
  return(long)
}