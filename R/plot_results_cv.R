plot_results_cv <- function(results_fitted_splits, cv_type, cv0_type)
  if (cv_type == 'cv0' & plot_PA) {
    if (cv0_type == 'leave-one-year-out') {
      if (length(unique(pheno$location)) == 1) {
        list_years <-
          sapply(splits, function(x)
            unique(as.character(x[[2]][, 'year'])))
        PA <-
          sapply(fitting_all_splits, function(x)
            as.numeric(x[['cor_pred_obs']]))
        df <- as.data.frame(cbind(list_years, PA))
        colnames(df) <- c('year', 'Prediction_accuracy')
        df$Prediction_accuracy <- as.numeric(df$Prediction_accuracy)
        
        p <-
          ggplot(df, aes(x = year, y = Prediction_accuracy)) + geom_bar(stat = 'identity') +
          xlab('Predicted year') + ylab(paste0('Prediction accuracy for the trait ', trait))
        ggsave(p, filename = paste0(path_plot_PA, 'cv0_leave1yearout.pdf'))
      } else {
        
      }
    }
    
    
    
    if (cv_type == 'cv1' & plot_PA) {
      list_years <-
        sapply(splits, function(x)
          unique(as.character(x[[2]][, 'year'])))
      PA <-
        sapply(fitting_all_splits, function(x)
          as.numeric(x[['cor_pred_obs']]))
      df <- as.data.frame(cbind(list_years, PA))
      colnames(df) <- c('year', 'Prediction_accuracy')
      df$Prediction_accuracy <-
        as.numeric(df$Prediction_accuracy)
      
      p <-
        ggplot(df, aes(x = year, y = Prediction_accuracy)) + geom_bar(stat = 'identity') +
        xlab('Predicted year') + ylab(paste0('Prediction accuracy for the trait ', trait))
      ggsave(p, filename = paste0(path_plot_PA, 'cv0_leave1yearout.pdf'))
      
    }
  }
