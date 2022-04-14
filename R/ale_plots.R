#' ALE plots feature-wise
#'
#' @description
#'   Makes an ALE plot for one single predictor variable of the model.
#'
#'
#' @name ALE_plot_split
#'
#' @param object an object of class `res_fitted_split`
#'
#' @param variable environmental or genomic variable (PC derived from SNPs
#'   matrix or single SNP) for which the ALE plot should be calculated
#'
#' @param path_plot path where the ALE plot should be saved
#'
#' @param min_interval_points minimum number of points contained within each
#'   interval (e.g. among the unique values taken by the feature). Values less
#'   than 4 should be avoided.
#'
#' @param nb_bins number of intervals. Note that if `nb_bins` is provided by the
#'   user, the `min_interval_points` argument is ignored.
#'
#'
#'
#' @author Cathy C. Westhues \email{cathy.jubin@@uni-goettingen.de}
#' @export




ALE_plot_split <-
  
  
  function(object,
           variable,
           path_plot,
           nb_bins = 10) {
    checkmate::assert_class(object, 'res_fitted_split')
    
    model <- object$fitted_model
    trait <-
      as.character(object$fitted_model$pre$actions$recipe$recipe$var_info[object$fitted_model$pre$actions$recipe$recipe$var_info$role ==
                                                                            'outcome', 'variable'])
    y_train <- as.numeric(as.data.frame(object$training %>%
                                          dplyr::select(all_of(trait)))[, 1])
    x_train <- object$training
    
    variable_df <- as.data.frame(x_train[,variable])
    colnames(variable_df) <- 'variable'
    
    # create custom predict function
    pred_wrapper <- function(model, newdata)  {
      results <- model %>% predict(new_data = newdata) %>%
        as.vector()
      
      return(as.numeric(as.vector(as.data.frame(results)[, '.pred'])))
    }
    
    # Create the intervals based on the number of bins given by the user
    
    if (round(length(unique(x_train$ClayProp)) / 2) < nb_bins) {
      nb_bins <- round(length(unique(x_train$ClayProp)) / 2)
      s <-
        split(x_train, Hmisc::cut2(x_train[, variable], g = nb_bins))
      
    } else{
      s <-
        split(x_train, Hmisc::cut2(x_train[, variable], g = nb_bins))
      
    }
    # Compute predictions for the upper and lower values of the defined bins
    min_max_interval <- function(s, v, variable) {
      if (v == 1) {
        return(c(min(s[[v]][, variable]), max(s[[v]][, variable])))
      }
      if (v > 1) {
        return(c(max(s[[v - 1]][, variable]), max(s[[v]][, variable])))
      }
      
    }
    min_max_by_interval <-
      lapply(1:length(s), function(x) {
        min_max_interval(s, v = x, variable = variable)
      })
    z = c(min(sapply(min_max_by_interval, "[[", 1)), sapply(min_max_by_interval, "[[", 2))
    
    pred_within_interval <- function(s, v, variable) {
      new_lower <- s[[v]]
      new_lower[, variable] <- min_max_by_interval[[v]][1]
      
      lower_pred = pred_wrapper(model, newdata = new_lower)
      
      new_upper <- s[[v]]
      new_upper[, variable] <- min_max_by_interval[[v]][2]
      
      upper_pred = pred_wrapper(model, newdata = new_upper)
      
      # Sum the effects and divide by the number of observations
      return(sum(upper_pred - lower_pred) / nrow(s[[v]]))
      
    }
    
    
    interval_effects <-
      unlist(lapply(1:length(s), function(x) {
        pred_within_interval(s, v = x, variable = variable)
      }))
    
    
    # Accumulate the effects
    
    interval_effects_acc <- c(0, cumsum(interval_effects))
    
    
    # Centering the effects by subtracting the avg effect (approx. per bin)
    fq <-
      unlist(lapply(s, function(x)
        nrow(x))) #frequency count of X_train[,variable] values falling into the respective bins
    
    
    centered_val = interval_effects_acc - sum((interval_effects_acc[1:nb_bins] +
                                                 interval_effects_acc[2:(nb_bins + 1)]) / 2 * fq) / sum(fq)
    
    
    # Plot function
    #plot(x = z, y = centered_val, type = 'l')
    df <- as.data.frame(cbind(z, centered_val))
    colnames(df) <- c('z', 'centered_val')
    
    ggplot(data = df) +
      geom_line(data = , aes(x = z,
                             y = centered_val), size = .35) +
      #theme_ipsum() +
      ylab(paste0('ALE of ', variable)) +
      xlab(variable) +
      theme_bw(base_size = 10) +
      theme(
        plot.title = element_text(size = 1),
        legend.position = 'none',
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14)
      )+
      geom_rug(data=variable_df, aes(x = variable), inherit.aes = F)
    
    
    ggsave(
      file = paste0(path_plot,
                    '/',
                    variable,
                    '_ale_plot.pdf'),
      width = 20,
      height = 10,
      units = "cm"
    )
  }
