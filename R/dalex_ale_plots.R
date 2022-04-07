#' ALE plots feature-wise with DALEX package
#' Wrapper functions
#'
#' @description
#'   Makes an ALE plot for a set of predictor variable of the model.
#'   Provides gateway to the package DALEX.
#'
#' @name ale_plot_split
#'
#' @param object an object of class `res_fitted_split`
#'
#' @param variable environmental or genomic variable (PC derived from SNPs
#'   matrix or single SNP) for which the ALE plot should be calculated
#'
#' @author Cathy C. Westhues \email{cathy.jubin@@uni-goettingen.de}
#' @references
#' \insertRef{dalex}{learnMET}
#' @export
ale_plot_split <- function(object, ...) {
  UseMethod("ale_plot_split")
}


#' @rdname ale_plot_split
#' @export
ale_plot_split.default <- function(object, ...) {
  stop('not implemented')
  
}

#' @rdname ale_plot_split
#' @export
ale_plot_split.fitted_xgb_reg_1 <-
  
  function(object, variable, path_plot) {
    model <- fitted_split$fitted_model
    trait <-
      as.character(fitted_split$fitted_model$pre$actions$recipe$recipe$var_info[fitted_split$fitted_model$pre$actions$recipe$recipe$var_info$role ==
                                                                                  'outcome', 'variable'])
    y_train <- as.numeric(as.data.frame(fitted_split$training %>%
                                          dplyr::select(all_of(trait)))[, 1])
    x_train <- fitted_split$training
    
    
    # create custom predict function
    pred_wrapper <- function(model, newdata)  {
      results <- model %>% predict(new_data = newdata) %>%
        as.vector()
      
      return(as.numeric(as.vector(as.data.frame(results)[, '.pred'])))
    }
    
    explainer <- DALEX::explain(
      model = model,
      data = x_train,
      y = y_train,
      predict_function = pred_wrapper,
      label = "xgb_reg_1",
      verbose = FALSE
    )
    
    ale <- DALEX::model_profile(explainer = explainer,
                                type       = "accumulated",
                                variables  = variable,
                                N = NULL)
    
    p <-  ggplot(data=as.data.frame(ale$agr_profiles))+
      geom_line(data=,aes(
        x = `_x_`,
        y = `_yhat_`
      ), size = .35) +
      #theme_ipsum() +
      ylab('Average prediction') + 
      xlab(variable) +
      theme_bw(base_size = 10) +
      theme(
        plot.title = element_text(size = 1),
        legend.position = 'none',
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14)
      )
    
    
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
