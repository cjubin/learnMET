plot_results_vip <-
  function(fitting_all_splits,
           info_environments,
           method_processing,
           splits,
           cv_type,
           cv0_type,
           filename_plot_PA,
           path_folder,
           nb_folds_cv1 = nb_folds_cv1,
           repeats_cv1 = repeats_cv1,
           nb_folds_cv2 = nb_folds_cv2,
           repeats_cv2 = nb_folds_cv2) {
    
    if (cv_type == 'cv0') {
      if (cv0_type == 'leave-one-environment-out') {
        list_envs <-
          sapply(fitting_all_splits, function(x)
            unique(as.character(x[["cor_pred_obs"]][, 'IDenv'])))
        
        PA <-
          sapply(fitting_all_splits, function(x)
            as.numeric(x[['cor_pred_obs']][, 'COR']))
        
        df <- as.data.frame(cbind(list_envs, PA))
        
        colnames(df) <- c('IDenv', 'Prediction_accuracy')
        df$Prediction_accuracy <- as.numeric(df$Prediction_accuracy)
        
        p <-
          ggplot(df,
                 mapping = aes(
                   x = reorder(IDenv, Prediction_accuracy),
                   y = Prediction_accuracy,
                   group = 1
                 )) + geom_line() +
          xlab('Predicted environment') + ylab(paste0('Prediction accuracy for the trait ', trait)) + ggtitle('Leave-one-environment-out CV scheme') +
          theme(axis.text.x = element_text(
            angle = 90,
            vjust = 0.5,
            hjust = 1
          ))
        ggsave(p, filename = paste0(path_folder, '/cv0_leave1yearout_',method_processing,'.pdf'))
      }
      
      
      if (cv0_type == 'leave-one-year-out') {
        
        
        list_predicted_years <-
          as.vector(sapply(fitting_all_splits, function(x)
            as.character(unique(x[['predictions_df']][,'year']))))
        
        VIP <-
          sapply(fitting_all_splits, function(x)
            x['ranking_vip'])
        
        
        for (j in 1:length(list_predicted_years)) {
          VIP[[j]]$year <- list_predicted_years[[j]]
          VIP[[j]] <- top_n(VIP[[j]],wt= Importance,n=40)
          
        }
        
        VIP <- do.call('rbind',VIP)
        
        p <- ggplot(VIP, aes(x=reorder_within(Variable,Importance,year), y=Importance)) + ylab('Relative importance') + xlab('Top 40 predictor variables\n for each training set')+
          geom_boxplot() + facet_wrap(~ year, ncol = length(list_predicted_years), scales = "free") + coord_flip()
        
        ggsave(
          p,
          filename = paste0(path_folder, '/cv0_', method_processing, '_Variable_Importance.pdf'),
          height = 8,
          width = 12
        )
        
        
        
        
      }
      
      if (cv0_type == 'forward_prediction') {
        list_envs <-
          as.vector(sapply(fitting_all_splits, function(x)
            as.character(unique(
              as.data.frame(x[["cor_pred_obs"]])[, 'IDenv']
            ))))
        
        
        list_years <-
          info_environments[match(list_envs, info_environments$IDenv), 'year']
        
        PA <-
          as.vector(sapply(fitting_all_splits, function(x)
            as.numeric(as.data.frame(x[['cor_pred_obs']])[, 'COR'])))
        
        df <- as.data.frame(cbind(list_years, PA))
        
        colnames(df) <- c('year', 'Prediction_accuracy')
        
        df$Prediction_accuracy <-
          as.numeric(df$Prediction_accuracy)
        
        df2 <- count(df,year)
        
        p <-
          ggplot(df,
                 mapping = aes(
                   x = reorder(year, Prediction_accuracy),
                   y = Prediction_accuracy
                 )) + geom_boxplot() +
          geom_text(data = df2, aes(x=year,y = 1, label = paste0(n,' environments'))) +
          xlab('Year to predict (average over all sites)') + ylab(paste0('Prediction accuracy for the trait ', trait)) + ggtitle('Forward prediction CV scheme') +
          theme(axis.text.x = element_text(
            angle = 90,
            vjust = 0.5,
            hjust = 1
          ))
        ggsave(p,
               filename = paste0(path_folder, '/cv0_forwardprediction_show_year_',method_processing,'.pdf'))
        
        
        PA <-
          as.vector(sapply(fitting_all_splits, function(x)
            as.character(unique(
              as.data.frame(x[["cor_pred_obs"]])[, 'COR']
            ))))
        
        
        df <- as.data.frame(cbind(list_envs, PA))
        
        colnames(df) <- c('IDenv', 'Prediction_accuracy')
        
        df$Prediction_accuracy <-
          as.numeric(df$Prediction_accuracy)
        
        p <-
          ggplot(df,
                 mapping = aes(
                   x = reorder(IDenv, Prediction_accuracy),
                   y = Prediction_accuracy,
                   group = 1
                 )) + geom_line() +
          xlab('Predicted environment') + ylab(paste0('Prediction accuracy for the trait ', trait)) + ggtitle('Forward prediction CV scheme') +
          theme(axis.text.x = element_text(
            angle = 90,
            vjust = 0.5,
            hjust = 1
          ))
        ggsave(p,
               filename = paste0(path_folder, '/cv0_forwardprediction_show_env_',method_processing,'.pdf'))
      }
      
      if (cv0_type == 'leave-one-site-out') {
        list_envs <-
          as.vector(sapply(fitting_all_splits, function(x)
            as.character(unique(
              as.data.frame(x[["cor_pred_obs"]])[, 'IDenv']
            ))))
        
        
        list_locations <-
          info_environments[match(list_envs, info_environments$IDenv), 'location']
        
        PA <-
          as.vector(sapply(fitting_all_splits, function(x)
            as.numeric(as.data.frame(x[['cor_pred_obs']])[, 'COR'])))
        
        df <- as.data.frame(cbind(list_locations, PA))
        
        colnames(df) <- c('location', 'Prediction_accuracy')
        
        df$Prediction_accuracy <-
          as.numeric(df$Prediction_accuracy)
        
        df2 <- count(df,location)
        
        p <-
          ggplot(df,
                 mapping = aes(
                   x = reorder(location, Prediction_accuracy),
                   y = Prediction_accuracy
                 )) + geom_boxplot() +
          geom_text(data = df2, aes(x=location,y = 1, label = paste0(n,' environments'))) +
          xlab('Location to predict (average over all years)') + ylab(paste0('Prediction accuracy for the trait ', trait)) + ggtitle('Leave-one-location-out CV scheme') +
          theme(axis.text.x = element_text(
            angle = 90,
            vjust = 0.5,
            hjust = 1
          ))
        ggsave(p,
               filename = paste0(
                 path_folder,
                 '/cv0_leave1locationout_show_location_',method_processing,'.pdf'
               ))
        
        
        PA <-
          as.vector(sapply(fitting_all_splits, function(x)
            as.character(unique(
              as.data.frame(x[["cor_pred_obs"]])[, 'COR']
            ))))
        
        
        df <- as.data.frame(cbind(list_envs, PA))
        
        colnames(df) <- c('IDenv', 'Prediction_accuracy')
        
        df$Prediction_accuracy <-
          as.numeric(df$Prediction_accuracy)
        
        p <-
          ggplot(df,
                 mapping = aes(
                   x = reorder(IDenv, Prediction_accuracy),
                   y = Prediction_accuracy,
                   group = 1
                 )) + geom_line() +
          xlab('Predicted environment') + ylab(paste0('Prediction accuracy for the trait ', trait)) + ggtitle('Leave-one-location-out CV scheme') +
          theme(axis.text.x = element_text(
            angle = 90,
            vjust = 0.5,
            hjust = 1
          ))
        ggsave(p,
               filename = paste0(path_folder, '/cv0_leave1locationout_show_env_',method_processing,'.pdf'))
      }
    }
    
    if (cv_type == 'cv1') {
      
      PA <-
        sapply(fitting_all_splits, function(x)
          as.numeric(x[['cor_pred_obs']]))
      
      df <- as.data.frame(PA)
      
      colnames(df) <- c('Prediction_accuracy')
      df$Prediction_accuracy <- as.numeric(df$Prediction_accuracy)
      
      p <-
        ggplot(df,
               mapping = aes(
                 x = 'CV1',
                 y = Prediction_accuracy,
                 group = 1
               )) + geom_boxplot() +
        xlab('CV1') + ylab(paste0('Average prediction accuracy for the trait ', trait,'\n using ',nb_folds_cv1,'-fold CV with ',repeats_cv1,'repeats')) + ggtitle('CV1 scheme') +
        theme(axis.text.x = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1
        ))
      ggsave(p, filename = paste0(path_folder, '/cv1_',method_processing,'.pdf'))
    }
    
    if (cv_type == 'cv2') {
      
      
      #VIP_all<- do.call('rbind',)
      # <- sapply(fitting_all_splits, function(x)
      #    as.numeric(x[['cor_pred_obs']]))
      
      
      
        
      p <-
        ggplot(df,
               mapping = aes(
                 x = 'CV2',
                 y = Prediction_accuracy,
                 group = 1
               )) + geom_boxplot() +
        xlab('CV2') + ylab(paste0('Average prediction accuracy for the trait ', trait,'\n using ',nb_folds_cv2,'-fold CV with ',repeats_cv2,'repeats')) + ggtitle('CV2 scheme') +
        theme(axis.text.x = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1
        ))
      ggsave(p, filename = paste0(path_folder, '/cv2_',method_processing,'.pdf'))
    }
    
  }

