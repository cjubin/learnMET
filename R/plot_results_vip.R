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
        
        predicted_years<- paste0('Predicted year: ',list_predicted_years)
        names(predicted_years) <- list_predicted_years
        
        
        
        p <- ggplot(VIP, aes(x=reorder_within(Variable,Importance,year), y=Importance)) + ylab('Relative importance') + xlab('Top 40 predictor variables\n for each training set')+
          geom_boxplot() + facet_wrap(~ year, ncol = length(list_predicted_years), scales = "free",labeller = as_labeller(predicted_years)) + coord_flip()
        
        

        ggsave(
          p,
          filename = paste0(path_folder, '/cv0_leave1yearout_', method_processing, '_Variable_Importance.pdf'),
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
        list_predicted_locations <-
          as.vector(sapply(fitting_all_splits, function(x)
            as.character(unique(x[['predictions_df']][,'location']))))
        
        VIP <-
          sapply(fitting_all_splits, function(x)
            x['ranking_vip'])
        
        
        for (j in 1:length(list_predicted_locations)) {
          VIP[[j]]$location <- list_predicted_locations[[j]]
          VIP[[j]] <- top_n(VIP[[j]],wt= Importance,n=40)
          
        }
        
        VIP <- do.call('rbind',VIP)
        
        predicted_loc <- paste0('Predicted location: ',list_predicted_locations)
        names(predicted_loc) <- list_predicted_locations
        
        
        
        p <- ggplot(VIP, aes(x=reorder_within(Variable,Importance,location), y=Importance)) + ylab('Relative importance') + xlab('Top 40 predictor variables\n for each training set')+
          geom_boxplot() + facet_wrap(~ location, ncol = length(list_predicted_locations), scales = "free",labeller = as_labeller(predicted_loc)) + coord_flip()
        
        ggsave(
          p,
          filename = paste0(path_folder, '/cv0_leave1locationout_', method_processing, '_Variable_Importance.pdf'),
          height = 8,
          width = 12
        )
        
        
      }
    }
    
    if (cv_type == 'cv1') {
      
      VIP <-
        sapply(fitting_all_splits, function(x)
          x['ranking_vip'])
      
      VIP <- do.call('rbind',VIP)
      
      for (j in unique(VIP$Variable)) {
        if (length(which(VIP$Variable==j))<6){
          print(j)
          m <- 6-length(which(VIP$Variable==j))
          supp <- matrix(c(j,0),nrow = m,ncol = 2,byrow = T)
          colnames(supp) <- colnames(VIP)
          VIP <- rbind(VIP,supp)
        }
        
      }
      
      VIP$Importance <- as.numeric(VIP$Importance)
      
      VIP <- as.data.frame(VIP %>%
        group_by(Variable) %>%
        dplyr::mutate(Mean = mean(Importance, na.rm=TRUE)))
      
      
        
      VIP_selected_var <- as.data.frame(unique(VIP[,c(1,3)])) %>% top_n(.,wt= Mean,n=40)
      
      VIP <- VIP[VIP$Variable%in%VIP_selected_var$Variable,]
      
      VIP$Mean <- as.numeric(VIP$Mean)
      
      p <- ggplot(VIP, aes(x=reorder(Variable,Importance), y=Importance)) + ylab('Average relative importance over models fitted on training sets from CV1') + xlab('Top 40 predictor variables\n (based on average relative importance)')+
        geom_boxplot()  + coord_flip()
      
      
      ggsave(
        p,
        filename = paste0(path_folder, '/cv1_', method_processing, '_Variable_Importance.pdf'),
        height = 8,
        width = 12
      )
      
      
      
    }
    
    if (cv_type == 'cv2') {
      
      VIP <-
        sapply(fitting_all_splits, function(x)
          x['ranking_vip'])
      
      VIP <- do.call('rbind',VIP)
      
      for (j in unique(VIP$Variable)) {
        if (length(which(VIP$Variable==j))<6){
          print(j)
          m <- 6-length(which(VIP$Variable==j))
          supp <- matrix(c(j,0),nrow = m,ncol = 2,byrow = T)
          colnames(supp) <- colnames(VIP)
          VIP <- rbind(VIP,supp)
        }
        
      }
      
      VIP$Importance <- as.numeric(VIP$Importance)
      
      VIP <- as.data.frame(VIP %>%
                             group_by(Variable) %>%
                             dplyr::mutate(Mean = mean(Importance, na.rm=TRUE)))
      
      
      
      VIP_selected_var <- as.data.frame(unique(VIP[,c(1,3)])) %>% top_n(.,wt= Mean,n=40)
      
      VIP <- VIP[VIP$Variable%in%VIP_selected_var$Variable,]
      
      VIP$Mean <- as.numeric(VIP$Mean)
      
      p <- ggplot(VIP, aes(x=reorder(Variable,Importance), y=Importance)) + ylab('Average relative importance over models fitted on training sets from CV2') + xlab('Top 40 predictor variables\n (based on average relative importance)')+
        geom_boxplot()  + coord_flip()
      
      
      ggsave(
        p,
        filename = paste0(path_folder, '/cv2_', method_processing, '_Variable_Importance.pdf'),
        height = 8,
        width = 12
      )
      
    }
    
  }

