data = readxl::read_xlsx(
  '~/results_learnmet_benchmarking/results_benchmarking_revisions_paper2.xlsx'
)
data[data$trait %in% c('GY', 'yld_bu_ac'), 'trait'] <- 'GY'
data[data$method %in% c('W-GW-GBLUP'), 'method']<- 'G-W-GW BLUP'
data[data$method == 'xgb_reg_1', 'method'] <- 'XGBoost (xgb_reg_1)'
data[data$method == 'stacking_reg_3', 'method'] <-
  'Stacked ensemble (stacked_reg_3)'
colnames(data)[2] <- 'Prediction method'
data$rmse<-as.numeric(data$rmse)

p2 <-
  ggplot(data = data[data$dataset %in% c('japonica', 'indica'),], aes(x = IDenv,
                                                                      y = rmse)) + geom_point(position = position_jitter(width = 0.35, height = 0),size = 3,
                                                                                             aes(colour = `Prediction method`,
                                                                                                 shape = dataset))  + theme(
                                                                                                   legend.title = element_text(size = 12),
                                                                                                   legend.text = element_text(size = 10),
                                                                                                   axis.title = element_text(size = 12),
                                                                                                   axis.text = element_text(size = 12),
                                                                                                   strip.text.x = element_text(size = 10),
                                                                                                   axis.text.x = element_text(
                                                                                                     angle = 90,
                                                                                                     vjust = 0.5,
                                                                                                     hjust = 1
                                                                                                   )
                                                                                                 ) + facet_wrap(~ trait,scales="free") + ylab('Root mean square error\n between predicted and\n observed values') 
ggsave("~/results_learnmet_benchmarking/rice_plot_rmse.PDF",width=25,height=12,units = 'cm')
p3 <-
  ggplot(data = data[data$dataset == 'G2F',], aes(x = IDenv,
                                                  y = rmse)) + geom_point(position = position_jitter(width = 0.2, height = 0),size = 3,
                                                                         aes(colour = `Prediction method`))  + theme(
                                                                           legend.title = element_text(size = 12),
                                                                           legend.text = element_text(size = 10),
                                                                           axis.title = element_text(size = 12),
                                                                           axis.text = element_text(size = 12),
                                                                           strip.text.x = element_text(size = 10),
                                                                           axis.text.x = element_text(
                                                                             angle = 90,
                                                                             vjust = 0.5,
                                                                             hjust = 1
                                                                           )
                                                                         ) + facet_wrap(~ trait, scales="free") + ylab('Root mean square error\n between predicted and\n observed values')

ggsave("~/results_learnmet_benchmarking/g2f_plot_rmse.PDF",width=25,height=12,units = 'cm')

gridExtra::grid.arrange(p2, p3, nrow = 3)



