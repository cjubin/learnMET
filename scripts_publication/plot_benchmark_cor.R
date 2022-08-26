data = readxl::read_xlsx(
  '~/results_learnmet_benchmarking/results_benchmarking_revisions_paper2.xlsx'
)
data[data$trait %in% c('GY', 'yld_bu_ac'), 'trait'] <- 'GY'
data[data$method %in% c('W-GW-GBLUP'), 'method']<- 'G-W-GW BLUP'
data[data$method == 'xgb_reg_1', 'method'] <- 'XGBoost (xgb_reg_1)'
data[data$method == 'stacking_reg_3', 'method'] <-
  'Stacked ensemble (stacked_reg_3)'
colnames(data)[2] <- 'Prediction method'


p2 <-
  ggplot(data = data[data$dataset %in% c('japonica', 'indica'),], aes(x = IDenv,
                                                                      y = cor)) + geom_point(position = position_jitter(width = 0.35, height = 0),size = 3,
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
                                                                                                 ) + facet_wrap(~ trait) + ylab('Pearson correlation\n between predicted and\n observed values') +
  ylim(c(-0.15, 1))

ggsave("~/results_learnmet_benchmarking/rice_plot_cor.PDF",width=25,height=12,units = 'cm')

p3 <-
  ggplot(data = data[data$dataset == 'G2F',], aes(x = IDenv,
                                                  y = cor)) + geom_point(position = position_jitter(width = 0.2, height = 0),size = 3,
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
                                                                         ) + facet_wrap(~ trait) + ylab('Pearson correlation\n between predicted and\n observed values') +
  ylim(c(-0.15, 1))


ggsave("~/results_learnmet_benchmarking/g2f_plot_cor.PDF",width=25,height=12,units = 'cm')


gridExtra::grid.arrange(p2, p3, nrow = 3)


data$benchmark_result <- NA
for (j in 1:nrow(data)) {
  data[j,'benchmark_result'] <- data[data$dataset==data$dataset[j] & data$`Prediction method` == 'G-W-GW BLUP' & data$IDenv == data$IDenv[j] & data$trait == data$trait[j],'cor']
}

data$comparison <-NA
data$comparison <- (data$cor-data$benchmark_result)/data$benchmark_result

