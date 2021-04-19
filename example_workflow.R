data(geno_indica)
data(map_indica)
data(pheno_indica)
data(info_environments_indica)
data(env_data_indica)
METdata_indica <- create_METdata(geno=geno_indica,pheno=pheno_indica,env_data = env_data_indica,unique_EC_by_geno = F,compute_ECs = F,info_environments = info_environments_indica,map = map_indica)
#METdata_indica$geno<-METdata_indica$geno[,1:2000]

METdata_indica2 <- select_markers(METdata_indica,
                                 trait='PH',
                                 method_marker_effects = 'BLINK',
                                 method_selection = c('effect_size_per_env'),
                                 size_subset_most_variable_markers = 200,
                                 size_top_markers_by_env = 400,
                                 plot_penalty_regression_coefficients = T,
                                 path_save_plot =  "/home/uni08/jubin1/Data/PackageMLpredictions/plots")