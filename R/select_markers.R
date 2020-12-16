#default elasticnet
#but add rqtl + ridge reg ,and RF/xgboost later


select_markers<-function(MET_Data,
                         method='elasticnet',
                         folds=5,
                         repeats=1)
