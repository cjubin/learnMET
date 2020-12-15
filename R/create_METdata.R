
## missing phenotypes for genotypes to predict should be indicated by NA


create_METdata<-function(geno=NULL,pheno=NULL,map=NULL,env_data=NULL,unique_EC_by_geno=FALSE,compute_ECs=FALSE,coordinates_locations=NULL){

  # check if one object is missing

  if(is.null(geno)){stop("genotypic data not provided")}

  if(is.null(pheno)){stop("phenotypic data not provided")}

  # test format of the genotypic data

  if(!is.matrix(geno)){ stop("genotypic data not provided as a matrix") }

  if(!is.numeric(geno)){ stop("genotypic data not provided as numeric") }

  # test that all genotypes present in the phenotypic data are also present in the genotypic data

  if(!identical(as.character(as.vector(unique(pheno[,1]))),row.names(geno))){

    stop("lines identified in the phenotypic data not identical to lines identified in the genotypic data")

  }

  # if marker matrix is given, test that the marker names are the same in the map and in the marker genotype matrices

  if(!is.null(map)&!identical(colnames(geno), marker[,1])){

    stop("marker names in genotypic data and in map are not the same")

  }

  # test phenotypic data
  # test correct class for the different columns of the phenotype data
  if(!is.data.frame(pheno)){}
  if(ncol(pheno)<4){stop('MET pheno data should contain at least 4 columns: genotype lines (col1), year (col2), location (col3) and phenotype values (from col4)')}

  if(!is.character(pheno[,1])){stop("the genotype names/IDs (first column of pheno) must be character")}
  if(!is.numeric(pheno[,2])){stop("the year (second column of pheno) must be numeric")}
  if(!is.numeric(pheno[,2])){stop("the year (second column of pheno) must be numeric")}


  # Assign names 3 first pheno columns and transform year + location to factor
  colnames(pheno)[1:3]<-c('geno_ID','year','location')
  pheno$year=as.factor(pheno$year)
  pheno$location=as.factor(pheno$location)


  # Give a numerical trait name if no name provided
  if(is.null(colnames(pheno)[4:ncol(pheno)])){

    trait_names <- paste0('trait', 4:dim(pheno)[2])
    colnames(pheno) <- trait_names

  }

  # if geographical coordinates data.frame provided, test that all locations in the pheno data have their geographical coordinates included
  # test that longitude and latitude numerically provided

  if(!is.null(coordinates_locations)&!identical(unique(coordinates_locations[,1]),unique(pheno$location))){

    stop("locations identified in the geographical coordinates data.frame are not identical to the locations present in the phenotypic data.")

  }

  if(!is.numeric(coordinates_locations[,3])){stop("Longitude is not numeric")}

  if(!is.numeric(coordinates_locations[,4])){stop("Latitude is not numeric")}

  ##CHANGE NAME COL COORDINAtes

  # Create unique ID environment based on the location x year combination
  pheno$IDenv<-paste0(pheno$location,'_',pheno$year)
  coordinates_locations$IDenv<paste0(coordinates_locations$location,'_',coordinates_locations$year)


  # if marker data.frame provided, test marker names + chromosome info + positions provided

  if(!is.null(map)){
    if(!is.character(map[, 1])){

      stop("the marker name (first column in map) must be character")

    }

  if(!is.numeric(map[, 2])){

    stop("the chromosome number (second column in map) must be numeric")

  }


  if(!is.numeric(map[, 3])){

    stop("the genetic position (third column in map) must be numeric")

  }



  } else {cat('No map provided')}

  # test environmental data

  if(unique_EC_by_geno==FALSE&&!is.null(env_data)&&nrow(env_data)!=length(unique(pheno$IDenv))){
    stop('The number of observations in the environmental data does not match the number of Year x Location combinations from the pheno file.')
  }


  if(unique_EC_by_geno==TRUE&&!is.null(env_data)&&nrow(env_data)!=nrow(pheno)){
    stop('The number of observations in the environmental data does not match the number of observations in the pheno file (Year x Location x Genotype).')
  }

  ## Test specific format for environmental data: if genotype based environmental covariates (taking phenology into account to compute ECs)

  if(unique_EC_by_geno==TRUE&&!is.null(env_data)&&!is.character(env_data[,1])){stop('The first column of environmental data should contain the genotype names/IDs as character.')}
  if(unique_EC_by_geno==TRUE&&!is.null(env_data)&&!is.numeric(env_data[,2])){stop('The second column of environmental data should contain the year as numeric.')}
  if(unique_EC_by_geno==TRUE&&!is.null(env_data)&&!is.character(env_data[,3])){stop('The third column of environmental data should contain the location as character.')}
  if(unique_EC_by_geno==TRUE&&!is.null(env_data)&&(!is.factor(env_data[,4:ncol(env_data)])||!is.numeric(env_data[,4:ncol(env_data)]))){stop('Col4+ of environmental data should contain the environmental variable as numeric or factor variable.')}

  if(unique_EC_by_geno==FALSE&&!is.null(env_data)&&!is.numeric(env_data[,1])){stop('The first column of environmental data should contain the year as numeric.')}
  if(unique_EC_by_geno==FALSE&&!is.null(env_data)&&!is.character(env_data[,2])){stop('The second column of environmental data should contain the location as character.')}
  if(unique_EC_by_geno==FALSE&&!is.null(env_data)&&(!is.factor(env_data[,3:ncol(env_data)])||!is.numeric(env_data[,3:ncol(env_data)]))){stop('Col3+ of environmental data should contain the environmental variable as numeric or factor variable.')}

  if(!is.null(env_data)){

  if(unique_EC_by_geno==TRUE){colnames(env_data)[1:3]<-c('geno_ID','year','location')}
  if(unique_EC_by_geno==FALSE){colnames(env_data)[1:2]<-c('year','location')}

  env_data$IDenv<paste0(env_data$location,'_',env_data$year)


  # Match the coordinates information with the environmental data

  env_data<-merge(env_data,coordinates_locations,by='IDenv',all.x = T)

  } else{
    env_data<-NULL}


  METpred_data<-list('geno'=geno,
                     'pheno'=pheno,
                     'compute_ECs'=compute_ECs,
                     'env_data'=env_data,
                     'map_markers'=map
                     )
  return(METpred_data)


}
