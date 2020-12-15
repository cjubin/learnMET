
## missing phenotypes for genotypes to predict should be indicated by NA


create_METdata<-function(geno=NULL,pheno=NULL,map=NULL,env_data=NULL,compute_ECs=FALSE,coordinates_locations=NULL){

  # check if one object is missing

  if(is.null(geno)){stop("genotypic data not provided")}

  if(is.null(pheno)){stop("phenotypic data not provided")}

  # test format of the genotypic data

  if(!is.matrix(geno)){ stop("genotypic data not provided as a matrix") }

  if(!is.numeric(geno)){ stop("genotypic data not provided as numeric") }

  # test that all genotypes present in the phenotypic data are also present in the genotypic data

  if(!identical(unique(phenos[,1]),row.names(geno))){

    stop("lines identified in the phenotypic data not identical to lines identified in th genotypic data")

  }

  # if marker matrix is given, test that the marker names are the same in the map and in the marker genotype matrices

  if(!is.null(map)&!identical(colnames(geno), marker[,1])){

    stop("marker names in genotypic data and in map are not the same")

  }

  # test phenotypic data
  # test correct class for the different columns of the phenotype data
  if(!is.data.frame(pheno)){}
  if(ncol(pheno)<4){stop('MET pheno data should contain at least 4 columns: genotype lines (col1), year (col2), location (col3) and phenotype values (from col4)')}

  if(!is.character(pheno[,1])){stop("the genotype names (first column of pheno) must be character")}
  if(!is.numeric(pheno[,2])){stop("the year (second column of pheno) must be numeric")}
  if(!is.numeric(pheno[,2])){stop("the year (second column of pheno) must be numeric")}


  # Assign names 3 first pheno columns and transform year + location to factor
  colnames(pheno)[1:3]<-c('geno_names','year','location')
  phenos$year=as.factor(phenos$year)
  phenos$location=as.factor(phenos$location)


  # Give a numerical trait name if no name provided
  if(is.null(colnames(pheno)[4:ncol(pheno)])){

    trait_names <- paste0('trait', 4:dim(pheno)[2])
    colnames(pheno) <- trait_names

  }

  # if geographical coordinates data.frame provided, test that all locations have their geographical coordinates included

  if(!is.null()  dentical(unique(phenos[,1]),row.names(geno))){

    stop("lines identified in the phenotypic data not identical to lines identified in th genotypic data")

  }


  # if marker data.frame provided, test marker names + chromosome info + position (bp) provided

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




}
