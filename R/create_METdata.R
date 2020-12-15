
## missing phenotypes for genotypes to predict should be indicated by NA


create_METdata<-function(geno=NULL,pheno=NULL,map=NULL,env_data=NULL,compute_ECs=FALSE){

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


  if(!is.matrix(pheno)){stop("phenotypic data is not a matrix")}






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

}
