ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
        install.packages(new.pkg, repos = "http://cran.us.r-project.org")
    sapply(pkg, require, character.only = TRUE, quietly = TRUE)
}

residuewise_shifts <- function(tmp, nuclei=c("C1'","C2'","C3'","C4'","C5'","C2","C5","C6","C8","H1'","H2'","H3'","H4'","H2","H5","H5'","H5''","H6","H8")){
  # helper function used to output residue-wise chemical shifts
  shifts <- NULL
  for (nucleus in nuclei){
    shifts <- c(shifts,ifelse(!is.null(tmp$cs[tmp$nucleus==nucleus]), tmp$cs[tmp$nucleus==nucleus], NA))
  }
  return(shifts)
}

load_cs_data <- function(cs_file, train = TRUE){
  # function to read in chemical shift data and convert the table into a matrix where each row corresponds to data for a single residue
  suppressPackageStartupMessages(require(plyr))
  suppressPackageStartupMessages(require(bcv))
  nuclei_names <- c("C1p","C2p","C3p","C4p","C5p","C2","C5","C6","C8","H1p","H2p","H3p","H4p","H2","H5","H5p","H5pp","H6","H8")
  nuclei_rna <- c("C1'","C2'","C3'","C4'","C5'","C2","C5","C6","C8","H1'","H2'","H3'","H4'","H2","H5","H5'","H5''","H6","H8")
  ifelse (train, names <- c("id", "resid", "resname", nuclei_names), names <- c("resid", "resname", nuclei_names))

  # read in chemical shift data
  ifelse(train, cs_names <- c("resname", "resid", "nucleus", "cs", "error", "id"), cs_names <- c("resname", "resid", "nucleus", "cs", "error"))
  cs <- read.table(cs_file, col.names = cs_names)
  types <- as.vector(unique(cs$nucleus))
  types <- types[types %in% nuclei_rna]
  cs$error <- 1
  ifelse(train, vars <- c("id", "resid", "resname"), vars <- c("resid", "resname"))
  if(train){
    cs <- ddply(.data = cs, .variables = vars, .fun = residuewise_shifts, .progress = "text", nuclei = nuclei_rna)
  } else {
    cs <- ddply(.data = cs, .variables = vars, .fun = residuewise_shifts, nuclei = nuclei_rna)
  }

  colnames(cs) <- names
  cs <- cs[cs$resname %in% c("ADE", "GUA", "CYT", "URA"),]
  tmp <- cs
  colnames(cs) <- names
  cs <- as.data.frame(cs)
  cs <- droplevels(cs)
  return(list(cs=cs, nuclei=types))
}

add_neighboring_residues <- function(data, rna, single = TRUE, train = TRUE, i = 1){
  # select cs features
  if(train){
    data <- data[order(data$id, data$resid),]
    names_1 <- c("id","resid","resname","C1p","C2p","C3p","C4p","C5p","C2","C5","C6","C8","H1p","H2p","H3p","H4p","H2","H5","H5p","H5pp","H6","H8","class")
    names_2 <- c("id","resid","class")
    cs_data <- data[, names_1]
    bp_data <- data[, names_2]
  } else {
    data <- data[order(data$resid),]
    #data$id <- rna
    names_1 <- c("id","resid","resname","C1p","C2p","C3p","C4p","C5p","C2","C5","C6","C8","H1p","H2p","H3p","H4p","H2","H5","H5p","H5pp","H6","H8")
    cs_data <- data[, names_1]
  }

  #combine features from neighboring residues
  if(single){
    cs_data <- neighbors_single_RNA(cs_data, i)
  } else {
    cs_data <- neighbors_multiple_RNAs(cs_data, i)
  }
  if(train){
    data <- merge(cs_data, bp_data, by=c("id","resid"), all=T)
    names = get_feature_names(i, train=TRUE)
    data <- data[,names]
    data$resname <- as.factor(data$resname)
  } else {
    data <- cs_data
  }

  data <- data[order(data$id, data$resid),]
  return(data)
}

add_neighboring_resnames <- function(data){
  data = data[order(data$id,data$resid),]
  data_resnames = matrix(NA,nrow(data),3)
  for(i in 1:nrow(data)){
    resid_i = data$resid[i]
    data_resnames[i,2] = as.character(data[i,"resname"])
    resid_i_minus_1 = resid_i - 1
    resid_i_plus_1 = resid_i + 1
    if(resid_i_minus_1 %in% data$resid){
      data_resnames[i,1] = as.character(data[data$resid==resid_i_minus_1,"resname"])
    } else {
      data_resnames[i,1] = "None"
    }
    if(resid_i_plus_1 %in% data$resid){
      data_resnames[i,3] = as.character(data[data$resid==resid_i_plus_1,"resname"])
    } else {
      data_resnames[i,3] = "None"
    }
  }
  data = data[,!(names(data) %in% "resname")]
  data = data.frame(data,data_resnames)
  nuclei = c("C1p","C2p","C3p","C4p","C5p","C2","C5","C6","C8","H1p","H2p","H3p","H4p","H2","H5","H5p","H5pp","H6","H8",
              "C1p.1","C2p.1","C3p.1","C4p.1","C5p.1","C2.1","C5.1","C6.1","C8.1","H1p.1","H2p.1","H3p.1","H4p.1","H2.1","H5.1","H5p.1","H5pp.1","H6.1","H8.1",
              "C1p.2","C2p.2","C3p.2","C4p.2","C5p.2","C2.2","C5.2","C6.2","C8.2","H1p.2","H2p.2","H3p.2","H4p.2","H2.2","H5.2","H5p.2","H5pp.2","H6.2","H8.2")
  cnames = c("id","resid",nuclei,"resname_before","resname","resname_after")
  colnames(data) = cnames
  cnames_edit = c("id","resid","resname_before","resname","resname_after",nuclei)
  return(data[,cnames_edit])
}

neighbors_single_RNA <- function(sub_rna, i){
  n_res <- nrow(sub_rna)
  # only non-exchangeable cs
  nuclei <- c("C1p","C2p","C3p","C4p","C5p","C2","C5","C6","C8","H1p","H2p","H3p","H4p","H2","H5","H5p","H5pp","H6","H8")
  colmean <- t(colMeans(sub_rna[,nuclei], na.rm = TRUE))
  if(i == 0){comb_res = combine_neighbors_0(sub_rna,nuclei)}
  if(i == 1){comb_res = combine_neighbors_1(sub_rna,nuclei)}
  if(i == 2){comb_res = combine_neighbors_2(sub_rna,nuclei)}
  if(i == 3){comb_res = combine_neighbors_3(sub_rna,nuclei)}
  colnames(comb_res) = c("C1p","C2p","C3p","C4p","C5p","C2","C5","C6","C8","H1p","H2p","H3p","H4p","H2","H5","H5p","H5pp","H6","H8",
                         "id","resid","resname",
                         "C1p.1","C2p.1","C3p.1","C4p.1","C5p.1","C2.1","C5.1","C6.1","C8.1","H1p.1","H2p.1","H3p.1","H4p.1","H2.1","H5.1","H5p.1","H5pp.1","H6.1","H8.1",
                         "C1p.2","C2p.2","C3p.2","C4p.2","C5p.2","C2.2","C5.2","C6.2","C8.2","H1p.2","H2p.2","H3p.2","H4p.2","H2.2","H5.2","H5p.2","H5pp.2","H6.2","H8.2"
                         )
  colid = which(names(comb_res) %in% c("class"))
  if(length(colid)!=0){comb_res = comb_res[,-which(names(comb_res) %in% c("class"))]}
  names = get_feature_names(i, train = FALSE)
  comb_res <- comb_res[,names]
  return(comb_res)
}

neighbors_multiple_RNAs <- function(cs_data,i){
  # This function combines chemical shifts from neighboring residues for train set
  rna_names <- as.character(unique(cs_data$id))

  # combine features from i-1, i, i+1 residues
  n_rna <- length(rna_names)
  comb_rna <- NULL
  for (m in 1:n_rna){
    sub_rna <- subset(cs_data, id==unique(cs_data$id)[m])
    comb_res <- neighbors_single_RNA(sub_rna, i)
    comb_rna <- rbind(comb_rna, comb_res)
  }

  return(comb_rna)
}

get_feature_names <- function(i, train=TRUE){
  if(i==0){
    names = c("id","resid","resname","C1p","C2p","C3p","C4p","C5p","C2","C5","C6","C8",
              "H1p","H2p","H3p","H4p","H2","H5","H5p","H5pp","H6","H8")
  }
  if(i==1){
    names = c("id","resid","resname","C1p","C2p","C3p","C4p","C5p","C2","C5","C6","C8","H1p","H2p","H3p","H4p","H2","H5","H5p","H5pp","H6","H8",
              "C1p.1","C2p.1","C3p.1","C4p.1","C5p.1","C2.1","C5.1","C6.1","C8.1","H1p.1","H2p.1","H3p.1","H4p.1","H2.1","H5.1","H5p.1","H5pp.1","H6.1","H8.1",
              "C1p.2","C2p.2","C3p.2","C4p.2","C5p.2","C2.2","C5.2","C6.2","C8.2","H1p.2","H2p.2","H3p.2","H4p.2","H2.2","H5.2","H5p.2","H5pp.2","H6.2","H8.2")
  }
  if(i==2){
    names = c("id","resid","resname","C1p","C2p","C3p","C4p","C5p","C2","C5","C6","C8","H1p","H2p","H3p","H4p","H2","H5","H5p","H5pp","H6","H8",
              "C1p.1","C2p.1","C3p.1","C4p.1","C5p.1","C2.1","C5.1","C6.1","C8.1","H1p.1","H2p.1","H3p.1","H4p.1","H2.1","H5.1","H5p.1","H5pp.1","H6.1","H8.1",
              "C1p.2","C2p.2","C3p.2","C4p.2","C5p.2","C2.2","C5.2","C6.2","C8.2","H1p.2","H2p.2","H3p.2","H4p.2","H2.2","H5.2","H5p.2","H5pp.2","H6.2","H8.2",
              "C1p.3","C2p.3","C3p.3","C4p.3","C5p.3","C2.3","C5.3","C6.3","C8.3","H1p.3","H2p.3","H3p.3","H4p.3","H2.3","H5.3","H5p.3","H5pp.3","H6.3","H8.3",
              "C1p.4","C2p.4","C3p.4","C4p.4","C5p.4","C2.4","C5.4","C6.4","C8.4","H1p.4","H2p.4","H3p.4","H4p.4","H2.4","H5.4","H5p.4","H5pp.4","H6.4","H8.4")
  }
  if(i==3){
    names = c("id","resid","resname","C1p","C2p","C3p","C4p","C5p","C2","C5","C6","C8","H1p","H2p","H3p","H4p","H2","H5","H5p","H5pp","H6","H8",
              "C1p.1","C2p.1","C3p.1","C4p.1","C5p.1","C2.1","C5.1","C6.1","C8.1","H1p.1","H2p.1","H3p.1","H4p.1","H2.1","H5.1","H5p.1","H5pp.1","H6.1","H8.1",
              "C1p.2","C2p.2","C3p.2","C4p.2","C5p.2","C2.2","C5.2","C6.2","C8.2","H1p.2","H2p.2","H3p.2","H4p.2","H2.2","H5.2","H5p.2","H5pp.2","H6.2","H8.2",
              "C1p.3","C2p.3","C3p.3","C4p.3","C5p.3","C2.3","C5.3","C6.3","C8.3","H1p.3","H2p.3","H3p.3","H4p.3","H2.3","H5.3","H5p.3","H5pp.3","H6.3","H8.3",
              "C1p.4","C2p.4","C3p.4","C4p.4","C5p.4","C2.4","C5.4","C6.4","C8.4","H1p.4","H2p.4","H3p.4","H4p.4","H2.4","H5.4","H5p.4","H5pp.4","H6.4","H8.4",
              "C1p.5","C2p.5","C3p.5","C4p.5","C5p.5","C2.5","C5.5","C6.5","C8.5","H1p.5","H2p.5","H3p.5","H4p.5","H2.5","H5.5","H5p.5","H5pp.5","H6.5","H8.5",
              "C1p.6","C2p.6","C3p.6","C4p.6","C5p.6","C2.6","C5.6","C6.6","C8.6","H1p.6","H2p.6","H3p.6","H4p.6","H2.6","H5.6","H5p.6","H5pp.6","H6.6","H8.6")
  }
  if(train){
    names = c(names, "class")
  }
  return(names)
}

combine_neighbors_1 <- function(sub_rna,nuclei){
  comb_res = NULL
  for(m in 1:max(sub_rna$resid)){
    if(m %in% sub_rna$resid){
      # one or more of m+1, m-1 are NA
      if_m_minus_1 = ifelse((m-1)%in%sub_rna$resid, TRUE, FALSE)
      if_m_plus_1 = ifelse((m+1)%in%sub_rna$resid, TRUE, FALSE)
      if(if_m_minus_1){
        row_m_minus_1 = as.vector(sub_rna[sub_rna$resid==m-1,nuclei])
      } else {
        row_m_minus_1 = as.vector(sub_rna[sub_rna$resid==m,nuclei])
      }
      if(if_m_plus_1){
        row_m_plus_1 = as.vector(sub_rna[sub_rna$resid==m+1,nuclei])
      } else {
        row_m_plus_1 = as.vector(sub_rna[sub_rna$resid==m,nuclei])
      }
      tmp = cbind(row_m_minus_1, sub_rna[sub_rna$resid==m,], row_m_plus_1)
      #tmp[tmp=="NaN"]=NA
      comb_res <- rbind(comb_res, tmp)
    }
  }

  return(comb_res)
}

impute_cs_data <- function(data, id ="test", droplist = "id", speed="slow", ss_table = ss_table) {
  if(speed=="slow"){
    m = 1
    maxit = 50
  } else {
    m = 1
    maxit = 10
  }
  cs_train = get(load(ss_table))
  cs_train = cs_train[, -c(ncol(cs_train))] # delete 'class' column
  data$id = id
  cnames = colnames(cs_train)
  data = data[,cnames]
  cs = rbind(data, cs_train)
  if (length(intersect(names(cs), droplist)) < length(droplist)) {
    stop("Droplist variables not found in data set")
  }
  predictorMatrix <- (1 - diag(1, ncol(cs)))
  for (term in droplist) {
    drop.index <- which(names(cs) == term)
    predictorMatrix[, drop.index] <- 0
  }
  mids.out <- mice(cs, m = m, maxit = maxit, method = "pmm", seed = 500,
                   predictorMatrix = predictorMatrix, print=FALSE)
  cs_test = complete(mids.out,1)  # 'tidyr' also contains a function called 'complete'

  return(cs_test[1:nrow(data),])
}

normalize_test <- function(test, ss_table_1){
  # normalize test
  cols = c("C1p","C2p","C3p","C4p","C5p","C2","C5","C6","C8","H1p","H2p","H3p","H4p","H2","H5","H5p","H5pp","H6","H8",
          "C1p.1","C2p.1","C3p.1","C4p.1","C5p.1","C2.1","C5.1","C6.1","C8.1","H1p.1","H2p.1","H3p.1","H4p.1","H2.1","H5.1","H5p.1","H5pp.1","H6.1","H8.1",
          "C1p.2","C2p.2","C3p.2","C4p.2","C5p.2","C2.2","C5.2","C6.2","C8.2","H1p.2","H2p.2","H3p.2","H4p.2","H2.2","H5.2","H5p.2","H5pp.2","H6.2","H8.2")

  train = read.table(ss_table_1)

  colnames(train) = c("id","resid","resname_before","resname","resname_after",cols,"class")
  colnames(test) = c("id","resid","resname_before","resname","resname_after",cols)
  resnames = c('resname_before','resname','resname_after')

  # perform one hot encoding on residue names
  dmy <- dummyVars(" ~ .", data = train[resnames])
  resnames_train <- data.frame(predict(dmy, newdata = train[resnames]))
  resnames_test <- data.frame(predict(dmy, newdata = test[resnames]))

  # scale data
  scalar = preProcess(train[cols])
  cs_test = predict(scalar, test[cols])
  cs_final = as.matrix(data.frame(cs_test, resnames_test))
  return(cs_final)
}

format_cs_file <- function(data, nuclei){
  total = NULL
  colnames(data) = c("id","resid","resname",nuclei)
  for(nucleus in nuclei){
    tmp = data[,c("resname","resid",nucleus)]
    tmp$nucleus = nucleus
    names(tmp) = c("resname","resid","cs","nucleus")
    tmp = tmp[,c("resname","resid","nucleus","cs")]
    if(is.null(total)){
      total = tmp
      } else {
      total = rbind(total,tmp)
      }
  }
  print(total)
  total$error = "."
  return(total)
}

load_and_predict <- function(resid = resid, test = test){
  data = as.data.frame(resid)
  for(i in 1:6){
    model = load_model_hdf5(paste0("models/nn_model_whole_",i,".h5"))
    pred = as.vector(predict_proba(model, test))
    data = cbind(data, pred)
  }
  colnames(data) = c("resid","r1","r2","r3","r4","r5","r6")
  data$mean = rowMeans(data[,2:7])
  data[,2:8]=format(round(data[,2:8], 4), nsmall = 4)
  return(data)
}
