my.vioplot <- function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL,
                        horizontal = FALSE, col = "magenta", border = "black", lty = 1,
                        lwd = 1, rectCol = "black", colMed = "white", pchMed = 19,
                        at, add = FALSE, wex = 1, drawRect = TRUE, las = 0, cex.axis = 1, aty = NULL, labely = NULL, lasy = 1) {
  require(sm)
  datas <- list(x, ...)
  n <- length(datas)
  if (missing(at))
    at <- 1:n
  upper <- vector(mode = "numeric", length = n)
  lower <- vector(mode = "numeric", length = n)
  q1 <- vector(mode = "numeric", length = n)
  q3 <- vector(mode = "numeric", length = n)
  med <- vector(mode = "numeric", length = n)
  base <- vector(mode = "list", length = n)
  height <- vector(mode = "list", length = n)
  baserange <- c(Inf, -Inf)
  args <- list(display = "none")
  if (!(is.null(h)))
    args <- c(args, h = h)
  for (i in 1:n) {
    data <- datas[[i]]
    data.min <- min(data)
    data.max <- max(data)
    q1[i] <- quantile(data, 0.25)
    q3[i] <- quantile(data, 0.75)
    med[i] <- median(data)
    iqd <- q3[i] - q1[i]
    upper[i] <- min(q3[i] + range * iqd, data.max)
    lower[i] <- max(q1[i] - range * iqd, data.min)
    est.xlim <- c(min(lower[i], data.min), max(upper[i],
                                               data.max))
    smout <- do.call("sm.density", c(list(data, xlim = est.xlim),
                                     args))
    hscale <- 0.4/max(smout$estimate) * wex
    base[[i]] <- smout$eval.points
    height[[i]] <- smout$estimate * hscale
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1], t[1])
    baserange[2] <- max(baserange[2], t[2])
  }
  if (!add) {
    xlim <- if (n == 1)
      at + c(-0.5, 0.5)
    else range(at) + min(diff(at))/2 * c(-1, 1)
    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  if (is.null(names)) {
    label <- 1:n
  }
  else {
    label <- names
  }
  boxwidth <- 0.05 * wex
  if (!add)
    plot.new()
  if (!horizontal) {
    if (!add) {
      plot.window(xlim = xlim, ylim = ylim)
      print(aty)
      axis(2, at = aty, label = labely, cex.axis = cex.axis, las = lasy )
      axis(1, at = at, label = label, cex.axis=cex.axis, las = las)
    }
    box()
    for (i in 1:n) {
      polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])),
              c(base[[i]], rev(base[[i]])), col = col[i], border = border,
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd,
              lty = lty)
        rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2,
             q3[i], col = rectCol)
        points(at[i], med[i], pch = pchMed, col = colMed)
      }
    }
  }
  else {
    if (!add) {
      plot.window(xlim = ylim, ylim = xlim)
      axis(1,cex.axis=1.5, at = aty, label = labely)
      axis(2, at = at, label = label, cex.axis=1.5)
    }
    box()
    for (i in 1:n) {
      polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]],
                                              rev(at[i] + height[[i]])), col = col[i], border = border,
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd,
              lty = lty)
        rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] +
               boxwidth/2, col = rectCol)
        points(med[i], at[i], pch = pchMed, col = colMed)
      }
    }
  }
  invisible(list(upper = upper, lower = lower, median = med,
                 q1 = q1, q3 = q3))
}


fraction_greater <- function(threshold, data){sum(data < threshold)}

plot_color_bar <- function(data, metric = "spec"){
  suppressPackageStartupMessages(require(tidyverse))
  data$group <- data$type
  data$group <- as.factor(data$group)
  data$value <- data[, metric]
  data$individual <- data$id
  data$id <- 1:nrow(data)
  data$colour <- "red"
  data$colour[data$group == "both"] <- "blue"
  # Set a number of 'empty bar' to add at the end of each group
  empty_bar=10
  to_add = data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
  colnames(to_add) = colnames(data)
  to_add$group=rep(levels(data$group), each=empty_bar)
  data=rbind(data, to_add)
  data = data %>% arrange(group, value)
  data$id=seq(1, nrow(data))
  # Get the name and the y position of each label
  label_data=data
  number_of_bar=nrow(label_data)
  angle= 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  label_data$hjust <- ifelse( angle < -90, 1, 0)
  label_data$angle<-ifelse(angle < -90, angle+180, angle)

  # prepare a data frame for base lines
  base_data=data %>%
    group_by(group) %>%
    summarize(start=min(id), end=max(id) - empty_bar) %>%
    rowwise() %>%
    mutate(title=mean(c(start, end)))
  base_data$mean <- mean(data$value[data$group == "both"], na.rm = TRUE)
  base_data$mean[base_data$group == "proton"] <- mean(data$value[data$group == "proton"], na.rm = TRUE)

  print(base_data)
  # prepare a data frame for grid (scales)
  grid_data = base_data
  grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
  grid_data$start = grid_data$start - 1
  grid_data=grid_data[-1,]

  # Make the plot
  p = ggplot(data, aes(x=as.factor(id), y=value, fill=data$colour)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar

    geom_bar(aes(x=as.factor(id), y=value, fill=data$colour), stat="identity", alpha=0.5) +

    # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.

    geom_segment(data=grid_data, aes(x = end, y = 0.80, xend = start, yend = 0.80), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 0.60, xend = start, yend = 0.60), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 0.40, xend = start, yend = 0.40), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 0.20, xend = start, yend = 0.20), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 0.00, xend = start, yend = 0.00), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +

    # Add text showing the value of each 100/75/50/25 lines
    annotate("text", x = rep(max(data$id), 6), y = c(0.00, 0.20, 0.40, 0.60, 0.80, 1.00), label = c("0.00", "0.20", "0.40", "0.60", "0.80", "1.00") , color="black", size=3.5 , angle=0, fontface="bold", hjust=1) +

    geom_bar(aes(x=as.factor(id), y=value, fill=group, colour = "black"), stat="identity", colour = "black") +
    ylim(-0.5,1.5) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(0,4), "cm")
    ) +
    coord_polar() +
    geom_text(data=label_data, aes(x=id, y=value+0.1, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +

    # Add base line information
    geom_segment(data=base_data, aes(x = start, y = 0.5, xend = end, yend = 0.5 ), colour = "red", alpha=0.8, size=0.5 , inherit.aes = FALSE, linetype= 2) +
    geom_segment(data=base_data, aes(x = start, y = mean, xend = end, yend = mean ), colour = "white", alpha=0.8, size=0.5 , inherit.aes = FALSE, linetype= 2) +
    geom_segment(data=base_data, aes(x = start, y = 0, xend = end, yend = 0 ), colour = "black", alpha=0.8, size=0.7 , inherit.aes = FALSE )
  p
}

draw_neural_network <- function(seed.val=2, num.vars=10, num.obs=100){
  suppressPackageStartupMessages(require(devtools))
  suppressPackageStartupMessages(require(clusterGeneration))
  suppressPackageStartupMessages(require(RSNNS))
  suppressPackageStartupMessages(require(nnet))
  
  suppressMessages(source_url('https://gist.githubusercontent.com/fawda123/7471137/raw/466c1474d0a505ff044412703516c34f1a4684a5/nnet_plot_update.r'))
  set.seed(seed.val)

  #input variables
  cov.mat<-genPositiveDefMat(num.vars,covMethod=c("unifcorrmat"))$Sigma
  rand.vars<-mvrnorm(num.obs,rep(0,num.vars),Sigma=cov.mat)

  #output variables
  parms<-runif(num.vars,-10,10)
  y1<-rand.vars %*% matrix(parms) + rnorm(num.obs,sd=20)
  parms2<-runif(num.vars,-10,10)
  y2<-rand.vars %*% matrix(parms2) + rnorm(num.obs,sd=20)

  #final datasets
  rand.vars<-data.frame(rand.vars)
  resp<-data.frame(y1,y2)
  names(resp)<-c('Y1','Y2')
  dat.in<-data.frame(resp,rand.vars)

  #neural net with three hidden layers, 9, 11, and 8 nodes in each
  mod<-mlp(rand.vars, resp, size=c(10,5))
  par(mar=numeric(4),family='serif')
  plot.nnet(mod, node.labs = FALSE, circle.col = "blue", pch=21, var.labs = FALSE, circle.cex = 4)
}

load_data <- function(file, colnames = NULL, header = FALSE){
  if(!is.null(colnames)){
    header <- FALSE
    return(read.table(file = file, col.names = colnames, header = header))
  } else {
    return(read.table(file = file, header = header))
  }

}

load_model_accuracy <- function(file = "data/NN_accuracy.txt", colnames = c("id", "sen", "spec", "overall"), header = FALSE){
  data_test <- load_data(file = file, colnames = colnames, header = header)
  load(file = "data/proton_only.RData")
  data_test$type <- "both"
  data_test$type[data_test$id %in% proton_only] <- "proton"
  data_test$id <- as.character(data_test$id)
  return(data_test)
}

print_model_summary <- function(data, metrics = c("sen", "spec", "overall")){
  # prints min, max, median
  for (metric in metrics){
    values <- data[, metric]
    cat(sprintf("metric: %7s min: %4.2f max: %4.2f mean: %4.2f median: %4.2f\n", metric, min(values), max(values), mean(values), median(values)))
  }
}

cs_filter <- function(data, filter1 = "consistency", filter2 = "energy", remove_native = FALSE){
  # does cs-filtering
  data <- data[complete.cases(data), ]
  if(remove_native){
    data <- data[!(data$PPV == 1 & data$sens == 1), ]
  }
  data <- data[order(data[, filter1], decreasing = TRUE), ]
  data <- data[data[, filter1] == max(data[, filter1]), ]
  data <- data[order(data[, filter2], decreasing = FALSE), ]
  return(data[1, ])
}


filter <- function(data, remove_native = FALSE, sortBy = "energy", decreasing = FALSE){
  # does cs-filtering
  data <- data[complete.cases(data), ]
  if(remove_native){
    data <- data[!(data$PPV == 1 & data$sens == 1), ]
  }
  data <- data[order(data[, sortBy], decreasing = decreasing), ]
  return(data[1, ])
}

get_results <- function(data, selection, selectBy, sortBy = NULL){
  data <- unique(data[(data[, selectBy] %in% selection), ])
  if(!is.null(sortBy)){data <- data[order(data[, sortBy]), ]}
  rownames(data) <- 1:nrow(data)
  return(data)
}

calculate_classification_metrics <- function(filename = "data/non_canonical/pred/thresholds/mean_1A60_0.5.txt", colnames = c("resid","true","pred"), metric = "sty"){
  data <- read.table(file = filename, col.names = colnames)
  data$true = factor(data$true, levels = sort(unique(data$true)))
  data$pred = factor(data$pred, levels = sort(unique(data$pred)))
  cfm <- confusionMatrix(data$pred, data$true, positive = '0')
  if(metric == "sty"){
    return(unname(cfm$byClass[1]))
  }
  if(metric == "spf"){
    return(unname(cfm$byClass[2]))
  }
  if(metric == "acc"){
    return(unname(cfm$overall[1]))
  }
}

thresholds_analysis <- function(thres = "0.5", rnas, metric = "sty", method="mean"){
  sum <- 0
  for(rna in rnas){
    sum = sum + calculate_classification_metrics(filename = paste0("data/non_canonical/pred/thresholds/",method,"_",rna,"_",thres,".txt"), metric = metric)
  }
  return(sum/length(rnas))
}

load_base_pairing_predictions <- function(rnas = NULL, rna = NULL, pred_cols = c("resid","bp_pred"), all_rnas = TRUE, method = "residue"){
  comb <- NULL
  if(!all_rnas){
    pred <- read.table(paste0("data/pred/",rna,"_pred_avg_binary_0.6.txt"), stringsAsFactors=F,
                       col.names=pred_cols)
    return(pred)
  } else {
    for(rna in rnas){
      ref <- read.table(paste0("info/2D_canonical/",rna,".ct"),skip=1,stringsAsFactors=F,
                       col.names=c("resid","resname","before","after","bp","resid2"))
      pred <- read.table(paste0("data/pred/",rna,"_pred_avg_binary_0.6.txt"), stringsAsFactors=F,
                        col.names=pred_cols)
      merged <- merge(ref, pred, by=c("resid"), all = T)
      merged$id <- rna
      merged <- merged[,c("id","resid","resname","bp",pred_cols[2])]
      if(method == "basepair"){
        # consider which residue it is base paired to
        for(i in 1:nrow(merged)){
          # if there is a paired residue, write its residue name and predicted base pairing status
          if(merged$bp[i] != 0){
            merged$resname_bp_ref[i] <- merged$resname[merged$resid == merged$bp[i]]
            merged$paired_bps_pred[i] <- merged[,pred_cols[2]][merged$resid == merged$bp[i]]
          } else {
            merged$resname_bp_ref[i] <- 0
            merged$paired_bps_pred[i] <- NA
          }
          if(merged$resname_bp_ref[i]!=0){
            # then write the base pair type
            merged$bp_type[i] <- paste(merged$resname[i],merged$resname_bp_ref[i], sep = "")  
          } else {
            merged$bp_type[i] <- 'none'
          }
          # GC and CG is the same thing
          if(merged$bp_type[i] == "CG"){merged$bp_type[i] <- "GC"}
          else if(merged$bp_type[i] == "UA"){merged$bp_type[i] <- "AU"}
          else if(merged$bp_type[i] == "UG"){merged$bp_type[i] <- "GU"}
        }
        merged$bp_type <- as.factor(merged$bp_type)
      } else {
        # only consider the base pairing status
        merged$class <- ifelse(merged$bp != 0, 0, 1)
      }
      if(is.null(comb)){
        comb <- merged
      } else {
        comb <- rbind(comb, merged)
      }
    }
    comb <- comb[order(comb$id,comb$resid),]
    return(comb)
  }
}

residue_wise_cs2bps_metrics <- function(residue, data, col_ref = "class", col_pred = "bp_pred", metric = "TPR"){
  data <- data[!is.na(data[,col_pred]),]
  x <- subset(data, resname == residue)
  if(metric == "number"){
    return(nrow(x))
  } else {
    pred <- as.factor(x[,col_pred])
    ref <- as.factor(x[,col_ref])
    conf <- caret::confusionMatrix(pred, ref, positive="0")
    if(metric == "TPR"){
      return(unname(round(conf$byClass["Sensitivity"],3)))
    }
    if(metric == "TNR"){
      return(unname(round(conf$byClass["Specificity"],3)))
    }
  }
}

basepair_wise_cs2bps_metrics <- function(basepair, data, col_pred = "paired_bps_pred", col_ref = "bp_pred", metric = "number"){
  # only consider residues with cs2bps prediction and paired residues with cs2bps prediction
  data <- data[!is.na(data[,col_ref]),] 
  data <- data[!is.na(data[,col_pred]),]
  x <- subset(data, bp_type == basepair)
  if(metric == "number"){
    return(nrow(x))
  }
  if(metric == "TPR"){
    return(round(sum(x$bp_pred == 0 & x$paired_bps_pred == 0)/nrow(x),2))
  }
}

structure_selection_using_cs2bps_and_energy <- function(rnas, from_path, to_path, methods, rewrite = F){
  # get energy of each predicted structure
  eng = matrix(NA, length(rnas), length(methods))
  for(i in 1:length(rnas)){
    rna = rnas[i]
    for(j in 1:length(methods)){
      FoldEnergy = read.table(paste0("data/ss_energy/energy_", methods[j], "_", rna, ".txt"))
      eng[i, j] = FoldEnergy[1,5]
    }
  }
  colnames(eng) = methods
  rownames(eng) = rnas
  
  # select structure based on consistency
  for(m in 1:108){
    rna = rnas[m]
    cs2bps = read.table(paste0("data/pred/",rna,"_pred_avg_binary_0.6.txt"),col.names = c("resid","bp_cs2bps"))
    x = matrix(NA,length(methods),3)
    
    for(i in 1:length(methods)){
      method = methods[i]
      csfold = read.table(paste0(from_path,method,"_",rna,".ct"),stringsAsFactors = F,
                           fill=T,skip=1,col.names = c("resid","resname","resid_before","resid_after","bp","resid2"))
      if(length(which(csfold$resname=="ENERGY"))!=0){
        # several structures, use lowest energy structure
        csfold <- csfold[1:(which(csfold$resname=="ENERGY")[1]-1),]
      }

      comp = merge(csfold,cs2bps,by="resid")
      comp$bp_csfold = ifelse(comp$bp==0, 1, 0)
      x[i,1] = method
      x[i,2] = sum(comp$bp_cs2bps==comp$bp_csfold)/nrow(comp)
      x[i,3] = eng[m,i]
    }
    x = as.data.frame(x)
    x$V2 = as.numeric(as.character(x$V2))
    x$V3 = as.numeric(as.character(x$V3))
    x = x[order(-x$V2,x$V3),]
    model = x[1,1]
    #print(paste0(rna, ": ", as.character(model)))
    
    # copy file to destination folder
    if(rewrite){
      file.copy(paste0(from_path,model,"_",rna,".ct"), paste0(to_path,rna,".ct"),overwrite = TRUE)
    }
  }
}

load_secondary_structure_predictions <- function(rnas, ss_path = "data/ss_selected/", ref_path = "info/2D_canonical/"){
  data <- NULL
  # look at all RNAs
  for(i in 1:length(rnas)){
    rna <- rnas[i]
    ref <- read.table(paste0(ref_path, rna, ".ct"),skip=1,stringsAsFactors=F,
                     col.names=c("resid","resname","before","after","bp","resid2"))
    ref[ref$resname == 'g',"resname"] <- 'G'
    ref[ref$resname == 'c',"resname"] <- 'C'
    ref[ref$resname == 'u',"resname"] <- 'U'
    ref[ref$resname == 'a',"resname"] <- 'A'
    
    pred <- read.table(paste0(ss_path, rna, ".ct"),skip=1,stringsAsFactors=F,
                      fill=T,col.names=c("resid","resname","before","after","bp_pred","resid2"))
    if(length(which(pred$resname=="ENERGY"))!=0){
      # several structures, use lowest energy structure
      pred <- pred[1:(which(pred$resname=="ENERGY")[1]-1),]
    }
    pred[pred$resname == 'g',"resname"] <- 'G'
    pred[pred$resname == 'c',"resname"] <- 'C'
    pred[pred$resname == 'u',"resname"] <- 'U'
    pred[pred$resname == 'a',"resname"] <- 'A'
    
    merged <- merge(ref, pred, by = c("resid","resname","before","after"), all = T)
    merged <- merged[, -which(colnames(merged)%in%c("resid2.x","resid2.y"))]
    merged <- merged[order(merged$resid),]
    for(j in 1:nrow(merged)){
      if(merged$bp[j]!=0){
        merged$resname_bp_ref[j] = merged$resname[merged$resid==merged$bp[j]]
      } else {
        merged$resname_bp_ref[j] = 0
      }
      if(merged$bp_pred[j]!=0){
        merged$resname_bp_pred[j] = merged$resname[merged$resid==merged$bp_pred[j]]  
      } else {
        merged$resname_bp_pred[j] = 0
      }
      # write basepair type
      if(merged$resname_bp_ref[j]!=0){
        # then write the base pair type
        merged$bp_type_ref[j] <- paste(merged$resname[j],merged$resname_bp_ref[j], sep = "")  
      } else {
        merged$bp_type_ref[j] <- 'none'
      }
      if(merged$resname_bp_pred[j]!=0){
        # then write the base pair type
        merged$bp_type_pred[j] <- paste(merged$resname[j],merged$resname_bp_pred[j], sep = "")  
      } else {
        merged$bp_type_pred[j] <- 'none'
      }
      
      # CG -> GC UA -> UA UG -> GU
      if(merged$bp_type_ref[j] == "CG"){merged$bp_type_ref[j] <- "GC"}
      else if(merged$bp_type_ref[j] == "UA"){merged$bp_type_ref[j] <- "AU"}
      else if(merged$bp_type_ref[j] == "UG"){merged$bp_type_ref[j] <- "GU"}
      if(merged$bp_type_pred[j] == "CG"){merged$bp_type_pred[j] <- "GC"}
      else if(merged$bp_type_pred[j] == "UA"){merged$bp_type_pred[j] <- "AU"}
      else if(merged$bp_type_pred[j] == "UG"){merged$bp_type_pred[j] <- "GU"}
    }
    
    
    merged$id <- rna
    if(is.null(data)){
      data <- merged
    } else {
      data <- rbind(data, merged)
    }
  } 
  data$bp_type_ref <- as.factor(data$bp_type_ref)
  data$bp_type_pred <- as.factor(data$bp_type_pred)
  return(data)
}

basepair_wise_csfold_metrics <- function(basepair, data, col_ref = "bp_type_ref", col_pred = "bp_type_pred", metric = "number"){
  # only consider residues with cs2bps prediction
  x <- data[data[,col_ref] == basepair,]
  if(metric == "number"){
    return(nrow(x))
  } 
  if(metric == "TPR"){
    return(round(sum(x[,col_pred] == basepair)/nrow(x),2))
  }
}

load_scorer_accuracy <- function(file, colnames = c("id","TPR","PPV")){
  data <- read.table(file, fill = T, col.names = c("id","TPR","PPV"))
  data <- data[data$id%in%rnas,] # delete sub-optimal structures (if they exist)
  data$TPR <- as.numeric(gsub("*%","",as.character(data$TPR)))/100
  data$PPV <- as.numeric(gsub("*%","",as.character(data$PPV)))/100
  return(data)
}