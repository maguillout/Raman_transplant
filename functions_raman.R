#### PACKAGES INSTALLATION #####

list.of.packages <- c('rmarkdown', "knitr", "mdatools", "dplyr", "RamanMP", "reshape2", "ggplot2", "plotly", "readxl","baseline", "randomForest")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

########### Packages ###########  
library("rmarkdown")
library("knitr")
library("mdatools")
library("dplyr")
library("RamanMP")
library("reshape2")
library("ggplot2")
library("plotly")
library("stringr")
library("readxl")
library("mixOmics")
library("baseline")
library("randomForest")
library("htmltools")

#### FUNCTIONS #####

#### a) DATA #### 

unshift_coord <- function(txt_file, freq = 1587){
  # 1157 cm is reference peak for liver
  wavelength_ref <- read.table("Z:/PROJETS RAMAN/Code/Code greffon hépatique/Raman_transplant/frequences.csv")[,1] 
  coord_x <- txt_file[,1]
  data <- txt_file[,2]
  ix1 <- findInterval(freq-10, coord_x, rightmost.closed = T)
  ix2 <- findInterval(freq+10, coord_x, rightmost.closed = T)
  tmp <- data[ix1:ix2]
  idx_freq <- which.max(tmp) + ix1 #indice du pic dans le txt
  idx_ref <- findInterval(freq, wavelength_ref, rightmost.closed = T) +1
  diff <- idx_freq - idx_ref
  return(diff:(1599+diff))
}

import_sp <- function(raw_path, freq=0, name="Raw spectra", sep="\t", dec=".") {
  wavelength_ref <- read.table("Z:/PROJETS RAMAN/Code/Code greffon hépatique/Raman_transplant/frequences.csv")[,1]
  if (endsWith(raw_path, '.txt')){
    #print(raw_path)
    txt_file <- read.table(raw_path, sep=sep, header=F, colClasses="numeric", dec=dec) 
    coord_x <- txt_file[,1]
    idx_min <- findInterval(301, coord_x, rightmost.closed = T)+1
    idx_max <- idx_min + 1599
    raman_data <- txt_file[idx_min:idx_max,2]
    rownames(raman_data) <- raw_path
  }
  else{
    # Import spectra of txt format stored in a given directory
    raw_files <- data.frame(filename = list.files(raw_path, 'txt'))
    df_filenames <- paste0(raw_path, raw_files$filename)
    raman_data <- c()
    
    for (i in 1:nrow(raw_files)){
      #print(df_filenames[i])
      txt_file <- read.table(df_filenames[i], sep=sep, header=F, colClasses="numeric", dec=dec)
      coord_x <- txt_file[,1]
      idx_min <- findInterval(301, coord_x, rightmost.closed = T)+1
      idx_max <- idx_min + 1599
      if (freq == 0){
        new_idx <- idx_min:idx_max
      }
      else{
        new_idx <- unshift_coord(txt_file, freq)
      }
      tmp_data <- txt_file[new_idx,2]
      #print(paste(idx_min, idx_max))
      raman_data <- rbind(raman_data, tmp_data)
    }
    rownames(raman_data) <- raw_files$filename
  }
  #attr(raman_data, "name") <- name
  attr(raman_data, "xaxis.values") <- wavelength_ref
  attr(raman_data, "xaxis.name") <- "Wavelength"
  attr(raman_data, "yaxis.name") <- "Intensity"
  
  colnames(raman_data) <- wavelength_ref
  return(raman_data)
}

get_samples_id <- function(data){
  # Extract identifier from the filename
  # Hidden code for confidentiality
}

regroup_by_patient <- function(data){
  # Split rownames by '_' and creates a new dataframe which where one row corresponds to the mean spectrum of each patient
  list_of_samples <- get_samples_id(data)
  
  unique_names <- unique(list_of_samples)
  patients <- data.frame()
  
  for (n in unique_names){
    idx <- which(list_of_samples==n)
    new_name <- paste("Average_", n)
    # new_name <- n
    moy <- average_sp(data[idx,], new_name)
    patients <- rbind(patients, moy)
  }
  return(patients)
}
get_varY <- function(dataset, info, var){
  # return qualitative variable (Y) associated to each spectra
  tab <- as.data.frame(get_samples_id(dataset))
  colnames(tab) <- c("ID2")
  merge <- merge(tab, info)[var]
  return(merge[,1])
}

select_spectra <- function(data, info, var, select){
  ### retourne les index des spectres d'une liste de modalités
  return(which(info[,var] %in% select, info))
}

filter_noised_spectra <- function(data){
  # remove noised spectra 
  deriv <- savitsky(normalize(data[,1300:1450]),2)
  avg_deriv <- apply(abs(deriv), 1, mean)  
  sd_deriv <- apply(abs(deriv), 1, sd)
  idx_to_keep <- which(avg_deriv < 0.01 & sd_deriv > 0.00000055 & get_values_peaks(1500, data ) > -50)
  removed <- rownames(data)[-idx_to_keep]
  cat("These", length(removed), "spectra have been removed from the dataset because too noised:\n")
  cat(removed, sep ="\n")
  return(idx_to_keep)
}

average_sp <- function(spectra, new_name=NULL) {
  # Compute the average spectrum of all spectra
  w <- colnames(spectra)
  avg <- t(as.data.frame(apply(spectra, 2, mean))) 
  colnames(avg) <- w
  attr(avg, "xaxis.values") <- w
  if (is.null(new_name)){
    new_name <- paste("Average spectra of", attr(spectra, "name")) 
  }
  attr(avg, "name") <- new_name
  rownames(avg) <- c(new_name)
  return(avg)
}

compute_average_factor <- function(data, groupes, select=1:nrow(data)){
  tab <- data.frame()
  data2 <- data[select,]
  groupes2 <- groupes
  for (u in unique(groupes2)){
    tmp <- data2[groupes2 == u,]
    tab <- rbind(tab, average_sp(tmp, paste("Average spectra of",u)))
  }  
  return(tab)
}

remove_silent_region <- function(data){
  coord_x <- as.double(as.character(colnames(data)))
  ix_700 <- findInterval(700, coord_x)
  ix_1700 <- findInterval(1700, coord_x)
  ix_2700 <- findInterval(2700, coord_x)
  ix_3100 <- findInterval(3100, coord_x)
  keep_col <- c(ix_700:ix_1700, ix_2700:ix_3100)
  return(as.data.frame(data[,keep_col]))
}


#### b) PEAKS #### 

ratio_2_pics <- function(data, c1, c2, pics=NULL, name=NULL){
  if (is.numeric(c1)){
    v1 <- c1
    v2 <- c2
  }
  else{
    v1 <- pics[pics$Correspondance == c1,'Pic']
    v2 <- pics[pics$Correspondance == c2,'Pic']    
  }
  coord_x <- as.double(as.character(colnames(data)))
  v1_ix <- findInterval(v1, coord_x)
  v2_ix <- findInterval(v2, coord_x)
  res <- as.data.frame(data[,v1_ix]/data[,v2_ix])
  if (is.null(name)){
    colnames(res) <- c(paste("Ratio", c1, 'on',c2))  
  }
  else{
    colnames(res) <- c(name)  
  }
  rownames(res) <- rownames(data)
  return(res)
}

ratio_m_pics <- function(data, val, peaks){
  table <- data.frame(row.names=rownames(data))
  for (p in peaks){
    table <- cbind(table, ratio_2_pics(data, val, p))
  }
  return(table)
}

mean_ratio_pics <- function(tab, new_name=NULL){
  # calcule la valeur moyenne (par colonne) d'un tableau de pics
  res <- t(as.data.frame(apply(tab, 2, mean)))
  rownames(res) <- c("Average ratio")
  return(res)
}

find_peak_value_around <- function(freq, data, max=F, tol=10){
  # return index of the peak for each spectra of the dataframe
  coord_x <- as.double(as.character(colnames(data)))
  if (max){
    ix1 <- findInterval(freq-tol, coord_x, rightmost.closed = T)
    ix2 <- findInterval(freq+tol, coord_x, rightmost.closed = T)
    tmp <- as.data.frame(data[,ix1:ix2])
    df <- as.data.frame(apply(tmp, 1, which.max))   
  }
  else{
    ix <- findInterval(freq, coord_x)
    df <- data[,ix]
    print(ix)
  }
  #rownames(df) <- rownames(data)
  return(df)
} 

get_intensities_peak <- function(idx_col, data){
  val_peaks <- c()
  for(i in 1:nrow(data)) {
    col <- idx_col[i,]
    val_peaks <- c(val_peaks, data[i,col])
  }
  return(val_peaks)
}

get_values_peaks <- function(peaks, data){
  tab <- data.frame(matrix(nrow=nrow(data), ncol=0))
  rownames(tab) <- rownames(data)
  for (p in peaks){
    val <- find_peak_value_around(p, data)
    tab <- cbind(tab, val)    
  }
  colnames(tab) <- peaks
  return(tab)
}


#### c) PRE-PROCESSING #### 
baseline_corr <- function(spectra) {
  # Correct baseline for all spectra
  w <- colnames(spectra)
  sp_name <- attr(spectra, "name")
  corr <- prep.alsbasecorr(as.matrix(spectra))
  colnames(corr) <- w
  attr(corr, "xaxis.values") <- w
  attr(corr, "name") <- paste("Baseline corrected spectra", "for", sp_name)
  rownames(corr) <- c(rownames(spectra))
  #corr[corr < 0] <- 0
  return(corr)
}

savitsky <- function(spectra, deg = 0, p=2) {
  # Compute the savitsky golay derivative at a specific degree
  w <- colnames(spectra)
  sp_name <- attr(spectra, "name")
  svg <- prep.savgol(as.matrix(spectra), width=15, porder=p, dorder=deg)
  colnames(svg) <- w
  attr(svg, "xaxis.values") <- w
  titre <- paste("Spectra at Savitsky-Golay's derivative", deg, "for", sp_name)
  attr(svg, "name") <- titre
  #rownames(svg) <- c(paste(rownames(spectra), "at SG's derivative", deg))
  return(svg)
}

normalize <- function(spectra, titre=NULL) {
  # Standard Normal Variate, makes all values from the same row to have zero mean and unit standard deviation 
  w <- colnames(spectra)
  sp_name <- attr(spectra, "name")
  norm <- prep.snv(as.matrix(spectra))
  colnames(norm) <- w
  attr(norm, "xaxis.values") <- w
  if (is.null(titre)){
    titre <- paste("Normalized", sp_name)
  }
  attr(norm, "name") <- titre
  #rownames(norm) <- c(paste(rownames(spectra), "normalized"))
  norm <- t(apply(norm, 1, function(s) s - min(s) + 1e-6))
  return(norm)
}

#### d) PLOTTING #### 

compare_groups_for_values <- function(vec, data, groups){
  # groups ne doivent pas être au format factor !!!
  coord_x <- as.double(as.character(colnames(data)))
  vec <- as.double(as.character(vec))
  for (freq in vec){
    ix <- findInterval(freq, coord_x, rightmost.closed = T)
    tab <- as.data.frame(cbind(groups,data[,ix]))
    tab$V2 <- as.numeric(tab$V2)
    titre <- paste("Comparaison of", freq, "intensities for groups:", paste(unique(groups), collapse = ", "))
    p <- ggplot(tab,aes(groups,V2))+geom_boxplot(aes(fill=groups))+  
      ggtitle(titre)
    plot(p)
  }
}

  
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

plot_spectra <- function(data, colors=gg_color_hue(8), titre=NULL){
  if (is.null(titre)){
    titre <- attr(data, "name")
  }
  coord_x <- as.double(as.character(colnames(data)))
  data_raman <- as.data.frame(cbind(coord_x, t(data)))
  data_long <- melt(data_raman, id = "coord_x")
  p <- plot_ly(data = data_long, x = ~coord_x, y = ~value, color= ~variable,
               name= ~variable, type = 'scatter', mode = 'lines',
               colors=colors) %>%
    layout(xaxis = list(title = "Raman Wavelength (cm -1)"),
           yaxis = list(title = "Intensity"),
           legend=list(title=list(text='<b> Spectres </b>')))
  
  if (! is.null(titre)){
    p <- layout(p, title = paste0('<b>', titre, '</b>'))
  }
  return(p)
}

