#' Filter features of data matrix (data) that are present in c (integer) consequtive fractions -1
#'
#' @param x Data frame containing elution profiles
#' @param c A required length of peak - number of fractions that each single peak has to span on
#' @return A data frame without elution profiles, which peaks span less than c fractions
#' @export
filterConsistentFeatures_ML <- function(x, c){
  conseqVec <- c() # create empty vector
  for (i in 1:dim(x)[1])
  {
    tmpNos <- x[i,] # copy intensities from row i of data x
    tmpCount <- which(tmpNos>0) # check which "cells" are higher than 0
    tmp<-c() # create empty vector
    for (p in 1:(length(tmpCount)-1))
    {
      tmp[p]<- tmpCount[p+1]-(tmpCount[p]) # check what is the distance between "non-zero" values
    }
    tmpLen <- length(which(tmp==1,arr.ind = TRUE)) # calculate number of consecutive fractions in which feature is present (e.g. two consecutive fractions will give you value = 1)
    if (tmpLen >=c) {
      conseqVec[i]<-TRUE
    }else conseqVec[i]<-FALSE # check if tmpLen is higher or equal than set "c" value
  }
  x <- x[conseqVec,]
}

#' Filter features of data matrix (data), which have a single appearance (e.g. vector of 1110101111 will be replaced by 1110001111)
#'
#' @param z Data frame containing elution profiles
#' @return A data frame containing elution profiles, which does not contain a single molecule appearances, which might be often caused by technical noise
#' @export
filterPeaksSpanFractions_ML <- function(z){
  dataConsPeaks <- as.data.frame(NULL)
  w <- 1
  while (w <= nrow(z)){
    tmpNos <- z[w,] # copy intensities from row w of data x
    tmpCount <- which(tmpNos>0) # check which "cells" are higher than 0
    tmp<-c() # create empty vector
    tmpr<-c() # create empty vector

    for (p in 1:length(tmpCount)){
      tmp[p]<- tmpCount[p+1]-(tmpCount[p]) # check what is the distance between "non-zero" values - 1st way
    }

    for (q in 1:length(tmpCount)){
      tmpr[q]<- tmpCount[q]-(tmpCount[q-1]) # check what is the distance between "non-zero" values - 2nd way
    }

    tmpb <- NULL
    for(i in 1:length(tmp)){
      tmpb[i] <- tmp[i] == 1 | tmpr[i] == 1 # check if at least one method gives 1 (positive value)
    }

    tmpb[is.na(tmpb)] <- FALSE # remove NA values from tmpb and replace by FALSE

    tmpCountReal <- tmpCount[tmpb] # check which cells are higher than 0 AND the distance between non zero values is 1
    #tmpNosReal <- tmpNos[,tmpCountReal] # probably no needed

    tmpm <- as.data.frame(matrix(0, nrow = nrow(tmpNos), ncol = ncol(tmpNos)))
    x <- 1
    y <- 1

    while(x <= length(tmpNos)){
      while(y <= length(tmpCountReal)){
        if(x == tmpCountReal[y]){
          tmpm[,x] <- tmpNos[x]
          y <- y+1
        }else{
          y <- y+1
        }
      }
      y <- 1
      x <- x+1
    }

    colnames(tmpm) <- colnames(tmpNos)
    rownames(tmpm) <- rownames(tmpNos)

    dataConsPeaks <- rbind(dataConsPeaks, tmpm)
    w <- w+1
  }
  dataConsPeaks # give the result of the function to the "clipboard"
}

#' Function combines filterConsistentFeatures_ML and filterPeaksSpanFractions_ML
#' Remove single appearances and short peaks
#' @param x Data frame containing elution profiles
#' @param nr_replicas Number of replicas used in the experiment
#' @param nr_fractions Number of fractions collected
#' @return A data frame containing elution profiles, which does not contain a single molecule appearances, which might be often caused by technical noise
#' @export
rmv_short_peaks <- function(x, nr_replicas, nr_fractions){
  nr_comparisons <- ((nr_replicas)*(nr_replicas-1))/2
  rownames(x) <- x$Name
  list_2frc <- as.data.frame(matrix(0, ncol = 0, nrow = nrow(x)))
  rnames <- x$Name
  list_2frc <- cbind(list_2frc, rnames)
  x2 <- x[,2:ncol(x)]

  i <- 1
  for(i in 1:nr_replicas){
    tmp <- paste("Rep_", i, sep = "")

    data_temp <- x2[,grep(tmp, colnames(x2))]

    names <- paste(tmp, "_Fraction_", c(1:nr_fractions), sep = "")

    colnames(data_temp) <- names

    temp_2frc <- as.data.frame(NULL)
    temp_2frc <- filterConsistentFeatures_ML(data_temp, 1) # filter for metabolites which peaks span on at least 2 fractions

    temp_long_peaks <- filterPeaksSpanFractions_ML(temp_2frc) # remove single appearance of metabolite in fractions 1111010111 then 1 in the middle will be remove since it doesn't create peak spanning on 2 fractions

    rnames <- rownames(temp_long_peaks)
    temp_long_peaks <- cbind(temp_long_peaks, rnames)

    list_2frc <- full_join(list_2frc, temp_long_peaks, by = "rnames")
  }
  colnames(list_2frc)[1] <- "Name"
  rownames(list_2frc) <- list_2frc$Name

  list_2frc <- list_2frc[,-1]

  tmp <- NULL
  list_2frc[list_2frc == 0] <- NA

  for(i in 1:nrow(list_2frc)){
    if(!all(is.na(list_2frc[i,]))){
      tmp <- rbind(tmp, list_2frc[i,])
    }
  }

  tmp[is.na(tmp)] <- 0

  tmp <- rownames_to_column(tmp, "Name")

  print(paste("Number of metabolites, which elution profiles spans at least 2 fractions: ", nrow(tmp), sep = ""))
  tmp
}

#' Filter for elution profiles of metabolites, which are reproducible
#' @param x Data frame containing elution profiles
#' @param req_repl A number that specify required number of replicas in which elution profiles were reproducible
#' @param repr_profile_threshold Minimal level of similarity between the elution profiles (in range from 0 to 1)
#' @return A data frame containing reproducible elution profiles.
#' @export
repr_profile <- function(x, req_repl, repr_profile_threshold){
  if(req_repl > 1){
    molecule_presence <- NULL
    i <- 1
    while(i <= nr_replicas){
      Repi_presence <- apply(x[,grep(paste0("Rep_", i), colnames(x))], 1, max) > 10000
      molecule_presence <- data.frame(cbind(molecule_presence, Repi_presence))
      i <- i+1
    }

    PresentInReplicas <- apply(molecule_presence, 1, sum)

    molecule_presence <- as.data.frame(cbind(molecule_presence, PresentInReplicas))

    x$PresentInReplicas <- molecule_presence$PresentInReplicas

    met_reproducible <- dplyr::filter(x, PresentInReplicas >= req_repl)

    i <- 1
    PCC_tmp <- NULL
    list_of_replicas <- c(1:nr_replicas)
    nr_comparisons <- ((nr_replicas)*(nr_replicas-1))/2
    PCC_All <- NULL

    while(i <= nrow(met_reproducible)){
      w <- 1
      tmp <- data.frame(matrix(NA, ncol = nr_comparisons, nrow = nr_fractions))
      PCC_tmp <- NULL
      PCC_collect <- NULL
      while(w <= nr_replicas){
        tmp[,w] <- t(met_reproducible[i,grep(paste0("Rep_", w), colnames(met_reproducible))])
        w <- w+1
      }

      a <- 1
      b <- a+1
      while(a <= (nr_replicas-1)){
        while(b <= nr_replicas){
          PCC_tmp <- cor(tmp[,a], tmp[,b])
          var_col_name <- paste0("Rep_", a, "_vs_", "Rep_", b)
          colnames(PCC_tmp) <- var_col_name
          PCC_collect <- cbind(PCC_collect, PCC_tmp)
          b <- b+1
        }
        a <- a+1
        b <- a+1
      }
      PCC_All <- rbind(PCC_All, PCC_collect)
      i <- i+1
    }


    PCC_All <- data.frame(PCC_All)
    list_of_comparisons <- as.list(colnames(PCC_All))
    Max_PCC <- apply(PCC_All, 1, max)
    PCC_All$Max_PCC <- Max_PCC

    met_reproducible <- cbind(met_reproducible, PCC_All)
    met_reproducible <- dplyr::filter(met_reproducible, Max_PCC >= repr_profile_threshold)
    tmp <- nrow(met_reproducible)

    print(paste("Number of metabolites having reproducible pattern: ", tmp, sep = ""))
    met_reproducible
  }else{
    return(x)
  }
}

#' Calculate single profile for specified molecule
#' @param x Data frame containing elution profiles
#' @param single_profile A function that will be used to determine the single elution profile ("sum" or "median")
#' @param nr_replicas A number of replicas used in the experiment
#' @return A data frame containing single elution profiles calculated based on defined number of replicates.
#' @export
calc_single_profile <- function(x, single_profile, nr_replicas){
  if(nr_replicas > 1){
    if(single_profile == "sum"){
      nr_comparisons <- ((nr_replicas)*(nr_replicas-1))/2
      get_col_names <- colnames(x)[(ncol(x)-nr_comparisons):(ncol(x)-1)]
      names <- paste("Fraction_", c(1:nr_fractions), sep = "")
      #x2 <- x[,-c((ncol(x)-nr_comparisons):ncol(x))]
      x2 <- x[,2:(ncol(x)-nr_comparisons-3)]
      rownames(x2) <- x$Name
      tmp_profile <- NULL
      w <- 1

      while(w <= nrow(x2)){
        tmp_all <- NULL
        i <- 1
        while(i <= length(get_col_names)){
          if(x[w,grep(get_col_names[i], colnames(x))] >= repr_profile_threshold){
            list_of_replicas <- as.character(unlist(strsplit(get_col_names[i], "_vs_")))
            tmp1 <- x2[w,grep(list_of_replicas[1], colnames(x2))]
            tmp2 <- x2[w,grep(list_of_replicas[2], colnames(x2))]
            colnames(tmp1) <- names
            colnames(tmp2) <- names

            tmp_all <- rbind(tmp_all, tmp1, tmp2)
            tmp_all <- unique(tmp_all)

            #(x[w,grep(list_of_replicas[1], colnames(x))], x[w,grep(list_of_replicas[2], colnames(x))])
            #tmp_all <- rbind()
          }
          #tmp_PCC <- x[w,grep(get_col_names[i], colnames(x))]
          i <- i+1
        }

        single_met_profile <- colSums(tmp_all, na.rm = TRUE)
        tmp_profile <- data.frame(rbind(tmp_profile, single_met_profile))
        w <- w+1
      }
      tmp_profile$Name <- x$Name
      rownames(tmp_profile) <- tmp_profile$Name
    }
    if(single_profile == "median"){
      nr_comparisons <- ((nr_replicas)*(nr_replicas-1))/2
      get_col_names <- colnames(x)[(ncol(x)-nr_comparisons):(ncol(x)-1)]
      names <- paste("Fraction_", c(1:nr_fractions), sep = "")
      x2 <- x[,2:(ncol(x)-nr_comparisons-2)]
      rownames(x2) <- x$Name
      tmp_profile <- NULL
      w <- 1

      while(w <= nrow(x2)){
        tmp_all <- NULL
        i <- 1
        while(i <= length(get_col_names)){
          if(x[w,grep(get_col_names[i], colnames(x))] >= repr_profile_threshold){
            list_of_replicas <- as.character(unlist(strsplit(get_col_names[i], "_vs_")))
            tmp1 <- x2[w,grep(list_of_replicas[1], colnames(x2))]
            tmp2 <- x2[w,grep(list_of_replicas[2], colnames(x2))]
            colnames(tmp1) <- names
            colnames(tmp2) <- names

            tmp_all <- rbind(tmp_all, tmp1, tmp2)
            tmp_all <- unique(tmp_all)

            #(x[w,grep(list_of_replicas[1], colnames(x))], x[w,grep(list_of_replicas[2], colnames(x))])
            #tmp_all <- rbind()
          }
          #tmp_PCC <- x[w,grep(get_col_names[i], colnames(x))]
          i <- i+1
        }
        tmp_all[tmp_all==0] <- NA
        single_met_profile <- apply(tmp_all, 2,  median, na.rm = TRUE)
        single_met_profile[is.na(single_met_profile)] <- 0

        tmp_profile <- data.frame(rbind(tmp_profile, single_met_profile))
        w <- w+1
      }
      tmp_profile$Name <- x$Name
      rownames(tmp_profile) <- tmp_profile$Name
    }
    return(tmp_profile)
  }else{
    tmp_profile <- x
    return(tmp_profile)
  }
}

#' Normalize all features to its max
#' @param x Data frame containing elution profiles
#' @return A data frame containing single elution normalized to maximum intensity in each row
#' @export
maxNormalize <- function(x){
  dataNew <- sweep(x, 1, apply(x, 1, max, na.rm = TRUE), FUN = "/")
  dataNew <- as.data.frame(dataNew)
}

#' Normalize all samples (columns) to its median factor (median/median of all)
#' @param x Data frame containing intensities
#' @return A data frame containing median normalized values
#' @export
medianNormalize <- function(x){
  # dataNew <- x
  x[x == 0] <- NA
  medFactor <- apply(x, 2, median, na.rm = TRUE)/median(apply(x, 2, median, na.rm = TRUE), na.rm = TRUE) ### 1) calculate median for each column; 2) calculate median for whole data set; 3) Calculate ratio of 1) and 2);
  dataNew <- sweep(x, 2, apply(x, 2, median, na.rm = TRUE)/median(as.matrix(x), na.rm = TRUE), "/") ### Take data set and DIVIDE each value in particular column by medFactor calculated above
  dataNew <- as.data.frame(dataNew) ### Save dataset and a data.frame
}

#' Deconvolute elution profiles
#' @param x Data frame containing single elution profiles
#' @param var_min_peak Determine minimum relative signal intensity accepted as peak e.g. 0.2 of max peak
#' @param var_limit Determine how broad will be base of the peak. What is the minimal intensity at which peak will be recorded
#' @param var_limit_down Set minimum drop of the intensity between the peaks
#' @param var_limit_up Set minimum rise of the intensity between the peaks
#' @param var_num_col Set number of fractions that were collected from SEC
#' @param var_size_range_name Define path to the file containing corresponding theoretical molecular weight of each fraction
#' @param var_exp_name Name of your output files
#' @param var_work_table Define whether your data will be normalized 0 - No normalization, 1 - log2 transformation, 2 - log10 transformation, 3 - normalize to maximum
#' @param var_plot_decon_peaks Define whether how plotting should be performed. 1 - without normalization, 2 - with normalization
#' @param dec_data_color Define coloro of plotted deconvoluted data
#' @param plot_singe_plot Define if single profiles of deconvoluted data should be created
#' @param tmp_dir Define temporary directory where data will be stored
#' @return A data frame containing single elution normalized to maximum intensity in each row
#' @export
deconvolution <- function(x, var_min_peak, var_limit, var_limit_down, var_limit_up, var_num_col,
                          var_size_range_name, var_exp_name, var_work_table, var_plot_decon_peaks,
                          dec_data_color, plot_single_plot, tmp_dir, fraction_names){
  ## Load output from MQ
  #tbl_data <- as.data.frame(read.delim(var_df_name, sep=var_sep_csv, header= TRUE, stringsAsFactors = FALSE, row.names = 1))
  var_date <- Sys.Date()  ###gets current date


  tbl_data <- x
  rownames(tbl_data) <- paste0(rownames(tbl_data), "_")
  fr_size_range <- read.delim(var_size_range_name)


  ## Create a matrix from data frame with row.names
  mat_data <- as.matrix(tbl_data)
  rownames(mat_data) <- row.names(tbl_data)
  mat_data <- mat_data[,1:var_num_col]

  ## Assign size ranges to column names and order them from the smallest to the largest
  #colnames(mat_data) <- round(fr_size_range[1:var_num_col,"Corresponding mass"],digits=1)
  colnames(mat_data) <- fr_size_range$Corresponding.mass
  #tbl_size_fract <- fr_size_range[,c("Fraction number","Corresponding mass")]
  tbl_size_fract <- fr_size_range$Corresponding.mass...value

  ## Assign mat_data2 file and remove rows which sum to 0.
  mat_data2 <- mat_data
  mat_data2[is.na(mat_data2)] <- 0
  mat_data2 <- mat_data2[rowSums(mat_data2)!=0, ]

  ###create peak matrix --> find peaks, remove peaks with low signal and peaks with too low MW for protein complex
  mat_data2 <- cbind(mat_data2,rep(0,length(mat_data2[,1]))) #!!!!!!!add one column of zeros at the end for peak picking
  mat_data2 <- cbind(rep(0,length(mat_data2[,1])),mat_data2) #!!!!!!!add one column of zeros at the FRONT for peak picking
  tbl_size_fract <- rbind(c(0.1,100),tbl_size_fract)#!!!!!!!add one row to size_fract table
  colnames(mat_data2)[ncol(mat_data2)] <- 0.1
  colnames(mat_data2)[1] <- 10000
  mat_peaks <- numeric() ##empty matrix for peak collection
  mat_peaks_factor <- numeric() ##empty matrix for peak factor collection
  #var_which_column <- sum(as.numeric(colnames(mat_data2))<var_size_cutoff)
  for(i in 1 :length(mat_data2[,1])){
    var_agi <- rownames(mat_data2)[i] ##get current AGI
    #if(var_agi=="AT4G28440"){stop()}
    #var_size_protein <- tbl_size_TAIR10[tbl_size_TAIR10[,1]==var_agi,3]/1000 ##get size of current protein
    var_max <- max(mat_data2[i,1:ncol(mat_data2)])##get maximum excluding signal in small size fractions
    #if(var_which_column==0 || var_size_cutoff==0 || var_size_protein<var_size_cutoff){var_max <- max(mat_data2[i,])} ##get maximum
    var_min_peak_accepted <- var_max * var_min_peak
    vec_peaks_factor <- rep(0,length(mat_data2[1,])) ##create vector with zeros --> no peaks detected yet
    vec_peaks <- extract(turnpoints(mat_data2[i,]), length(mat_data2[1,]), peak=1, pit=0)
    for(j in 1:(length(vec_peaks)-1)){ #delete all peaks which are too low in signal AND monomeric peaks
      #var_size_current_fract <- tbl_size_fract[j,2]
      if(mat_data2[i,j] < var_min_peak_accepted){vec_peaks[j]<- 0} ##if signal to low set vec_peaks position to 0 = delete "peak"
      #if(var_size_current_fract < var_size_protein * var_peak_factor){vec_peaks[j]<- 0} ##if current peak position too small for protein complex --> delete peak
      #vec_peaks_factor[j] <- var_size_current_fract/var_size_protein
    }
    mat_peaks <- rbind(mat_peaks, vec_peaks) ##collect peak data
    #mat_peaks_factor <- rbind(mat_peaks_factor, vec_peaks_factor) ##collect peak data
  }
  rownames(mat_peaks) <- rownames(mat_data2)
  mat_peaks <- mat_peaks[rowSums(mat_peaks)!=0, ]

  mat_peaks_factor <- matrix(rep(3), nrow=nrow(mat_peaks),ncol=ncol((mat_peaks)))
  rownames(mat_peaks_factor) <- rownames(mat_data2)

  ####deconvolution of peaks
  mat_data2_decon <- numeric()
  vec_decon_peakfactor <- numeric()
  vec_decon_peakwidth <- numeric()
  for(i in 1:length(mat_data2[,1])){
    var_record <- 0
    var_count <- 1 ##variable for counting peaks for each gene
    var_agi <- rownames(mat_data2)[i] ##get name of protein
    #if(var_agi=="AT1G17100"){stop()}
    vec_data <- mat_data2[i,] ##get abundance data for protein
    vec_data_peakfactor <- mat_peaks_factor[i,]
    vec_data_peaks <- mat_peaks[i,] ##get peak data for protein
    vec_data_record <- rep(0,length(vec_data)) ##create empty vector for abundance data
    vec_data_record_peaks <- rep(0,length(vec_data))##create empty vector for peak data
    for(j in 1:length(vec_data_record)){
      var_max_current_peak <- max(vec_data_record) ###gets the current maximum in the recorded profile (can also be zero)
      if(var_record == 0 & vec_data[j]/max(vec_data)>var_limit){ ###if abundance value > minimum peak size
        var_record <-1
        if(j>1){vec_data_record[(j-1)]<-vec_data[(j-1)]} ##record previous abundance value (the last one below threshhold)
      }
      if(var_record == 1 & vec_data[j]/max(vec_data)>var_limit){ ###keep recording abundance and peak data
        vec_data_record[j]<-vec_data[j]
        vec_data_record_peaks[j]<-vec_data_peaks[j]
        if(vec_data[j]<var_max_current_peak*var_limit_down){var_record <- 2} #if abundance value is falling more than x percent --> set var_record to 2 --> entering next "phase" of scanning
      }
      if(var_record == 1 & vec_data[j]/max(vec_data)<=var_limit){#### ...until abundance falls below threshold
        vec_data_record[j]<-vec_data[j] ##take last abundance value (which is already below threshold)
        vec_data_record_peaks[j]<-vec_data_peaks[j]
        if(sum(vec_data_record_peaks)>0){  #only transfer data to final table if valid peak part of recorded abundance
          mat_data2_decon <- rbind(mat_data2_decon,vec_data_record)  ###transfer data to new table
          rownames(mat_data2_decon)[length(mat_data2_decon[,1])]<- paste(var_agi,var_count,sep="_") ###give tranfered data name of protein + peak index
          vec_decon_peakfactor <- c(vec_decon_peakfactor,max(vec_data_peakfactor[vec_data_record_peaks>0]))
          names(vec_decon_peakfactor)[length(vec_decon_peakfactor)]<- paste(var_agi,var_count,sep="_") ###give tranfered data name of protein + peak index
          var_count <- var_count + 1 #increase peak count
        }
        vec_data_record <- rep(0,length(vec_data)) ##create new empty vector
        vec_data_record_peaks <- rep(0,length(vec_data))
        var_record <-0

      }
      if(var_record == 2){
        if(vec_data[j]/max(vec_data)>var_limit & vec_data[j]*var_limit_up<=vec_data[(j-1)]){ ###keep recording abundance and peak data if abundance is still falling
          vec_data_record[j]<-vec_data[j]
          vec_data_record_peaks[j]<-vec_data_peaks[j]
          if(vec_data[j]<var_max_current_peak*var_limit_down){var_record <- 2} #if abundance value is falling more than x percent --> set var_record to 2 --> entering next "phase" of scanning
        }
        if(vec_data[j]/max(vec_data)<=var_limit){#### ...until abundance falls below threshold
          vec_data_record[j]<-vec_data[j] ##take last abundance value (which is already below threshold)
          vec_data_record_peaks[j]<-vec_data_peaks[j]
          if(sum(vec_data_record_peaks)>0){  #only transfer data to final table if valid peak part of recorded abundance
            mat_data2_decon <- rbind(mat_data2_decon,vec_data_record)  ###transfer data to new table
            rownames(mat_data2_decon)[length(mat_data2_decon[,1])]<- paste(var_agi,var_count,sep="_") ###give tranfered data name of protein + peak index
            vec_decon_peakfactor <- c(vec_decon_peakfactor,max(vec_data_peakfactor[vec_data_record_peaks>0]))
            names(vec_decon_peakfactor)[length(vec_decon_peakfactor)]<- paste(var_agi,var_count,sep="_") ###give tranfered data name of protein + peak index
            var_count <- var_count + 1 #increase peak count
          }
          vec_data_record <- rep(0,length(vec_data)) ##create new empty vector
          vec_data_record_peaks <- rep(0,length(vec_data))
          var_record <-0
        }
        if(vec_data[j]*var_limit_up>vec_data[(j-1)]){####stop recording if abundance value increases
          ###do not record current abundance HERE as already par tof new peak
          if(sum(vec_data_record_peaks)>0){  #only transfer data to final table if valid peak part of recorded abundance
            mat_data2_decon <- rbind(mat_data2_decon,vec_data_record)  ###transfer data to new table
            rownames(mat_data2_decon)[length(mat_data2_decon[,1])]<- paste(var_agi,var_count,sep="_") ###give tranfered data name of protein + peak index
            vec_decon_peakfactor <- c(vec_decon_peakfactor,max(vec_data_peakfactor[vec_data_record_peaks>0]))
            names(vec_decon_peakfactor)[length(vec_decon_peakfactor)]<- paste(var_agi,var_count,sep="_") ###give tranfered data name of protein + peak index
            var_count <- var_count + 1 #increase peak count
          }
          vec_data_record <- rep(0,length(vec_data)) ##create new empty vector
          vec_data_record_peaks <- rep(0,length(vec_data))
          var_record <-1 ##start recording immediately
          vec_data_record[j]<-vec_data[j] ##take last abundance value (which is already below threshold)
          vec_data_record_peaks[j]<-vec_data_peaks[j]
          vec_data_record[(j-1)]<-vec_data[(j-1)] ##also tranfer previous point
        }
      }

    }
  }

  ###remove additional column of zeros from data tables
  mat_data2 <- mat_data2[,-(length(mat_data2[1,]))] ###remove additional column of zeros
  mat_peaks <- mat_peaks[,-(length(mat_peaks[1,]))]  ###remove additional column of zeros
  mat_data2_decon <- mat_data2_decon[,-(length(mat_data2_decon[1,]))] ###remove additional column of zeros
  mat_data2_decon <- mat_data2_decon[,-1]
  mat_peaks <- mat_peaks[,-1]
  mat_data2 <- mat_data2[,-1]



  colnames(mat_data2_decon) <- colnames(mat_data)
  write.table(mat_data2_decon, file = paste(wd, var_exp_name, ".txt", sep = ""), sep="\t")


  ### Normalization ###

  ##########################
  ##create data table to work with - normalise all deconvoluted peaks (normalise according to selection)
  if(var_work_table == 0){mat_data2_decon_normalized <- mat_data2_decon} ##no normalisation
  if(var_work_table == 1){
    mat_data2_decon_normalized <- mat_data2_decon ##copy data to new table
    for(i in 1:length(mat_data2_decon_normalized[,1])){ ##loop through rows
      for(j in 1:length(mat_data2_decon_normalized[1,])){ ###loop through columns
        if(mat_data2_decon_normalized[i,j]>0){mat_data2_decon_normalized[i,j] <- log10(mat_data2_decon_normalized[i,j])} ##create log value if abundance > 0 (otherwise there would be an error)
      }
    }
  }
  if(var_work_table == 2){ ##log2 normalisation
    mat_data2_decon_normalized <- mat_data2_decon ##copy data to new table
    for(i in 1:length(mat_data2_decon_normalized[,1])){ ##loop through rows
      for(j in 1:length(mat_data2_decon_normalized[1,])){ ###loop through columns
        if(mat_data2_decon_normalized[i,j]>0){mat_data2_decon_normalized[i,j] <- log2(mat_data2_decon_normalized[i,j])} ##create log value if abundance > 0 (otherwise there would be an error)
      }
    }
  }
  if(var_work_table == 3){ #normalisation to 1 (recommended)
    mat_data2_decon_normalized <- mat_data2_decon
    for(i in 1:length(mat_data2_decon[,1])){
      mat_data2_decon_normalized[i,] <- mat_data2_decon_normalized[i,]/max(mat_data2_decon_normalized[i,])
    }
  }


  ########
  ####### Normalization of mat_data2
  ########

  mat_data2_norm <- mat_data2
  for(i in 1:length(mat_data2[,1])){
    mat_data2_norm[i,] <- mat_data2_norm[i,]/max(mat_data2_norm[i,])
  }

  #########################################################################################################
  #########################################################################################################
  ################################SINGLE MOLECULE PLOTTING ################################################
  #########################################################################################################
  #########################################################################################################
  if(plot_single_plot == TRUE){
    print(paste0("Deconvoluted profiles are now stored in ", tmp_dir))
    i <- 1
    while(i <= nrow(mat_data2_norm)){
      var_met <- rownames(mat_data2_norm)[i]

      var_orig_profile <- matrix(mat_data2_norm[grep(paste0("^", var_met), rownames(mat_data2_norm), fixed = FALSE),], ncol = nr_fractions)
      var_decon_profile <- matrix(mat_data2_decon_normalized[grep(paste0("^", var_met), rownames(mat_data2_decon_normalized), fixed = FALSE),], ncol = nr_fractions)

      plot_name <- paste(tmp_dir, "/", rownames(x)[i], ".jpeg", sep = "")
      jpeg(plot_name, width = 8, height = 8, units = 'cm', res = 600)
      par(mfrow = c(1,1), xpd = TRUE, mar=c(4,2.5,3,1.5))
      plot(1:nr_fractions, var_orig_profile, col=dec_data_color, type = "p", ylim = c(0,1), lwd = 1, pch = 21, axes = FALSE, ann=FALSE, lty = 1, cex = 0.3)
      w <- 1
      while(w <= nrow(var_decon_profile)){
        lines(1:nr_fractions, var_decon_profile[w,], col=w+1, type = "l", lwd = 1, lty = 1)
        w <- w+1
      }
      axis(1, at=1:nr_fractions, lab=fraction_names, las = 2, cex.axis = 0.3)
      axis(2, las=1, cex.axis = 0.4)
      box()
      title(main = list(rownames(x)[i], cex = 0.5))
      title(ylab = list("Relative intensity (Max)", cex = 0.5), line = 1.5)
      title(xlab = list("Fraction number", cex = 0.5))
      dev.off()
      i <- i+1
    }

  }

  #########################################################################################################
  #########################################################################################################
  ########################################## PLOTTING #####################################################
  #########################################################################################################
  #########################################################################################################

  ######################################################
  if(var_plot_decon_peaks == 1){

    ##plot all protein profiles incl. decon

    colnames(mat_data2) <- 1:ncol(mat_data2)
    colnames(mat_data2_decon) <- 1:ncol(mat_data2)
    vec_unique_AGIs <- unique(rownames(mat_data2))



    #create PDF to plot all profiles
    pdf(file = sprintf(paste(tmp_dir, "/", var_exp_name, ".pdf", sep=""), var_exp_name), paper='A4r')
    #set x graphs per page per page
    par(mfrow=c(3,1)) ###six per page
    for(t in 1:length(vec_unique_AGIs)){
      var_current_agi <- vec_unique_AGIs[t]
      vec_plot_original <- mat_data2[var_current_agi,] ##get original profile from data
      mat_plot <- mat_data2_decon[substr(rownames(mat_data2_decon),1,nchar(var_current_agi))==var_current_agi,,drop=FALSE] ###get matrix of deconvoluted profiles
      plot(names(vec_plot_original),vec_plot_original)
      #axis(1,at=1:ncol(mat_plot))
      lines(names(vec_plot_original),vec_plot_original)
      title(var_current_agi)
      if(length(mat_plot)>0){
        for(c in 1:length(mat_plot[,1])){
          lines(names(vec_plot_original),mat_plot[c,],col=(c+1), lwd=2)
        }
      }
    }
    dev.off()

  }


  ######################################################
  #setwd("H:/Software, devices, machines/R/PROMIS package/Results/Deconvoluted profiles")

  if(var_plot_decon_peaks == 2){

    ##plot all protein profiles incl. decon

    colnames(mat_data2_norm) <- 1:ncol(mat_data2_norm)
    colnames(mat_data2_decon_normalized) <- 1:ncol(mat_data2_norm)
    vec_unique_AGIs <- unique(rownames(mat_data2_norm))



    #create PDF to plot all profiles
    pdf(file = sprintf(paste(tmp_dir, "/", var_exp_name, "_normalized", ".pdf", sep=""), var_exp_name), paper='A4r')
    #set x graphs per page per page
    par(mfrow=c(3,1)) ###six per page
    for(t in 1:length(vec_unique_AGIs)){
      var_current_agi <- vec_unique_AGIs[t]
      vec_plot_original <- mat_data2_norm[var_current_agi,] ##get original profile from data
      mat_plot <- mat_data2_decon_normalized[substr(rownames(mat_data2_decon_normalized),1,nchar(var_current_agi))==var_current_agi,,drop=FALSE] ###get matrix of deconvoluted profiles
      plot(names(vec_plot_original),vec_plot_original)
      #axis(1,at=1:ncol(mat_plot))
      lines(names(vec_plot_original),vec_plot_original)
      title(var_current_agi)
      if(length(mat_plot)>0){
        for(c in 1:length(mat_plot[,1])){
          lines(names(vec_plot_original),mat_plot[c,],col=(c+1), lwd=2)
        }
      }
    }
    dev.off()

  }
  mat_data2_decon_normalized <- data.frame(mat_data2_decon_normalized, check.names = FALSE)
  return(mat_data2_decon_normalized)
}

#' Perform adduct detection
#' @param x Data frame containing basic information about the clusters (Cluster name, m/z, RT, Charge etc)
#' @param y Data frame containing list of allowed adducts and defined m/z differences
#' @param mass_deviation Maximum allowed mass_deviation in Da
#' @param RT_deviation Maximum allowed RT deviation
#' @return A data frame containing theoretical adduct assignment
#' @export
adduct_detection <- function(x, y, mass_deviation, RT_deviation){
  my_data <- x

  adduct_pos <- y

  my_data[, "Adduct"] <- "NA"
  my_data[, "Group"] <- "NA"

  RTdiff_set <- RT_deviation ############ SET MAXIMUM RT DIFFERENCE IN MINUTES - RECOMMENDED IS BETWEEN 0.015 AND 0.2##############
  Massdiff_set <- mass_deviation ############ SET MAXIMUM MASS DIFFERENCE IN DALTONS - RECOMMENDED 0.005 ##############

  i<-1
  x<-1
  y<-1
  w<-1

  while (x <= nrow(my_data)) {

    while (y <= nrow(my_data)) {

      #Collect information about rt and m.z

      r.t_data_x <- my_data$RT[x]
      r.t_data_y <- my_data$RT[y]
      m.z_data_x <- my_data$m.z[x]
      m.z_data_y <- my_data$m.z[y]

      Massdiff <- my_data$m.z[x] - my_data$m.z[y] ############ CALCULATE MASS DIFFERENCE BETWEEN X AND Y ##############

      RTdiff <- abs(my_data$RT[x] - my_data$RT[y]) ############ CALCULATE RT DIFFERENCE BETWEEN X AND Y ##############

      charge_x <- my_data$Charge[x]
      charge_y <- my_data$Charge[y]

      if(RTdiff < RTdiff_set & charge_x == charge_y) {
        if(i <= nrow(adduct_pos)){
          if(abs(Massdiff - adduct_pos$m.z_diff[i]) < Massdiff_set){
            if(my_data$Group[x]=="NA"){
              my_data$Group[x] <- w
              my_data$Group[y] <- w
              w <- w+1
            }else{
              my_data$Group[y] <- my_data$Group[x]
              w <- w+1
            }
            if(my_data$m.z[x] > my_data$m.z[y]){
              my_data$Adduct[x] <- as.character(adduct_pos$Adduct_2[i])
              my_data$Adduct[y] <- as.character(adduct_pos$Adduct_1[i])
              y<-y+1
            } else {
              my_data$Adduct[x] <- as.character(adduct_pos$Adduct_1[i])
              my_data$Adduct[y] <- as.character(adduct_pos$Adduct_2[i])
              my_data$Group[x] <- x
              my_data$Group[y] <- x
              y<-y+1
            }
          }else{
            i<-i+1
          }
        }else{
          i <- 1
          y<-y+1
        }
      }else{
        i <- 1
        y<-y+1
      }


    }
    x <- x+1
    y <- 1
    i <- 1
    print(x)

  }
  return(my_data)
}

#' Prepare correlation matrix based on two input matrices
#' @param x 1st data frame containing (deconvoluted) elution profiles
#' @param y 2nd data frame containing (deconvoluted) elution profiles
#' @return A data frame containing calculated Pearson correlation coefficients for each peaks pair
#' @export
calc_correlation <- function(x, y){

  mat1 <- as.matrix(x) # format to matrix, without first row
  mat1 <- t(mat1) # transpose data set

  mat2 <- as.matrix(y) # format to matrix, without first row
  mat2 <- t(mat2) # transposte data set

  column_names <- as.matrix(colnames(mat1)) # make df as matrix
  row_names <- as.matrix(colnames(mat2)) # make df as matrix

  ###
  x <- 1 # for loop
  y <- 1 # for loop

  correlation_mat <- numeric() # create empty matrix to store correlation data
  correlation_final <- numeric() # create empty matrix to store correlation data


  for(y in 1:ncol(mat2)) {

    for(i in 1:ncol(mat1)) {
      correlation <- as.matrix(cor(mat1[,x],mat2[,y]))
      colnames(correlation) <- column_names[x,]
      rownames(correlation) <- row_names[y,]

      correlation_mat <- as.matrix(cbind(correlation_mat, correlation))

      x <- x+1
    }

    correlation_final <- as.matrix(rbind(correlation_final, correlation_mat))
    correlation_mat <- numeric()

    y <- y+1
    x <- 1
    #i <- 1

  }

  correlation_final_t <- t(correlation_final)
  return(correlation_final_t)
}

#' Filter for unique features (it take roughly 5-10minutes)
#' @param FUF_df Data frame containing elution profiles and cluster identifiers
#' @param FUF_info Data frame containing info regarding the clusters
#' @param nr_replicas Define how many replicas were used in the experiment
#' @param nr_fraction Define how many fractions were collected in the experiment
#' @param start_column Define position of the column containing molecule intensity measured for the first collected fraction
#' @param end_column Define position of the column containing molecule intensity measured for the last collected fraction
#' @return A data frame containing unique mass features
#' @export
find_unique_features <- function(FUF_df, FUF_info, nr_replicas, nr_fractions, start_column, end_column){

  my_data <- FUF_df
  rownames(my_data) <- my_data$Peak.ID
  my_data <- my_data[,-1]
  my_data_info <- FUF_info

  list_2frc <- as.data.frame(matrix(0, ncol = 0, nrow = nrow(my_data)))
  rnames <- rownames(my_data)
  list_2frc <- cbind(list_2frc, rnames)

  for(i in 1:nr_replicas){
    tmp <- paste("Rep_", i, sep = "")

    data_temp <- my_data[,grep(tmp, colnames(my_data))]

    names <- paste(tmp, "_Fraction_", c(1:nr_fractions), sep = "")

    colnames(data_temp) <- names

    temp_2frc <- as.data.frame(NULL)
    temp_2frc <- filterConsistentFeatures_ML(data_temp, 1) # filter for metabolites which peaks span on at least 2 fractions

    temp_long_peaks <- filterPeaksSpanFractions_ML(temp_2frc) # remove single appearance of metabolite in fractions 1111010111 then 1 in the middle will be remove since it doesn't create peak spanning on 2 fractions

    rnames <- rownames(temp_long_peaks)
    temp_long_peaks <- cbind(temp_long_peaks, rnames)

    list_2frc <- full_join(list_2frc, temp_long_peaks, by = "rnames")
  }
  ###
  rownames(list_2frc) <- list_2frc$rnames

  list_2frc <- list_2frc[,-1]

  tmp <- NULL

  for(i in 1:nrow(list_2frc)){
    if(!all(is.na(list_2frc[i,]))){
      tmp <- rbind(tmp, list_2frc[i,])
    }
  }

  tmp[is.na(tmp)] <- 0
  tmp$Sum <- apply(tmp, 1, sum)

  # Merge file with intensities and with cluster information

  tmp$Peak.ID <- rownames(tmp)

  PP1_data <- inner_join(tmp, my_data_info, by = "Peak.ID")

  #write.table(PP1_data, "PP1_PCF_longpeak2020.txt", sep = "\t")

  ################################################
  ### Remove redundancy in pos or neg mode
  ################################################

  x <- 1 ############ STARTING POINT FOR LOOP - DO NOT CHANGE ##############
  y <- 1 ############ STARTING POINT FOR LOOP - DO NOT CHANGE ##############
  w <- start_column ############ STARTING POINT OF YOUR DATASET - CHANGE IT ##############
  z <- end_column ############ ENDING POINT OF YOUR DATASET - CHANGE IT ##############

  ########################### YOU DO NOT HAVE TO CHANGE ANYTHING HERE ###########################

  PP1_Group0 <- filter(PP1_data, Group == 0)
  PP1_Group <- filter(PP1_data, Group != 0)

  pb <- winProgressBar(title="Example progress bar", label="0% done", min=0, max=100, initial=0)

  while (x <= nrow(PP1_Group)) {

    while (y <= nrow(PP1_Group)) {


      #Collect information about to which group x or y were assigned
      groupx <- as.matrix(PP1_Group$Group[x])
      groupy <- as.matrix(PP1_Group$Group[y])

      sumx <- PP1_Group$Sum[x] ############ TAKE SUM OF THE INTENSITY IN ROW X ##############
      sumy <- PP1_Group$Sum[y] ############ TAKE SUM OF THE INTENSITY IN ROW Y ##############

      identifierx <- as.matrix(PP1_Group$Peak.ID[x])
      identifiery <- as.matrix(PP1_Group$Peak.ID[y])

      PCC <- cor(t(PP1_Group[x,w:z]), t(PP1_Group[y,w:z]))

      if(identifierx != identifiery & groupx == groupy & PCC >= 0.9) {

        if(sumx > sumy){
          PP1_Group <- PP1_Group[-y,]
          y <- y+1
        } else {
          PP1_Group <- PP1_Group[-x,]
        }

      } else {
        y <- y+1
      }


    }
    x <- x+1
    y <- 1

    info <- sprintf("%d%% done", round((x/nrow(PP1_Group))*100))
    setWinProgressBar(pb, x/(nrow(PP1_Group))*100, label=info)
  }
  PP1_c <- rbind(PP1_Group, PP1_Group0)

  # Filter for features present in 3/3 replicas and having reproducible elution profile in at least 2 replicas
  PP1_c_copy <- PP1_c

  i <- 1
  while(i <= nrow(PP1_c)){
    present_in_replicas <- cbind(max(PP1_c[i, grep("Rep_1", colnames(PP1_c))]),
                                 max(PP1_c[i, grep("Rep_2", colnames(PP1_c))]),
                                 max(PP1_c[i, grep("Rep_3", colnames(PP1_c))]))
    PP1_c_copy$freq[i] <- sum(present_in_replicas > 0)
    i <- i+1
  }

  PP1_c_frequent <- filter(PP1_c_copy, freq == 3)

  cor_collect <- as.data.frame(matrix(0, ncol = 4, nrow = nrow(PP1_c_frequent)))
  colnames(cor_collect) <- c("1vs2", "1vs3", "2vs3", "Max")
  i <- 1
  while(i <= nrow(PP1_c_frequent)){
    data_row <- cbind(t(PP1_c_frequent[i, grep("Rep_1", colnames(PP1_c_frequent))]),
                      t(PP1_c_frequent[i, grep("Rep_2", colnames(PP1_c_frequent))]),
                      t(PP1_c_frequent[i, grep("Rep_3", colnames(PP1_c_frequent))]))
    colnames(data_row) <- c("Rep_1", "Rep_2", "Rep_3")
    cor_tmp <- cor(data_row, method = "pearson")
    cor_tmp <- unique(as.numeric(cor_tmp))
    cor_tmp <- cor_tmp[cor_tmp != 1]
    cor_tmp <- c(cor_tmp, max(cor_tmp))
    cor_collect[i,] <- cor_tmp
    i <- i+1
  }

  unique_features <- filter(cor_collect, Max >= 0.9)

  return(unique_features)

}

#' Filter for protein peaks with apparent mass that exceeds theoretical molecular weight of a protein, which indicates formation of a complex
#' @param deconv_profiles A path to file containing deconvoluted protein profiles
#' @param prot_theor_mass A path to file containing theoretical mass of a proteins
#' @param column_calibration A path to file containing column calibration curve based on mass of reference proteins
#' @param oligomeric_state_ratio The threshold to be used to determine proteins present in molecular complexes
#' @return A data frame with proteins, which oligomeric state ratio above the set threshold
#' @export
find_proteins_in_complex <- function(deconv_profiles, prot_theor_mass, column_calibration, oligomeric_state_ratio){

  mass_calibration <- read.delim(column_calibration, header = TRUE, check.names = FALSE, row.names = 1)
  protein_yeast_mine <- read.delim(prot_theor_mass, header = TRUE, check.names = FALSE)
  prot_prof <- read.delim(deconv_profiles, header = TRUE, check.names = FALSE, row.names = 1)

  prot_prof_subset <- prot_prof[,1:38]
  max_fraction <- mass_calibration$Mass[max.col(prot_prof_subset)]

  prot_prof <- cbind(prot_prof, max_fraction)

  SEC_mine <- inner_join(prot_prof, protein_yeast_mine, by = c("Gene" = "Name"))
  ratio <- SEC_mine$max_fraction/SEC_mine$MW
  SEC_mine <- cbind(SEC_mine, ratio)
  SEC_mine_complex <- filter(SEC_mine, ratio >= oligomeric_state_ratio)
  rownames(SEC_mine_complex) <- SEC_mine_complex$Code
  return(SEC_mine_complex)
}

#' Find elution profiles of known protein-protein complexes
#' @param deconv_profiles A path to file containing deconvoluted protein profiles
#' @param prot_theor_mass A path to file containing theoretical mass of a proteins
#' @param column_calibration A path to file containing column calibration curve based on mass of reference proteins
#' @param oligomeric_state_ratio The threshold to be used to determine proteins present in molecular complexes
#' @param ref_database A path to string database file
#' @param def_heatmap_name A path to save heatmap showing elution profiles of preselected protein-protein complexes
#' @return Heatmap and a data frame containing elution profiles of preselected known protein-protein complexes
#' @export
find_known_PPI <- function(deconv_profiles, prot_theor_mass, column_calibration,
                           oligomeric_state_ratio, ref_database, def_heatmap_name){

  mass_calibration <- read.delim(column_calibration, header = TRUE, check.names = FALSE, row.names = 1)
  protein_yeast_mine <- read.delim(prot_theor_mass, header = TRUE, check.names = FALSE)
  prot_prof <- read.delim(deconv_profiles, header = TRUE, check.names = FALSE, row.names = 1)
  protein_database <- read.table(ref_database, header = TRUE, check.names = FALSE)

  prot_prof_subset <- prot_prof[,1:38]
  max_fraction <- mass_calibration$Mass[max.col(prot_prof_subset)]

  prot_prof <- cbind(prot_prof, max_fraction)

  SEC_mine <- inner_join(prot_prof, protein_yeast_mine, by = c("Gene" = "Name"))
  ratio <- SEC_mine$max_fraction/SEC_mine$MW
  SEC_mine <- cbind(SEC_mine, ratio)
  SEC_mine_complex <- filter(SEC_mine, ratio >= oligomeric_state_ratio)
  rownames(SEC_mine_complex) <- SEC_mine_complex$Code

  #take subset of the combined profiles and yeast mine and prepare correlation matrix
  complexes_subset <- SEC_mine_complex[,1:38]
  complexes_subset_t <- t(complexes_subset)
  cor_complexes <- cor(complexes_subset_t)

  #prepare list of unique protein names appearing in the dataset
  tmp1 <- data.frame(SEC_mine_complex$Gene)
  tmp2 <- data.frame(SEC_mine_complex$protein)
  protein_names <- cbind(tmp1, tmp2)
  colnames(protein_names) <- c("Name", "Protein")
  protein_names_unique <- unique(protein_names)

  #number of unique proteins creating PPI#

  unique_prot_PPI <- unique(SEC_mine_complex$Gene)

  #determine overlap between SEC and STRING

  overlap_SEC_string <- inner_join(protein_names_unique, protein_database, by = c("Protein" = "protein1")) # protein x
  overlap_SEC_string <- inner_join(overlap_SEC_string, protein_names_unique, by = c("protein2" = "Protein")) # protein y

  overlap_SEC_string_subset <- filter(overlap_SEC_string, experiments >= 800) # subset of string database

  #####remove redundancy from the data PER1 - PER3 vs PER3 - PER1
  i <- 1
  while(i <= nrow(overlap_SEC_string_subset)){
    protein_x <- as.character(overlap_SEC_string_subset$Name.x[i])
    protein_y <- as.character(overlap_SEC_string_subset$Name.y[i])
    prot_delete <- which(overlap_SEC_string_subset$Name.x == protein_y & overlap_SEC_string_subset$Name.y == protein_x)
    if(length(prot_delete) == 0){
      i <- i+1
    }else{
      overlap_SEC_string_subset <- overlap_SEC_string_subset[-prot_delete,]
      i <- i+1
    }
    print(paste((i/nrow(overlap_SEC_string))*100, "%", sep = ""))
  }

  ##Collect correlation coefficient for all proteins that are present in String and SEC database
  cor_value_collect <- NULL
  i <- 1
  while(i <= nrow(overlap_SEC_string_subset)){
    prot_x <- paste(overlap_SEC_string_subset$Name.x[i], "_", sep = "")
    prot_y <- paste(overlap_SEC_string_subset$Name.y[i], "_", sep = "")
    cor_value <- data.frame(cor_complexes[,grep(prot_x, colnames(cor_complexes), fixed = TRUE)])
    cor_value_exact <- data.frame(cor_value[grep(prot_y, rownames(cor_value), fixed = TRUE),])
    cor_value_exact_max <- max(cor_value_exact)

    cor_value_collect <- data.frame(rbind(cor_value_collect, cor_value_exact_max))
    i <- i+1
    print((i/nrow(overlap_SEC_string_subset))*100)
  }

  colnames(cor_value_collect) <- "Correlation"
  overlap_SEC_string_subset <- cbind(overlap_SEC_string_subset, cor_value_collect)

  known_captured_sec <- filter(overlap_SEC_string_subset, Correlation >= 0.8)

  ################ Analysis of representative yeast complexes

  rep_complexes <- select(known_captured_sec, Name.x, Protein, Name.y, protein2, Correlation, combined_score)
  colnames(rep_complexes) <- c("Protein_1", "Gene_1", "Protein_2", "Gene_2", "Correlation", "Combined_score")

  #Proteasome 19S

  proteasome_26S <- SEC_mine_complex[grep("26S proteasome", SEC_mine_complex$Description, ignore.case = TRUE),]
  proteasome_26S_intensity <- proteasome_26S[,1:38]
  proteasome_26S_max <- maxNormalize(proteasome_26S_intensity)
  max_temp <- max.col(proteasome_26S_intensity)
  proteasome_26S_max_order <- proteasome_26S_max[order(max_temp),]
  proteasome_26S_HM <- proteasome_26S_max_order[which(max.col(proteasome_26S_max_order) == 5),]
  proteasome_26S_HM_2 <- proteasome_26S_max_order[which(max.col(proteasome_26S_max_order) == 5 | max.col(proteasome_26S_max_order) == 6),]
  tmp_class <- as.data.frame(matrix(0, ncol = 1, nrow = nrow(proteasome_26S_HM)))
  rownames(tmp_class) <- rownames(proteasome_26S_HM)
  tmp_class[,1] <- "19S proteasome"
  colnames(tmp_class) <- "Complex"

  hm_all <- proteasome_26S_HM
  class_all <- tmp_class

  #Proteasome 20S

  proteasome_20S <- SEC_mine_complex[grep("20S proteasome", SEC_mine_complex$Description, ignore.case = TRUE),]
  proteasome_20S_intensity <- proteasome_20S[,1:38]
  proteasome_20S_max <- maxNormalize(proteasome_20S_intensity)
  max_temp <- max.col(proteasome_20S_intensity)
  proteasome_20S_max_order <- proteasome_20S_max[order(max_temp),]
  proteasome_20S_HM <- proteasome_20S_max_order[which(max.col(proteasome_20S_max_order) == 9),]
  tmp_class <- as.data.frame(matrix(0, ncol = 1, nrow = nrow(proteasome_20S_HM)))
  rownames(tmp_class) <- rownames(proteasome_20S_HM)
  tmp_class[,1] <- "20S proteasome"
  colnames(tmp_class) <- "Complex"

  hm_all <- rbind(hm_all, proteasome_20S_HM)
  class_all <- rbind(class_all, tmp_class)

  #Anaphase promoting complex

  APC <- SEC_mine_complex[grep("APC", SEC_mine_complex$Description, ignore.case = TRUE),]
  APC_intensity <- APC[,1:38]
  APC_max <- maxNormalize(APC_intensity)
  max_temp <- max.col(APC_intensity)
  APC_max_order <- APC_max[order(max_temp),]
  APC_HM <- APC_max_order[which(max.col(APC_max_order) == 5) ,]
  tmp_class <- as.data.frame(matrix(0, ncol = 1, nrow = nrow(APC_HM)))
  rownames(tmp_class) <- rownames(APC_HM)
  tmp_class[,1] <- "Anaphase promoting complex"
  colnames(tmp_class) <- "Complex"

  hm_all <- rbind(hm_all, APC_HM)
  class_all <- rbind(class_all, tmp_class)

  #Arp2/3 protein complex

  ARP2 <- SEC_mine_complex[grep("ARP2/3", SEC_mine_complex$Description, ignore.case = TRUE),]
  ARP2_intensity <- ARP2[,1:38]
  ARP2_max <- maxNormalize(ARP2_intensity)
  max_temp <- max.col(ARP2_intensity)
  ARP2_max_order <- ARP2_max[order(max_temp),]
  ARP2_HM <- ARP2_max_order[which(max.col(ARP2_max_order) == 19) ,]
  tmp_class <- as.data.frame(matrix(0, ncol = 1, nrow = nrow(ARP2_HM)))
  rownames(tmp_class) <- rownames(ARP2_HM)
  tmp_class[,1] <- "ARP2/3 protein complex"
  colnames(tmp_class) <- "Complex"

  hm_all <- rbind(hm_all, ARP2_HM)
  class_all <- rbind(class_all, tmp_class)

  #CCT

  CCT <- SEC_mine_complex[grep("CCT", SEC_mine_complex$Alias, ignore.case = TRUE),]
  CCT_intensity <- CCT[,1:38]
  CCT_max <- maxNormalize(CCT_intensity)
  max_temp <- max.col(CCT_intensity)
  CCT_max_order <- CCT_max[order(max_temp),]
  CCT <- CCT[order(max_temp),]
  CCT_HM <- CCT_max_order[which(max.col(CCT_max_order) < 7) ,]
  CCT_HM <- CCT_HM[which(max.col(CCT_HM) > 3) ,]
  tmp_class <- as.data.frame(matrix(0, ncol = 1, nrow = nrow(CCT_HM)))
  rownames(tmp_class) <- rownames(CCT_HM)
  tmp_class[,1] <- "CCT complex"
  colnames(tmp_class) <- "Complex"

  hm_all <- rbind(hm_all, CCT_HM)
  class_all <- rbind(class_all, tmp_class)

  #Coatomer

  Coatomer <- SEC_mine_complex[grep("Coatomer", SEC_mine_complex$Alias, ignore.case = TRUE),]
  Coatomer_intensity <- Coatomer[,1:38]
  Coatomer_max <- maxNormalize(Coatomer_intensity)
  max_temp <- max.col(Coatomer_intensity)
  Coatomer_max_order <- Coatomer_max[order(max_temp),]
  Coatomer_HM <- Coatomer_max_order[which(max.col(Coatomer_max_order) < 8) ,]
  tmp_class <- as.data.frame(matrix(0, ncol = 1, nrow = nrow(Coatomer_HM)))
  rownames(tmp_class) <- rownames(Coatomer_HM)
  tmp_class[,1] <- "Coatomer complex"
  colnames(tmp_class) <- "Complex"

  hm_all <- rbind(hm_all, Coatomer_HM)
  class_all <- rbind(class_all, tmp_class)

  #COG - 1

  COG <- SEC_mine_complex[grep("oligomeric Golgi", SEC_mine_complex$Description, ignore.case = TRUE),]
  COG_intensity <- COG[,1:38]
  COG_max <- maxNormalize(COG_intensity)
  max_temp <- max.col(COG_intensity)
  COG_max_order <- COG_max[order(max_temp),]
  COG <- COG[order(max_temp),]
  COG_HM <- COG_max_order[which(max.col(COG_max_order) == 5) ,]
  tmp_class <- as.data.frame(matrix(0, ncol = 1, nrow = nrow(COG_HM)))
  rownames(tmp_class) <- rownames(COG_HM)
  tmp_class[,1] <- "COG complex - 1"
  colnames(tmp_class) <- "Complex"

  hm_all <- rbind(hm_all, COG_HM)
  class_all <- rbind(class_all, tmp_class)

  #COG - 2

  COG2 <- SEC_mine_complex[grep("oligomeric Golgi", SEC_mine_complex$Description, ignore.case = TRUE),]
  COG2_intensity <- COG2[,1:38]
  COG2_max <- maxNormalize(COG2_intensity)
  max_temp <- max.col(COG2_intensity)
  COG2_max_order <- COG2_max[order(max_temp),]
  COG2 <- COG2[order(max_temp),]
  COG2_HM <- COG2_max_order[which(max.col(COG2_max_order) == 9) ,]
  tmp_class <- as.data.frame(matrix(0, ncol = 1, nrow = nrow(COG2_HM)))
  rownames(tmp_class) <- rownames(COG2_HM)
  tmp_class[,1] <- "COG complex - 2"
  colnames(tmp_class) <- "Complex"

  hm_all <- rbind(hm_all, COG2_HM)
  class_all <- rbind(class_all, tmp_class)

  #Cohesin

  Cohesin <- SEC_mine_complex[grep("Cohesin", SEC_mine_complex$Description, ignore.case = TRUE),]
  Cohesin_intensity <- Cohesin[,1:38]
  Cohesin_max <- maxNormalize(Cohesin_intensity)
  max_temp <- max.col(Cohesin_intensity)
  Cohesin_max_order <- Cohesin_max[order(max_temp),]
  Cohesin_HM <- Cohesin_max_order[which(max.col(Cohesin_max_order) == 14) ,]
  tmp_class <- as.data.frame(matrix(0, ncol = 1, nrow = nrow(Cohesin_HM)))
  rownames(tmp_class) <- rownames(Cohesin_HM)
  tmp_class[,1] <- "Cohesin complex"
  colnames(tmp_class) <- "Complex"

  hm_all <- rbind(hm_all, Cohesin_HM)
  class_all <- rbind(class_all, tmp_class)

  #Condensin

  Condensin <- SEC_mine_complex[grep("Condensin", SEC_mine_complex$Description, ignore.case = TRUE),]
  Condensin_intensity <- Condensin[,1:38]
  Condensin_max <- maxNormalize(Condensin_intensity)
  max_temp <- max.col(Condensin_intensity)
  Condensin_max_order <- Condensin_max[order(max_temp),]
  Condensin_HM <- Condensin_max_order[which(max.col(Condensin_max_order) < 7) ,]
  tmp_class <- as.data.frame(matrix(0, ncol = 1, nrow = nrow(Condensin_HM)))
  rownames(tmp_class) <- rownames(Condensin_HM)
  tmp_class[,1] <- "Condensin complex"
  colnames(tmp_class) <- "Complex"

  hm_all <- rbind(hm_all, Condensin_HM)
  class_all <- rbind(class_all, tmp_class)

  #Exocyst

  Exocyst <- SEC_mine_complex[grep("Exocyst", SEC_mine_complex$Description, ignore.case = TRUE),]
  Exocyst_intensity <- Exocyst[,1:38]
  Exocyst_max <- maxNormalize(Exocyst_intensity)
  max_temp <- max.col(Exocyst_intensity)
  Exocyst_max_order <- Exocyst_max[order(max_temp),]
  Exocyst <- Exocyst[order(max_temp),]
  Exocyst_HM <- Exocyst_max_order[which(max.col(Exocyst_max_order) == 5) ,]
  tmp_class <- as.data.frame(matrix(0, ncol = 1, nrow = nrow(Exocyst_HM)))
  rownames(tmp_class) <- rownames(Exocyst_HM)
  tmp_class[,1] <- "Exocyst complex"
  colnames(tmp_class) <- "Complex"

  hm_all <- rbind(hm_all, Exocyst_HM)
  class_all <- rbind(class_all, tmp_class)

  #Exosome - 1

  Exosome1 <- SEC_mine_complex[grep("Exosome", SEC_mine_complex$Description, ignore.case = TRUE),]
  Exosome1_intensity <- Exosome1[,1:38]
  Exosome1_max <- maxNormalize(Exosome1_intensity)
  max_temp <- max.col(Exosome1_intensity)
  Exosome1_max_order <- Exosome1_max[order(max_temp),]
  Exosome1 <- Exosome1[order(max_temp),]
  Exosome1_HM <- Exosome1_max_order[which(max.col(Exosome1_max_order) < 7) ,]
  Exosome1_HM <- Exosome1_HM[which(max.col(Exosome1_HM) > 3) ,]
  tmp_class <- as.data.frame(matrix(0, ncol = 1, nrow = nrow(Exosome1_HM)))
  rownames(tmp_class) <- rownames(Exosome1_HM)
  tmp_class[,1] <- "Exosome complex - 1"
  colnames(tmp_class) <- "Complex"

  hm_all <- rbind(hm_all, Exosome1_HM)
  class_all <- rbind(class_all, tmp_class)


  #Exosome - 2

  Exosome2 <- SEC_mine_complex[grep("Exosome", SEC_mine_complex$Description, ignore.case = TRUE),]
  Exosome2_intensity <- Exosome2[,1:38]
  Exosome2_max <- maxNormalize(Exosome2_intensity)
  max_temp <- max.col(Exosome2_intensity)
  Exosome2_max_order <- Exosome2_max[order(max_temp),]
  Exosome2 <- Exosome2[order(max_temp),]
  Exosome2_HM <- Exosome2_max_order[which(max.col(Exosome2_max_order) < 13) ,]
  Exosome2_HM <- Exosome2_HM[which(max.col(Exosome2_HM) > 9) ,]
  tmp_class <- as.data.frame(matrix(0, ncol = 1, nrow = nrow(Exosome2_HM)))
  rownames(tmp_class) <- rownames(Exosome2_HM)
  tmp_class[,1] <- "Exosome complex - 2"
  colnames(tmp_class) <- "Complex"

  hm_all <- rbind(hm_all, Exosome2_HM)
  class_all <- rbind(class_all, tmp_class)

  #Lsm2-8

  LSM <- SEC_mine_complex[grep("LSM", SEC_mine_complex$Description, ignore.case = TRUE),]
  LSM_intensity <- LSM[,1:38]
  LSM_max <- maxNormalize(LSM_intensity)
  max_temp <- max.col(LSM_intensity)
  LSM_max_order <- LSM_max[order(max_temp),]
  LSM <- LSM[order(max_temp),]
  LSM_HM <- LSM_max_order[which(max.col(LSM_max_order) == 12) ,]
  tmp_class <- as.data.frame(matrix(0, ncol = 1, nrow = nrow(LSM_HM)))
  rownames(tmp_class) <- rownames(LSM_HM)
  tmp_class[,1] <- "Lsm2-8 complex"
  colnames(tmp_class) <- "Complex"

  hm_all <- rbind(hm_all, LSM_HM)
  class_all <- rbind(class_all, tmp_class)

  #Prefoldin

  prefoldin <- SEC_mine_complex[grep("Prefoldin", SEC_mine_complex$Description, ignore.case = TRUE),]
  prefoldin_intensity <- prefoldin[,1:38]
  prefoldin_max <- maxNormalize(prefoldin_intensity)
  max_temp <- max.col(prefoldin_intensity)
  prefoldin_max_order <- prefoldin_max[order(max_temp),]
  prefoldin <- prefoldin[order(max_temp),]
  prefoldin_HM <- prefoldin_max_order[which(max.col(prefoldin_max_order) == 23) ,]
  tmp_class <- as.data.frame(matrix(0, ncol = 1, nrow = nrow(prefoldin_HM)))
  rownames(tmp_class) <- rownames(prefoldin_HM)
  tmp_class[,1] <- "Prefoldin complex"
  colnames(tmp_class) <- "Complex"

  hm_all <- rbind(hm_all, prefoldin_HM)
  class_all <- rbind(class_all, tmp_class)

  #septin

  septin <- SEC_mine_complex[grep("Septin", SEC_mine_complex$Description, ignore.case = TRUE),]
  septin_intensity <- septin[,1:38]
  septin_max <- maxNormalize(septin_intensity)
  max_temp <- max.col(septin_intensity)
  septin_max_order <- septin_max[order(max_temp),]
  septin <- septin[order(max_temp),]
  septin_HM <- septin_max_order[which(max.col(septin_max_order) < 7) ,]
  septin_HM <- septin_HM[which(max.col(septin_HM) > 4) ,]
  tmp_class <- as.data.frame(matrix(0, ncol = 1, nrow = nrow(septin_HM)))
  rownames(tmp_class) <- rownames(septin_HM)
  tmp_class[,1] <- "Septin complex"
  colnames(tmp_class) <- "Complex"

  hm_all <- rbind(hm_all, septin_HM)
  class_all <- rbind(class_all, tmp_class)

  pheatmap(hm_all,
           annotation_row = class_all,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           clustering_distance_rows = "correlation",
           clustering_method = "complete",
           show_rownames = FALSE,
           show_colnames = FALSE,
           fontsize = 6,
           fontsize_row = 6,
           fontsize_col = 6,
           main = "",
           #cellwidth = 6,
           filename = def_heatmap_name,
           height = 4.3,
           na_col = "grey",
           width = 3.64,
           color = colorpanel(256,low = "white", high = "#3366cc"))

  return(hm_all)

}

#' Find known protein-metabolite complexes in PROMIS dataset and prepare ROC curve
#' @param ref_library A path to library containing known or predicted protein-metabolite complexes
#' @param prot_theor_mass A path to file containing theoretical mass of a proteins
#' @param prot_db A path to file containing protein data
#' @param annotated_met A path to file containing deconvoluted and annotated metabolites
#' @param met_chem_id A path to file containing pubchem ids of metabolites
#' @param cor_matrix A path to correlation matrix showing similarity between protein and metabolite elution profiles
#' @return A data frame with known protein-metabolite complexes
#' @export
find_known_PMI <- function(ref_library, prot_theor_mass, prot_db, annotated_met, met_chem_id, cor_matrix){
  stitch_lib <- read.table(ref_library, header = TRUE)
  yeast_mine <- read.delim(prot_theor_mass, fill = TRUE)
  protein_data <- read.delim(prot_db)
  protein_names <- as.data.frame(protein_data$Gene.names)
  colnames(protein_names) <- "Name"
  protein_code <- inner_join(protein_names, yeast_mine, by = "Name")

  protein_interacting <- inner_join(protein_code, stitch_lib, by = "protein") # gives a data set containing all protein that are available in SEC data and their potential interaction partners

  # add metabolite data
  met_PCF <- read.delim(annotated_met)
  lib_met <- read.delim(met_chem_id)
  met_PCF <- unique(met_PCF)
  met_PCF_CID <- inner_join(met_PCF, lib_met, by = "Metabolite")
  met_PCF_CID <- unique(met_PCF_CID)

  # All possible interactions found when reference library was compared with PROMIS experiment
  pred_interaction <- inner_join(protein_interacting, met_PCF_CID, by = c("chemical" = "CID"))
  # List of predicted interactions
  pred_interaction_exp <- filter(pred_interaction, experimental >= 800)

  # Add correlation matrix
  correlation_matrix <- read.delim(cor_matrix, row.names = 1, header = TRUE, check.names = FALSE)

  # Collect correlation coefficient for all proteins that are present in created subset
  tmp_4_cor <- pred_interaction_exp
  cor_value_collect <- NULL
  i <- 1
  while(i <= nrow(pred_interaction_exp)){
    prot_x <- paste(pred_interaction_exp$Name[i], "_", sep = "")
    prot_y <- paste(pred_interaction_exp$Metabolite_nick[i], sep = "")
    cor_value <- select(correlation_matrix, starts_with(prot_y))
    cor_value_exact <- data.frame(cor_value[grep(prot_x, rownames(cor_value), fixed = TRUE),])
    if(nrow(cor_value_exact) == 0){
      tmp_4_cor$Correlation[i] <- NA
      i <- i+1
      print((i/nrow(pred_interaction))*100)
    }else{
      cor_value_exact_max <- max(cor_value_exact)
      tmp_4_cor$Correlation[i] <- cor_value_exact_max
      i <- i+1
      print((i/nrow(pred_interaction))*100)
    }
  }

  tmp_4_cor <- tmp_4_cor[is.finite(tmp_4_cor$Correlation),]


  #Calculate correlation by chance
  cor_by_chance <- NULL
  for(i in 1:100){
    tmp <- correlation_matrix[, sample(1:ncol(correlation_matrix), 1)]
    tmp <- sample(tmp, nrow(tmp_4_cor))
    tmp <- sort(tmp)
    cor_by_chance <- cbind(cor_by_chance, tmp)
  }
  cor_by_chance_median <- apply(cor_by_chance, 1, median)

  nr_list <- c(1:length(cor_by_chance_median))
  tmp_cor <- tmp_4_cor$Correlation
  tmp_cor <- tmp_cor[order(tmp_cor, decreasing = TRUE)]
  tmp_chance <- cor_by_chance_median[order(cor_by_chance_median, decreasing = TRUE)]

  #Prepare ROC curve
  tmp_chance[tmp_chance < 0] <- 0
  tmp_cor[tmp_cor < 0] <- 0
  tpr <- NULL
  fpr <- NULL
  threshold <- 0
  i <- 1
  curve <- data.frame(matrix(ncol = 4, nrow = 0))

  while(threshold <= 1) {
    SEC_captured_tpr <- tmp_cor[tmp_cor >= threshold]
    tpr <- (length(SEC_captured_tpr)/length(tmp_cor))*100

    SEC_captured_fpr <- tmp_chance[tmp_chance >= threshold]
    fpr <- (length(SEC_captured_fpr)/length(tmp_chance))*100

    curve_tmp <- cbind(i, threshold, tpr, fpr)
    curve <- rbind(curve, curve_tmp)
    threshold <- threshold + 0.01
    i <- i+1
  }

  colnames(curve) <- c("Number", "Correlation", "TPR", "FPR")
  curve$Enrichment <- curve$TPR/curve$FPR

  # Create ROC curve
  middlex <- seq(0, 100, by = 10)
  middley <- seq(0, 100, by = 10)

  jpeg("Prot-met_ROC.jpeg", width = 8, height = 8, units = 'cm', res = 600)
  par(mfrow = c(1,1), xpd = FALSE, mar=c(5.1, 4.1, 2.1, 2.1))
  plot(curve$FPR,
       curve$TPR,
       type = "b",
       xlab = "False positive rate",
       ylab = "True positive Rate",
       cex.lab = 1,
       cex.axis = 1,
       col = "blue",
       pch = 20,
       main = "",
       xaxt = 'n',
       yaxt = 'n',
       cex = 0.4)
  axis(1, at=seq(0, 100, 25), cex.axis = 1)
  axis(2, at=seq(0, 100, 25), cex.axis = 1, las = 1)
  #points(3.45, 16.09,
  #col = "red",
  #pch = 19,
  #cex = 0.6)
  par(new = TRUE)
  plot(middlex,
       middley,
       type = "l",
       lty = 2,
       axes = FALSE,
       ann = FALSE)
  dev.off()

  return(tmp_4_cor)
}

#' Find predicted protein-metabolite complexes in PROMIS dataset
#' @param ref_library A path to library containing known or predicted protein-metabolite complexes
#' @param prot_theor_mass A path to file containing theoretical mass of a proteins
#' @param prot_db A path to file containing protein data
#' @param annotated_met A path to file containing deconvoluted and annotated metabolites
#' @param met_chem_id A path to file containing pubchem ids of metabolites
#' @param cor_matrix A path to correlation matrix showing similarity between protein and metabolite elution profiles
#' @return A data frame with predicted protein-metabolite complexes, which find experimental evidence in PROMIS dataset
#' @export
find_predicted_PMI <- function(ref_library, prot_theor_mass, prot_db, annotated_met, met_chem_id, cor_matrix){
  stitch_lib <- read.table(ref_library, header = TRUE)
  yeast_mine <- read.delim(prot_theor_mass, fill = TRUE)
  protein_data <- read.delim(prot_db)
  protein_names <- as.data.frame(protein_data$Gene.names)
  colnames(protein_names) <- "Name"
  protein_code <- inner_join(protein_names, yeast_mine, by = "Name")

  protein_interacting <- inner_join(protein_code, stitch_lib, by = "protein") # gives a data set containing all protein that are available in SEC data and their potential interaction partners

  # add metabolite data
  met_PCF <- read.delim(annotated_met)
  lib_met <- read.delim(met_chem_id)
  met_PCF <- unique(met_PCF)
  met_PCF_CID <- inner_join(met_PCF, lib_met, by = "Metabolite")
  met_PCF_CID <- unique(met_PCF_CID)

  # All possible interactions found when reference library was compared with PROMIS experiment
  pred_interaction <- inner_join(protein_interacting, met_PCF_CID, by = c("chemical" = "CID"))
  # List of predicted interactions
  pred_interaction_exp <- filter(pred_interaction, combined_score >= 400)
  pred_interaction_exp <- filter(pred_interaction_exp, experimental >= 150)
  pred_interaction_exp <- filter(pred_interaction_exp, experimental < 800)

  # Add correlation matrix
  correlation_matrix <- read.delim(cor_matrix, row.names = 1, header = TRUE, check.names = FALSE)

  # Collect correlation coefficient for all proteins that are present in created subset
  tmp_4_cor <- pred_interaction_exp
  cor_value_collect <- NULL
  i <- 1
  while(i <= nrow(pred_interaction_exp)){
    prot_x <- paste(pred_interaction_exp$Name[i], "_", sep = "")
    prot_y <- paste(pred_interaction_exp$Metabolite_nick[i], sep = "")
    cor_value <- select(correlation_matrix, starts_with(prot_y))
    cor_value_exact <- data.frame(cor_value[grep(prot_x, rownames(cor_value), fixed = TRUE),])
    if(nrow(cor_value_exact) == 0){
      tmp_4_cor$Correlation[i] <- NA
      i <- i+1
      print((i/nrow(pred_interaction))*100)
    }else{
      cor_value_exact_max <- max(cor_value_exact)
      tmp_4_cor$Correlation[i] <- cor_value_exact_max
      i <- i+1
      print((i/nrow(pred_interaction))*100)
    }
  }

  tmp_4_cor <- tmp_4_cor[is.finite(tmp_4_cor$Correlation),]

  return(tmp_4_cor)
}

#' Determine accumulation pattern of metabolites during yeast growth
#' @param metabolites_df A path to file containing collected metabolite intensity
#' @param annotation_df A path to file containing annotation information
#' @param s_order A path to file containing correctly ordered and groupped sample names
#' @param tmp_file_name A directory in which a heatmap showing accumulation pattern of metabolites will be saved
#' @return A heatmap showin accumulation pattern of annotated metabolites during growth
#' @export
metabolites_accumulation <- function(metabolites_df, annotation_df, s_order, tmp_file_name){

  ### Load Data

  d1 <- read.delim(metabolites_df, check.names = FALSE)
  annotation <- read.delim(annotation_df)

  rownames(d1) <- d1$Peak.ID

  d1_order <- as.data.frame(read.delim(s_order))

  for_order <- as.character(d1_order$Name)

  d1_ordered <- d1[for_order]


  ### Filtering

  #filter intensity 10000 and save only features which appear in 31 samples
  freq_d1 <- rowSums(d1_ordered > 10000)

  df1 <- cbind(d1_ordered, freq_d1)

  df1 <- df1[df1$freq_d1 > 30,]
  df1 <- df1[,1:90]

  ### Normalization

  #Normalization to median
  df1_n <- medianNormalize(df1)

  #Pick annotated metabolites

  df1_a <- rownames_to_column(df1, "Peak.ID")
  df1_a <- inner_join(df1_a, annotation, by = "Peak.ID")

  df1_n_a <- rownames_to_column(df1_n, "Peak.ID")
  df1_n_a <- inner_join(df1_n_a, annotation, by = "Peak.ID")

  rownames(df1_a) <- df1_a$Entry_ID
  rownames(df1_n_a) <- df1_n_a$Entry_ID

  df1_a <- df1_a[,2:91]
  df1_n_a <- df1_n_a[,2:91]

  #### CALCULATE MEAN FOR EACH GROUP

  df1_n_a_mean <- df1_n_a

  samp <- d1_order
  groups <- factor(samp$Sample, levels = unique(samp$Sample))

  tmpR <- apply(df1_n_a_mean, 1, function(x) {tapply(unlist(x), groups, mean, na.rm=TRUE)})
  tmpR <- as.data.frame(t(tmpR))

  df1_n_a_mean <- tmpR

  ### Fold change to time 0 - control conditions
  df1_control <- df1_n_a_mean[,grep("30C", colnames(df1_n_a_mean))]

  control_time0_fold <- df1_control
  control_time0_fold[is.na(control_time0_fold)] <- 0

  control_time0_fold <- control_time0_fold[control_time0_fold[,1] > 0,]

  i <- 1
  control_time0_fold_collect <- data.frame(matrix(ncol = 0, nrow = nrow(control_time0_fold)))
  rownames(control_time0_fold_collect) <- rownames(control_time0_fold)
  while(i <= ncol(control_time0_fold)){
    control_time0_fold_tmp <- control_time0_fold[,i] / control_time0_fold[,1]
    control_time0_fold_collect <- cbind(control_time0_fold_collect, control_time0_fold_tmp)
    i <- i+1
  }

  control_time0_fold_collect_log2 <- log2(control_time0_fold_collect)
  control_time0_fold_collect_log2[control_time0_fold_collect_log2 == -Inf] <- 0

  colnames(control_time0_fold_collect_log2) <- colnames(control_time0_fold)
  control_time0_fold_collect_log2 <- control_time0_fold_collect_log2[order(rownames(control_time0_fold_collect_log2)),]

  response <- control_time0_fold_collect_log2

  paletteLength <- 100
  myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)

  myBreaks <- c(seq(min(response), 0, length.out=ceiling(paletteLength/2 + 1)),
                seq(max(response)/paletteLength, max(response), length.out=floor(paletteLength/2)))

  pheatmap(response,
           #annotation_col = class_all,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           clustering_distance_rows = "correlation",
           clustering_method = "complete",
           show_rownames = TRUE,
           fontsize_row = 7,
           fontsize_col = 8,
           main = "",
           width = 3.15,
           filename = tmp_file_name,
           height = 7.08,
           #cellheight = 8,
           #cellwidth = 12,
           na_col = "grey",
           #breaks = breaksList,
           #color = colorRampPalette(colorpanel(256, "blue", "white", "red"))(length(breaksList)),
           breaks = myBreaks,
           color = myColor,
           treeheight_row = 0,
           labels_col = c("0min", "15min", "30min", "45min", "60min", "120min", "180min", "240min", "360min", "1440min"))
  return(response)
}

#' Function to calculate the mean and the standard deviation (of sample) for each group
#' @param data A data frame
#' @param varname The name of a column containing the variable to be summarised
#' @param groupnames A vector of column names to be used as groupping variables
#' @return A data frame containing the mean and the standard deviation (of sample) for each group
#' @export
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

#' Determine 13C fraction intensity
#' @param enrichment_df A path to file with calculated enrichment
#' @param lab_annotation A path to file containing annotation information
#' @param s_order A path to file containing correctly ordered and groupped sample names
#' @param metabolite_poolsize A path to file containing calculated metabolite poolsizes
#' @param tmp_barplot_name A name of a pdf file with all the plots
#' @param tmp_dir_plot A directory in which single plots showing 13C fraction intensity for each metabolite will be saved
#' @return Generate barplots showing 13C fraction intensity of annotated metabolites
#' @export
analyze_13C_data <- function(enrichment_df, s_order, metabolite_poolsize, lab_annotation, tmp_barplot_name, tmp_dir_plot){
  #################prepare enrichment dataset, order it etc
  df_enrichment <- read.delim(enrichment_df)
  df_order <- read.delim(s_order)
  df_enrichment$NA_count <- rowSums(is.na(df_enrichment))
  df_enrichment <- filter(df_enrichment, NA_count != 72) # remove all rows which contain only NA

  df_enrichment <- dplyr::select(df_enrichment, Name_Analyte, as.character(df_order$Name))

  colnames(df_enrichment)[2:73] <- as.character(df_order$Sample_complex)

  groups <- factor(df_order$Sample2, levels = unique(df_order$Sample2))

  #################prepare poolsize dataset
  df_pool <- read.delim(metabolite_poolsize)

  #################prepare a list of peaks that are present in both datasets and calculate level of 13C metabolite
  met_list <- inner_join(df_enrichment, df_pool, by = "Name_Analyte")
  met_list <- as.character(met_list$Name_Analyte)

  df_enrichment[is.na(df_enrichment)] <- 0
  df_pool[is.na(df_pool)] <- 0

  i <- 1
  met_13_collect <- data.frame(matrix(NA, ncol = 0, nrow = 0))
  while(i <= length(met_list)){
    met_tmp_enr <- dplyr::filter(df_enrichment, Name_Analyte == met_list[i])
    met_tmp_pool <- dplyr::filter(df_pool, Name_Analyte == met_list[i])
    met_tmp_13C <- met_tmp_enr[,2:73]*met_tmp_pool[,2:73]
    met_tmp_13C$Name_Analyte <- met_list[i]
    met_13_collect <- rbind(met_13_collect, met_tmp_13C)
    i <- i+1
  }

  ### Select only these metabolites which are quite frequent in the experiment

  freq_13 <- rowSums(met_13_collect > 0)
  met_13_collect$Freq <- freq_13
  met_13_collect <- filter(met_13_collect, Freq >= 14)

  ### Annotate peaks
  df_annotation <- read.delim(lab_annotation)

  met_13_collect <- inner_join(met_13_collect, df_annotation, by = c("Name_Analyte" = "Peak.ID"))
  met_13_collect$Entry_ID <- as.character(met_13_collect$Entry_ID)

  names_met <- make.names(met_13_collect$Entry_ID, unique = TRUE)
  rownames(met_13_collect) <- names_met
  met_13_collect <- met_13_collect[,1:72]
  met_13_collect[is.na(met_13_collect)] <- 0

  ### Log transform data
  met_13_collect[met_13_collect == 0] <- NA
  met_13_collect <- log(met_13_collect, 2)
  met_13_collect[is.na(met_13_collect)] <- 0

  ##prepare functions need in loop##
  anova.2 <- function(x,y) anova(lm(x ~ y))$Pr[1]
  tukey.2 <- function(x,y) TukeyHSD(aov(x ~ y))

  getPairs <- function(g) {
    z <- levels(g)
    out <- character(length(z)*(length(z)-1)/2)
    k <- 1
    for(i in 1:(length(z)-1))
      for(j in (i+1):length(z)) {
        out[k] <- paste(z[j],z[i],sep='-')
        k <- k+1
      }
    out
  }

  ### PREPARE A DATAFRAME CONTAINING TUKEYHSD P-VALUES BUT ONLY FOR COMARISON BETWEEN THE SAMPLES FROM SPECIFIC TIME-POINTS AND NOT ALL VS ALL

  tuk_p_values_collect <- matrix(NA, ncol = nrow(met_13_collect), nrow = 0) # matrix to collect TukeyHSD
  condition <- c("_5min", "_15min", "_30min", "_45min", "_60min", "_120min", "_180min", "_240min")
  i <- 1

  while(i <= length(condition)){
    metab_tmp <- met_13_collect[,grep(condition[i], colnames(met_13_collect))] # subset df for a specific condition
    groups_tmp <- factor(groups[grep(condition[i], colnames(met_13_collect))], levels = unique(groups[grep(condition[i], colnames(met_13_collect))])) # subset groups for specific condition

    tuk <- apply(metab_tmp,1, tukey.2, groups_tmp) # perform tukey for a subset (specific condition)

    # Function to get a list of all group-pairwise combinations
    gp <- getPairs(groups_tmp)

    # extract the TukeyHSD p-values
    tuk.p.values <- sapply(tuk, function(x) {
      z <- x$y[,"p adj"]
      z <- z[gp]
      names(z) <- gp
      z
    })

    tuk_p_values_collect <- rbind(tuk_p_values_collect, tuk.p.values) # collect tukey p-values for a specific group
    i <- i+1
  }

  tuk_p_values_all <- tuk_p_values_collect # save tukey results as p-values

  ##############transform tuk_p_values into signif label
  tuk_p_values_signif <- tuk_p_values_collect
  tuk_p_values_signif[tuk_p_values_collect > 0.05] <- "NS"
  tuk_p_values_signif[tuk_p_values_collect <= 0.05] <- "*"
  tuk_p_values_signif[tuk_p_values_collect <= 0.01] <- "**"
  tuk_p_values_signif[tuk_p_values_collect == "NaN"] <- "NS"


  ###PREPARE BARPLOTS FOR ENRICHMENT
  tmp <- met_13_collect
  my_ggplot_list_barplot <- vector('list', nrow(tmp)) ##prepare empty list for barplots
  my_ggplot_list_barplot2 <- vector('list', nrow(tmp)) ##prepare empty list for barplots
  tmp <- rownames_to_column(tmp, var = "Name_Analyte")
  positions <- c("5", "15", "30", "45", "60", "120", "180", "240")

  i <- 1
  while(i <= nrow(tmp)){
    tmp_long <- gather(tmp[i,], key = "Sample", value = "Ratio", 2:73)
    tmp_long$Group1 <- df_order$Sample3
    tmp_long$Group2 <- df_order$Condition
    tmp_avg <- data_summary(tmp_long, varname = "Ratio", groupnames = "Group1")
    colnames(tmp_avg) <- c("Sample", "Ratio_avg", "SD")

    tmp_desc <- t(as.data.frame(strsplit(as.character(tmp_avg$Sample), split = "_")))
    colnames(tmp_desc) <- c("Group2", "Time")
    tmp_avg <- cbind(tmp_avg, tmp_desc[,1:2])
    tmp_avg$Time <- as.character(tmp_avg$Time)
    rownames(tmp_avg) <- NULL

    tmp_avg <- tmp_avg %>%
      arrange(factor(Group2, levels = c("C", "D", "A"))) %>%
      arrange(factor(Time, levels = c("5", "15", "30", "45", "60", "120", "180", "240")))

    tmp_long$Time <- groups
    #signif_test <- compare_means(Ratio ~ Group1, tmp_long, method = "t.test", group.by = "Time")
    signif_test <- data.frame(tuk_p_values_signif[,i])
    #signif_test <- data.frame(tuk_p_values_all[,i])
    colnames(signif_test) <- "p.signif"
    groups_max <- aggregate(tmp_avg$Ratio_avg, by = list(tmp_avg$Time), max)
    groups_max <- groups_max %>%
      arrange(factor(Group.1, levels = c("5", "15", "30", "45", "60", "120", "180", "240")))
    groups_SD_max <- aggregate(tmp_avg$SD, by = list(tmp_avg$Time), max)
    groups_SD_max <- groups_SD_max %>%
      arrange(factor(Group.1, levels = c("5", "15", "30", "45", "60", "120", "180", "240")))
    groups_max <- groups_max$x + groups_SD_max$x
    y_max_scale <- max(groups_max)+4

    df_y <- vector()
    w <- 1
    while(w <= length(groups_max)){
      df_y <- append(df_y, c(groups_max[w]+0.5, groups_max[w]+2, groups_max[w]+3.5))
      w <- w+1
    }

    #########define x position for line
    nr_time_points <- 8
    w <- 1
    x <- NULL
    while(w <= nr_time_points){
      x_tmp <- c(w-0.3, w-0.3, w)
      x <- c(x, x_tmp)
      w <- w+1
    }

    #########define xend position for line
    nr_time_points <- 8
    w <- 1
    xend <- NULL
    while(w <= nr_time_points){
      xend_tmp <- c(w, w+0.3, w+0.3)
      xend <- c(xend, xend_tmp)
      w <- w+1
    }

    ###########define position for asterisk
    nr_time_points <- 8
    w <- 1
    p_x_position <- NULL
    while(w <= nr_time_points){
      p_x_position_tmp <- c(w-0.15, w, w+0.15)
      p_x_position <- c(p_x_position, p_x_position_tmp)
      w <- w+1
    }
    aster_position = data.frame(x=x,
                                xend=xend,
                                p_x_position=p_x_position,
                                y=df_y,
                                p_y_position=df_y+0.05) # define how high on y scale asterisk is

    signif_test <- cbind(signif_test, aster_position)
    signif_test_tmp <- filter(signif_test, p.signif != "NS")

    tmp_avg$Group2 <- factor(tmp_avg$Group2, levels = c("C", "D", "A"))
    tmp_avg$Time <- factor(tmp_avg$Time, levels = c("5", "15", "30", "45", "60", "120", "180", "240"))

    p2 <- ggplot(data = tmp_avg, aes(x = Time, y=Ratio_avg, fill = Group2))+
      geom_bar(stat="identity", color="black",
               position = position_dodge())+
      labs(title=tmp$Name_Analyte[i], x="Time [min]", y = "13C fraction intensity [Log2]", fill = "Treatment")+
      theme_classic()+
      scale_fill_manual(values=c("grey90", "bisque3", "azure4"),
                        labels = c("Control", "Ser-Leu", "Ser + Leu"))+
      geom_errorbar(aes(ymin=Ratio_avg-SD, ymax=Ratio_avg+SD), width=.2,
                    position=position_dodge(.9))+
      annotate("segment", x = signif_test_tmp$x,
               xend = signif_test_tmp$xend,
               y = signif_test_tmp$y,
               yend = signif_test_tmp$y,
               colour = "black")+
      annotate("text", x = signif_test_tmp$p_x_position,
               y = signif_test_tmp$p_y_position,
               label = signif_test_tmp$p.signif)


    #print(p2)

    my_ggplot_list_barplot[[i]] <- p2

    p3 <- ggplot(tmp_avg, aes(x = Time, y=Ratio_avg, fill = Group2))+
      geom_bar(stat="identity", color="black",
               position = position_dodge(),
               size = 0.001)+
      #labs(x="Time [min]", y = "Intensity [Log2]")+
      labs(y = "13C fraction intensity [Log2]")+
      theme_classic()+
      scale_fill_manual(values=c("grey90", "bisque3", "azure4"),
                        labels = c("Control", "Ser-Leu", "Ser + Leu"))+
      geom_errorbar(aes(ymin=Ratio_avg-SD, ymax=Ratio_avg+SD), width=.2,
                    position=position_dodge(.9),
                    size = 0.2)+
      annotate("segment", x = signif_test_tmp$x,
               xend = signif_test_tmp$xend,
               y = signif_test_tmp$y,
               yend = signif_test_tmp$y,
               colour = "black",
               size = 0.1)+
      annotate("text", x = signif_test_tmp$p_x_position,
               y = signif_test_tmp$p_y_position,
               label = signif_test_tmp$p.signif,
               size = 2)+
      theme(axis.text.x = element_text(color = "grey20", size = 5, angle = 90, hjust = .5, vjust = .5, face = "plain"),
            axis.text.y = element_text(color = "grey20", size = 5, angle = 0, hjust = 1, vjust = 0, face = "plain"),
            #axis.title.x = element_text(color = "grey20", size = 6, angle = 0, hjust = .5, vjust = 0, face = "plain"),
            axis.title.y = element_text(color = "grey20", size = 4.5, angle = 90, hjust = .5, vjust = .5, face = "plain"),
            legend.title = element_blank(),
            legend.position = "none",
            #legend.key.size = unit(0.5, "cm"),
            axis.title = element_blank())
    #print(p3)

    my_ggplot_list_barplot2[[i]] <- p3


    i <- i+1

  }

  pdf(tmp_barplot_name, onefile = TRUE)
  my_ggplot_list_barplot
  dev.off()

  for(i in 1:nrow(tmp)){
    ggsave(filename = paste0(tmp_dir_plot, as.character(tmp$Name_Analyte[i]), ".jpeg"), plot = my_ggplot_list_barplot2[[i]], dpi = 1200,
           width = 5, height = 2.5, units = "cm")
  }

  return(tuk_p_values_signif)

}
