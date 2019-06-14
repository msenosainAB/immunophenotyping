#---------------------------------------------------------------------------------------
# SET UP ENVIRONMENT
#---------------------------------------------------------------------------------------

# Installing/loading dependencies
env_load <- function(){
    require(dplyr)
    require(psych)
    require(flowCore)
    require(tools)
    require(tidyverse)
    require(purrr)
    require(tibble)
    require(FlowSOM)
    require(cytofkit)
    require(parallel)
    require(MASS)
    require(magrittr)
    require(gplots)
    require(RColorBrewer)
}


#---------------------------------------------------------------------------------------
# READ DATA
#---------------------------------------------------------------------------------------

# Function to read files into a list of data.frames
extract_data <- function (file_type = '.fcs',
                          sampling = F,
                          sample_size = 10) {
    files_list <- list.files(pattern=file_type)
    if (sampling == T){
        files_list <- sample(files_list, size = sample_size)
    }
    exprs_data <- list()
    if (file_type == '.fcs'){
        for (i in 1:length(files_list)){
            exprs_data[[i]] <-as.data.frame(flowCore::read.FCS(files_list[i], transformation = FALSE)@exprs)
        }}
    if (file_type == '.csv'){
        for (i in 1:length(files_list)){
            exprs_data[[i]] <-as.data.frame(read.csv(files_list[i]))
        }}
    if (file_type == '.txt'){
        for (i in 1:length(files_list)){
            exprs_data[[i]] <-as.data.frame(read.csv(files_list[i], sep = '\t'))
        }
    } else {
        warning("Data type not supported", call. = F)
    }
    names(exprs_data) <- tools::file_path_sans_ext(basename(files_list))
    return(exprs_data)
}

# Arcsinh transformation
t_asinh <- function(df, cofactor = 5) {
    nums <- vapply(df, is.numeric, FUN.VALUE = logical(1)) # only operates in numeric
    df[,nums] <- asinh(df[,nums]/cofactor)
    (df)
}


#---------------------------------------------------------------------------------------
# Determine flowSOM k
#---------------------------------------------------------------------------------------

# Metaclustering
#' MetaClustering
#'
#' Cluster data with automatic number of cluster determination for
#' several algorithms
#'
#' @param data   Matrix containing the data to cluster
#' @param method Clustering method to use
#' @param max    Maximum number of clusters to try out
#' @param ...    Extra parameters to pass along
#'
#' @return Numeric array indicating cluster for each datapoint
#' @seealso   \code{\link{metaClustering_consensus}}
#'
#' @examples
#'    # Read from file, build self-organizing map and minimal spanning tree
#'    fileName <- system.file("extdata","lymphocytes.fcs",package="FlowSOM")
#'    flowSOM.res <- ReadInput(fileName, compensate=TRUE,transform=TRUE,
#'                             scale=TRUE)
#'    flowSOM.res <- BuildSOM(flowSOM.res,colsToUse=c(9,12,14:18))
#'    flowSOM.res <- BuildMST(flowSOM.res)
#'
#'    # Apply metaclustering
#'    metacl <- MetaClustering(flowSOM.res$map$codes,
#'                             "metaClustering_consensus",
#'                             max=10)
#'
#'    # Get metaclustering per cell
#'    flowSOM.clustering <- metacl[flowSOM.res$map$mapping[,1]]
#'
#' @export
MetaClustering <- function(data,method,max=20,...){
    res <- DetermineNumberOfClusters(data,max,method,...)
    method <- get(method)
    method(data,k=res)
}

DetermineNumberOfClusters <- function(data,max,method,plot=FALSE,smooth=0.2,
                                      ...){
    # Try out a clustering algorithm for several numbers of clusters and
    # select optimal
    #
    # Args:
    #     data:     Matrix containing the data to cluster
    #     max:        Maximum number of clusters to try
    #     method: Clustering method to use
    #     plot:     Whether to plot the results for different k
    #     smooth: Smoothing option to find elbow:
    #             0: no smoothing, 1: maximal smoothing
    #
    # Returns:
    #     Optimal number of clusters
    if(method ==    "metaClustering_consensus"){
        results <- consensus(data,max,...)
        res <- rep(0,max)
        res[1] <- SSE(data,rep(1,nrow(data)))
        for(i in 2:max){
            c <- results[[i]]$consensusClass
            res[i] <- SSE(data,c)
        }
    } else {
        method <- get(method)
        res <- rep(0,max)
        for(i in 1:max){
            c <- method(data, k=i,...)
            res[i] <- SSE(data,c)
        }
    }

    for(i in 2:(max-1)){
        res[i] <- (1-smooth)*res[i]+(smooth/2)*res[i-1]+(smooth/2)*res[i+1]
    }


    if(plot) plot(1:max, res, type="b", xlab="Number of Clusters",
                  ylab="Within-cluster sum of squares", las=1)

    #abline(v = x, col = 'red')

    findElbow(res,max)
}

findElbow <- function(data, max){
    n <- length(data)
    data <- as.data.frame(cbind(1:n,data))
    colnames(data) <- c("X","Y")
    #r_plt <- c() #MF
    min_r <- Inf
    optimal <- 1
    for(i in 2:(n-1)){
        f1 <- stats::lm(Y~X,data[1:(i-1),])
        f2 <- stats::lm(Y~X,data[i:n,])
        r <- sum(abs(c(f1$residuals,f2$residuals)))
        #r_plt <- c(r_plt, r) #MF
        if(r < min_r){
            min_r <- r
            optimal <-i
        }
    }
    #plot(2:(max-1), r_plt, type="b", xlab="Number of Clusters",
    #ylab="Residuals") #MF
    #abline(v = optimal, col = 'red')
    optimal
}

#' MetaClustering
#'
#' Cluster data using hierarchical consensus clustering with k clusters
#'
#' @param data Matrix containing the data to cluster
#' @param k    Number of clusters
#' @param seed Seed to pass to consensusClusterPlus
#'
#' @return  Numeric array indicating cluster for each datapoint
#' @seealso \code{\link{MetaClustering}}
#' @examples
#'    # Read from file, build self-organizing map and minimal spanning tree
#'    fileName <- system.file("extdata","lymphocytes.fcs",package="FlowSOM")
#'    flowSOM.res <- ReadInput(fileName, compensate=TRUE,transform=TRUE,
#'                             scale=TRUE)
#'    flowSOM.res <- BuildSOM(flowSOM.res,colsToUse=c(9,12,14:18))
#'    flowSOM.res <- BuildMST(flowSOM.res)
#'
#'    # Apply consensus metaclustering
#'    metacl <- metaClustering_consensus(flowSOM.res$map$codes,k=10)
#'
#' @export
metaClustering_consensus <- function(data, k=7,seed=NULL){
    results <- suppressMessages(ConsensusClusterPlus::ConsensusClusterPlus(
        t(data),
        maxK=k, reps=100, pItem=0.9, pFeature=1,
        title=tempdir(), plot="pdf", verbose=FALSE,
        clusterAlg="hc", # "hc","km","kmdist","pam"
        distance="euclidean" ,
        #"euclidean","pearson","spearman","binary","maximum","canberra","minkowski"
        seed=seed
    ))

    results[[k]]$consensusClass
}

consensus <- function(data,max,...){
    results <- suppressMessages(ConsensusClusterPlus::ConsensusClusterPlus(
        t(data),
        maxK=max, reps=100, pItem=0.9, pFeature=1,
        title=tempdir(), plot="pdf", verbose=FALSE,
        clusterAlg="hc", # "hc","km","kmdist","pam"
        distance="euclidean"
        #"euclidean","pearson","spearman","binary","maximum","canberra","minkowski"
    ))
}

metaClustering_hclust <- function(data, k=7){
    d <- stats::dist(data, method = "minkowski")
    fit <- stats::hclust(d, method="ward.D2")
    stats::cutree(fit, k=k)
}

metaClustering_kmeans <- function(data, k=7){
    stats::kmeans(data, centers=k)$cluster
}

metaClustering_som <- function(data, k=7){
    s <- SOM(data,xdim=k,ydim=1,rlen=100)
    s$unit.classif
}

SSE <- function(data,clustering){
    if(class(clustering)!= "numeric")
        clustering <- as.numeric(as.factor(clustering))
    c_wss <- 0
    for(j in seq_along(clustering)){
        if(sum(clustering==j) > 1){
            c_wss <- c_wss + (nrow(data[clustering==j,,drop=FALSE])-1)*
                sum(apply(data[clustering==j,,drop=FALSE],2,stats::var))
        }
    }
    c_wss
}




#' F measure
#'
#' Compute the F measure between two clustering results
#'
#' @param realClusters Array containing real cluster labels for each sample
#' @param predictedClusters Array containing predicted cluster labels for each
#'                          sample
#' @param silent    Logical, if FALSE (default), print some information about
#'                  precision and recall
#'
#' @return  F measure score
#' @examples
#' # Generate some random data as an example
#' realClusters <- sample(1:5,100,replace = TRUE)
#' predictedClusters <- sample(1:6, 100, replace = TRUE)
#'
#' # Calculate the FMeasure
#' FMeasure(realClusters,predictedClusters)
#' @export
FMeasure <- function(realClusters, predictedClusters,silent=FALSE){
    if (sum(predictedClusters)==0)
        return(0);
    a <- table(realClusters, predictedClusters);
    p <- t(apply(a,1,function(x)x/colSums(a)))
    r <- apply(a,2,function(x)x/rowSums(a))
    if(!silent) message("Precision: ",
                        sum(apply(p,1,max) * (rowSums(a)/sum(a))),
                        "\nRecall: ",sum(apply(r,1,max) * (rowSums(a)/sum(a))),"\n")
    f <- 2*r*p / (r+p)
    f[is.na(f)] <- 0
    sum(apply(f,1,max) * (rowSums(a)/sum(a)))
}

#' MetaclusterMFIs
#'
#' Compute the median fluorescence intensities for the metaclusters
#'
#' @param fsom Result of calling the FlowSOM function
#' @return  Metacluster MFIs
#' @examples
#' fileName <- system.file("extdata","lymphocytes.fcs",package="FlowSOM")
#' ff <- flowCore::read.FCS(fileName)
#' ff <- flowCore::compensate(ff,ff@@description$SPILL)
#' ff <- flowCore::transform(ff,
#'          flowCore::transformList(colnames(ff@@description$SPILL),
#'                                 flowCore::logicleTransform()))
#' flowSOM.res <- FlowSOM(ff,scale=TRUE,colsToUse=c(9,12,14:18),maxMeta=10)
#' mfis <- MetaclusterMFIs(flowSOM.res)
#' @export
MetaclusterMFIs <- function(fsom){
    MFIs <- t(sapply(seq_along(levels(fsom$metaclustering)),
                     function(i) {
                         apply(subset(fsom$FlowSOM$data,
                                      fsom$metaclustering[fsom$FlowSOM$map$mapping[,1]] == i),
                               2,
                               stats::median)
                     }))
    rownames(MFIs) <- seq_len(nrow(MFIs))
    return(MFIs)
}

#' MetaclusterCVs
#'
#' Compute the coefficient of variation for the metaclusters
#'
#' @param fsom Result of calling the FlowSOM function
#' @return  Metacluster CVs
#' @examples
#' fileName <- system.file("extdata","lymphocytes.fcs",package="FlowSOM")
#' ff <- flowCore::read.FCS(fileName)
#' ff <- flowCore::compensate(ff,ff@@description$SPILL)
#' ff <- flowCore::transform(ff,
#'          flowCore::transformList(colnames(ff@@description$SPILL),
#'                                 flowCore::logicleTransform()))
#' flowSOM.res <- FlowSOM(ff,scale=TRUE,colsToUse=c(9,12,14:18), nClus=10)
#' cvs <- MetaclusterCVs(flowSOM.res)
#' @export
MetaclusterCVs <- function(fsom){
    CVs <- t(sapply(seq_along(levels(fsom$metaclustering)),
                    function(i) {
                        apply(subset(fsom$FlowSOM$data,
                                     fsom$metaclustering[fsom$FlowSOM$map$mapping[,1]] == i),
                              2,
                              function(y){
                                  if(length(y) > 0 && mean(y) != 0){
                                      stats::sd(y)/mean(y)
                                  } else {
                                      NA
                                  }})
                    }))
    return(CVs)
}

#---------------------------------------------------------------------------------------
# FlowSOM in functions
#---------------------------------------------------------------------------------------

fsom_k <- function(data_set, cols_n, ...){
    data_set <- dir(pattern = ".fcs")
    fSOM <- FlowSOM::ReadInput(data_set, pattern = '.fcs', transform = TRUE,
                  toTransform = cols_n,
                  transformFunction = flowCore::arcsinhTransform(a=0, b=0.0067))
    fSOM <- FlowSOM::BuildSOM(fSOM, colsToUse = cols_n)
    fSOM <- FlowSOM::BuildMST(fSOM, tSNE=FALSE)
    k <- DetermineNumberOfClusters(fSOM$map$codes,max=50,'metaClustering_consensus',
        plot=T,smooth=0.2, seed=42)
    return(k)
}

fsom_wp <- function(data_set, cols_n, k, ...){
    x = FlowSOM::FlowSOM(data_set,
                # Input options:
                compensate = F,transform = TRUE, toTransform=cols_n,
                transformFunction = flowCore::arcsinhTransform(a=0, b=0.0067), scale = TRUE,
                # SOM options:
                colsToUse = cols_n, xdim = 7, ydim = 7,
                # Metaclustering options:
                nClus = k,
                # Seed for reproducible results:
                seed = 101)
    return(x)
}

#---------------------------------------------------------------------------------------
# Functions to work with fSOM output
#---------------------------------------------------------------------------------------

# Extract data and generate a df with sample and cluster ID cols
fsom_data <- function(fSOM, data_set){
    pts_fS <- as.data.frame(fSOM$FlowSOM$data)
    data_set_n <- gsub(".fcs", "", data_set, perl = TRUE)
    
    # Add sample and cluster id columns
    pts_fS['TP_ID'] <- 1
    for (i in 1:length(data_set)){
    pts_fS[fSOM$FlowSOM$metaData[[i]][1]:fSOM$FlowSOM$metaData[[i]][2], 'TP_ID'] <- data_set_n[i]
    }
    pts_fS['cluster_ID'] <- fSOM[[2]][fSOM[[1]]$map$mapping[,1]]
    return(pts_fS)
}


# Generate a summary table
cluster_summary <- function(pts_fS, sm_stat = median, write_CSV = TRUE){
    flowSOM_median <- pts_fS %>%
                        group_by(cluster_ID) %>%
                        summarise_if(is.numeric, sm_stat, na.rm = TRUE)
    if (write_CSV == TRUE) {
        write.csv(flowSOM_median, file ='flowSOM_median.csv', row.names = F)
    }
    return(flowSOM_median)
}

# Get sample size table
smp_size <- function(pts_fS, data_set){
    sample_size <- c()
    tp_smp <- unique(pts_fS$TP_ID)
    data_set_n <- gsub(".fcs", "", data_set, perl = TRUE)
    for (i in 1:length(data_set)){
        sample_size <- c(sample_size, length(which(pts_fS$TP_ID==tp_smp[i])))
    }
    df <- data.frame(matrix(NA, nrow = length(data_set_n), ncol = 2))
    df[,1] <- data_set_n
    df[,2] <- sample_size
    colnames(df) <- c('sample_ID', 'size')
    return(df)
}

# Get cluster size table
cl_size <- function(pts_fS){
    cl_size <- c()
    for (i in 1:length(unique(pts_fS$cluster_ID))){
        cl_size <- c(cl_size, length(which(pts_fS$cluster_ID==i)))
    }
    df <- data.frame(matrix(NA, nrow = length(cl_size), ncol = 2))
    df[,1] <- c(1:length(cl_size))
    df[,2] <- cl_size
    colnames(df) <- c('cluster_ID', 'size')
    return(df)
}

# Generate a table with quantities of clusters and samples
qnts <- function(pts_fS, data_set, write_CSV = TRUE){
    data_set_n <- gsub(".fcs", "", data_set, perl = TRUE)
    k <- length(unique(pts_fS$cluster_ID))
    tp <- length(unique(pts_fS$TP_ID))
    df_q <- data.frame(matrix(NA, nrow = k, ncol = tp))
    for(i in 1:k){
        for (ii in 1:tp){
            df_q[i,ii] <- length(which(pts_fS$cluster_ID == i & pts_fS$TP_ID == tp_smp[ii]))
        }   
    }
    colnames(df_q) <- data_set_n
    rownames(df_q) <- paste('cluster_',1:k, sep = '')
    if (write_CSV == TRUE) {
        write.csv(df_q, file ='summary_clusters.csv', row.names = F)
    return(df_q)
    }
}

# Generate a table with percentages of cluster per sample
clust_p_smp <- function(pts_fS, df_q, write_CSV = TRUE){
    k <- length(unique(pts_fS$cluster_ID))
    tp <- length(unique(pts_fS$TP_ID))
    df_per_tp <- df_q
    for(i in 1:tp){
        sum_tp <- sum(df_q[,i])
        for (ii in 1:k){
            df_per_tp[ii,i] <- df_q[ii,i]/sum_tp
        }
    }
    if (write_CSV == TRUE) {
        write.csv(df_per_tp, file ='summary_clust_p_smp.csv', row.names = F)
    return(df_per_tp)
    }
}


# Generate a table with percentages of sample per cluster
smp_p_clust <- function(pts_fS, df_q, write_CSV = TRUE){
    k <- length(unique(pts_fS$cluster_ID))
    tp <- length(unique(pts_fS$TP_ID))
    df_per_k <- df_q
    for(i in 1:k){
        sum_k <- sum(df_q[i,])
        for (ii in 1:tp){
            df_per_k[i,ii] <- df_q[i,ii]/sum_k
        }
    }
    if (write_CSV == TRUE) {
        write.csv(df_per_k, file ='summary_smp_p_clust.csv', row.names = F)
    return(df_per_k)
    }
}
