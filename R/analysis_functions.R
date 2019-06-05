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
}


#---------------------------------------------------------------------------------------
# READ DATA
#---------------------------------------------------------------------------------------

# Function to read files into a list of data.frames
extract_data <- function (file_type = '.fcs',
                          sampling = F,
                          sample_size = 10) {
    #file_type <- as.factor(file_type)
    #setwd(path)
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
