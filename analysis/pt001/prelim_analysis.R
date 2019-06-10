# set environment
source("https://raw.githubusercontent.com/msenosainAB/immunophenotyping/master/R/analysis_functions.R")
env_load()
# Read data (only to explore as df)
#pt001 <- extract_data()

#--------------------------------------------
# Get k for flowSOM
#--------------------------------------------
data_set <- dir(pattern = ".fcs")
cols_n <- c("TIGIT", "CD45-RA", "HLA-DR","CD226","CD8","CD3","CD56","LAG3","FoxP3","PD-1","Ki-67", "CD62L")
fSOM <- ReadInput(data_set, pattern = '.fcs', transform = TRUE,
                  toTransform = cols_n,
                  transformFunction = flowCore::arcsinhTransform(a=0, b=0.0067))

# BuildSOM
fSOM2 <- BuildSOM(fSOM, colsToUse = cols_n)

# BuildMST

fSOM3 <- BuildMST(fSOM2, tSNE=FALSE)
PlotStars(fSOM3)
PlotStars(fSOM3, view = 'grid')
PlotStars(fSOM3, view = "tSNE")
PlotMarker(fSOM3)

# MetaClustering

par(mar=c(6,6,4,4), lwd=2) # space of the plot
par(mgp=c(3.5,1,0), cex.lab=1.2) # distance of the axis and axis labels
DetermineNumberOfClusters(fSOM3$map$codes,max=80,'metaClustering_consensus',plot=F,smooth=0.2, seed=42)
#abline(v = 16, col = 'red')
# 16

#--------------------------------------------
# Run Cytofkit
#--------------------------------------------
cols_n <- c("TIGIT<TIGIT>", "CD45-RA<CD45-RA>", "HLA-DR<HLA-DR>","CD226<CD226>","CD8<CD8>","CD3<CD3>","CD56<CD56>","LAG3<LAG3>","FoxP3<FoxP3>","PD-1<PD-1>","Ki-67<Ki-67>", "CD62L<CD62L>")
cytofkit(fcsFiles = getwd(), markers = cols_n,
         projectName = "PT001_ck", ifCompensation = FALSE,
         transformMethod = "arcsinh",
         mergeMethod = "all", fixedNum = NULL,
         dimReductionMethod = "tsne",
         clusterMethods = "FlowSOM",
         visualizationMethods = "tsne",
         progressionMethod = "isomap",
         FlowSOM_k = 16, seed = 101,
         resultDir = getwd(), saveResults = TRUE,
         saveObject = TRUE, openShinyAPP = FALSE, a_a = 0, a_b = 0.0067, a_c=0)

#--------------------------------------------
# Run FlowSOM independently
#--------------------------------------------
data_set <- dir(pattern = ".fcs")
cols_n <- c("TIGIT", "CD45-RA", "HLA-DR","CD226","CD8","CD3","CD56","LAG3","FoxP3","PD-1","Ki-67", "CD62L")
fSOM <- FlowSOM(data_set,
                # Input options:
                compensate = F,transform = TRUE, toTransform=cols_n,
                transformFunction = flowCore::arcsinhTransform(a=0, b=0.0067), scale = TRUE,
                # SOM options:
                colsToUse = cols_n, xdim = 7, ydim = 7,
                # Metaclustering options:
                nClus = 16,
                # Seed for reproducible results:
                seed = 101)
PlotStars(fSOM$FlowSOM, backgroundValues = as.factor(fSOM$metaclustering))

# Plot by groups/samples
groupRes <- CountGroups(fSOM[[1]],
                        groups=data_set,
                        backgroundValues = as.factor(fSOM$metaclustering))

# Generate a dataframe with all info (sample ID, cluster ID)
# data
pt001_fS<- as.data.frame(fSOM$FlowSOM$data)
# sample id
pt001_fS['TP_ID'] <- 1
for (i in 1:length(data_set)){
    pt001_fS[fSOM$FlowSOM$metaData[[i]][1]:fSOM$FlowSOM$metaData[[i]][2], 'TP_ID'] <- i
}
# clusters id
pt001_fS['cluster_ID'] <- fSOM[[2]][fSOM[[1]]$map$mapping[,1]]

# generate summary table (median expression per cluster per marker)
flowSOM_median <- aggregate(. ~ cluster_ID, pt001_fS, median)
flowSOM_median[,'TP_ID'] <- NULL
flowSOM_median[,2:7] <- NULL
flowSOM_median[,'Time'] <- NULL
write.csv(flowSOM_median, file ='flowSOM_median.csv', row.names = F)

library(gplots)
heatmap.2(as.matrix(flowSOM_median[,2:13]), scale = 'column')

# tsne by cluster, by sample

tsne <- Rtsne::Rtsne(pt001_fS[,7:18], dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)

#--------------------------------------------
# Run Monocle
#--------------------------------------------

