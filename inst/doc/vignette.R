## ----setup, include = FALSE,eval=FALSE-------------------------------------
#  knitr::opts_chunk$set(
#    collapse = TRUE,
#    comment = "#>"
#  )

## ----style, echo = FALSE, results = 'asis'---------------------------------
BiocStyle::markdown()

## ----packages1,eval=TRUE,message=FALSE,warning=FALSE-----------------------
library(MLML2R)
library(minfi)
library(GEOquery)
library(IlluminaHumanMethylation450kmanifest)

## ----getData1,eval=FALSE---------------------------------------------------
#  getGEOSuppFiles("GSE63179")
#  untar("GSE63179/GSE63179_RAW.tar", exdir = "GSE63179/idat")
#  
#  list.files("GSE63179/idat", pattern = "idat")
#  files <- list.files("GSE63179/idat", pattern = "idat.gz$", full = TRUE)
#  sapply(files, gunzip, overwrite = TRUE)

## ----readData1,eval=FALSE--------------------------------------------------
#  rgSet <- read.metharray.exp("GSE63179/idat")

## ----eval=FALSE------------------------------------------------------------
#  pData(rgSet)

## ----getPheno1,eval=FALSE--------------------------------------------------
#  if (!file.exists("GSE63179/GSE63179_series_matrix.txt.gz"))
#  download.file(
#  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63179/matrix/GSE63179_series_matrix.txt.gz",
#  "GSE63179/GSE63179_series_matrix.txt.gz")
#  
#  geoMat <- getGEO(filename="GSE63179/GSE63179_series_matrix.txt.gz",getGPL=FALSE)
#  pD.all <- pData(geoMat)
#  
#  #Another option
#  #geoMat <- getGEO("GSE63179")
#  #pD.all <- pData(geoMat[[1]])
#  
#  pD <- pD.all[, c("title", "geo_accession", "characteristics_ch1.1",
#                   "characteristics_ch1.2","characteristics_ch1.3")]
#  pD

## ----eval=FALSE------------------------------------------------------------
#  sampleNames(rgSet) <- sapply(sampleNames(rgSet),function(x)
#    strsplit(x,"_")[[1]][1])
#  rownames(pD) <- pD$geo_accession
#  pD <- pD[sampleNames(rgSet),]
#  pData(rgSet) <- as(pD,"DataFrame")
#  rgSet

## ----Preprocess1,eval=FALSE------------------------------------------------
#  MSet.noob<- preprocessNoob(rgSet)

## ----eval=FALSE------------------------------------------------------------
#  MethylatedBS <- getMeth(MSet.noob)[,c(1,3,5,6)]
#  UnMethylatedBS <- getUnmeth(MSet.noob)[,c(1,3,5,6)]
#  MethylatedOxBS <- getMeth(MSet.noob)[,c(7,8,2,4)]
#  UnMethylatedOxBS <- getUnmeth(MSet.noob)[,c(7,8,2,4)]

## ----MLML2Rexact1,eval=FALSE-----------------------------------------------
#  results_exact <- MLML(T.matrix = MethylatedBS , U.matrix = UnMethylatedBS,
#                        L.matrix = UnMethylatedOxBS, M.matrix = MethylatedOxBS)

## ----MLML2REM1,eval=FALSE--------------------------------------------------
#  results_em <- MLML(T.matrix = MethylatedBS , U.matrix = UnMethylatedBS,
#                     L.matrix = UnMethylatedOxBS, M.matrix = MethylatedOxBS,
#                     iterative = TRUE)

## ----echo=TRUE,eval=FALSE--------------------------------------------------
#  all.equal(results_exact$hmC,results_em$hmC,scale=1)

## ----plot,echo=FALSE,eval=FALSE--------------------------------------------
#  pdf(file="Real1_estimates.pdf",width=15,height=5)
#  par(mfrow =c(1,3))
#  densityPlot(results_exact$hmC,main= "5-hmC estimates - PAVA",cex.axis=1.4,cex.main=1.5,cex.lab=1.4,xlab="Proportion")
#  densityPlot(results_exact$mC,main= "5-mC estimates - PAVA",cex.axis=1.4,cex.main=1.5,cex.lab=1.4,xlab="Proportion")
#  densityPlot(results_exact$C,main= "uC estimates - PAVA",cex.axis=1.4,cex.main=1.5,cex.lab=1.4,xlab="Proportion")
#  dev.off()

## ----echo=FALSE,fig.width=15,fig.height=5,fig.cap="Estimated proportions of hydroxymethylation, methylation and unmethylation for the CpGs in the dataset using the MLML function with default options."----
knitr::include_graphics("Real1_estimates.pdf") 

## ----packages2,eval=TRUE,message=FALSE,warning=FALSE-----------------------
library(MLML2R)
library(minfi)
library(GEOquery)
library(IlluminaHumanMethylation450kmanifest)
library(wateRmelon)

## ----getData2,eval=FALSE---------------------------------------------------
#  getGEOSuppFiles("GSE71398")
#  untar("GSE71398/GSE71398_RAW.tar", exdir = "GSE71398/idat")
#  
#  list.files("GSE71398/idat", pattern = "idat")
#  files <- list.files("GSE71398/idat", pattern = "idat.gz$", full = TRUE)
#  sapply(files, gunzip, overwrite = TRUE)

## ----readData2,eval=FALSE--------------------------------------------------
#  rgSet <- read.metharray.exp("GSE71398/idat")

## ----eval=FALSE------------------------------------------------------------
#  pData(rgSet)

## ----getPheno2,eval=FALSE--------------------------------------------------
#  if (!file.exists("GSE71398/GSE71398_series_matrix.txt.gz"))
#  download.file(
#  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE71nnn/GSE71398/matrix/GSE71398_series_matrix.txt.gz",
#  "GSE71398/GSE71398_series_matrix.txt.gz")
#  
#  geoMat <- getGEO(filename="GSE71398/GSE71398_series_matrix.txt.gz",getGPL=FALSE)
#  pD.all <- pData(geoMat)
#  
#  #Another option
#  #geoMat <- getGEO("GSE71398")
#  #pD.all <- pData(geoMat[[1]])
#  
#  pD <- pD.all[, c("title", "geo_accession", "source_name_ch1",
#                   "tabchip or bschip:ch1","hypoxia status:ch1",
#                   "tumor name:ch1","batch:ch1","platform_id")]
#  pD$method <- pD$`tabchip or bschip:ch1`
#  pD$group <- pD$`hypoxia status:ch1`
#  pD$sample <- pD$`tumor name:ch1`
#  pD$batch <- pD$`batch:ch1`

## ----eval=FALSE------------------------------------------------------------
#  sampleNames(rgSet) <- sapply(sampleNames(rgSet),function(x)
#    strsplit(x,"_")[[1]][1])
#  rownames(pD) <- as.character(pD$geo_accession)
#  pD <- pD[sampleNames(rgSet),]
#  pData(rgSet) <- as(pD,"DataFrame")
#  rgSet

## ----preprocess2,eval=FALSE------------------------------------------------
#  ## Noob
#  MSet.noob<- preprocessNoob(rgSet)
#  
#  BSindex <- which(pD$method=="BSchip")
#  TABindex <- which(pD$method=="TABchip")
#  
#  ## BMIQ
#  anno <- getAnnotation(MSet.noob)
#  beta.b <- getBeta(MSet.noob, type = "Illumina")
#  design.v <- as.vector(anno$Type)
#  design.v[design.v == "I"] = 1
#  design.v[design.v == "II"] = 2
#  design.v <- as.numeric(design.v)
#  coln = colnames(beta.b)
#  
#  beta.noob.bmiq <- BMIQ(beta.b, design.v = design.v,sampleID = 1:48)
#  
#  beta_BS <- beta.noob.bmiq$nbeta[,BSindex]
#  beta_TAB <- beta.noob.bmiq$nbeta[,TABindex]
#  
#  # Total Signal = methylated + unmethylated
#  TotalBS <- getMeth(MSet.noob[,BSindex]) + getUnmeth(MSet.noob[,BSindex])
#  TotalTAB <- getMeth(MSet.noob[,TABindex]) + getUnmeth(MSet.noob[,TABindex])
#  
#  MethylatedBS <- beta_BS*TotalBS
#  UnMethylatedBS <- (1-beta_BS)*TotalBS
#  MethylatedTAB <- beta_TAB*TotalTAB
#  UnMethylatedTAB <- (1-beta_TAB)*TotalTAB

## ----eval=FALSE------------------------------------------------------------
#  results_exact <- MLML(T.matrix = MethylatedBS , U.matrix = UnMethylatedBS,
#                        G.matrix = UnMethylatedTAB, H.matrix = MethylatedTAB)

## ----MLML2REM2,eval=FALSE--------------------------------------------------
#  results_em <- MLML(T.matrix = MethylatedBS , U.matrix = UnMethylatedBS,
#                     G.matrix = UnMethylatedTAB, H.matrix = MethylatedTAB,
#                     iterative = TRUE)

## ----echo=TRUE,eval=FALSE--------------------------------------------------
#  all.equal(results_exact$hmC,results_em$hmC,scale=1)

## ----echo=TRUE,eval=FALSE--------------------------------------------------
#  all.equal(results_exact$mC,results_em$mC,scale=1)

## ----plot3,echo=FALSE,eval=FALSE-------------------------------------------
#  pdf(file="Real2a_estimates.pdf",width=15,height=5)
#  par(mfrow =c(1,3))
#  densityPlot(results_exact$hmC,main= "5-hmC estimates - PAVA",cex.axis=1.4,cex.main=1.5,cex.lab=1.4,xlab="Proportion",sampGroups=pD$group[BSindex])
#  densityPlot(results_exact$mC,main= "5-mC estimates - PAVA",cex.axis=1.4,cex.main=1.5,cex.lab=1.4,xlab="Proportion",sampGroups=pD$group[BSindex])
#  densityPlot(results_exact$C,main= "uC estimates - PAVA",cex.axis=1.4,cex.main=1.5,cex.lab=1.4,xlab="Proportion",sampGroups=pD$group[BSindex])
#  dev.off()

## ----echo=FALSE,fig.width=15,fig.height=5,fig.cap="Estimated proportions of hydroxymethylation, methylation and unmethylation for the CpGs in the dataset using the MLML function with default options."----
knitr::include_graphics("Real2a_estimates.pdf") 

## ----simulation,echo=TRUE,eval=FALSE---------------------------------------
#  set.seed(112017)
#  
#  index <- sample(1:dim(results_exact$mC)[1],1000,replace=FALSE) # 1000 CpGs
#  
#  Coverage <- round(MethylatedBS+UnMethylatedBS)[index,1:2] # considering 2 samples
#  
#  temp1 <- data.frame(n=as.vector(Coverage),
#                      p_m=c(results_exact$mC[index,1],
#                            results_exact$mC[index,1]),
#                      p_h=c(results_exact$hmC[index,1],
#                            results_exact$hmC[index,1]))
#  
#  MethylatedBS_temp <- c()
#  for (i in 1:dim(temp1)[1])
#  {
#    MethylatedBS_temp[i] <- rbinom(n=1, size=temp1$n[i],
#                                   prob=(temp1$p_m[i]+temp1$p_h[i]))
#  }
#  
#  
#  UnMethylatedBS_sim2 <- matrix(Coverage - MethylatedBS_temp,ncol=2)
#  MethylatedBS_sim2 <- matrix(MethylatedBS_temp,ncol=2)
#  
#  
#  MethylatedOxBS_temp <- c()
#  for (i in 1:dim(temp1)[1])
#  {
#    MethylatedOxBS_temp[i] <- rbinom(n=1, size=temp1$n[i], prob=temp1$p_m[i])
#  }
#  
#  UnMethylatedOxBS_sim2 <- matrix(Coverage - MethylatedOxBS_temp,ncol=2)
#  MethylatedOxBS_sim2 <- matrix(MethylatedOxBS_temp,ncol=2)
#  
#  
#  MethylatedTAB_temp <- c()
#  for (i in 1:dim(temp1)[1])
#  {
#    MethylatedTAB_temp[i] <- rbinom(n=1, size=temp1$n[i], prob=temp1$p_h[i])
#  }
#  
#  
#  UnMethylatedTAB_sim2 <- matrix(Coverage - MethylatedTAB_temp,ncol=2)
#  MethylatedTAB_sim2 <- matrix(MethylatedTAB_temp,ncol=2)
#  
#  true_parameters_sim2 <- data.frame(p_m=results_exact$mC[index,1],
#                                     p_h=results_exact$hmC[index,1])
#  true_parameters_sim2$p_u <- 1-true_parameters_sim2$p_m-true_parameters_sim2$p_h

## ----eval=FALSE,echo=FALSE-------------------------------------------------
#  save(true_parameters_sim2,MethylatedBS_sim2,UnMethylatedBS_sim2,MethylatedOxBS_sim2,UnMethylatedOxBS_sim2,MethylatedTAB_sim2,UnMethylatedTAB_sim2,file="Data_sim2.rds")

## ----plot2,echo=FALSE,eval=FALSE-------------------------------------------
#  pdf(file="True_parameters.pdf",width=15,height=5)
#  par(mfrow =c(1,3))
#  plot(density(results_exact$hmC[index,1]),main= "True 5-hmC",xlab="Proportions",xlim=c(0,1),ylim=c(0,10),cex.axis=1.5,cex.main=1.5,cex.lab=1.5)
#  plot(density(results_exact$mC[index,1]),main= "True 5-mC",xlab="Proportions",xlim=c(0,1),ylim=c(0,10),cex.axis=1.5,cex.main=1.5,cex.lab=1.5)
#  plot(density(results_exact$C[index,1]),main= "True uC",ylim=c(0,10),xlab="Proportions",xlim=c(0,1),cex.axis=1.5,cex.main=1.5,cex.lab=1.5)
#  dev.off()

## ----echo=FALSE,fig.width=15,fig.height=5,fig.cap="True proportions of hydroxymethylation, methylation and unmethylation for the CpGs used to generate the datasets."----
knitr::include_graphics("True_parameters.pdf") 

## ----echo=FALSE,eval=TRUE--------------------------------------------------
load("Data_sim2.rds")

## --------------------------------------------------------------------------
library(MLML2R)
 results_exactBO1 <- MLML(T.matrix = MethylatedBS_sim2 , 
                          U.matrix = UnMethylatedBS_sim2,
                          L.matrix = UnMethylatedOxBS_sim2, 
                          M.matrix = MethylatedOxBS_sim2)

## --------------------------------------------------------------------------
 results_emBO1 <- MLML(T.matrix = MethylatedBS_sim2, 
                       U.matrix = UnMethylatedBS_sim2,
                       L.matrix = UnMethylatedOxBS_sim2, 
                       M.matrix = MethylatedOxBS_sim2,
                       iterative=TRUE)

## --------------------------------------------------------------------------
 all.equal(results_emBO1$hmC,results_exactBO1$hmC,scale=1)

## --------------------------------------------------------------------------
 library(microbenchmark)
 mbmBO1 = microbenchmark(
    EXACT = MLML(T.matrix = MethylatedBS_sim2, 
                 U.matrix = UnMethylatedBS_sim2,
                 L.matrix = UnMethylatedOxBS_sim2, 
                 M.matrix = MethylatedOxBS_sim2),
    EM =    MLML(T.matrix = MethylatedBS_sim2, 
                 U.matrix = UnMethylatedBS_sim2,
                 L.matrix = UnMethylatedOxBS_sim2, 
                 M.matrix = MethylatedOxBS_sim2,
                 iterative=TRUE),
    times=10)
 mbmBO1

## --------------------------------------------------------------------------
all.equal(true_parameters_sim2$p_h,results_exactBO1$hmC[,1],scale=1)

## --------------------------------------------------------------------------
all.equal(true_parameters_sim2$p_h,results_emBO1$hmC[,1],scale=1)

## --------------------------------------------------------------------------
results_exactBT1 <- MLML(T.matrix = MethylatedBS_sim2, 
                         U.matrix = UnMethylatedBS_sim2,
                         G.matrix = UnMethylatedTAB_sim2, 
                         H.matrix = MethylatedTAB_sim2)

## --------------------------------------------------------------------------
 results_emBT1 <- MLML(T.matrix = MethylatedBS_sim2, 
                       U.matrix = UnMethylatedBS_sim2,
                       G.matrix = UnMethylatedTAB_sim2, 
                       H.matrix = MethylatedTAB_sim2,
                       iterative=TRUE)

## --------------------------------------------------------------------------
 all.equal(results_emBT1$hmC,results_exactBT1$hmC,scale=1)

## --------------------------------------------------------------------------
 mbmBT1 = microbenchmark(
    EXACT = MLML(T.matrix = MethylatedBS_sim2, 
                 U.matrix = UnMethylatedBS_sim2,
                 G.matrix = UnMethylatedTAB_sim2, 
                 H.matrix = MethylatedTAB_sim2),
    EM =    MLML(T.matrix = MethylatedBS_sim2, 
                 U.matrix = UnMethylatedBS_sim2,
                 G.matrix = UnMethylatedTAB_sim2, 
                 H.matrix = MethylatedTAB_sim2,
                 iterative=TRUE),
    times=10)
 mbmBT1

## --------------------------------------------------------------------------
all.equal(true_parameters_sim2$p_h,results_exactBT1$hmC[,1],scale=1)

## --------------------------------------------------------------------------
all.equal(true_parameters_sim2$p_h,results_emBT1$hmC[,1],scale=1)

## --------------------------------------------------------------------------
 results_exactOT1 <- MLML(L.matrix = UnMethylatedOxBS_sim2, 
                          M.matrix = MethylatedOxBS_sim2,
                          G.matrix = UnMethylatedTAB_sim2, 
                          H.matrix = MethylatedTAB_sim2)

## --------------------------------------------------------------------------
 results_emOT1 <- MLML(L.matrix = UnMethylatedOxBS_sim2, 
                       M.matrix = MethylatedOxBS_sim2,
                       G.matrix = UnMethylatedTAB_sim2, 
                       H.matrix = MethylatedTAB_sim2,
                       iterative=TRUE)

## --------------------------------------------------------------------------
 all.equal(results_emOT1$hmC,results_exactOT1$hmC,scale=1)

## --------------------------------------------------------------------------
 mbmOT1 = microbenchmark(
    EXACT = MLML(L.matrix = UnMethylatedOxBS_sim2, 
                 M.matrix = MethylatedOxBS_sim2,
                 G.matrix = UnMethylatedTAB_sim2, 
                 H.matrix = MethylatedTAB_sim2),
    EM =    MLML(L.matrix = UnMethylatedOxBS_sim2, 
                 M.matrix = MethylatedOxBS_sim2,
                 G.matrix = UnMethylatedTAB_sim2, 
                 H.matrix = MethylatedTAB_sim2,
                 iterative=TRUE),
    times=10)
 mbmOT1

## --------------------------------------------------------------------------
all.equal(true_parameters_sim2$p_h,results_exactOT1$hmC[,1],scale=1)

## --------------------------------------------------------------------------
all.equal(true_parameters_sim2$p_h,results_emOT1$hmC[,1],scale=1)

## --------------------------------------------------------------------------

results_exactBOT1 <- MLML(T.matrix = MethylatedBS_sim2, 
                          U.matrix = UnMethylatedBS_sim2,
                          L.matrix = UnMethylatedOxBS_sim2, 
                          M.matrix = MethylatedOxBS_sim2,
                          G.matrix = UnMethylatedTAB_sim2, 
                          H.matrix = MethylatedTAB_sim2)

## --------------------------------------------------------------------------
 results_emBOT1 <- MLML(T.matrix = MethylatedBS_sim2, 
                        U.matrix = UnMethylatedBS_sim2,
                        L.matrix = UnMethylatedOxBS_sim2, 
                        M.matrix = MethylatedOxBS_sim2,
                        G.matrix = UnMethylatedTAB_sim2, 
                        H.matrix = MethylatedTAB_sim2,iterative=TRUE)

## --------------------------------------------------------------------------
 all.equal(results_emBOT1$hmC,results_exactBOT1$hmC,scale=1)

## ----computationCost-------------------------------------------------------
 mbmBOT1 = microbenchmark(
    EXACT = MLML(T.matrix = MethylatedBS_sim2, 
                 U.matrix = UnMethylatedBS_sim2,
                 L.matrix = UnMethylatedOxBS_sim2, 
                 M.matrix = MethylatedOxBS_sim2,
                 G.matrix = UnMethylatedTAB_sim2, 
                 H.matrix = MethylatedTAB_sim2),
    EM =    MLML(T.matrix = MethylatedBS_sim2, 
                 U.matrix = UnMethylatedBS_sim2,
                 L.matrix = UnMethylatedOxBS_sim2, 
                 M.matrix = MethylatedOxBS_sim2,
                 G.matrix = UnMethylatedTAB_sim2, 
                 H.matrix = MethylatedTAB_sim2,
                 iterative=TRUE),
    times=10)
 mbmBOT1

## --------------------------------------------------------------------------
all.equal(true_parameters_sim2$p_h,results_exactBOT1$hmC[,1],scale=1)

## --------------------------------------------------------------------------
all.equal(true_parameters_sim2$p_h,results_emBOT1$hmC[,1],scale=1)

