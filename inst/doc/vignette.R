## ----setup, include = FALSE,eval=FALSE-----------------------------------
#  knitr::opts_chunk$set(
#    collapse = TRUE,
#    comment = "#>"
#  )

## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()

## ----eval=FALSE,message=FALSE--------------------------------------------
#  library(MLML2R)
#  library(minfi)
#  library(GEOquery)

## ----eval=FALSE----------------------------------------------------------
#  getGEOSuppFiles("GSE63179")
#  untar("GSE63179/GSE63179_RAW.tar", exdir = "GSE63179/idat")
#  
#  list.files("GSE63179/idat", pattern = "idat")
#  files <- list.files("GSE63179/idat", pattern = "idat.gz$", full = TRUE)
#  sapply(files, gunzip, overwrite = TRUE)

## ----eval=FALSE----------------------------------------------------------
#  rgSet <- read.metharray.exp("GSE63179/idat")

## ----eval=FALSE,echo=FALSE-----------------------------------------------
#  rgSet <- read.metharray.exp("../data-raw/example1/GSE63179/idat")

## ----eval=FALSE----------------------------------------------------------
#  pData(rgSet)

## ----eval=FALSE----------------------------------------------------------
#  geoMat <- getGEO("GSE63179")
#  pD.all <- pData(geoMat[[1]])
#  pD <- pD.all[, c("title", "geo_accession", "characteristics_ch1.1",
#                   "characteristics_ch1.2","characteristics_ch1.3")]
#  pD

## ----eval=FALSE----------------------------------------------------------
#  sampleNames(rgSet) <- sapply(sampleNames(rgSet),function(x)
#    strsplit(x,"_")[[1]][1])
#  rownames(pD) <- pD$geo_accession
#  pD <- pD[sampleNames(rgSet),]
#  pData(rgSet) <- as(pD,"DataFrame")
#  rgSet

## ----eval=FALSE----------------------------------------------------------
#  MSet.noob<- preprocessNoob(rgSet)

## ----eval=FALSE----------------------------------------------------------
#  MethylatedBS <- getMeth(MSet.noob)[,c(1,3,5,6)]
#  UnMethylatedBS <- getUnmeth(MSet.noob)[,c(1,3,5,6)]
#  MethylatedOxBS <- getMeth(MSet.noob)[,c(7,8,2,4)]
#  UnMethylatedOxBS <- getUnmeth(MSet.noob)[,c(7,8,2,4)]

## ----eval=FALSE----------------------------------------------------------
#  results_exact <- MLML(T.matrix = MethylatedBS , U.matrix = UnMethylatedBS,
#                        L.matrix = UnMethylatedOxBS, M.matrix = MethylatedOxBS)

## ----eval=FALSE----------------------------------------------------------
#  results_em <- MLML(T.matrix = MethylatedBS , U.matrix = UnMethylatedBS,
#                     L.matrix = UnMethylatedOxBS, M.matrix = MethylatedOxBS,
#                     iterative = TRUE)

## ----echo=TRUE,eval=FALSE------------------------------------------------
#  all.equal(results_exact$hmC,results_em$hmC,scale=1)

## ----echo=FALSE,eval=FALSE-----------------------------------------------
#  load("../data-raw/example1/results_exact.rds")

## ----plot,echo=FALSE,eval=FALSE------------------------------------------
#  pdf(file="Real1_estimates.pdf",width=15,height=5)
#  par(mfrow =c(1,3))
#  densityPlot(results_exact$hmC,main= "5-hmC using PAVA",cex.axis=1.4,cex.main=1.5,cex.lab=1.4,xlab="Proportion")
#  densityPlot(results_exact$mC,main= "5-mC using PAVA",cex.axis=1.4,cex.main=1.5,cex.lab=1.4,xlab="Proportion")
#  densityPlot(results_exact$C,main= "5-C using PAVA",cex.axis=1.4,cex.main=1.5,cex.lab=1.4,xlab="Proportion")
#  dev.off()

## ----echo=FALSE,fig.width=15,fig.height=5,fig.cap="Estimated proportions of hydroxymethylation, methylation and unmethylation for the CpGs in the dataset using the MLML function with default options."----
knitr::include_graphics("Real1_estimates.pdf") 

## ----echo=FALSE,eval=FALSE-----------------------------------------------
#  load("../data-raw/example1/results_exact.rds")
#  load("../data-raw/Input.rds")

## ----echo=TRUE,eval=FALSE------------------------------------------------
#  set.seed(112017)
#  
#  index <- sample(1:dim(results_exact$mC)[1],1000,replace=FALSE) # 1000 CpGs
#  
#  Coverage <- round(MethylatedBS+UnMethylatedBS)[index,1:2] # considering 2 samples
#  
#  temp1 <- data.frame(n=as.vector(Coverage),
#                      p_m=c(results_exact$mC[index,1],results_exact$mC[index,1]),
#                      p_h=c(results_exact$hmC[index,1],results_exact$hmC[index,1]))
#  
#  MethylatedBS_temp <- c()
#  for (i in 1:dim(temp1)[1])
#  {
#    MethylatedBS_temp[i] <- rbinom(n=1, size=temp1$n[i], prob=(temp1$p_m[i]+temp1$p_h[i]))
#  }
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
#  true_parameters_sim2 <- data.frame(p_m=results_exact$mC[index,1],p_h=results_exact$hmC[index,1])
#  true_parameters_sim2$p_u <- 1-true_parameters_sim2$p_m-true_parameters_sim2$p_h

## ----eval=FALSE,echo=FALSE-----------------------------------------------
#  save(true_parameters_sim2,MethylatedBS_sim2,UnMethylatedBS_sim2,MethylatedOxBS_sim2,UnMethylatedOxBS_sim2,MethylatedTAB_sim2,UnMethylatedTAB_sim2,file="Data_sim2.rds")

## ----plot2,echo=FALSE,eval=FALSE-----------------------------------------
#  pdf(file="True_parameters.pdf",width=15,height=5)
#  par(mfrow =c(1,3))
#  plot(density(results_exact$hmC[index,1]),main= "True 5-hmC",xlab="Proportions",xlim=c(0,1),ylim=c(0,10),cex.axis=1.5,cex.main=1.5,cex.lab=1.5)
#  plot(density(results_exact$mC[index,1]),main= "True 5-mC",xlab="Proportions",xlim=c(0,1),ylim=c(0,10),cex.axis=1.5,cex.main=1.5,cex.lab=1.5)
#  plot(density(results_exact$C[index,1]),main= "True 5-C",ylim=c(0,10),xlab="Proportions",xlim=c(0,1),cex.axis=1.5,cex.main=1.5,cex.lab=1.5)
#  dev.off()

## ----echo=FALSE,fig.width=15,fig.height=5,fig.cap="True proportions of hydroxymethylation, methylation and unmethylation for the CpGs used to generate the datasets."----
knitr::include_graphics("True_parameters.pdf") 

## ----echo=FALSE,eval=TRUE------------------------------------------------
load("Data_sim2.rds")

## ------------------------------------------------------------------------
library(MLML2R)
 results_exactBO1 <- MLML(T.matrix = MethylatedBS_sim2 , U.matrix = UnMethylatedBS_sim2,
 L.matrix = UnMethylatedOxBS_sim2, M.matrix = MethylatedOxBS_sim2)

## ------------------------------------------------------------------------
 results_emBO1 <- MLML(T.matrix = MethylatedBS_sim2 , U.matrix = UnMethylatedBS_sim2,
 L.matrix = UnMethylatedOxBS_sim2, M.matrix = MethylatedOxBS_sim2,iterative=TRUE)

## ------------------------------------------------------------------------
 all.equal(results_emBO1$hmC,results_exactBO1$hmC,scale=1)

## ------------------------------------------------------------------------
 library(microbenchmark)
 mbmBO1 = microbenchmark(
    EXACT = MLML(T.matrix = MethylatedBS_sim2 , U.matrix = UnMethylatedBS_sim2,
                 L.matrix = UnMethylatedOxBS_sim2, M.matrix = MethylatedOxBS_sim2),
    EM =    MLML(T.matrix = MethylatedBS_sim2, U.matrix = UnMethylatedBS_sim2,
                 L.matrix = UnMethylatedOxBS_sim2, M.matrix = MethylatedOxBS_sim2,
                 iterative=TRUE),
    times=10)
 mbmBO1

## ------------------------------------------------------------------------
all.equal(true_parameters_sim2$p_h,results_exactBO1$hmC[,1],scale=1)

## ------------------------------------------------------------------------
all.equal(true_parameters_sim2$p_h,results_emBO1$hmC[,1],scale=1)

## ------------------------------------------------------------------------
results_exactBT1 <- MLML(T.matrix = MethylatedBS_sim2 , U.matrix = UnMethylatedBS_sim2,
G.matrix = UnMethylatedTAB_sim2, H.matrix = MethylatedTAB_sim2)

## ------------------------------------------------------------------------
 results_emBT1 <- MLML(T.matrix = MethylatedBS_sim2 , U.matrix = UnMethylatedBS_sim2,
 G.matrix = UnMethylatedTAB_sim2, H.matrix = MethylatedTAB_sim2,iterative=TRUE)

## ------------------------------------------------------------------------
 all.equal(results_emBT1$hmC,results_exactBT1$hmC,scale=1)

## ------------------------------------------------------------------------
 mbmBT1 = microbenchmark(
    EXACT = MLML(T.matrix = MethylatedBS_sim2, U.matrix = UnMethylatedBS_sim2,
                 G.matrix = UnMethylatedTAB_sim2, H.matrix = MethylatedTAB_sim2),
    EM =    MLML(T.matrix = MethylatedBS_sim2, U.matrix = UnMethylatedBS_sim2,
                 G.matrix = UnMethylatedTAB_sim2, H.matrix = MethylatedTAB_sim2,
                 iterative=TRUE),
    times=10)
 mbmBT1

## ------------------------------------------------------------------------
all.equal(true_parameters_sim2$p_h,results_exactBT1$hmC[,1],scale=1)

## ------------------------------------------------------------------------
all.equal(true_parameters_sim2$p_h,results_emBT1$hmC[,1],scale=1)

## ------------------------------------------------------------------------
 results_exactOT1 <- MLML(L.matrix = UnMethylatedOxBS_sim2, M.matrix = MethylatedOxBS_sim2,
 G.matrix = UnMethylatedTAB_sim2, H.matrix = MethylatedTAB_sim2)

## ------------------------------------------------------------------------
 results_emOT1 <- MLML(L.matrix = UnMethylatedOxBS_sim2, M.matrix = MethylatedOxBS_sim2,
 G.matrix = UnMethylatedTAB_sim2, H.matrix = MethylatedTAB_sim2,iterative=TRUE)

## ------------------------------------------------------------------------
 all.equal(results_emOT1$hmC,results_exactOT1$hmC,scale=1)

## ------------------------------------------------------------------------
 mbmOT1 = microbenchmark(
    EXACT = MLML(L.matrix = UnMethylatedOxBS_sim2, M.matrix = MethylatedOxBS_sim2,
                 G.matrix = UnMethylatedTAB_sim2, H.matrix = MethylatedTAB_sim2),
    EM =    MLML(L.matrix = UnMethylatedOxBS_sim2, M.matrix = MethylatedOxBS_sim2,
                 G.matrix = UnMethylatedTAB_sim2, H.matrix = MethylatedTAB_sim2,
                 iterative=TRUE),
    times=10)
 mbmOT1

## ------------------------------------------------------------------------
all.equal(true_parameters_sim2$p_h,results_exactOT1$hmC[,1],scale=1)

## ------------------------------------------------------------------------
all.equal(true_parameters_sim2$p_h,results_emOT1$hmC[,1],scale=1)

## ------------------------------------------------------------------------

results_exactBOT1 <- MLML(T.matrix = MethylatedBS_sim2 , U.matrix = UnMethylatedBS_sim2,
L.matrix = UnMethylatedOxBS_sim2, M.matrix = MethylatedOxBS_sim2,
G.matrix = UnMethylatedTAB_sim2, H.matrix = MethylatedTAB_sim2)

## ------------------------------------------------------------------------
 results_emBOT1 <- MLML(T.matrix = MethylatedBS_sim2 , U.matrix = UnMethylatedBS_sim2,
 L.matrix = UnMethylatedOxBS_sim2, M.matrix = MethylatedOxBS_sim2,
 G.matrix = UnMethylatedTAB_sim2, H.matrix = MethylatedTAB_sim2,iterative=TRUE)

## ------------------------------------------------------------------------
 all.equal(results_emBOT1$hmC,results_exactBOT1$hmC,scale=1)

## ------------------------------------------------------------------------
 mbmBOT1 = microbenchmark(
    EXACT = MLML(T.matrix = MethylatedBS_sim2, U.matrix = UnMethylatedBS_sim2,
                 L.matrix = UnMethylatedOxBS_sim2, M.matrix = MethylatedOxBS_sim2,
                 G.matrix = UnMethylatedTAB_sim2, H.matrix = MethylatedTAB_sim2),
    EM =    MLML(T.matrix = MethylatedBS_sim2, U.matrix = UnMethylatedBS_sim2,
                 L.matrix = UnMethylatedOxBS_sim2, M.matrix = MethylatedOxBS_sim2,
                 G.matrix = UnMethylatedTAB_sim2, H.matrix = MethylatedTAB_sim2,
                 iterative=TRUE),
    times=10)
 mbmBOT1

## ------------------------------------------------------------------------
all.equal(true_parameters_sim2$p_h,results_exactBOT1$hmC[,1],scale=1)

## ------------------------------------------------------------------------
all.equal(true_parameters_sim2$p_h,results_emBOT1$hmC[,1],scale=1)

