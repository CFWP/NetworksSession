###############################################################################
###############################################################################
## 
## Network Exercises
## On the basis of rags2ridges 2.2.1
## Course: Statistics for Omics: Networks
## Place:  Amsterdam, 27/06/2019
##
## Code:   Carel F.W. Peeters
##         Department of Epidemiology & Biostatistics
##         VU University medical center Amsterdam
##         Amsterdam, the Netherlands
##         cf.peeters@vumc.nl
## Date:   13/06/2019, Amsterdam, VUmc
##
###############################################################################
###############################################################################



#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Preliminaries**
#' **------------------------------------------------------------------------**

## Set working directory
setwd("")  ## Set to convenience

## Needed libraries
library(rags2ridges)



#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 1: Get acquainted with the data**
#' **------------------------------------------------------------------------**

## Data packaged as 3 data objects in ADdata object
## Will invoke basic functions for looking at data

## Invoke data, get to know objects
data(ADdata)
objects()

## Look at metabolic abundancies/expression
head(ADmetabolites)

## Look at Clinical (sample) information
head(sampleInfo)

## Look at Feature information
colnames(variableInfo)
table(variableInfo)




###############################################################################
###############################################################################
## 
## Section 1: Single Network
##
###############################################################################
###############################################################################

#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 2: Obtain Precision Matrix for AD Class 2 Data**
#' **------------------------------------------------------------------------**

## Extract data for AD Class 2
ADclass2 <- ADmetabolites[,sampleInfo$ApoEClass == "Class 2"]

## Check dimensions
dim(ADclass2)

## Transpose and scale data
ADclass2 <- scale(t(ADclass2))

## Use kCVauto to find optimal regularized precision matrix
OPT <- optPenalty.kCVauto(ADclass2, 
                          lambdaMin = 1e-07, 
                          lambdaMax = 20,
                          target = default.target(covML(ADclass2), type = "DUPV"))

## Have a look at optimal penalty
OPT$optLambda

## Have a look at optimal precision
head(OPT$optPrec)

## Heatmap of precision matrix
edgeHeat(OPT$optPrec, diag = FALSE, textsize = 3)



#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 3: Assess Conditioning**
#' **------------------------------------------------------------------------**

## Assess conditioning of optimal precision using CN plot
CNplot(covML(ADclass2), 
       lambdaMin = 1e-07, 
       lambdaMax = 20, 
       step = 5000,
       target = default.target(covML(ADclass2), type = "DUPV"),
       Iaids = TRUE,
       vertical = TRUE,
       value = OPT$optLambda)



#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 4: Extract Network**
#' **------------------------------------------------------------------------**

## Threshold the optimal precision matrix
## Retain those elements whose posterior probability of being present >= .999
P0 <- sparsify(OPT$optPrec,
               threshold = "localFDR",
               FDRcut = .999,
               verbose = FALSE)

## Heatmap of sparsified partial correlation matrix
edgeHeat(P0$sparseParCor, diag = FALSE, textsize = 3)



#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 5: Visualize Network**
#' **------------------------------------------------------------------------**

## Simple visualization
## Circular layout, no pruning
Ugraph(P0$sparseParCor, 
       type = "fancy", 
       Vsize = 3, 
       Vcex = .1)


## Pruning can give more clear idea of structure
Ugraph(P0$sparseParCor, 
       type = "fancy", 
       Vsize = 3, 
       Vcex = .1, 
       prune = TRUE)


## Type of layout is important
## Layout with FR algorithm
Ugraph(P0$sparseParCor, 
       type = "fancy", 
       lay = "layout_with_fr", 
       prune = TRUE, 
       Vsize = 7, 
       Vcex = .3)


## Colorings may add additional information
## Creating vector of colors (simple way) to indicate compound family
PcorP  <- pruneMatrix(P0$sparseParCor)
Colors <- rownames(PcorP)
Colors[grep("Amine", rownames(PcorP))]     <- "lightblue"
Colors[grep("Org.Acid", rownames(PcorP))]  <- "orange"
Colors[grep("Lip", rownames(PcorP))]       <- "yellow"
Colors[grep("Ox.Stress", rownames(PcorP))] <- "purple"

Ugraph(PcorP, 
       type = "fancy", 
       lay = "layout_with_fr", 
       Vcolor = Colors, 
       prune = TRUE, 
       Vsize = 7, 
       Vcex = .3)


## One can also add legends for additional information
Compounds <- c("Amines","Lipids", "Organic Acids", "Oxidative Stress")
Colors2   <- c("lightblue", "yellow", "orange", "purple")
legend(x=-1.3, y=-.8, Compounds, pch=21,
       col = "#777777", pt.bg = Colors2, 
       pt.cex = 2, cex = .6, bty = "n", ncol = 1)
legend(x=.7, y=-.8, 
       c("positive partial cor.", "negative partial cor."), 
       lwd = 2, lty = c(1,2),
       col = c("black", "black") , pt.cex = 2, 
       cex = .6, bty = "n", ncol = 1)


## Export to pdf for inspection
pdf(file = "ADnetworkClass2.pdf", width = 11, height = 11)
Ugraph(PcorP, 
       type = "fancy", 
       lay = "layout_with_fr", 
       Vcolor = Colors, 
       prune = TRUE, 
       Vsize = 7, 
       Vcex = .3)
dev.off()



#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 6: Simple Network Analysis: Node Metrics**
#' **------------------------------------------------------------------------**

## Obtain simple network statistics
## GGMnetworkStats(P0$sparseParCor) or
NwkSTATS <- GGMnetworkStats(PcorP)

## Look at degree
NwkSTATS$degree
DEGREE <- NwkSTATS$degree[order(NwkSTATS$degree, decreasing = TRUE)]
head(DEGREE)

## Look at Betweennes Centrality
NwkSTATS$betweenness
BETWEEN <- NwkSTATS$betweenness[order(NwkSTATS$betweenness, decreasing = TRUE)]
head(BETWEEN)



#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 7: Simple Network Analysis: Path Metrics**
#' **------------------------------------------------------------------------**

## Find 2 strongest paths between Amines 1 and 2
## Determine if these are mediating or moderating paths
Paths <- GGMpathStats(P0$sparsePrecision, node1 =  1, node2 = 2, 
                      nrPaths = 2, Vsize = 3, Vcex = .2, prune = TRUE)



#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 8: Simple Network Analysis: Community Search**
#' **------------------------------------------------------------------------**

## Find and visualize the communities for the extracted network
## Finding communities on the basis of betweenness-centrality
Commy <- Communities(PcorP, Vcolor = Colors, Vsize = 7,  Vcex = .3)

## Peeking at memberships
Commy$membership




###############################################################################
###############################################################################
## 
## Section 2: Multiple Networks
##
###############################################################################
###############################################################################

#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 9: Lists of Class Data**
#' **------------------------------------------------------------------------**

## Subset
ADclass1 <- ADmetabolites[, sampleInfo$ApoEClass == "Class 1"]
ADclass2 <- ADmetabolites[, sampleInfo$ApoEClass == "Class 2"]

## Transpose and scale data
ADclass1 <- scale(t(ADclass1))
ADclass2 <- scale(t(ADclass2))

## Correlations for subsets
rAD1 <- cor(ADclass1)
rAD2 <- cor(ADclass2)

## Constructing list of correlation matrices
Rlist = list(rAD1 = rAD1, rAD2 = rAD2)
samps = c(dim(ADclass1)[1], dim(ADclass2)[1])

## Constructing list of target matrices and data
Tlist <- default.target.fused(Slist = Rlist, ns = samps, type = "DUPV")
Ylist <- list(AD1data = ADclass1, AD2data = ADclass2)



#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 10: Class Precision Matrices**
#' **------------------------------------------------------------------------**

##########################
## If not enough time during session, we load the OPTf result from 'FusedOptimal.Rdata'
## Otherwise, we let the 'optPenalty.fused' function below run its course
#load("FusedOptimal.Rdata")
##########################

## Finding optimal penalty parameters and precision matrices
OPTf <- optPenalty.fused(Ylist = Ylist, 
                         Tlist = Tlist, 
                         cv.method = "LOOCV", 
                         verbose = TRUE)

## Have a look at optimal penalties
OPTf$lambda.unique

## Heatmaps of optimal precision matrices
edgeHeat(OPTf$Plist$AD1data, diag = FALSE, textsize = 3)
edgeHeat(OPTf$Plist$AD2data, diag = FALSE, textsize = 3)



#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 11: Extract Networks**
#' **------------------------------------------------------------------------**

## Threshold the optimal precision matrices
## Retain those elements whose posterior probability of being present >= .999
P0s <- sparsify.fused(OPTf$Plist, 
                      threshold = "localFDR", 
                      FDRcut = .999,
                      verbose = FALSE)

## Heatmaps of sparsified partial correlation matrices
edgeHeat(P0s$AD1data$sparseParCor, diag = FALSE, textsize = 3)
edgeHeat(P0s$AD2data$sparseParCor, diag = FALSE, textsize = 3)



#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 12: Visualize Networks**
#' **------------------------------------------------------------------------**

## Visualizing classes with same coordinates helps (visual) comparison
## First, need to retain union of features over class-networks
TST <- Union(P0s$AD1data$sparseParCor, P0s$AD2data$sparseParCor)

## Convenience, rename sparsified partial correlation matrices over union
PCclass1 <- TST$M1subset
PCclass2 <- TST$M2subset

## Set coloring
Colors <- rownames(PCclass2)
Colors[grep("Amine", rownames(PCclass2))]     <- "lightblue"
Colors[grep("Org.Acid", rownames(PCclass2))]  <- "orange"
Colors[grep("Lip", rownames(PCclass2))]       <- "yellow"
Colors[grep("Ox.Stress", rownames(PCclass2))] <- "purple"

## Visualize first network, retain coordinates
Coords <- Ugraph(PCclass2, 
                 type = "fancy", 
                 lay = "layout_with_fr", 
                 Vcolor = Colors, 
                 prune = FALSE, 
                 Vsize = 7, 
                 Vcex = .3)

## Visualize second network with retained coordinates
Ugraph(PCclass1, 
       type = "fancy",
       lay = NULL,
       coords = Coords,
       Vcolor = Colors, 
       prune = FALSE, 
       Vsize = 7, 
       Vcex = .3)

## Export
pdf(file = "ADnetworks.pdf", width = 11, height = 11)
Coords <- Ugraph(PCclass2, 
                 type = "fancy", 
                 lay = "layout_with_fr", 
                 Vcolor = Colors, 
                 prune = FALSE, 
                 Vsize = 7, 
                 Vcex = .3)
Ugraph(PCclass1,
       type = "fancy",
       lay = NULL,
       coords = Coords,
       Vcolor = Colors, 
       prune = FALSE, 
       Vsize = 7, 
       Vcex = .3)
dev.off()



#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 13: Simple Network Analysis: Node metrics**
#' **------------------------------------------------------------------------**

## Make list of sparse Partial correlation (or precision) matrices
PC0list = list(PCclass1 = PCclass1, 
               PCclass2 = PCclass2)

## GGMnetworks.fused
NwkSTATSList <- GGMnetworkStats.fused(PC0list)

## Quick look at various metrics
head(NwkSTATSList)

## Compare top degree centralities
DegreesAD1  <- data.frame(rownames(NwkSTATSList), NwkSTATSList$PCclass1.degree)
DegreesAD2  <- data.frame(rownames(NwkSTATSList), NwkSTATSList$PCclass2.degree)
DegreesAD1o <- DegreesAD1[order(DegreesAD1[,2], decreasing = TRUE),]
DegreesAD2o <- DegreesAD2[order(DegreesAD2[,2], decreasing = TRUE),]
head(DegreesAD1o, 7)
head(DegreesAD2o, 7)



#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 14: Simple Network Analysis: Communities**
#' **------------------------------------------------------------------------**

## Finding feature-communities on the basis of betweenness-centrality
CommyC1 <- Communities(PCclass1, Vcolor = Colors, Vsize = 7,  Vcex = .3)
CommyC2 <- Communities(PCclass2, Vcolor = Colors, Vsize = 7,  Vcex = .3)

## Peeking at memberships
CommyC1$membership
CommyC2$membership

## Exporting for convenience
pdf(file = "Communities.pdf", width = 11, height = 11)
CommyC1 <- Communities(PCclass1, 
                       Vcolor = Colors, 
                       Vsize = 7,  
                       Vcex = .3,
                       main = "Modules AD Class 1")
CommyC2 <- Communities(PCclass2, 
                       Vcolor = Colors, 
                       Vsize = 7, 
                       Vcex = .3,
                       main = "Modules AD Class 2")
dev.off()



#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 15: Simple Network Analysis: Testing**
#' **------------------------------------------------------------------------**

## Testing closely corresponds to assessing global network properties
## Testing degree distributions with dependent 2-group Wilcoxon Signed Rank Test
## Plot degree distributions
plot(density(DegreesAD1[,2]), col = "red", main = "Degree distributions")
lines(density(DegreesAD2[,2]), col = "blue")
legenda <- c("AD class 1", "AD class 2")
legend(10, .2, legend = legenda, lwd = rep(1,2), lty = rep(1,2), 
       col = c("red", "blue"), cex = .7)

## Test
wilcox.test(DegreesAD1[,2],DegreesAD2[,2], paired = TRUE, alternative = "less")


## Testing Entropy
## If sigaR is not installed, then uncomment and run the following:
#source("https://bioconductor.org/biocLite.R")
#biocLite("sigaR")
library(sigaR)
DATA  <- rbind(ADclass1, ADclass2)
ID    <- c(rep(1,dim(ADclass1)[1]),rep(0,dim(ADclass2)[1]))
Etest <- entropyTest(DATA, ID, 
                     nPerm = 20000, 
                     method = "knn",
                     lowCiThres = 0.1)
summary(Etest)



#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Additional Exercise: Simple Network Analysis: Differential Graphs**
#' **------------------------------------------------------------------------**

## Visualize the edges unique to each class
DiffGraph(PCclass1, PCclass2,
          lay = NULL,
          coords = Coords,
          Vcolor = Colors, 
          Vsize = 7, 
          Vcex = .3)

## Exporting
pdf(file = "DifferentialNetworks.pdf", width = 11, height = 11)
Coords <- Ugraph(PCclass2, 
                 type = "fancy", 
                 lay = "layout_with_fr", 
                 Vcolor = Colors, 
                 prune = FALSE, 
                 Vsize = 7, 
                 Vcex = .3,
                 main = "AD Class 2")
Ugraph(PCclass1,
       type = "fancy",
       lay = NULL,
       coords = Coords,
       Vcolor = Colors, 
       prune = FALSE, 
       Vsize = 7, 
       Vcex = .3,
       main = "AD Class 1")
DiffGraph(PCclass1, PCclass2,
          lay = NULL,
          coords = Coords,
          Vcolor = Colors, 
          Vsize = 7,
          Vcex = .3,
          main = "Differential Network")
dev.off()




###############################################################################
###############################################################################
## 
## Section 3: Hidden Gems
##
###############################################################################
###############################################################################

## Get yourself some unsolicited advice
rags2ridges:::.TwoCents()

## No alternative fact
rags2ridges:::.JayZScore()

## Thank you for flying rags2ridges
rags2ridges:::.rags2logo()
