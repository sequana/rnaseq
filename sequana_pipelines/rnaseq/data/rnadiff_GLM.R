###################################################
### DESeq2_GLM parameters: to be modified by the user
###################################################
rm(list=ls())                                        # remove all the objects of the R session

workspace <- "."                                     # workspace for the R session

projectName <- "BXXXX"                               # name of the project (cannot contain any ".")
analysisVersion <- "vN"                              # name of the analysis version (cannot contain any ".")

author <- "FILLME (Biomics platform - Institut Pasteur)"   # author of the statistical report
researcher <- "FILLME"                               #  name of the researcher
chief <- ""                                          # name of the head of unit

varInt1 <- "varInt1"                                 # first factor of interest
varInt2 <- "varInt2"                                 # second factor of interest
condRef1 <- "condRef1"         			             # reference biological condition for varInt1
condRef2 <- "condRef2"                               # reference biological condition for varInt2
design <- ~ varInt1 + varInt2 + varInt1:varInt2      # design du mod?le statistique

outfile <- TRUE                                      # TRUE to export figures, FALSE to display them in R
colors <- c("#f3c300", "#875692", "#f38400", "#a1caf1", "#be0032", # vector of colors of each group on the plots
            "#c2b280", "#848482", "#008856", "#e68fac", "#0067a5", 
            "#f99379", "#604e97", "#f6a600", "#b3446c", "#dcd300", 
            "#882d17", "#8db600", "#654522", "#e25822", "#2b3d26")

cooksCutoff <- NULL                                  # outliers detection threshold (NULL to leave DESeq2 choosing it, Inf to keep outliers)
independentFiltering <- TRUE                         # FALSE to turn off the independent filtering (default is TRUE)
alpha <- 0.05                                        # threshold of statistical significance
adjMethod <- "BH"                                    # p-value adjustment method: "BH" (default) or "BY"

type.trans <- "VST"                                  # transformation for exploratory analysis: "VST" ou "rlog" (if size factors vary very widely)
locfunc <- "median"                                  # "median" (default) or "shorth" with library(genefilter) (to estimate the size factors)
interestingFeatures <- NULL                          # vector of features for which to plot the expression
featuresToRemove <- c("alignment_not_unique",        # names of the features to be removed (default is the HTSeq-count specific lines)
                      "ambiguous", "no_feature",
                      "not_aligned", "too_low_aQual") 

fitType <- "parametric"				                 # mean-variance relationship: "parametric" (default) or "local"

#####################################
# INPUT FILES
#####################################
geneLengthFile <- "input_gene_lengths.tsv"      # path to the genes lenghts file (default is NULL)
targetFile <- "target.txt"                      # path to the design/target file
infoFile <- "input_info.tsv"                    # path to the annotation file (needed if 0 counts not in counts files)
rawDir <- "feature_counts"                   # path to the directory containing raw counts files

###################################################
### code chunk number 1: construction autres parametres et divers chargements
###################################################
setwd(workspace)
library(RNADiff)
library(knitr)
if (locfunc=="shorth") library(genefilter)

versionName <- paste(projectName, analysisVersion, sep="-")
ncol <- NULL                                         # largeur des tableaux dans le rapport

cat("Creation des dossiers d'exports\n") 
dir.create("figures", showWarnings=FALSE)
dir.create("tables", showWarnings=FALSE)            

###################################################
### code chunk number 2: loadData
###################################################
cat("Chargement des annotations et longueurs des genes si besoin\n")
if (!is.null(infoFile)) print(head(info <- read.delim(infoFile, sep="\t", header=TRUE, stringsAsFactors=FALSE))) else info <- NULL
if (!is.null(geneLengthFile)) print(head(glength <- read.table(geneLengthFile, sep="\t", header=TRUE, stringsAsFactors=FALSE))) else glength <- NULL

cat("Chargement du target file\n")
print(target <- loadTargetFile(targetFile, varInt=c(varInt1,varInt2), condRef=c(condRef1,condRef2)))

cat("Chargement des donnees\n")
counts <- loadCountData(target, rawDir=rawDir, versionName=versionName, featuresToRemove=featuresToRemove)

cat("Verifier que les echantillons de counts sont dans le meme ordre que le target\n")
print(cbind(target=as.character(target[,1]),counts=colnames(counts)))

cat("Verifier que les identifiants dans info et glength sont les memes que dans les comptages\n")
checkInfoGlength(counts=counts, info=info, glength=glength)

####################################################
#### code chunk number 3: description of raw data
####################################################
cat("\nFigure : nombre de reads par echantillon\n")
barplotTC(counts=counts, group=target[,c(varInt1,varInt2)], col=colors, out=outfile, versionName=versionName)

cat("Figure : nombre de comptages nuls par echantillon\n")
barplotNul(counts=counts, group=target[,c(varInt1,varInt2)], col=colors, out=outfile, versionName=versionName)
N <- nrow(counts) - nrow(removeNul(counts))
cat("\nNombre de genes avec que des comptages nuls :", N,"\n")

cat("\nFigure : estimation de la densite des comptages de chaque echantillon\n")
densityPlot(counts=counts, group=target[,c(varInt1,varInt2)], col=colors, out=outfile, versionName=versionName)

cat("\nFigure + tableau : sequences majoritaires pour chaque echantillon\n")
majSequences <- majSequences(counts=counts, group=target[,c(varInt1,varInt2)], versionName=versionName, col=colors, out=outfile)

cat("\nCalcul des SERE\n")
print(sere <- pairwiseSERE(counts, versionName=versionName))

cat("\nFigure : pairwise scatterplots of samples\n")
pairwiseScatterPlots(counts=counts, group=target[,c(varInt1,varInt2)], out=outfile, versionName=versionName)

###################################################
### code chunk number 4: creating DESeqDataSet object, normalization and estimateDispersion
###################################################
dds <- DESeqDataSetFromMatrix(countData=counts, colData=target, design=design)
print(design(dds))

cat("Estimation des size factors\n")
dds <- estimateSizeFactors(dds, locfunc=eval(as.name(locfunc)))
print(sf <- sizeFactors(dds))
cat("\nFigure : diagnostic des size factors\n")
diagSizeFactors(dds=dds, group=target[,c(varInt1,varInt2)], col=colors, out=outfile, versionName=versionName)

cat("\nCalcul des dispersions et graph relation mean-dispersion\n")
dds <- estimateDispersions(dds, fitType=fitType)
plotDispEstimates(dds=dds, out=outfile, versionName=versionName)
cat("\nFigure : diagnostic de log-normalite des dispersions\n")
diagLogNormalityDisp(dds=dds, out=outfile, versionName=versionName)

####################################################
### code chunk number 5: Boxplot avant et apres normalisation
####################################################
cat("Figure : boxplots sur comptages bruts et normalises\n")
boxplotCounts(counts=counts(dds), group=target[,c(varInt1,varInt2)], col=colors, out=outfile, versionName=versionName)
boxplotCounts(counts=counts(dds, normalized=TRUE), group=target[,c(varInt1,varInt2)], col=colors, type="norm", out=outfile, versionName=versionName)

###################################################
### code chunk number 6: clustering + PCA of samples
###################################################
cat("Figure : dendrogramme de la classification sur comptages transformes\n")
if (type.trans == "VST") counts.trans <- assay(varianceStabilizingTransformation(dds))
if (type.trans == "rlog") counts.trans <- assay(rlogTransformation(dds))
clusterPlot(counts=counts.trans, out=outfile, versionName=versionName)

cat("Figure : premier plan de l'ACP sur les comptages transformes\n")
PCAPlot(dds=dds, group=target[,c(varInt1,varInt2)], col=colors, type.trans=type.trans, out=outfile, versionName=versionName)

###################################################
### code chunk number 7: analyse differentielle
###################################################
cat("Tests statistiques\n")
dds <- nbinomWaldTest(dds)

resultsNames(dds)
#  [1] "Intercept"                 "soucheSEG"                 "soucheB6"                 
#  [4] "infectionNI"               "infectionImoins"           "infectionIplus"           
#  [7] "soucheSEG.infectionNI"     "soucheB6.infectionNI"      "soucheSEG.infectionImoins"
# [10] "soucheB6.infectionImoins"  "soucheSEG.infectionIplus"  "soucheB6.infectionIplus" 

to_test <- list("B6-NI_vs_SEG-NI"=c(0,-1,1,0,0,0,-1,1,0,0,0,0),
                "B6-Imoins_vs_SEG-Imoins"=c(0,-1,1,0,0,0,0,0,-1,1,0,0),
                "(SEG-Iplus_vs_SEG-Imoins)_vs_(B6-Iplus_vs_B6-Imoins)"=c(0,0,0,0,0,0,0,0,-1,1,1,-1))

checkContrasts(coefs=resultsNames(dds),contrasts=to_test,versionName=versionName)
				
results <- vector("list",length(to_test)); names(results) <- names(to_test);
for (name in names(to_test)){
  results[[name]] <- results(dds, contrast=to_test[[name]], pAdjustMethod=adjMethod,
                             cooksCutoff=ifelse(!is.null(cooksCutoff), cooksCutoff, TRUE),
                             independentFiltering=independentFiltering, alpha=alpha)
}

###################################################
### code chunk number 8: results of the independent filtering
###################################################
if(independentFiltering){
  cat("Tableau : independent filtering\n")
  print(tabIndepFiltering <- tabIndepFiltering(results, versionName=versionName), quote=FALSE)
}

###################################################
### code chunk number 9: export tables
###################################################
cat("Export des resultats\n")
complete <- exportComplete.DESeq2(dds=dds, results=results, alpha=alpha, cooksCutoff=cooksCutoff,
                                  group=paste(target[,varInt1], target[,varInt2], sep="-"),
                                  conds=unique(paste(target[,varInt1], target[,varInt2], sep="-")),
                                  versionName=versionName, info=info, export=TRUE)

cat("# genes up, down et total par comparaison\n")
print(nDiffTotal <- nDiffTotal(complete, alpha=alpha, versionName=versionName), quote=FALSE)

cat("Figure : nb de genes DE selon seuil FDR\n")
nbDiffSeuil(complete=complete, out=outfile, versionName=versionName)

if (!is.null(geneLengthFile)){
  cat("Export : comptages normalises par la longueur des genes\n")
  normGeneLength(counts=counts(dds, normalized=TRUE), glength=glength, versionName=versionName)
  geneLengthEffect(counts, complete, glength, out=outfile, versionName=versionName)
}

###################################################
### code chunk number 10: distribution of raw p-values and MA-plot
###################################################
cat("Figure : distribution des log2(Fold-Changes)\n")
diagLogFC(complete=complete, out=outfile, versionName=versionName)

cat("Figure : histogramme des p-valeurs brutes\n")
histoRawp(complete=complete, out=outfile, versionName=versionName)

cat("\nFigure : MA-plot\n")
MAplotDE(complete=complete, pvalCutoff=alpha, out=outfile, versionName=versionName)

cat("\nFigure : volcano-plot\n")
volcanoPlotDE(complete=complete, pvalCutoff=alpha, out=outfile, versionName=versionName)

cat("\nFigure : Venn diagram\n")
vennDiagramDE(complete=complete, alpha=alpha, out=outfile, versionName=versionName)

cat("\nFigure : heatmap\n")
heatmapDE(counts.trans=counts.trans, complete=complete, alpha=alpha, out=outfile,
          key.xlab=paste0(type.trans, "-centered data"), versionName=versionName)

cat("\nFigure : interesting features\n")
if (!is.null(interestingFeatures)){
  plotEvolution(mat=log2(counts(dds, normalized=TRUE)+1), features=interestingFeatures,
                target=target, varInt1=varInt2, varInt2=varInt1, colors=colors,
                ylab=expression(log[2] ~ norm ~ counts + 1), out=outfile, versionName=versionName)
}

###################################################
### code chunk number 11: sessionInfo and saving
###################################################
cat("Sauvegarde des resultats\n")
sessionInfo <- sessionInfo()
pckVersionRNADiff <- packageVersion("RNADiff")
pckVersionDESeq2 <- packageVersion("DESeq2")
save.image(file=paste0(versionName, ".RData"))
# export RData for PF2heatmaps
results <- lapply(results, as.data.frame)
pf2heatmaps_objects <- c("varInt1", "varInt2", "target", "type.trans", "counts.trans", "results", "info")
save(list=pf2heatmaps_objects, file=paste0(versionName, "_PF2heatmaps.RData"), version=2)
# export RData for PF2toolsFilter
extract_col <- function(comp, info=NULL){
  if (is.null(info)){
    comp[, c("Id","baseMean", "log2FoldChange","padj")]
  } else{
    comp[, c(1:ncol(info), which(names(comp) %in%  c("baseMean", "log2FoldChange","padj")))]
  }
}
complete <- lapply(complete, extract_col, info=info)
save(complete, file=paste0(versionName, "_PF2toolsFilter.RData"), version=2)

###################################################
### code chunk number 12: knitr compilation
###################################################
if (!outfile){
  cat("outfile is FALSE: report and slides cannot be generated\n")
} else{
  cat("Creation du rapport et des slides\n")
  knit(system.file("reportGLM.Rnw", package="RNADiff"), paste0("report-", versionName, ".tex"), quiet=TRUE)
  knit(system.file("slidesGLM.Rnw", package="RNADiff"), paste0("slides-", versionName, ".tex"), quiet=TRUE)
  cat("Compilation du rapport\n")
  system(paste0("pdflatex report-", versionName, ".tex"))
  system(paste0("bibtex report-", versionName, ".aux"))
  system(paste0("pdflatex report-", versionName, ".tex"))
  system(paste0("pdflatex report-", versionName, ".tex"))
}
