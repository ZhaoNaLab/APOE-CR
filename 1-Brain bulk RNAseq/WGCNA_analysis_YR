########
## Before opening R section, do the following to set /tmp so that anRichment has enough disk space to load

export SQLITE_TMPDIR=/tmp/
sample.infor <- read.delim("covar.txt")
dim(sample.infor)		
names(sample.infor)

datExpr0 = as.data.frame(sample.infor[, -1]);
rownames(datExpr0) <- sample.infor$ID
dim(datExpr0)


sample.infor <- datExpr0
#########
## Load Expression Values: need to make sure that the rows are samples, and columns are genes. 


exp <- read.delim("residual.txt")

dim(exp)

datExpr0 = as.data.frame(exp[, -c(1:8)]) 


datExpr0[1:5,1:4]


rownames(datExpr0) <-exp$ID

rownames(datExpr0)

exp<-datExpr0
dim(exp)

######
save(exp, sample.infor, file = "Input.RData")

## Load WGCNA library

library(WGCNA);
options(stringsAsFactors = FALSE);
allowWGCNAThreads()

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(exp, powerVector = powers, networkType="signed hybrid",verbose = 5)


sft

#[1] 3 for all datasets

> # Scale-free topology fit index as a function of the soft-thresholding power
  pdf("soft_thresholding_powers.pdf", width=9, height=5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h

abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# After viewing the plots, 3 is a good value for both datasets

net = blockwiseModules(exp, power = 3, networkType="signed hybrid", 
                       TOMType = "signed", minModuleSize = 60,
                       reassignThreshold = 0, mergeCutHeight = 0.4,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "TOM", 
                       verbose = 3)


table(net$colors) 
## min module size of 50 and merge of 0.3 yielded more than 30 modules
##20, 23 and 22 modules, many significant for genotype, group and gender



pdf("module_colors.pdf", width=12,height=9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "networkConstruction-auto.RData")

# Define numbers of genes and samples
nGenes = ncol(exp);
nSamples = nrow(exp);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(exp, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, sample.infor, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


pdf("ME_trait_correlation.pdf", width=18,height=10)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\t(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(sample.infor),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()




#### calculate MEs for each module
MEs = moduleEigengenes(exp, moduleColors)$eigengenes
write.table(MEs,"MEs.txt")

###to make the plot without any text in the modules

pdf("ME_trait_correlation_no_value.pdf", width=10,height=10)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(sample.infor),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

######
## Quantify Gene relationship to our trait of interest, ranking for group
# Define variable Genotype containing the Genotype column of sample.infor 
group = as.data.frame(sample.infor$group);
names(group) = "group"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(exp, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(exp, group, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(group), sep="");
names(GSPvalue) = paste("p.GS.", names(group), sep="");

#######
##  Use biomaRt to annotate Ensembl IDs with HUGO gene symbol and entrez ID. The HUGO gene symbol is for the PI to get the first impression of interesting genes, and the entrez ID is for using the GO analysis in WGCNA

library("biomaRt")
ensembl  <-  useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
bm <-getBM(attributes=c("chromosome_name","start_position","end_position","external_gene_name","entrezgene_id", "ensembl_gene_id"), mart=ensembl)
dim(bm)
#[1]53980      6
#write.csv(bm, file = "gene_anno.csv")
save(bm,file="bm.RData")


#match Ensembl IDs from the experiment to the bm file which contains #Ensembl ID, Entrez ID, Hugo gene symbol etc.

probes = names(exp)
probes2bm = match(probes, bm$ensembl_gene_id)
sum(is.na(probes2bm))

# Create the starting data frame
geneInfo0 = data.frame(Ensemble_ID = probes,
                       geneSymbol = bm$external_gene_name[probes2bm],
                       Chr = bm$chromosome_name[probes2bm],
                       Start= bm$start_position[probes2bm],
                       Stop= bm$end_position[probes2bm],
                       LocusLinkID = bm$entrezgene[probes2bm],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

# Order modules by their significance for Chr
modOrder = order(-abs(cor(MEs, group, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by #geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.group));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "geneInfo.csv")




## Annotate modules using anRichment
###match ensemble ID to entrez ID 

LocusLinkID = bm$entrezgene[probes2bm]

##Run enrichment analysis

library("anRichment")
GOcollection = buildGOcollection(organism = "mouse")


GOenrichment = enrichmentAnalysis( classLabels = moduleColors, identifiers = LocusLinkID, refCollection = GOcollection, useBackground = "given", threshold = 1e-4, thresholdType = c("Bonferroni", "FDR", "nominal"),nBestDataSets=10, getOverlapEntrez = FALSE, getOverlapSymbols = TRUE, maxReportedOverlapGenes=100, ignoreLabels = "grey", getFDR=TRUE, getBonferroniCorrection=TRUE);

collectGarbage();
write.csv(GOenrichment$enrichmentTable, file = "AnRichment-enrichmentTable.csv", row.names = FALSE);






