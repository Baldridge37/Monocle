# R script

```
setwd("C:/Documents/SingleCell/")

source("https://bioconductor.org/biocLite.R")

lib_loc="C:/Rpackages/"

biocLite("monocle")
library(monocle)

#install.packages("plyr")
#library(plyr)

#Cell numbers need to be in the same order in all files
ann<- read.table("C:/Documents/SingleCell/Combined_annotation2.txt", header=TRUE, sep="\t")
#rownames=sample names
counts<-read.table("C:/Documents/SingleCell/new_counts.txt", header = TRUE)
#rownames=genes header=sample names
feat<-read.table("C:/Documents/SingleCell/gene_feature.txt", header=TRUE, sep="\t")
#rownames=genes gene_short_name =genes #only have gene atm, add feature data

pd <- new("AnnotatedDataFrame", data = ann)
fd <- new("AnnotatedDataFrame", data = feat)
Tap1 <- newCellDataSet(as.matrix(counts), phenoData = pd, featureData = fd)

#It  is  often  convenient  to  know  how  many  express  a  particular  gene,  or  how  many  genes  are  expressed  by  a  given
#cell.  Monocle provides a simple function to compute those statistics:

Tap1 <- detectGenes(Tap1, min_expr = 0.1)
print(head(fData(Tap1)))

#Create a vector of genes that are expressed in >= 50 cells
expressed_genes <- row.names(subset(fData(Tap1), num_cells_expressed >= 50))

#filter cells based on pdata 
#valid_cells <- row.names(subset(pData(HSMM), Cells.in.Well == 1 & Control == FALSE & Clump == FALSE & Debris == FALSE & Mapped.Fragments > 1000000))
#HSMM <- HSMM[,valid_cells]
valid_cells <- row.names(subset(pData(Tap1), cell_number == 1 ))
Tap1 <- Tap1[,valid_cells]

# Log-transform each value in the expression matrix.
L <- log(exprs(Tap1[expressed_genes,])+1)
# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily"
#install.packages("plyr", lib_loc)
#library(plyr, lib.loc=lib_loc)

#using reshape2 instead of plyr
melted_dens_df <- melt(t(scale(t(L))))

qplot(value, geom="density", data=melted_dens_df) +  stat_function(fun = dnorm, size=0.5, color='red') +
  xlab("Standardized log(FPKM)") +
  ylab("Density")

#Differential expression
#marker_genes <- row.names(subset(fData(Tap1), gene_short_name %in% c("AT1G01020", "AT5G07230", "AT4G11130", "AT4G20420")))
#marker_genes <- intersect(marker_genes, expressed_genes)
                                
#diff_test_res <- differentialGeneTest(Tap1[marker_genes, ], fullModelFormulaStr="expression~plate")

#diff_test_res <- differentialGeneTest(Tap1, fullModelFormulaStr="expression~plate")

# Select genes that are significant at an FDR < 10%
#sig_genes <- subset(diff_test_res, qval < 0.1)
# Attach the HUGO symbols and other featureData for these genes#sig_genes <- merge(fData(HSMM), sig_genes, by="row.names")
#sig_genes[,c("gene_short_name", "pval", "qval")]

#Plot jitterplot for cells expressing RDR2 and A9 by plate
#RDR2_A91 <- Tap1[row.names(subset(fData(Tap1),gene_short_name %in% c("AT4G11130", "AT5G07230", "AT4G20420"))),]
#plot_genes_jitter(RDR2_A91, grouping="Cell_type", ncol=3)


#diff_test_res <- differentialGeneTest(Tap1[expressed_genes, ], reducedModelFormulaStr="~1")
#ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

Tap1<- setOrderingFilter(Tap1, expressed_genes)
Tap1@phenoData@data$Size_Factor<-1
Tap1 <- reduceDimension(Tap1, use_irlba=FALSE)

Tap1 <- orderCells(Tap1, reverse=TRUE)
plot_spanning_tree(Tap1, color_by = "plate")
plot_spanning_tree(Tap1, color_by = "Cell_type")
plot_spanning_tree(Tap1, color_by = "num_gene_expressed")
plot_spanning_tree(Tap1, color_by = "Pseudotime")
plot_spanning_tree(Tap1, color_by = "cell_number")

Tap1_filtered <- Tap1[expressed_genes, pData(Tap1)$State != 1]
my_genes <- row.names(subset(fData(Tap1_filtered), gene_short_name %in% c("AT5G07230")))
cds_subset <- Tap1_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by="plate")
```
