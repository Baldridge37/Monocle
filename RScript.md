# R script for Monocle
```
setwd("C:/Documents/SingleCell/")

source("https://bioconductor.org/biocLite.R")

#lib_loc="C:/Rpackages/"

#biocLite("monocle")
library(monocle)

#install.packages("plyr")
#library(plyr)

#Cell numbers need to be in the same order in all files
ann<- read.table("C:/Documents/normalise_workflow/Combined_annotation2.txt", header=TRUE, sep="\t")
#rownames=sample names
counts<-read.table("C:/Documents/SingleCell/new_counts.txt", header = TRUE)
#rownames=genes header=sample names
feat<-read.table("C:/Documents/SingleCell/gene_feature.txt", header=TRUE, sep="\t")
#rownames=genes gene_short_name =genes #only have gene atm, add feature data
tpms<-read.table("protoplast_tpm.txt")

pd <- new("AnnotatedDataFrame", data = ann)
fd <- new("AnnotatedDataFrame", data = feat)
#Norm_tpms from the normalisation script
Tap1 <- newCellDataSet(as.matrix(tpms), phenoData = pd, featureData = fd)


#It  is  often  convenient  to  know  how  many  express  a  particular  gene,  
#or  how  many  genes  are  expressed  by  a  given cell.  
#Monocle provides a simple function to compute those statistics:

Tap1 <- detectGenes(Tap1, min_expr = 0.1)
print(head(fData(Tap1)))

#Create a vector of genes that are expressed in >= 50 cells
expressed_genes <- row.names(subset(fData(Tap1), num_cells_expressed >= 10))

#filter cells based on pdata 
#valid_cells <- row.names(subset(pData(HSMM), Cells.in.Well == 1 & Control == FALSE & Clump == FALSE & Debris == FALSE & Mapped.Fragments > 1000000))
#HSMM <- HSMM[,valid_cells]
#Can only take two arguments at a time
valid_cells1 <- row.names(subset(pData(Tap1),cell_number == 1))
valid_cells2 <- row.names(subset(pData(Tap1), feature_number >= 1500, read_count >= 250000)) 

valid_cells <- intersect(valid_cells1, valid_cells2)

Tap1 <- Tap1[,valid_cells]

## Log-transform each value in the expression matrix.
L <- log(exprs(Tap1[expressed_genes,])+1)
## Standardize each gene, so that they are all on the same scale,
## Then melt the data with plyr so we can plot it easily"
#install.packages("plyr", lib_loc)
#library(plyr, lib.loc=lib_loc)

#install.packages("reshape2")
library(reshape2)

#using reshape2 instead of plyr
melted_dens_df <- melt(t(scale(t(L))))

qplot(value, geom="density", data=melted_dens_df) +  stat_function(fun = dnorm, size=0.5, color='red') +
  xlab("Standardized log(FPKM)") +
  ylab("Density")


Tap1<- setOrderingFilter(Tap1, expressed_genes)
Tap1@phenoData@data$Size_Factor<-1
Tap1 <- reduceDimension(Tap1, use_irlba=FALSE)

Tap1 <- orderCells(Tap1, reverse=TRUE)
p1<-plot_spanning_tree(Tap1, color_by = "plate")
p2<-plot_spanning_tree(Tap1, color_by = "Cell_type")
p3<-plot_spanning_tree(Tap1, color_by = "num_gene_expressed")
p4<-plot_spanning_tree(Tap1, color_by = "Pseudotime")
p5<-plot_spanning_tree(Tap1, color_by = "cell_number")
p6<-plot_spanning_tree(Tap1, color_by = "read_count")
p7<-plot_spanning_tree(Tap1, color_by = "feature_number")
p8<-plot_spanning_tree(Tap1, color_by = "State")


Tap1_filtered <- Tap1[expressed_genes,]
my_genes <- row.names(subset(fData(Tap1_filtered), gene_short_name %in% c("AT4G27330","AT1G02050")))
cds_subset <- Tap1_filtered[my_genes,]
Pseduo_plot<-plot_genes_in_pseudotime(cds_subset, nrow=3,ncol=2, color_by="State")
```
