## SeuratDevelop
# Description
Some extensibility function applied to the Seurat package
# Installation
You should install Seurat package first,and then install SeuratDevelop package.
~~~
devtools::install_github("ShengFang-0829/SeuratDevelop")
~~~
# Usage
Here, we use the pbmc3k dataset that comes with SeuratData.

# Prepare Data
~~~
library(SeuratData)
InstallData("pbmc3k")
data("pbmc3k")
rds <- pbmc3k.final
~~~

# Some Function
## subsetSeurat
~~~
library(SeuratDevelop)
# You can filter based on a column in Seurat metadata,meta parameter can choose colname in Seurat's metadata
# For example,you can screen out "Naive CD4 T" and "Memory CD4 T".
rds1 <- subsetSeurat(rds,meta="seurat_annotations",value=c("Naive CD4 T","Memory CD4 T"))

# You can alsoo filter based on a gene expression in Seurat assays.meta parameter can choose a gene symbol.
# You can use compair parameter,it can choose ">","<","=",">=","<=".
# When a gene is selected by meta, fill in the threshold of gene expression by using value parameter.
# When meta choose a gene,you can use slot parameter,it can choose "counts","data".
# For example,you can screen cells with CD3d expression greater than 0.

rds2 <- subsetSeurat(rds,meta="CD3d",value=0,compair=">")
~~~

## Seurat_diff
~~~
# You can select any combination based on any column of Seurat metadata to compare two cell populations and screen out
# For example,we select a column in metadata such as "Naive CD4 T" and "Naive CD4 T" in "seurat_annotations" and compare them.
data <- Seurat_diff(rds,ident.1="Naive CD4 T",ident.2="Memory CD4 T",meta="seurat_annotations")
head(data)
~~~