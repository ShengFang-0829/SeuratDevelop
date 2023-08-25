#' subsetSeurat
#' @param rds seurat object
#' @param meta choose meta sub(cluster,sample,group,genename .etc)
#' @param value You can select the type to be filtered in the meta.When a gene is selected by meta, fill in the threshold of gene expression here
#' @param slot when meta choose a gene,you can use slot parameter,it can choose "counts","data"
#' @param compair when meta choose a gene,you can use compair parameter,it can choose ">","<","=",">=","<="
#' @return Return a seurat object.
#' @export
#' @examples
subsetSeurat <- function(rds,
                         meta,
                         value = "all",
                         slot = "data",
                         compair = ">",
                         verbose = TRUE) {
  metadata <- rds@meta.data
  if (!"all" %in% value) {
    if (meta == "cluster") {
      tmp <- Idents(rds)
      barcode <- names(tmp[which(tmp %in% value)])
      rds <- rds[, barcode]
      levels(rds) <- value
      if (verbose) {
        print(paste0(
          "cluster list after subset: ",
          paste(levels(rds), collapse = ', ')
        ))
      }
    } else if (meta %in% colnames(metadata)) {
      barcode <- rownames(metadata[which(metadata[, meta] %in% value), ])
      rds <- rds[, barcode]
      rds@meta.data[, meta] <-
        factor(rds@meta.data[, meta], levels = value)
      if (verbose) {
        print(paste0(
          meta,
          " list after subset: ",
          paste(unique(rds@meta.data[, meta]), collapse = ', ')
        ))
      }
    } else if (meta %in% rownames(rds)) {
      metadata <- FetchData(rds,
                            vars = c("orig.ident", meta),
                            slot = slot)
      if (compair == ">") {
        barcode <- rownames(metadata[which(metadata[, meta] > value), ])
      } else if (compair == ">=") {
        barcode <- rownames(metadata[which(metadata[, meta] >= value), ])
      } else if (compair == "=") {
        barcode <- rownames(metadata[which(metadata[, meta] == value), ])
      } else if (compair == "<") {
        barcode <- rownames(metadata[which(metadata[, meta] < value), ])
      } else if (compair == "<=") {
        barcode <- rownames(metadata[which(metadata[, meta] <= value), ])
      }
      if (length(barcode) != 0) {
        rds <- rds[, barcode]
      } else {
        stop(
          paste0(
            "error: No cell satisfies this regulation ",
            meta,
            " gene ",
            compair,
            " ",
            value
          )
        )
      }
    } else {
      stop ("error: meta parameter is not in seurat object,please check meta parameter!")
    }

  } else {
    rds <- rds
  }
  return(rds)
}

#' subsetRDS
#' @param rds seurat object
#' @param group subset group,default:all
#' @param sample subset sample,default:all
#' @param cluster subset cluster,default:all
#' @return Return a seurat object.
#' @export
#' @examples
subsetRDS <-
  function(rds,
           group = "all",
           sample = "all",
           cluster = "all",
           verbose = TRUE) {
    if (!identical(group, 'all')) {
      rds <- subsetSeurat(rds,
                          meta = "group",
                          value = group,
                          verbose = FALSE)
    }
    if (!identical(sample, 'all')) {
      rds <- subsetSeurat(rds,
                          meta = "sample",
                          value = sample,
                          verbose = FALSE)
    }
    if (!identical(cluster, 'all')) {
      rds <- subsetSeurat(rds,
                          meta = "cluster",
                          value = cluster,
                          verbose = FALSE)
    }
    if (verbose) {
      print(paste0("cluster list after subset: ", paste(levels(rds), collapse =
                                                          ', ')))
      if ("sample" %in% rds@meta.data) {
        print(paste0("sample list after subset: ", paste(unique(
          rds@meta.data[, "sample"]
        ), collapse = ', ')))
      }
      if ("group" %in% rds@meta.data) {
        print(paste0("group list after subset: ", paste(unique(
          rds@meta.data[, "group"]
        ), collapse = ', ')))
      }

    }
    return(rds)
  }



#' Seurat_diff
#' @param rds seurat object
#' @param meta choose meta colname(cluster,sample,group,genename .etc)
#' @param ident.1 Identity class to define markers for
#' @param ident.2 A second identity class for comparison. If NULL (default) use all other cells for comparison.
#' @param rename.ident.1 Rename ident.1.If NULL(default),do not rename it.
#' @param rename.ident.2 Rename ident.2.If NULL(default),rename "Others".
#' @param assay of assay to fetch data for (default is RNA).
#' @param slot Slot to pull data from; note that if \code{test.use} is "negbinom", "poisson", or "DESeq2",
#' \code{slot} will be set to "counts"
#' @param logfc.threshold Limit testing to genes which show, on average, at least
#' X-fold difference (log-scale) between the two groups of cells. Default is 0.25
#' Increasing logfc.threshold speeds up the function, but can miss weaker signals.
#' @param test.use Denotes which test to use. Available options are:
#' \itemize{
#'  \item{"wilcox"} : Identifies differentially expressed genes between two
#'  groups of cells using a Wilcoxon Rank Sum test (default)
#'  \item{"bimod"} : Likelihood-ratio test for single cell gene expression,
#'  (McDavid et al., Bioinformatics, 2013)
#'  \item{"roc"} : Identifies 'markers' of gene expression using ROC analysis.
#'  For each gene, evaluates (using AUC) a classifier built on that gene alone,
#'  to classify between two groups of cells. An AUC value of 1 means that
#'  expression values for this gene alone can perfectly classify the two
#'  groupings (i.e. Each of the cells in cells.1 exhibit a higher level than
#'  each of the cells in cells.2). An AUC value of 0 also means there is perfect
#'  classification, but in the other direction. A value of 0.5 implies that
#'  the gene has no predictive power to classify the two groups. Returns a
#'  'predictive power' (abs(AUC-0.5) * 2) ranked matrix of putative differentially
#'  expressed genes.
#'  \item{"t"} : Identify differentially expressed genes between two groups of
#'  cells using the Student's t-test.
#'  \item{"negbinom"} : Identifies differentially expressed genes between two
#'   groups of cells using a negative binomial generalized linear model.
#'   Use only for UMI-based datasets
#'  \item{"poisson"} : Identifies differentially expressed genes between two
#'   groups of cells using a poisson generalized linear model.
#'   Use only for UMI-based datasets
#'  \item{"LR"} : Uses a logistic regression framework to determine differentially
#'  expressed genes. Constructs a logistic regression model predicting group
#'  membership based on each feature individually and compares this to a null
#'  model with a likelihood ratio test.
#'  \item{"MAST"} : Identifies differentially expressed genes between two groups
#'  of cells using a hurdle model tailored to scRNA-seq data. Utilizes the MAST
#'  package to run the DE testing.
#'  \item{"DESeq2"} : Identifies differentially expressed genes between two groups
#'  of cells based on a model using DESeq2 which uses a negative binomial
#'  distribution (Love et al, Genome Biology, 2014).This test does not support
#'  pre-filtering of genes based on average difference (or percent detection rate)
#'  between cell groups. However, genes may be pre-filtered based on their
#'  minimum detection rate (min.pct) across both cell groups. To use this method,
#'  please install DESeq2, using the instructions at
#'  https://bioconductor.org/packages/release/bioc/html/DESeq2.html
#' }
#' @param min.pct  only test genes that are detected in a minimum fraction of
#' min.pct cells in either of the two populations. Meant to speed up the function
#' by not testing genes that are very infrequently expressed. Default is 0.1.
#' @param pseudocount.use Pseudocount to add to averaged expression values when
#' calculating logFC. 1 by default.
#' @return Return a seurat object.
#' @export
#' @examples
Seurat_diff <-
  function(rds,
           meta,
           ident.1 = NULL,
           ident.2 = NULL,
           rename.ident.1 = NULL,
           rename.ident.2 = NULL,
           assay = "RNA",
           slot = "data",
           logfc.threshold = 0.25,
           test.use = "wilcox",
           min.pct = 0.1,
           pseudocount.use = 1) {
    print("group compair")
    if (!meta %in% colnames(rds@meta.data)) {
      stop(paste0(meta, " is not in seurat metadata, please check it!"))
    }
    if (is.null(rename.ident.1)){
      rename.ident.1 <- ident.1
    }
    if (is.null(rename.ident.2)){
      rename.ident.2 <- ident.2
    }
    # print(table(rds@meta.data[,meta]))
    before_time <- Sys.time()
    oldgroup1 <- ident.1
    newgroup1 <- rep(rename.ident.1, length(oldgroup1))
    if (!is.null(ident.2)) {
      oldgroup2 <- ident.2
      newgroup2 <- rep(rename.ident.2, length(oldgroup2))
    } else{
      others <-
        setdiff(unique(as.character(rds@meta.data[, meta])), oldgroup1)
      ident.2 <- others
      oldgroup2 <- others
      newgroup2 <- rep("Others", length(oldgroup2))
      rename.ident.2 <- "Others"
    }
    oldgroup <- c(oldgroup1, oldgroup2)
    newgroup <- c(newgroup1, newgroup2)
    print(paste0(
      rename.ident.1,
      "(",
      paste(ident.1, collapse = ","),
      ") vs ",
      rename.ident.2,
      "(",
      paste(ident.2, collapse = ","),
      ")"
    ))
    #change ident
    level_orign <- levels(factor(rds@meta.data[, meta]))
    df <- data.frame(group = rds@meta.data[, meta])
    df$group <- as.character(df$group)
    for (i in 1:length(oldgroup)) {
      level_orign[level_orign == oldgroup[i]] <- newgroup[i]
      df$group[df$group == oldgroup[i]] <- newgroup[i]
    }
    rds@meta.data[, meta] <-
      factor(df$group, levels = unique(level_orign))
    rds <- SetIdent(rds, value = rds@meta.data[, meta])
    rds <- subsetRDS(rds, cluster = c(rename.ident.1, rename.ident.2))
    print(table(rds@active.ident))
    diffdata <-
      FindMarkers(
        object = rds,
        assay = assay,
        slot = slot ,
        ident.1 = rename.ident.1,
        ident.2 = rename.ident.2,
        test.use = test.use,
        min.pct = min.pct,
        logfc.threshold = logfc.threshold,
        pseudocount.use = pseudocount.use,
        verbose = F
      )
    diffdata$ident.1 <- rename.ident.1
    diffdata$ident.2 <- rename.ident.2
    diffdata$gene_id <- rownames(diffdata)
    average <- AverageExpression(rds,  verbose = F)$RNA
    average <- as.data.frame(average)
    average <- average[, c(rename.ident.1, rename.ident.2)]
    colnames(average) <- c("ave.1", "ave.2")
    average$gene_id <- rownames(average)
    average <- average[rownames(diffdata), ]
    diffdata_final <- merge(diffdata, average, by = "gene_id")
    diffdata_final <-
      diffdata_final[order(diffdata_final$avg_log2FC, decreasing = T), ]
    after_time <- Sys.time()
    time_use <- difftime(after_time , before_time, units = "mins")
    print(paste0("group compair use :", round(time_use, 3), " mins"))
    return(diffdata_final)
  }
