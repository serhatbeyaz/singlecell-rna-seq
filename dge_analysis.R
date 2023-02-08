## Title: Single Cell RNA-Seq analysis
## Author: serhat.beyaz

## Load required packages ----
library("glmGamPoi")
library("SingleCellExperiment")
library("scater")
library("ggplot2")
library("cowplot")
library("tibble")
library("DEGreport")
library("DESeq2")
library("tidyr")
library("dplyr")


today = gsub("-", "", Sys.Date())

## GLOBAL VARIABLES TO DEFINE PATHS ----

# Path to data
parent_dir = "/home/serhat.beyaz/count_matrices_celltype"

# Path to export files and tables
output_dir = "/home/serhat.beyaz/DGE_results_sb"
dge_results  = "/DGE results"
QC_path = "/QC and exploratory figures"
tables = "/tables"



#' Create SingleCellExperiment (sce) object
#' 
#' @description This function serves for creating sce object from csv count matrix and metadata
#' 
#' @param [parent_path] Path for parent dictionary <- character
#'        [cell_type] name of the cell type <- character
#' 
#' @return S4 type SingleCellExperiment object
#' 
#' @author serhat.beyaz
#' 
create_sce = function(parent_path, cell_type){
  
  # Join parent_dir and cell_type to go inside the cell specific folder
  cell_path = paste(parent_path, cell_type, sep = "/")
  
  # Get full_path for count matrix and metadata file
  matrix_path = list.files(cell_path, pattern = "raw", full.names = TRUE)
  metadata_path = list.files(cell_path, pattern = "metadata", full.names = TRUE)
  
  # Read in the files
  count_matrix = read.csv(matrix_path, row.names = "index")
  metadata = data.frame(read.csv(metadata_path, row.names = "index", stringsAsFactors = TRUE))
  
  # Discard these samples from the analysis (failed during quality check)
  discard_idx = metadata$sample_id == "MUC29230" | metadata$sample_id == "MUC29233" 
  
  count_matrix = count_matrix[,!discard_idx]
  
  metadata = metadata[!discard_idx,]
  
  metadata = droplevels(metadata)
  
  # Create SingleCellExperiment object
  sce = SingleCellExperiment(assays = list(counts = as.matrix(count_matrix)),
                             colData = DataFrame(metadata))
  
  # Add log transformation of the counts 
  assay(sce, "logcounts") = log2(counts(sce) + 1)
  
  sce$Condition = gsub(" ", "_", sce$Condition, fixed = TRUE)
  
  sce$Condition = factor(sce$Condition)
  
  # Remove variables holding csv content from the memory
  rm(count_matrix, metadata)
  
  return(sce)
}


#' Create QC metric plots
#' 
#' @param [sce] SingleCellExperiment object
#' 
#' @return it returns the plot just for demonstration, but before returning it renders
#' the plot content into a pdf, and save it to output_dir
#' 
#' @author serhat.beyaz
QC_plotter = function(sce){
  
  cell_info <- as.data.frame(colData(sce))
  
  cell_type = unique(sce$cell.type.per.sample)
  
  p1 = ggplot(data = cell_info, aes(x = sample_id, y = log_counts))+
    geom_violin(fill = 'brown') + theme_bw() + 
    geom_jitter(position = position_jitter(seed = 1, width = 0.2), size = 0.2) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
  
  p2 = ggplot(data = cell_info, aes(x = Condition, y = n_counts))+
    geom_violin(fill = 'brown') + theme_bw() + 
    geom_jitter(position = position_jitter(seed = 1, width = 0.2), size = 0.2) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
    
  p3 = ggplot(data = cell_info, aes(x = Condition, y = n_genes))+
    geom_violin(fill = 'brown') + theme_bw() + 
    geom_jitter(position = position_jitter(seed = 1, width = 0.2), size = 0.2) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
  
  p4 = ggplot(data = cell_info, aes(x = Condition, y = mt_fraction))+
    geom_violin(fill = 'brown') + theme_bw() + 
    geom_jitter(position = position_jitter(seed = 1, width = 0.2), size = 0.2) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
  
  plots = plot_grid(p1,p2,p3,p4, labels = list("A", "B", "C", "D"))
  
  title = ggdraw() + draw_label(
    paste0("QC plots for: ", cell_type),
    x = 0,
    hjust = 0
  )+theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )
  
  whole_plot = plot_grid(title, plots, rel_heights = c(0.1, 1), ncol=1)
  
  fig_name = paste(today, "QCplots", cell_type, ".pdf", sep = "_")
  fig_path = paste0(output_dir, "/", cell_type, QC_path, "/", fig_name)
  
  pdf(fig_path)
  print(whole_plot)
  dev.off()
  
  return(whole_plot)
  
}


#' Function to find highly variable genes
#' 
#' @param [sce] SingleCellExperiment object
#'        [top_n] Total gene number to be found
#'        
#' @description This function finds highly variable genes and
#' subsequently clusters them and plots PCA, which is saved as PDF
#' 
#' @return NULL
#' 
#' @author serhat.beyaz

find_hvg = function(sce, top_n){

  cell_type = unique(sce$cell.type.per.sample)
  
  model_gene_var = modelGeneVar(sce)
  
  top1000 = getTopHVGs(stats = model_gene_var, n = top_n)
  
  top1000.qc = sce[top1000,]
  
  top_clust = quickCluster(top1000.qc, min.size = 30)
  
  top1000.qc = computeSumFactors(top1000.qc, clusters = top_clust)
  
  top1000.qc = logNormCounts(top1000.qc)
  
  top1000.qc = runPCA(top1000.qc)
  
  pca = plotPCA(top1000.qc, colour_by = "Condition")
  
  pca2 = plotPCA(top1000.qc, colour_by = "sample_id")
  
  plots = plot_grid(pca, pca2, align = "v")
  
  fig_name = paste(today, "PCA", cell_type, ".pdf", sep = "_")
  fig_path = paste0(output_dir, "/", cell_type, QC_path, "/", fig_name)
  
  pdf(fig_path)
  print(plots)
  dev.off()
  
}

#' Likelihood Ratio Test (LRT) for DGE analysis
#' 
#' @description Whole LRT pipeline optimized for scRNA-Seq data
#' Specifically designed to make comparison within cell types, between Conditions
#' 
#' @param [whole_sce] SingleCellExperiment Object
#' 
#' @return resulting dds object, but just for manual testing after the analysis
#' 
#' @author serhat.beyaz
#' 
LRT_test = function(whole_sce){
  
  # Get the cell type
  cell_type = unique(sce$cell.type.per.sample)
  
  # Create unique combinations between Conditions for subsequent DGE analysis
  levels = levels(sce$Condition)
  combs = combn(unique(sce$Condition), 2)
  
  # Start loop for each comparison
  for (i in 1:length(levels)){
    
    # Get the comparison name
    comparison = paste0(combs[,i], collapse = "_vs_")
    
    cat("DGE analysis started for: ", comparison , "\n" )
    
    # Find which condition should be excluded
    exclude = levels[which(levels %in% combs[,i] == FALSE)]
    
    # Exclude that condition by creating a copy of the whole S4 object
    # but not taking that condition
    sce = whole_sce[, whole_sce$Condition != exclude]
    
    # Drop levels 
    sce$Condition = droplevels(sce$Condition)
    
    # Create DESeq object
    dds = DESeqDataSet(sce, design = ~ Condition )
    
    # Use size factors that are calculated by Maren Buttner, which seems more convenient to use
    sizeFactors(dds) = dds$size_factors
    
    # Start DGE analysis. The parameters are optimized for scRNA-seq data
    # Details can be found in DESeq2 documentation
    dds = DESeq(dds, useT=TRUE, minmu=1e-6, minReplicatesForReplace=Inf, 
                fitType = "glmGamPoi", test = "LRT", reduced = ~ 1)
    
    # Get the result coefficients
    coefs = resultsNames(dds)

    # Get the results table
    res <- results(dds, name = coefs[2])

    # Order genes according to their p-values, from smaller to bigger
    resOrdered <- res[order(res$pvalue),]
    
    # Set 0 p-values to 2e-300. The reason why we have 0 p-values is that
    # Some p-values are smaller than the minimum representable number in R
    # thus converges them to 0
    resOrdered$padj[resOrdered$padj == 0] = 2e-300
    
    # Same for p-values, above was for padj
    resOrdered$pvalue[resOrdered$pvalue == 0] = (2e-300)/nrow(resOrdered)
    
    # Drop NAs, and also drop some genes that have extreme log2FoldChange due unbalanced distribution across Conditions
    sig_DGE = as.data.frame(resOrdered) %>% drop_na(log2FoldChange) %>% filter(abs(log2FoldChange) < 1e+07,
                                                                            padj < 0.05)
    
    # Create volcano plot
    volcano_plot = EnhancedVolcano::EnhancedVolcano(sig_DGE,
                                     lab = rownames(sig_DGE),
                                     x = 'log2FoldChange',
                                     y = 'padj',
                                     title = "DGE analysis LRT test with glmGamPoi",
                                     subtitle = paste0(cell_type, " ", coefs[2])
                                        )
    
    # Create tibble for a nice clustering plot
    res_LRT_tb <- resOrdered %>%
      data.frame() %>%
      rownames_to_column(var="gene") %>% 
      as_tibble()
    
    # Subset to return genes with padj < 0.05
    sigLRT_genes <- res_LRT_tb %>% 
      filter(padj < 0.05,
             abs(log2FoldChange) < 1e+07)
    
    
    clustering_sig_genes <- sigLRT_genes %>%
      arrange(padj) %>%
      head(n=1000)
    
    
    # Obtain rlog values for those significant genes
    cluster_rlog <- assay(sce, "logcounts")[clustering_sig_genes$gene, ]
    
    metadata = as.data.frame(colData(sce))
    
    cells = colnames(cluster_rlog)
    
    new_metadata = metadata[cells, ]
    
    # Create the clustering plot
    cluster_plot <- degPatterns(cluster_rlog, metadata = new_metadata, time = "sample_id", col=NULL)
    
    # Rest is for saving csv tables and plots
    
    file_name = paste0("DEG_list", "_", comparison, "_", coefs[2], ".csv")

    path = paste0(output_dir, "/", cell_type, dge_results)

    full_path = paste0(path, "/", file_name)

    write.csv(as.data.frame(sig_DGE), file= full_path)
    
    fig_name = paste(today, "Volcano", cell_type, comparison, ".pdf", sep = "_")
    fig_path = paste0(output_dir, "/", cell_type, dge_results, "/", fig_name)
    
    pdf(fig_path)
    print(volcano_plot)
    dev.off()
    
    fig_name = paste(today, "Gene_clusters", cell_type, comparison, ".pdf", sep = "_")
    fig_path = paste0(output_dir, "/", cell_type, dge_results, "/", fig_name)
    
    pdf(fig_path)
    print(cluster_plot)
    dev.off()
    
    
  }
  
  return(dds)
  
}



#' Wald Test for DGE analysis
#' 
#' @description Whole Wald test pipeline 
#' Specifically designed to make comparison within cell types, between Conditions
#' 
#' @param [whole_sce] SingleCellExperiment Object
#' 
#' @return resulting dds object, but just for manual testing after the analysis
#' 
#' @author serhat.beyaz
#' 
#' @seealso LRT_test(), most of the steps are same, except for DESeq function
#' so you can take look at the explanation above function  
#'
Wald_test = function(whole_sce){
  
  cell_type = unique(sce$cell.type.per.sample)
  
  levels = levels(sce$Condition)
  combs = combn(unique(sce$Condition), 2)
  for (i in 1:length(levels)){
    
    comparison = paste0(combs[,i], collapse = "_vs_")
    
    cat("DGE analysis started for: ", comparison , "\n" )
    
    exclude = levels[which(levels %in% combs[,i] == FALSE)]
    
    sce = whole_sce[, whole_sce$Condition != exclude]
    
    sce$Condition = droplevels(sce$Condition)
    
    dds = DESeqDataSet(sce, design = ~ Condition )
    
    sizeFactors(dds) = dds$size_factors
    
    dds = DESeq(dds, test = "Wald", useT=TRUE, minmu=1e-6, minReplicatesForReplace=Inf)
    
    coefs = resultsNames(dds)
    
    res <- results(dds, name = coefs[2], lfcThreshold = 0.58, altHypothesis = "greaterAbs")
    
    resOrdered <- res[order(res$pvalue),]
    
    resOrdered$padj[resOrdered$padj == 0] = 2e-300
    
    resOrdered$pvalue[resOrdered$pvalue == 0] = (2e-300)/nrow(resOrdered)
    
    sig_DGE = as.data.frame(resOrdered) %>% drop_na(log2FoldChange) %>% filter(abs(log2FoldChange) < 1e+07,
                                                                               padj < 0.05)
    
    volcano_plot = EnhancedVolcano::EnhancedVolcano(sig_DGE,
                                                    lab = rownames(sig_DGE),
                                                    x = 'log2FoldChange',
                                                    y = 'padj',
                                                    title = "DGE analysis Wald test",
                                                    subtitle = paste0(cell_type, " ", coefs[2])
    )
    
    
    res_LRT_tb <- resOrdered %>%
      data.frame() %>%
      rownames_to_column(var="gene") %>% 
      as_tibble()
    
    # Subset to return genes with padj < 0.05
    sigLRT_genes <- res_LRT_tb %>% 
      filter(padj < 0.05,
             abs(log2FoldChange) < 1e+07)
    
    
    clustering_sig_genes <- sigLRT_genes %>%
      arrange(padj) %>%
      head(n=1000)
    
    
    # Obtain rlog values for those significant genes
    cluster_rlog <- assay(sce, "logcounts")[clustering_sig_genes$gene, ]
    
    metadata = as.data.frame(colData(sce))
    
    cells = colnames(cluster_rlog)
    
    new_metadata = metadata[cells, ]
    
    cluster_plot <- degPatterns(cluster_rlog, metadata = new_metadata, time = "sample_id", col=NULL)
    
    file_name = paste0("DEG_list", "_", comparison, "_", coefs[2], ".csv")
    
    path = paste0(output_dir, "/", cell_type, dge_results, "/Wald")
    
    full_path = paste0(path, "/", file_name)
    
    write.csv(as.data.frame(sig_DGE), file= full_path)
    
    fig_name = paste(today, "Volcano", cell_type, comparison, ".pdf", sep = "_")
    fig_path = paste0(output_dir, "/", cell_type, dge_results, "/Wald/", fig_name)
    
    pdf(fig_path)
    print(volcano_plot)
    dev.off()
    
    fig_name = paste(today, "Gene_clusters", cell_type, comparison, ".pdf", sep = "_")
    fig_path = paste0(output_dir, "/", cell_type, dge_results, "/Wald/", fig_name)
    
    pdf(fig_path)
    print(cluster_plot)
    dev.off()
    
    
  }
  
  return(dds)
  
}



#' Function that establishes the main pipeline for DGE analysis
#' 
#' @description Main pipeline starting from calling count matrices/metadata to the DGE test results
#' 
#' @param [data_path] path to data, where count matrices and metadata is expected
#'        [test_type] Type of statistical test to be used for finding DEG
#'        
#' @author serhat.beyaz 
#' 
main_function = function(data_path, test_type = "Wald"){
  
  # Get the files in the data directory to iterate whole pipeline for each cell type
  files_list = list.files(parent_dir, recursive = FALSE, full.names = FALSE)
  
  tryCatch(
    {
      for (i in files_list){
        
        cat("DGE analysis for ", i, " is started!\n")
        sce = create_sce(parent_dir, i)
        
        cat("Plotting QC metrics\n")
        QC_plotter(sce)
        
        cat("Finding highly variable genes\n")
        find_hvg(sce)
        
        cat(test_type, "Test started!\n")
        if (test_type == "Wald"){
  
          Wald_test(sce)
          
        }else if (test_type == "LRT"){
          
          LRT_test(sce)
          
        }
        
        rm(sce)
        cat("DGE analysis for ", i, " is finished!\n")
        
      }
    }
  )
  
  
}


main_function(parent_dir, "Wald")

