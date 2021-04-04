#!/usr/bin/env Rscript

..last2Deseq <- function(in_dir, df_sampleInfo, design, gene_isoform_mapping=NULL, countsFromAbundance="no"){
  library("DESeq2", quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE)
  library("tximport", quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE)
  
  samplenames <- df_sampleInfo[,1]
  file_counts <- file.path(in_dir, samplenames, paste0(samplenames, ".counts"))
  
  if(is.null(gene_isoform_mapping)){
    txi <- tximport(file_counts, type = "salmon", txIn = TRUE, txOut = TRUE)
  }else{
    txi <- tximport(file_counts, type = "salmon", txIn = TRUE, txOut = FALSE, tx2gene = gene_isoform_mapping, countsFromAbundance = countsFromAbundance)
  }
  
  deseq <- DESeq2::DESeqDataSetFromTximport(txi, colData = df_sampleInfo, design = design)
  
  return(deseq)
}


last2deseq <- function(in_dir, out_dir, sample_info, factors, design, gene_isoform_mapping = NULL){
  sample_info <- rjson::fromJSON(sample_info)
  factors <- rjson::fromJSON(factors)
  design <- as.formula(paste0("~", design))
  
  # ------- Data Validation ------- # 
  # Check if in_dir exists
  if(!dir.exists(in_dir)) { 
    stop("File directory `in_dir` does not exist.")
  }
  
  # Check if gene_isoform_mapping exists
  if(!is.null(gene_isoform_mapping) && !file.exists(gene_isoform_mapping)) { 
    stop("File `gene_isoform_mapping` does not exist.")
  }
  
  list_files <- list.files(in_dir, recursive = TRUE)
  n_factors <- length(factors)
  for(i in names(sample_info)){
    # Check if length of factors is same as length of sample_info rows
    if(length(sample_info[[i]]) != n_factors){ 
      stop(paste0("Sample info of `", i, "` does not have the same length as factors.")) 
    }
    
    # Check if file counts exists
    file_count <- paste0(i, "/", i, ".counts")
    if(!file_count %in% list_files){ 
      stop(paste0("File count `", file_count, "` does not exist in ", in_dir, ".")) 
    }
    
  }
  
  # Check if formula contains the factors declared
  for(f in all.vars(design)){
    if(!f %in% factors){
      stop(paste0("Design contains undeclared factor `", i, "`."))
    }
  }
  
  # ------- Sample Info Preparation ------- # 
  df_sampleInfo <- data.frame(samples = names(sample_info), do.call(rbind,sample_info), row.names = NULL, stringsAsFactors = TRUE)
  colnames(df_sampleInfo)[-1] <- factors
  for(i in 2:ncol(df_sampleInfo)){ 
    # Re-level first sample as reference group 
    df_sampleInfo[,i] <- relevel(df_sampleInfo[,i], ref = sample_info[[1]][[(i-1)]])
  }
  
  # ------- DESeq2 Analysis ------- # 
  deseq <- ..last2Deseq(in_dir = in_dir, df_sampleInfo = df_sampleInfo, design = design, gene_isoform_mapping = gene_isoform_mapping)
  deseq <- DESeq2::DESeq(deseq)
  
  
  # ------- Output ------- # 
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(deseq, file = file.path(out_dir, "DESeq2_fit.RDS"))
  sink(file.path(out_dir, "DE_count_summary.txt"))
  for(r in resultsNames(deseq)[-1]){
    res <- DESeq2::results(deseq, name = r)
    
    cat(paste0("Significance Testing: ", r))
    summary(res)
    cat("\n")
    
    write.csv(res, file.path(out_dir, paste0("Test_result-",r,".csv")))
  }
  sink()
  
  message("Last counts to DESeq2 analysis done!")
}

if(sys.nframe() == 0L) {  # mimics `__name__ == "main"` of python
  library("optparse", quietly = TRUE, warn.conflicts = FALSE)
  option_list = list(
    make_option("--in_dir",  help="directory of LAST counts",
                type="character", default=NULL),
    
    make_option("--out_dir",  help="",
                type="character", default="r output"),
    
    make_option("--sample_info",  help="json array of sample info. same format as the one in SnakeMake",
                type="character", default=NULL),
    
    make_option("--factors",  help="",
                type="character", default=NULL),
    
    make_option("--design",  help="formula which expresses how the counts depend on the variables in sample info",
                type="character", default=NULL),
    
    make_option("--gene_isoform_mapping", help="Tab delimited file with columns `Gene.stable.ID` and `Transcript.stable.ID`", 
                type="character", default=NULL)
  )
  
  opt <- parse_args(OptionParser(option_list=option_list))
  attach(opt)
  
  # print("--------Passed Arguments---------")
  # for(i in names(opt)){
  #   print(paste0(i,": ", opt[[i]]))
  # }
  
  last2deseq(in_dir = in_dir,
             out_dir = out_dir,
             sample_info = sample_info,
             factors = factors,
             design = design,
             gene_isoform_mapping = gene_isoform_mapping)
}