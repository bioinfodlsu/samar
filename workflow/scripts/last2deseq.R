#!/usr/bin/env Rscript

last2Deseq <- function(in_dir, df_sampleInfo, design, df_protein2gene=NULL, countsFromAbundance="no"){
  library("DESeq2", quietly = TRUE, warn.conflicts = FALSE)
  library("tximport", quietly = TRUE, warn.conflicts = FALSE)
  
  samplenames <- df_sampleInfo[,1]
  filenames_counts <- file.path(in_dir, samplenames, paste0(samplenames, ".counts"))
  
  if(is.null(df_protein2gene)){
    txi <- tximport(filenames_counts, type = "salmon", txIn = TRUE, txOut = TRUE)
  }else{
    txi <- tximport(filenames_counts, type = "salmon", txIn = TRUE, txOut = FALSE, tx2gene = df_protein2gene, countsFromAbundance = countsFromAbundance)
  }
  
  deseq <- DESeq2::DESeqDataSetFromTximport(txi, colData = df_sampleInfo, design = design)
  
  return(deseq)
}

if(sys.nframe() == 0L) {  # mimics __name__ == "main" of python
  library("optparse", quietly = TRUE, warn.conflicts = FALSE)
  option_list = list(
    make_option(c("--in_dir"),  help="directory of LAST counts",
                type="character", default=NULL),

    make_option(c("--out_dir"),  help="",
                type="character", default="r output"),
    
    make_option(c("--sample_info"),  help="",
                type="character", default=NULL),
    
    make_option(c("--factors"),  help="",
                type="character", default=NULL),

    make_option(c("--design"),  help="formula which expresses how the counts depend on the variables in sample info",
                type="character", default=NULL)
  )

  opt_parser <- OptionParser(option_list=option_list)
  opt <- parse_args(opt_parser)

  # print("--------Passed Arguments---------")
  # for(i in names(opt)){
  #   print(paste0(i,": ", opt[[i]]))
  # }
  
  in_dir <- opt$in_dir
  out_dir <- opt$out_dir # TODO: ask if kyle creates this or I create this
  sample_info <- rjson::fromJSON(opt$sample_info)
  factors <- rjson::fromJSON(opt$factors)
  design <- as.formula(paste0("~", opt$design))
  
  # ------- Data Validation ------- # 
  # Check if in_dir exists
  if(!dir.exists(in_dir)) { 
    stop("File directory `in_dir` does not exist.")
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
  str_design <- strsplit(opt$design, "\\+|\\*")[[1]]
  str_design <- gsub(" ", "", str_design)
  for(i in str_design){
    if(!i %in% factors){
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
  deseq <- last2Deseq(in_dir = in_dir, df_sampleInfo = df_sampleInfo, design = design)
  deseq <- DESeq2::DESeq(deseq)
  
  
  # ------- Output ------- # 
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(deseq, file = file.path(out_dir, "DESeq2_fit.RDS"))
  write.csv(rowData(deseq), file.path(out_dir, "DESeq2_fit.csv"))
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
