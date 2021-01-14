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

  print("--------Passed Arguments---------")
  for(i in names(opt)){
    print(paste0(i,": ", opt[[i]]))
  }
  
  in_dir <- opt$in_dir
  out_dir <- opt$out_dir # TODO: ask if kyle creates this or I create this
  sample_info <- rjson::fromJSON(opt$sample_info)[[1]]
  factors <- rjson::fromJSON(opt$factors)
  design <- as.formula(paste0("~", opt$design))
  
  # TODO : data validation
  if(is.null(opt$in_dir)){
    stop("No `in_dir` argument passed. Call --help to see options.")
  }
  # Check if in_dir exists
  # Check if length of factors is same as length of sample_info rows
  # Check if all sample file counts exists
  
  df_sampleInfo <- data.frame(samples = names(sample_info), do.call(rbind,sample_info), row.names = NULL, stringsAsFactors = TRUE)
  colnames(df_sampleInfo)[-1] <- factors
  
  deseq <- last2Deseq(in_dir = in_dir, df_sampleInfo = df_sampleInfo, design = design)
  deseq <- DESeq2::DESeq(deseq)
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  write.csv(rowData(deseq), file.path(out_dir, "DESeq2 result.csv"))
  sink(file.path(out_dir, "DESeq2 result summary.txt"))
  for(r in resultsNames(deseq)[-1]){
    res <- DESeq2::results(deseq, name = r)
    
    cat(paste0("Significance Testing: ", r))
    summary(res)
    cat("\n")
    
    write.csv(res, file.path(out_dir, paste0("Test result - ",r,".csv")))
  }
  sink()
  
  message("Last counts to DESeq2 analysis done!")
}
