# user functions
source("functions.r")
packages <- c("optparse", "mice", "caret", "reticulate", "tensorflow", "keras", "bcv", "plyr")
suppressPackageStartupMessages(ipak(packages))

option_list = list(
    make_option(c("-i","--id"), type = "character", default = "test",
                help = "ID tag of the test RNA"),
    make_option(c("-s","--speed"), type = "character", default = "slow",
                help = "whether to perform fast or slow imputation"),
    make_option(c("--ss_table_0"), type = "character", default = "data/ss_table_0.txt",
                help = "training set after imputation"),
    make_option(c("--ss_table_1"), type = "character", default = "data/ss_table_1.txt",
                help = "training set after adding neighboring information"),
    make_option(c("-p","--program"), type = "character", default = "impute",
                help = "impute or cs2bps"),
    make_option(c("-o","--output"), type = "character", default = "output/",
                help = "name of output file"),
    make_option(c("-w","--whole_set_prediction"), type = "logical", default = "FALSE",
                help = "whether to output predictions from all classifiers")
)

parser = OptionParser(usage = "%prog [options] path_to_chemical_shift_file",
                      option_list = option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

if(length(arguments$args) != 1) {
  cat("Incorrect number of required positional arguments\n\n")
  print_help(parser)
  stop()
} else {
  cat("cs2bps -- Chemical Shift to Base Pairing Status predictions\n")
  cat("Author: Kexin Zhang\n")
  cat("Author: Aaron T. Frank\n")

  # get arguments
  cs_file_path = arguments$args[1]

  # get options
  speed = opt$speed
  output = opt$output
  id = opt$id
  program = opt$program
  whole_set = opt$whole_set_prediction
  ss_table_0 = opt$ss_table_0
  ss_table_1 = opt$ss_table_1
  fasta = opt$fasta

   # load and impute cs
  cat("Running CS2BPS for RNA:",id,"\n")
  cs = load_cs_data(cs_file_path, train=F)
  nuclei = cs$nuclei
  cs = cs$cs
  cat("Imputation...\n")
  cs_imp = impute_cs_data(cs = cs, id = id, speed = speed, ss_table = ss_table_0)

  # whether predict or impute chemical shifts
  # add time stamp
  currentDate <- Sys.Date()
  if(program == "impute"){
    cs_final = format_cs_file(cs_imp, nuclei)
    write.table(cs_final, file = paste0(output,currentDate, "_", id, "_impute.dat"), row.names = F, col.names = F, quote = F)
  } else {
    # add neighboring information and normalize based on train data
    test_add_neighbor_cs = add_neighboring_residues(cs_imp, rna=id, single=T, train=F, i=1)
    test_add_neighbor_cs_and_resnames = add_neighboring_resnames(test_add_neighbor_cs)
    test_normalized = normalize_test(test_add_neighbor_cs_and_resnames, ss_table_1 = ss_table_1)
    
    # load model and make predictions
    cat("-------------------------------------------------------------------------------------\n")
    cat("Load model and predict...\n")
    pred = load_and_predict(cs_imp['resid'], test_normalized, rna=id) # pred is ensemble of predictions
    
    if(!whole_set){
      pred = pred[,c(1,ncol(pred))]
    }
    write.table(pred, file = paste0(output,currentDate, "_", id, "_cs2bps.txt"), row.names = F, col.names = F, quote = F)
    cat("-------------------------------------------------------------------------------------\n")
    cat("CS2BPS prediction done!\n")
    cat("-------------------------------------------------------------------------------------\n")
    cat("Instruction on folding secondary structure using RNAstructure with CS2BPS as restraints:\n")
    cat("-------------------------------------------------------------------------------------\n")
    cat("Fold test-sequence.fasta test-structure.ct -sh test-cs2bps.txt\n\n")
    cat("OR\n\n")
    cat("partition test-sequence.fasta test-parition.pfs -sh test-cs2bps.txt\n")
    cat("MaxExpect test-partition.pfs test-structure.ct\n\n")
    cat("OR\n\n")
    cat("partition test-sequence.fasta test-parition.pfs -sh test-cs2bps.txt\n")
    cat("ProbKnot test-partition.pfs test-structure.ct\n")
    cat("-------------------------------------------------------------------------------------\n")
  }
}

