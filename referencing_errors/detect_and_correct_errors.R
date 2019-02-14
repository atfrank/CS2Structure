library(nmR)
# Goto referencing error fold
location <- "/Users/atfrank/Documents/GitHub/CS2Structure/referencing_errors/"
setwd(location)

# If needed create fold to store corrected shift files
system("mkdir -p corrected_shifts")

# loop over RNAs and correct individual dataset (if necessary)
expnames <- c("resname","resid","nucleus","expCS","error")
prednames <- c("model","resid","resname","nucleus","predCS","id")  
rnas <- unlist(strsplit("28SR 1KKA 1L1W 1LC6 1LDZ 1NC0 1OW9 1PJY 1R7W 1SCL 1UUU 1XHP 1YSV 1Z2J 1ZC5 2FDT 2JWV 2K66 2KOC 2L1V 2L3E 2LBJ 2LBL 2LDL 2LDT 2LHP 2LI4 2LK3 2LP9 2LPA 2LPS 2LQZ 2LU0 2LUB 2LUN 2LV0 2M4W 2M5U 2M8K 2M12 2M21 2M22 2M24 2MEQ 2MFD 2MHI 2MIY 2MNC 2MXL 2N2O 2N2P 2N4L 2N6S 2N6T 2N6W 2N6X 2NCI 2QH2 2QH4 2RVO 2Y95 4A4S 4A4T 4A4U 5A17 5A18 5IEM 5KH8 5KMZ 5KQE 5LSN 5LWJ 5N5C 5UF3 5UZT 5V16 5V17 5WQ1 6EZ0"," "))
for (rna in rnas){
  # Detect errors
  # Read in data
  expcs <- read.table(sprintf("%s/%s/chemical_shifts.txt.gz", location, rna), col.names = expnames, stringsAsFactors = FALSE)
  expcs <- subset(expcs, expCS != 0.0)
  predcs_larmord <- read.table(sprintf("%s/%s/chemical_shifts.larmord.txt.gz", location, rna), col.names = prednames, stringsAsFactors = FALSE)

  # LARMORD
  reference <- detect_referencing_error(observed_chemical_shifts = expcs, computed_chemical_shifts = predcs_larmord,  verbose = FALSE)  
  accuracy <- accuracy_estimate(observed_chemical_shifts = expcs, computed_chemical_shifts = predcs_larmord, verbose = FALSE)  
  save(reference, file=sprintf("%s/%s/referencing_errors_larmord.RData", location, rna))
  save(accuracy, file=sprintf("%s/%s/accuracy_larmord.RData", location, rna))
  
  # Corrected Shifts
  # LARMORD
  errors <- correct_referencing_error(observed_chemical_shifts = expcs, computed_chemical_shifts = predcs_larmord, ratio_mean_sd = 5, threshold = -2, verbose = FALSE)  
  expcs_corrected <- errors$shifts
  errors <- errors$errors
  save(expcs_corrected, file=sprintf("%s/%s/chemical_shifts_corrected_larmord.RData", location, rna))
  save(errors, file=sprintf("%s/%s/detected_errors_larmord.RData", location, rna))
  write.table(expcs_corrected, file=sprintf("%s/%s/chemical_shifts_corrected_larmord.txt", location, rna) , quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(errors, file=sprintf("%s/detected_errors_larmord.txt", location, rna) , quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(expcs_corrected, file=sprintf("%s/corrected_shifts/observed_shifts_corrected_larmord_%s.txt", location, rna) , quote = FALSE, row.names = FALSE, col.names = FALSE)
  cat(sprintf("done correcting referencing errors %s\n",rna))
}


