GATK_hom_call <- paste("callability/", list.files("callability/")[grep("homo_callability_tableGQ", list.files("callability/"))], sep = "")
GATK_het_call <- paste("callability/", list.files("callability/")[grep("het_callability_tableGQ", list.files("callability/"))], sep = "")

homc_list <- lapply(GATK_hom_call, function(hom_table) {
  GATK_homc <- read.table(hom_table, header = TRUE)
  return(GATK_homc[["bam_filtered"]] / GATK_homc[["denominator"]])
})
homc_mat <- matrix(unlist(homc_list), nrow = 4, byrow = TRUE)

hetc_list <- lapply(GATK_het_call, function(het_table) {
  GATK_hetc <- read.table(het_table, header = TRUE)
  return(GATK_hetc[["ab_filtered"]] / GATK_hetc[["denominator"]])
})
hetc_mat <- matrix(unlist(hetc_list), nrow = 4, byrow = TRUE)

callability <- hetc_mat * homc_mat
colnames(callability) <- paste("trio", 1, sep = "")
rownames(callability) <- paste("GQ", seq(20,50,10), sep = "")

pile_files <- list.files("pilecount/")[grep("chr", list.files("pilecount/"))]

hapsize_table <- sapply(pile_files, function(file) {
  ctable <- read.table(paste("pilecount/", file, sep = ""))
  return(ctable$V2[order(ctable$V1)])
})
