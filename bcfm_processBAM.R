# Must be in macaque/bcfmendel_raw

rawbam_bcfcalls <- read.table("bcfpileup_mcompetition.decap.vcf", comment.char = "", header = TRUE)
names(rawbam_bcfcalls) <- gsub("X\\.", "", names(rawbam_bcfcalls))
names(rawbam_bcfcalls) <- gsub("NHP.", "", names(rawbam_bcfcalls))
names(rawbam_bcfcalls) <- gsub("X", "", names(rawbam_bcfcalls))
names(rawbam_bcfcalls) <- gsub("A$", "", names(rawbam_bcfcalls))

rb_mvcf <- rawbam_bcfcalls[,c(1,2,10:ncol(rawbam_bcfcalls))]
rb_mvcf_gt <- rb_mvcf[,3:ncol(rb_mvcf)]

rb_mvcf <- rawbam_bcfcalls[,c(1,2,10:ncol(rawbam_bcfcalls))]
rb_mvcf_gt <- as.matrix(rb_mvcf[,3:ncol(rb_mvcf)])
rb_mvcf_AD <- as.data.frame(apply(rb_mvcf_gt, 2, function(col) { return(strsplit(col, ":") %>% sapply(., function(entry) { return(entry[length(entry)]) })) }))
rb_mvcf_ADF <- as.data.frame(apply(rb_mvcf_gt, 2, function(col) { return(strsplit(col, ":") %>% sapply(., function(entry) { return(entry[length(entry)-2]) })) }))
rb_mvcf_ADR <- as.data.frame(apply(rb_mvcf_gt, 2, function(col) { return(strsplit(col, ":") %>% sapply(., function(entry) { return(entry[length(entry)-1]) })) }))

rb_mvcf_homoref <- as.data.frame(apply(rb_mvcf_AD, 2, function(col) {
  splitcol <- sapply(col, strsplit, ",")
  single_allele <- sapply(splitcol, length) < 2                      # No ALT entry, e.g. AD = 23
  homo_allele <- as.numeric(sapply(splitcol, tail, 1)) == 0          # No read is alternate, e.g. AD = 23,0
  
  return(single_allele | homo_allele)
}))

rb_mvcf_homorefAD1 <- as.data.frame(apply(rb_mvcf_AD, 2, function(col) {
  splitcol <- sapply(col, strsplit, ",")
  single_allele <- sapply(splitcol, length) < 2                      # No ALT entry, e.g. AD = 23
  homo_allele <- as.numeric(sapply(splitcol, tail, 1)) < 2           # No more than 1 read is alternate, e.g. AD = 22,1
  
  return(single_allele | homo_allele)
}))

rb_mvcf_hetF <- as.data.frame(apply(rb_mvcf_ADF, 2, function(col) {
  splitcol <- sapply(col, strsplit, ",")
  single_allele <- sapply(splitcol, length) < 2                     # Must have an ALT entry
  het_allele <- as.numeric(sapply(splitcol, tail, 1)) > 0           # ALT entry must be > 0
  
  return(!single_allele & het_allele)
}))

rb_mvcf_hetR <- as.data.frame(apply(rb_mvcf_ADR, 2, function(col) {
  splitcol <- sapply(col, strsplit, ",")
  single_allele <- sapply(splitcol, length) < 2
  het_allele <- as.numeric(sapply(splitcol, tail, 1)) > 0
  
  return(!single_allele & het_allele)
}))

rb_mvcf_homoref_df <- rb_mvcf
rb_mvcf_homoref_df[,3:ncol(rb_mvcf_homoref_df)] <- rb_mvcf_homoref
rb_mvcf_homoref_df$upos <- paste(rb_mvcf_homoref_df$CHROM, rb_mvcf_homoref_df$POS)

rb_mvcf_homorefAD1_df <- rb_mvcf
rb_mvcf_homorefAD1_df[,3:ncol(rb_mvcf_homorefAD1_df)] <- rb_mvcf_homorefAD1
rb_mvcf_homorefAD1_df$upos <- paste(rb_mvcf_homorefAD1_df$CHROM, rb_mvcf_homorefAD1_df$POS)

rb_mvcf_hetF_df <- rb_mvcf
rb_mvcf_hetF_df[,3:ncol(rb_mvcf_hetF_df)] <- rb_mvcf_hetF
rb_mvcf_hetF_df$upos <- paste(rb_mvcf_hetF_df$CHROM, rb_mvcf_hetF_df$POS)

rb_mvcf_hetR_df <- rb_mvcf
rb_mvcf_hetR_df[,3:ncol(rb_mvcf_hetR_df)] <- rb_mvcf_hetR
rb_mvcf_hetR_df$upos <- paste(rb_mvcf_hetR_df$CHROM, rb_mvcf_hetR_df$POS)

momlookup <- as.character(trio_table$Mother)
names(momlookup) <- trio_table$trio
dadlookup <- as.character(trio_table$Father)
names(dadlookup) <- trio_table$trio
kidlookup <- as.character(trio_table$Child)
names(kidlookup) <- trio_table$trio

get_bamparents_homobool <- function(mv_frame, trio) {
  rb_bool_table <- rb_mvcf_homoref_df[match(paste(mv_frame$CHROM, mv_frame$POS), rb_mvcf_homoref_df$upos),]
  rb_bool_matrix <- as.matrix(rb_bool_table[3:6])
  rownames(rb_bool_matrix) <- 1:nrow(rb_bool_matrix)
  
  include_bamparents_bool <-
    mapply(function(r,c) { return(rb_bool_matrix[r,c]) }, 1:nrow(rb_bool_matrix), momlookup[trio]) &
    mapply(function(r,c) { return(rb_bool_matrix[r,c]) }, 1:nrow(rb_bool_matrix), dadlookup[trio])
  
  return(include_bamparents_bool)
}

get_bamparents_homoboolAD1 <- function(mv_frame, trio) {
  rb_bool_table <- rb_mvcf_homorefAD1_df[match(paste(mv_frame$CHROM, mv_frame$POS), rb_mvcf_homorefAD1_df$upos),]
  rb_bool_matrix <- as.matrix(rb_bool_table[3:6])
  rownames(rb_bool_matrix) <- 1:nrow(rb_bool_matrix)
  
  include_bamparents_bool <-
    mapply(function(r,c) { return(rb_bool_matrix[r,c]) }, 1:nrow(rb_bool_matrix), momlookup[trio]) &
    mapply(function(r,c) { return(rb_bool_matrix[r,c]) }, 1:nrow(rb_bool_matrix), dadlookup[trio])
  
  return(include_bamparents_bool)
}

get_bamchild_hetbool <- function(mv_frame, trio) {
  rb_bool_tableF <- rb_mvcf_hetF_df[match(paste(mv_frame$CHROM, mv_frame$POS), rb_mvcf_hetF_df$upos),]
  rb_bool_tableR <- rb_mvcf_hetR_df[match(paste(mv_frame$CHROM, mv_frame$POS), rb_mvcf_hetR_df$upos),]
  
  rb_bool_matrixF <- as.matrix(rb_bool_tableF[3:6])
  rb_bool_matrixR <- as.matrix(rb_bool_tableR[3:6])
  
  rownames(rb_bool_matrixF) <- 1:nrow(rb_bool_matrixF)
  rownames(rb_bool_matrixR) <- 1:nrow(rb_bool_matrixR)
  
  include_bamchild_bool <-
    mapply(function(r,c) { return(rb_bool_matrixF[r,c] & rb_bool_matrixR[r,c]) }, 1:nrow(rb_bool_matrixF), kidlookup[trio])
  
  return(include_bamchild_bool)
}

filter_bam <- function(vcf_tables) {
  parentbools <- mapply(get_bamparents_homobool, vcf_tables, names(vcf_tables), SIMPLIFY = FALSE)
  childbools <- mapply(get_bamchild_hetbool, vcf_tables, names(vcf_tables), SIMPLIFY = FALSE)
  
  return(mapply(function(candidate_frame, pbool, cbool) { return(candidate_frame[(pbool & cbool),]) },
                vcf_tables, parentbools, childbools, SIMPLIFY = FALSE))
}

filter_bamAD1 <- function(vcf_tables) {
  parentbools <- mapply(get_bamparents_homoboolAD1, vcf_tables, names(vcf_tables), SIMPLIFY = FALSE)
  childbools <- mapply(get_bamchild_hetbool, vcf_tables, names(vcf_tables), SIMPLIFY = FALSE)
  
  return(mapply(function(candidate_frame, pbool, cbool) { return(candidate_frame[(pbool & cbool),]) },
                vcf_tables, parentbools, childbools, SIMPLIFY = FALSE))
}
