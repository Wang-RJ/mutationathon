# $ #!/bin/bash
# $ bcftools annotate -x INFO/CSQ,INFO/InbreedingCoeff,INFO/BaseQRankSum,INFO/ExcessHet,INFO/ClippingRankSum,INFO/FS,INFO/MLEAC,INFO/MLEAF,INFO/MQRankSum,INFO/ReadPosRankSum,INFO/SOR /N/dcwan/projects/hahnlab-phi/macaque/vcf-files/1/snps.trio1.phaseByTransmission.alleleBalance.vcf.gz -R cpos_DP_GQ.txt > gatk_candidates_DPGQ.vcf
# $ bcftools annotate -x INFO/CSQ,INFO/InbreedingCoeff,INFO/BaseQRankSum,INFO/ExcessHet,INFO/ClippingRankSum,INFO/FS,INFO/MLEAC,INFO/MLEAF,INFO/MQRankSum,INFO/ReadPosRankSum,INFO/SOR /N/dcwan/projects/hahnlab-phi/macaque/vcf-files/14/snps.trio14.phaseByTransmission.alleleBalance.vcf.gz -R cpos_DP_GQ.txt > gatk_candidates_DPGQ.trio14.vcf
# $ grep -v '##' gatk_candidates_DPGQ.vcf > gatk_candidates_DPGQ.decap.vcf
# $ grep -v '##' gatk_candidates_DPGQ.trio14.vcf > gatk_candidates_DPGQ.trio14.decap.vcf
#
# $ bcftools view -Oz -s ^NHP-39238,NHP-39214,NHP-39237,NHP-39243,39317 gatk_candidates_DPGQ.vcf -o gatk_candidates_DPGQ_unphased_3-14.vcf.gz
# $ bcftools view -Oz -s NHP-39238,NHP-39214,NHP-39237,NHP-39243,39317 gatk_candidates_DPGQ.trio14.vcf -o gatk_candidates_DPGQ_unphased_1-2.vcf.gz
# $ bcftools merge -m none gatk_candidates_DPGQ_unphased_3-14.vcf.gz gatk_candidates_DPGQ_unphased_1-2.vcf.gz -o gatk_candidates_DPGQ_unphased.vcf
# $ grep -v '##' gatk_candidates_DPGQ_unphased.vcf > gatk_candidates_DPGQ_unphased.decap.vcf

# fix some of the sample names
gatkgt_vcf <- read.table("gatk_candidates_DPGQ.decap.vcf", header = TRUE, comment.char = "")
names(gatkgt_vcf) <- gsub("X\\.", "", names(gatkgt_vcf))
names(gatkgt_vcf) <- gsub("NHP.", "", names(gatkgt_vcf))
names(gatkgt_vcf) <- gsub("X", "", names(gatkgt_vcf))
names(gatkgt_vcf) <- gsub("A$", "", names(gatkgt_vcf))

# reduce to just the genotype
gatkgt_brief <- gatkgt_vcf[,c(1,2,4,5,10:length(gatkgt_vcf))]
gatkgt_brief[,5:8] <- sapply(gatkgt_brief[,5:8], function(x) { return(substr(x, 1, 3)) })
gatkgt_brief[,5:8] <- apply(gatkgt_brief[,5:8], MARGIN = c(1,2), function(x) { gsub("\\|", "/", x) })


# downstream_descendants pulls all related descendants of that trio (includes descendants of half-sibs to children of the trio, e.g. potential germline mosaicism)
downstream_descendants <- function(trio, trio_table) {
  trio_idx <- match(trio, trio_table$trio)
  trio_entry <- trio_table[trio_idx,]
  
  if(is.na(trio_entry)[4])
    return(NULL)
  
  else {
    descendant_trios <- unique(c(which(trio_entry$Father == trio_table$Father),
                                 which(trio_entry$Mother == trio_table$Mother),
                                 which(trio_entry$Child == trio_table$Mother),
                                 which(trio_entry$Child == trio_table$Father)))
    return(unique(
      c(as.numeric(trio_entry[3]), unlist(sapply(trio_table$trio[descendant_trios], downstream_descendants, trio_table[-trio_idx,])))))
  }
}

# make a list of downstream descendants and parents, which should both be excluded for each trio
dd_memo <- sapply(trio_table$trio, downstream_descendants, trio_table)
parents_memo <- split(trio_table[,1:2], trio_table$trio)
exclude_memo <- sapply(trio_table$trio, function(trio) { return(c(dd_memo[[trio]], as.numeric(parents_memo[[trio]]))) })

# Cross-reference entries in mv_autosomes to make sure candidates do not appear in unrelated individuals
# exclude parents and downstream descendants from unrelated individuals

# create a unique position ID to match up mvcf and mv_autosome dataframes
gatkgt_brief_upos <- paste(gatkgt_brief$CHROM, gatkgt_brief$POS)

# is there a nonref allele in an unrelated individual?
crossref_mv_exclude <- function(upos, trio) {
  entry <- gatkgt_brief[match(upos, gatkgt_brief_upos),][,-(1:4)]
  entry <- entry[!(names(entry) %in% exclude_memo[[trio]])]
  
  return("0/1" %in% entry | "1/1" %in% entry)
}

filter_rel <- function(vcf_tables) {
  rel_exclude_bools <- mapply(function(vcf_frame, trio) { return(sapply(vcf_frame, crossref_mv_exclude, trio)) },
                              lapply(vcf_tables, function(candidate_frame) { paste(candidate_frame$CHROM, candidate_frame$POS) }),
                              names(vcf_tables), SIMPLIFY = FALSE)
  
  return(mapply(function(candidate_frame, excl_bool) { return(candidate_frame[!excl_bool,]) },
                vcf_tables, rel_exclude_bools, SIMPLIFY = FALSE))
}

# Was an allele at a given position transmitted to a trio's grandchild?
transmission_check <- function(upos, trio) {
  trio_child <- trio_table$Child
  names(trio_child) <- trio_table$trio

  # Hacked for macaque competition with only 1 grandchild
  trio_grandchild <- "Hoegaarde"
  names(trio_grandchild) <- trio_table$trio
  
  grandchild <- trio_grandchild
  
  entry <- gatkgt_brief[match(upos, gatkgt_brief_upos),][,-(1:4)]
  
  # Check if grandchild has same genotype, only have to worry about heterozygotes since we are excluding non-reference in P generation
  return(entry[as.character(trio_child[trio])] == entry[as.character(trio_grandchild[trio])])
}

extract_transflag <- function(vcf_tables) {
  mapply(function(candidate_frame, trio) {
    return(mapply(transmission_check, paste(candidate_frame$CHROM, candidate_frame$POS), trio))
  },
  vcf_tables, names(vcf_tables))
}

