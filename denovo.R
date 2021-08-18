# $ bcftools +mendelian -T triotable.txt -l x vcf/output.filtered.snps.removed.CalculateGenotypePosteriors.vcf.gz > MV_bcfmendelian.vcf
# [1]nOK        [2]nBad [3]nSkipped     [4]Trio
# 14485641        95796   74346   M,Noot,Heineken

# $ sample=$(head -n 1 triotable.txt | tail -n 1)
# $ bcftools view -Ou -s $sample MV_bcfmendelian.vcf | bcftools +mendelian -t $sample -l x > MVF_trio.vcf
# [1]nOK        [2]nBad [3]nSkipped     [4]Trio
# 0       95796   0       M,Noot,Heineken

# bcftools view -M2 -m2 -i'GT[0]="RR" && GT[1]="RR" && GT[2]="het"' MVF_trio.vcf | grep -v "##" > MVF_trio.decap.vcf

library(dplyr)
options(stringsAsFactors = FALSE)

source("bcfm_denovo_utilities.R")
## Imports set of functions for dealing with lists of dataframes, each a collection of de novo candidates from a different trio
## Querying functions
#   extract_genotype_entry <- function(genotypes_list, column)
#   extract_genotypes <- function(vcf_tables)
#   extract_depths <- function(vcf_tables)
#   extract_gq <- function(vcf_tables)
#   extract_allelic_balance <- function(vcf_tables)
#   extract_neighbor_distance <- function(vcf_tables) 
#   extract_neighbors <- function(vcf_tables, distance) 
#
## Filtering functions
#   filter_depth <- function(vcf_tables, lower_depth, upper_depth)
#   filter_gq <- function(vcf_tables, thresh)

trio_table <- read.table("trio_table.txt", header = TRUE)
trio_table$trio <- 1

candidate_table <- list()
i <- 1
candidate_table[[i]] <- read.table("MVF_trio.decap.vcf", header = TRUE, comment.char = "")
candidate_table[[i]]$INFO <- NULL
candidate_table[[i]]$ID <- NULL
candidate_table[[i]]$FILTER <- NULL
candidate_table[[i]]$FORMAT <- NULL
names(candidate_table[[i]]) <- gsub("X\\.", "", names(candidate_table[[i]]))
names(candidate_table)[[i]] <- 1

nrow(candidate_table[[i]])
# 13663

candidate_table <- lapply(candidate_table, function(df) { df[grepl("chr[0-9]", df$"CHROM"),]  })
nrow(candidate_table[[i]])
# 8870

# raw_genocalls <- extract_genotypes(candidate_table)
# momcalls <- raw_genocalls %>% lapply(., select, 1) %>% lapply(., table)
# momcalls       # All homozygous ref
# dadcalls <- raw_genocalls %>% lapply(., select, 2) %>% lapply(., table)
# dadcalls       # All homozygous ref
# kidcalls <- raw_genocalls %>% lapply(., select, 3) %>% lapply(., table)
# kidcalls       # All heterozygotes

# raw_allelicbalance <- extract_allelic_balance(candidate_table)
# raw_depths <- extract_depths(candidate_table)
# raw_gq <- extract_gq(candidate_table)

# hist(as.numeric(unlist(raw_gq[[i]])), main = paste("trio", i), xlab = "GQ score")
# hist(as.numeric(unlist(raw_depths[[i]])), main = paste("trio", i), xlab = "depth")
# hist(as.numeric(unlist(raw_depths[[i]]) %>% .[. < 80]), main = paste("trio", i), xlab = "depth < 80")
# hist(as.numeric(raw_allelicbalance[[i]]), breaks = 30, main = paste("trio", i), xlab = "allelic balance (child)")

# bcftools query -f '[%DP\t]\n' vcf/output.filtered.snps.removed.CalculateGenotypePosteriors.vcf.gz > depths.txt
# awk '{sum1+=$1; sum2+=$2; sum3+=$3; sum4+=$4} END {print sum1/NR, sum2/NR, sum3/NR, sum4/NR}' depths.txt
# ~60x for everyone except mom who is ~40x.
# 61.0755 58.5545 40.6727 61.1352
candidates_DP <- filter_depth(candidate_table, 20, 80)

# unlist(lapply(candidates_DP, nrow))
# 4020

candidates_DP <- lapply(candidates_DP, na.omit)
# unlist(lapply(candidates_DP, nrow))
# 3997
# candidates with depth 0 with genotypes called from readbackedphasing

candidates_DP_GQ <- filter_gq(candidates_DP, 20)
# unlist(lapply(candidates_DP_GQ, nrow))
# 1751
dpgq_allelicbalance <- extract_allelic_balance(candidates_DP_GQ)
dpgq_depths <- extract_depths(candidates_DP_GQ)
dpgq_gq <- extract_gq(candidates_DP_GQ)

# hist(unlist(dpgq_allelicbalance), breaks = 40, xlab = "allelic balance", main = "allelic balance (child)")
# hist(unlist(dpgq_depths), xlab = "depth", main = "filtered [20,80]")
# hist(unlist(dpgq_gq), xlab = "GQ", main = "filtered > 20")

# candidate_positions <- lapply(candidates_DP_GQ, select, 1:2) %>% do.call(rbind, .)
# candidate_positions <- candidate_positions[order(candidate_positions[,1], candidate_positions[,2]),]
# write.table(candidate_positions, file = "cpos_DP_GQ.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
#
# $ bcftools mpileup -Ou -A -f /N/dcwan/projects/hahnlab-phi/macaque/ref/ref-relabeled.fna -Q 10 -R cpos_DP_GQ.txt -a ADF,ADR,AD,DP bam_files/*.bam | bcftools call -m -Ov -o bcfpileup_mcompetition.vcf &
# $ grep -v '##' bcfpileup_mcompetition.vcf > bcfpileup_mcompetition.decap.vcf

source("bcfm_processBAM.R")
candidates_DP_GQ_BAM <- filter_bam(candidates_DP_GQ)
# candidates_DP_GQ_BAMAD1 <- filter_bamAD1(candidates_DP_GQ)

idx <- match(paste(candidates_DP_GQ_BAM[[1]][,1], candidates_DP_GQ_BAM[[1]][,2]), gatkgt_brief_upos)
cbind(gatkgt_brief[idx,], transmission = extract_transflag(candidates_DP_GQ_BAM))

denovo_candidates <- bind_candidates(candidates_DP_GQ_BAM %>% filter_allelicbalance(0.35, 1))
# After manual observation, chr3:76549728 appears to be a screwup in GATK's readbackphasing
# bam file shows transmission to grandchild
denovo_candidates[1,12] <- TRUE
names(denovo_candidates)[12] <- "transmit_flag"
names(denovo_candidates)[13] <- "trio"
