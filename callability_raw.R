GATK_hom_call <- paste("callability/", list.files("callability/")[grep("homo_callabilityraw_tableGQ", list.files("callability/"))], sep = "")
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

auto_size <- rowSums(hapsize_table)
names(auto_size) <- paste("trio", 1:14, sep = "")

denovo_candidates <- list()
for(GQ in seq(20,80,10)) {
  denovo_candidates[[paste("GQ", GQ, sep = "")]] <- filter_gq(candidates_DP, GQ) %>% filter_bam %>% filter_rel %>% filter_allelicbalance(0.35, 1) %>% bind_candidates
}
GATK_ratebytrio <- lapply(names(denovo_candidates), function(dnf_idx) {
  dnf <- denovo_candidates[[dnf_idx]]
  
  n_mut <- sapply(split(dnf, dnf$trio), nrow)
  t_vec <- names(n_mut)
  
  rate <- n_mut / (2 * auto_size[t_vec] * callability[dnf_idx,t_vec])
  return(rate)
})
GATK_rates <- sapply(GATK_ratebytrio, mean)

callability_data <- data.frame(GQ = seq(20,80,10),
                               mutations = unlist(lapply(denovo_candidatesAD1, nrow)),
                               meanCallability = rowSums(callability_AD1) / 14,
                               meanCallableSites = (callability_AD1 %*% auto_size) / (14 * 1e9),
                               meanMutationRate = GATK_ratesAD1 * 1e9)
p1 <- ggplot(callability_data, aes(x = GQ, y = mutations)) + geom_line() + geom_point() + theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
p2 <- ggplot(callability_data, aes(x = GQ, y = meanCallability)) + ylab("Callability") + geom_line() + geom_point() + theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
p3 <- ggplot(callability_data, aes(x = GQ, y = meanCallableSites)) + ylab("Callable Sites") + geom_line() + geom_point() + theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
p4 <- ggplot(callability_data, aes(x = GQ, y = meanMutationRate)) + ylab("Mutation Rate") + geom_line() + geom_point() + theme_bw() + ylim(5, 6) +
  geom_hline(yintercept = mean(rate_FBAD1) * 1e9, col = 'blue')
GGp1 <- ggplotGrob(p1)
GGp2 <- ggplotGrob(p2)
GGp3 <- ggplotGrob(p3)
GGp4 <- ggplotGrob(p4)
grid::grid.newpage()
grid::grid.draw(rbind(GGp1, GGp2, GGp3, GGp4))

plot(unlist(lapply(denovo_candidates, nrow)), type = 'o', xaxt = 'n', xlab = "", ylab = "number of mutations", pch = 21, bg = 'firebrick', lwd = 0.5)
plot(rowSums(callability) / 14, type = 'o', xaxt = 'n', xlab = "", ylab = "mean callability", pch = 21, bg = 'grey', lwd = 0.5)
plot(as.vector(callability %*% auto_size) / (14 * 1e9), type = 'o', xaxt = 'n', xlab = "", ylab = expression("callable sites × 10"^9*" (callability * pileup count)"), pch = 21, bg = 'darkblue', lwd = 0.5)
plot(GATK_rates * 1e9, type = 'o', xaxt = 'n', xlab = "", ylab = expression("mutation rate (× 10"^9*") per bp per generation"), ylim = c(6.5, 7.5), pch = 21, bg = 'firebrick', lwd = 0.5)

plot(unlist(lapply(denovo_candidatesAD1, nrow)), type = 'o', xaxt = 'n', xlab = "", ylab = "number of mutations", pch = 21, bg = 'firebrick', lwd = 0.5)
abline(h = sum(n_mut_FBAD1), col = 'blue')
plot(rowSums(callability_AD1) / 14, type = 'o', xaxt = 'n', xlab = "", ylab = "mean callability", pch = 21, bg = 'grey', lwd = 0.5)
abline(h = mean(callability_FB_AD1))
plot(as.vector(callability_AD1 %*% auto_size) / (14 * 1e9), type = 'o', xaxt = 'n', xlab = "", ylab = expression("callable sites × 10"^9*" (callability * pileup count)"), pch = 21, bg = 'darkblue', lwd = 0.5)
plot(GATK_ratesAD1 * 1e9, type = 'o', xaxt = 'n', xlab = "", ylab = expression("mutation rate (× 10"^9*") per bp per generation"), pch = 21, ylim = c(5, 6.5), bg = 'firebrick', lwd = 0.5)
abline(h = mean(rate_FBAD1) * 1e9, lwd = 2, col = 'blue')

GATK_ratebytrio_mat <- matrix(unlist(GATK_ratebytrio), nrow = 7, byrow = TRUE)
plot(0, xlim = c(1,7), ylim = c(1e-9,1.7e-8), xlab = "", ylab = "")
for(trio in 1:14) {
  lines(as.vector(GATK_ratebytrio_mat[,trio]), type = 'o', pch = 21, bg = colorRampPalette(c("red", "blue"))(14)[trio])
  text(7, GATK_ratebytrio_mat[7,trio], labels = paste("trio", trio), pos = 1, cex = 0.6)
}  


mutation_count <- aggregate(denovo_candidates[[4]]$trans_flag, by = list(denovo_candidates[[4]]$trio), FUN = length)
names(mutation_count) <- c("trio", "mutations")
mutation_count <- mutation_count[order(as.numeric(substring(mutation_count$trio, 5))),]

mutation_count[["het_callability"]] <- as.vector(hetc_mat[4,])
mutation_count[["homo_callability"]] <- as.vector(homAD1c_mat[4,])
mutation_count[["pile_size"]] <- auto_size

attach(mutation_count)
mutation_count[["rate"]] <- mutations / (2 * het_callability * homo_callability * pile_size)
mutation_count[["rate_SN"]] <- mutation_count$rate * 1e8
detach(mutation_count)

DOB_table <- read.table("../DOB.txt")
names(DOB_table) <- c("Individual", "Date")
DOB_table$Date <- as.Date(DOB_table$Date, format = '%m/%d/%Y')
DOB_vector <- DOB_table$Date
names(DOB_vector) <- DOB_table$Individual
gestation_length <- 166

DOB_idx <- match(mutation_count$trio, trio_table$trio)
mutation_count$GTMale <- DOB_vector[as.character(trio_table[DOB_idx,]$Child)] - DOB_vector[as.character(trio_table[DOB_idx,]$Father)] - gestation_length
mutation_count$GTFemale <- DOB_vector[as.character(trio_table[DOB_idx,]$Child)] - DOB_vector[as.character(trio_table[DOB_idx,]$Mother)] - gestation_length
mutation_count$GTMale <- as.numeric(mutation_count$GTMale / 365)
mutation_count$GTFemale <- as.numeric(mutation_count$GTFemale / 365)

col1 <- "#d95f02"
col2 <- "#7570b3"

ggplot(mutation_count, aes(x = GTMale, y = rate_SN)) + coord_cartesian(xlim = c(1, 26), ylim = c(0,2)) +
  scale_x_continuous(limits = c(1.5,40), breaks = seq(0,30,5)) +
  scale_y_continuous(limits = c(0,2), breaks = seq(0,2,0.25)) +
  geom_point(bg = col1, size = 4, shape = 21) + geom_smooth(fullrange = TRUE, method = lm, color = col1, fill = col1) + theme_classic()

ggplot(mutation_count, aes(x = GTFemale, y = rate_SN)) + coord_cartesian(xlim = c(1, 26), ylim = c(0,2)) +
  scale_x_continuous(limits = c(1.5,40), breaks = seq(0,30,5)) +
  scale_y_continuous(limits = c(0,2), breaks = seq(0,2,0.25)) +
  geom_point(bg = col2, size = 4, shape = 21) + geom_smooth(fullrange = TRUE, method = lm, color = col2, fill = col2) + theme_classic()

###

denovo_candidates_GATK <- candidates_DP %>% filter_gq(., 70) %>% filter_bam %>% filter_rel %>% filter_allelicbalance(0.3, 1) %>% bind_candidates

isect_table <- c()
for(GQ in seq(20,80,10)) {
  denovo_candidates_GATK <- candidates_DP %>% filter_gq(., GQ) %>% filter_bam %>% filter_rel %>% filter_allelicbalance(0.3, 1) %>% bind_candidates
  
  dnc_GATK_idx <- denovo_candidates_GATK[,1:2]
  dnc_FBidx <- denovo_candidates_FB[,1:2]
  
  isect_table <- rbind(isect_table, c(nrow(dnc_GATK_idx),
                                      nrow(dnc_FBidx),
                                      nrow(inner_join(dnc_GATK_idx, dnc_FBidx)), # GATK & FB
                                      nrow(anti_join(dnc_GATK_idx, dnc_FBidx)),  # GATK & !FB
                                      nrow(anti_join(dnc_FBidx, dnc_GATK_idx))))  # FB & !GATK
}

c("GQ", "GATK", "FB", "GATK & FB", "GATK & !FB", "FB & !GATK")
