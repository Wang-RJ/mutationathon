gqs <- extract_gq(candidates_DP_GQ_BAM)[[1]]
ggpf <- rbind(data.frame(id = "M", GQ = gqs[,1]),
              data.frame(id = "Noot", GQ = gqs[,2]),
              data.frame(id = "Heineken", GQ = gqs[,3]))
ggab <- data.frame(allelic_balance = extract_allelic_balance(candidates_DP_GQ_BAM)[[1]])

ggdpgq_dp <- rbind(data.frame(id = "M", DP = dpgq_depths[[1]][,1]),
                   data.frame(id = "Noot", DP = dpgq_depths[[1]][,2]),
                   data.frame(id = "Heineken", DP = dpgq_depths[[1]][,3]))

ggdpgq_gq <- rbind(data.frame(id = "M", GQ = dpgq_gq[[1]][,1]),
                   data.frame(id = "Noot", GQ = dpgq_gq[[1]][,2]),
                   data.frame(id = "Heineken", GQ = dpgq_gq[[1]][,3]))

ggdpgq_ab <- data.frame(allelic_balance = extract_allelic_balance(candidates_DP_GQ)[[1]])

ggplot(ggdpgq_dp, aes(DP, fill = id)) + geom_histogram(position = "dodge") + theme_bw() +
  theme(
    text = element_text(size = 20),
    legend.text = element_text(size = 10),
    axis.text = element_text(size = 20),
    legend.justification=c(0,0), legend.position=c(0.05,0.8))

ggplot(ggdpgq_gq, aes(GQ, fill = id)) + geom_histogram(position = "dodge") + theme_bw() +
  theme(
    text = element_text(size = 20),
    legend.text = element_text(size = 10),
    axis.text = element_text(size = 20),
    legend.justification=c(0,0), legend.position=c(0.05,0.8))

ggplot(ggdpgq_ab, aes(allelic_balance)) + geom_histogram(binwidth = 0.05, fill = "white", color = "black") + theme_bw() +
  theme(
    text = element_text(size = 20),
    legend.text = element_text(size = 10),
    axis.text = element_text(size = 20),
    legend.justification=c(0,0), legend.position=c(0.05,0.8))

ggplot(ggab, aes(allelic_balance)) + geom_histogram(binwidth = 0.05, fill = "white", color = "black") + theme_bw() +
  theme(
    text = element_text(size = 20),
    legend.text = element_text(size = 10),
    axis.text = element_text(size = 20),
    legend.justification=c(0,0), legend.position=c(0.05,0.8))

ggplot(ggpf, aes (GQ, fill = id)) + geom_histogram(position = "dodge", binwidth = 4) + theme_bw() +
  theme(
    text = element_text(size = 20),
    legend.text = element_text(size = 10),
    axis.text = element_text(size = 20),
    legend.justification=c(0,0), legend.position=c(0.05,0.8))