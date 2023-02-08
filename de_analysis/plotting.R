pdf("Raw_expression_counts_plots.pdf", width = 10, height = )

ggplot(data) +
  geom_histogram(aes(x = AIPL1_6_1), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")
#jpeg('rplot.jpg')

ggplot(data) +
  geom_histogram(aes(x = AIPL1_6_1), stat = "bin", bins = 200) + 
  xlim(-5, 500)  +
  xlab("Raw expression counts") +
  ylab("Number of genes")
dev.off() 

## So
pdf("mean_counts_plots.pdf")
mean_counts_AIPL3 <- apply(data[, 5:7], 1, mean)
variance_counts <- apply(data[, 5:7], 1, var)
df <- data.frame(mean_counts_AIPL3, variance_counts)

ggplot(df) +
  geom_point(aes(x=mean_counts_AIPL3, y=variance_counts)) + 
  geom_line(aes(x=mean_counts_AIPL3, y=mean_counts_AIPL3, color="red")) +
  scale_y_log10() +
  scale_x_log10()


mean_counts_AIPL16 <- apply(data[, 8:10], 1, mean)
variance_counts <- apply(data[, 8:10], 1, var)
df <- data.frame(mean_counts_AIPL16, variance_counts)

ggplot(df) +
  geom_point(aes(x=mean_counts_AIPL16, y=variance_counts)) + 
  geom_line(aes(x=mean_counts_AIPL16, y=mean_counts_AIPL16, color="red")) +
  scale_y_log10() +
  scale_x_log10()


mean_counts_CC <- apply(data[, 11:13], 1, mean)
variance_counts <- apply(data[, 11:13], 1, var)
df <- data.frame(mean_counts_CC, variance_counts)

ggplot(df) +
  geom_point(aes(x=mean_counts_CC, y=variance_counts)) + 
  geom_line(aes(x=mean_counts_CC, y=mean_counts_CC, color="red")) +
  scale_y_log10() +
  scale_x_log10()


mean_counts_GFP <- apply(data[, 14:16], 1, mean)
variance_counts <- apply(data[, 14:16], 1, var)
df <- data.frame(mean_counts_GFP, variance_counts)

ggplot(df) +
  geom_point(aes(x=mean_counts_GFP, y=variance_counts)) + 
  geom_line(aes(x=mean_counts_GFP, y=mean_counts_GFP, color="red")) +
  scale_y_log10() +
  scale_x_log10()
dev.off() ## Close writing in pdf-file


