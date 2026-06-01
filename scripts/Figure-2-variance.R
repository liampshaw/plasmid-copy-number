source('scripts/load-data.R')

# RELATIONSHIP WITH VARIANCE
# Bin by length
df <- maddamsetti %>%
  mutate(length_bin = cut(Size, 
                          breaks = seq(0, max(Size) + 1000, by = 1000),
                          right = FALSE))
variance_by_bin <- df %>%
  group_by(length_bin) %>%
  summarise(
    bin_midpoint = mean(Size),  # for plotting
    pcn_variance = var(PCN),
    pcn_sd = sd(PCN),
    n = n(),
    mean_pcn = mean(PCN),
    .groups = 'drop'
  )

# Data frame for Becker et al. dataset
becker = read.csv('data/Becker-2016.csv', header=T)

becker.long = becker %>% pivot_longer(cols=starts_with("X"))
becker.long$size = as.numeric(gsub("X", "", gsub("kb.plasmid", "", becker.long$name)))*1000
becker.long$PCN = as.numeric(becker.long$value)/becker.long$Chromosome

becker.long$plasmid = ordered(gsub(".plasmid", "", gsub("X", "", becker.long$name)), 
                              levels=c("3.8kb", "4.8kb", "362kb"))
becker.df = becker.long %>% group_by(size) %>% summarise(var=sd(PCN))
becker.df$bin_midpoint = becker.df$size
becker.df$pcn_variance = becker.df$var
pdf("figures/Figure-2-copy-number-variance.pdf", width=6, height=4)
p.variance = ggplot(variance_by_bin, aes(bin_midpoint, pcn_variance))+
  geom_point(aes(size=n), colour="grey")+
  theme_bw()+
  scale_x_log10(limits=c(0.7e3, 3e6))+
  scale_y_log10()+  xlab("Plasmid size (1kb bins)")+
  ylab("Copy number variance")+
  stat_smooth(method="lm", se=FALSE, colour="black", linetype='dashed')
p.variance + 
  stat_smooth(method="lm", data=becker.df, se=FALSE,colour="red", linetype='dashed')+
  geom_point(data=becker.df, colour="red", shape=17, size=4)
dev.off()