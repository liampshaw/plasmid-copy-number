library(dplyr)
library(ggplot2)

source('scripts/load-data.R')

# Subset data by AMR for each size bin
maddamsetti %>% group_by(Size)

summary(lm(amr.N ~ Size+PCN, data=maddamsetti))

df = maddamsetti
df$amr_present <- as.numeric(df$amr.N > 0)

p <- ggplot(df, aes(x = Size, y = PCN, z = amr_present)) +
  stat_summary_hex(fun = mean, bins = 25) +
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_viridis_c(name = "Proportion\nwith AMR", limits = c(0, 1)) +
  labs(x = "Plasmid size (kb, log scale)", y = "Plasmid copy number (log scale)",
       title = "Proportion of plasmids carrying AMR genes, by size and copy number") +
  theme_minimal()

ggplot(maddamsetti, aes(Size, amr.N))+
  geom_point()+
  scale_x_log10()

ggplot(maddamsetti, aes(PCN, amr.N))+
  geom_point()+
  scale_x_log10()

amr.high.copy.plasmids = maddamsetti[which(maddamsetti$amr=="yes" & maddamsetti$PCN>10),"SeqID"]
sort(table(amrfinder[which(amrfinder$Contig.id %in% amr.high.copy.plasmids),"Class"]))
# How many of the plasmids carry a beta-lactamase gene?
length(amr.high.copy.plasmids)
length(unique(sort(amrfinder[which(amrfinder$Contig.id %in% amr.high.copy.plasmids & amrfinder$Class=="BETA-LACTAM"),"Contig.id"])))
      

amr.low.copy.plasmids = maddamsetti[which(maddamsetti$amr=="yes" & maddamsetti$PCN<10),"SeqID"]
length(amr.low.copy.plasmids)
length(unique(sort(amrfinder[which(amrfinder$Contig.id %in% amr.low.copy.plasmids & amrfinder$Class=="BETA-LACTAM"),"Contig.id"])))
1781/2716

amr.1.low.copy.plasmids = maddamsetti[which(maddamsetti$amr.N==1 & maddamsetti$PCN<3),"SeqID"]
length(amr.1.low.copy.plasmids)
length(unique(sort(amrfinder[which(amrfinder$Contig.id %in% amr.1.low.copy.plasmids & amrfinder$Class=="BETA-LACTAM"),"Contig.id"])))
271/541
