# Analyse Clowes (1972) Table 8 on PCN and size
library(ggplot2)
library(scales)
clowes = read.csv('data/Clowes-1972.csv', header=T)

base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}

p.clowes = ggplot(clowes, aes(Length, PCN))+
  stat_smooth(method="lm")+
  geom_point()+
  ggrepel::geom_text_repel(aes(label=Plasmid))+
  scale_y_continuous(trans = log_trans(), breaks = base_breaks(5)) +
  scale_x_continuous(trans = log_trans(), breaks = base_breaks(), labels=prettyNum) +
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(axis.text=element_text(colour="black", size=8))
pdf('figures/Figure-1-Clowes.pdf', width=4, height=4)
p.clowes
dev.off()

summary(lm(log10(PCN)~log10(Length), data=clowes))
