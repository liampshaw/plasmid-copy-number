
source('load-data.R')


IncF.plasmids = maddamsetti[grep("IncF", maddamsetti$rep.types),]
those.with.sopB = read.csv('data/IncF-sopB-hits-ids.txt', header=F)$V1
IncF.plasmids$with.sopB = ifelse(IncF.plasmids$SeqID %in% those.with.sopB,
                                 "sopB hit", "No sopB hit")
IncFIA.plasmids = IncF.plasmids[which(IncF.plasmids$rep.types=="IncFIA"),]

log_labels <- function(x) {
  parse(text = ifelse(
    x < 10,
    as.character(x),
    paste0(mantissa <- x / 10^floor(log10(x)),
           ifelse(mantissa == 1,
                  paste0("10^", floor(log10(x))),
                  paste0(x / 10^floor(log10(x)), " %*% 10^", floor(log10(x)))))
  ))
}

pdf('figures/Fig4-IncFIA.pdf', width=8, height=4)
ggplot(IncFIA.plasmids, aes(Size, PCN, group=with.sopB))+
  scale_y_log10()+
  stat_smooth(method="lm", colour="grey")+
  geom_point()+
  facet_wrap(~ with.sopB) + 
  theme_bw()+
  scale_x_log10(breaks = c(1e4, 3e4, 1e5,3e5),
labels = scales::parse_format()(c("10^4", "3 %*% 10^4", "10^5", "3 %*% 10^5")))+
  ylab("Copy number")
dev.off()


# Single vs. segmented model
incF.lm = lm(log10_PCN ~ log10_Size, data=IncFIA.plasmids)
incF.segmented.lm = segmented(
  incF.lm,
  seg.Z = ~log10_Size,
  psi = list(log10_Size = 4.25))
AIC(incF.lm, incF.segmented.lm)

# Instead, look at with/without sopB
summary(lm(log10_PCN ~ log10_Size, data=IncFIA.plasmids))
summary(lm(log10_PCN ~ log10_Size, data=IncFIA.plasmids[which(IncFIA.plasmids$with.sopB=="yes"),]))
summary(lm(log10_PCN ~ log10_Size, data=IncFIA.plasmids[which(IncFIA.plasmids$with.sopB=="no"),]))


# With sopB
incF.sopB.lm = lm(log10_PCN ~ log10_Size, data=IncF.plasmids[which(IncF.plasmids$with.sopB=="yes"),])
incF.sopB.segmented.lm = segmented(
  incF.sopB.lm,
  seg.Z = ~log10_Size,
  psi = list(log10_Size = 5))
AIC(incF.sopB.lm, incF.sopB.segmented.lm)
