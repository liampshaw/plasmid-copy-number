# Fitting biphasic law

library(segmented)
library(ggplot2)
library(dplyr)
library(ggExtra)
#library(tidyclust)
library(scales)
library(purrr)
library(tidyr)
library(stringr)
library(lme4)
library(lmerTest)  # gives p-values for fixed effects


source('load-data.R')

# Ramiro-Martinez
# Summarise PTUs
ramiro.martinez.PTU.summary = ramiro.martinez %>% filter(Predicted_PTU!="-" & Predicted_PTU!="PTU-?") %>% group_by(Predicted_PTU) %>% 
  summarise(median_size=median(Size),
            min_size=min(Size),
            iqr_upper_size=quantile(Size, prob=0.75),
            iqr_lower_size=quantile(Size, prob=0.25),
            max_size=max(Size),
            median_PCN=median(PCN),
            min_PCN=min(PCN),
            max_PCN=max(PCN),
            iqr_upper_PCN=quantile(PCN, prob=0.75),
            iqr_lower_PCN=quantile(PCN, prob=0.25),
            n=length(Size))

p.ramiro.martinez.PTUs = ggplot(ramiro.martinez.PTU.summary, aes(median_size, median_PCN))+
  geom_point(aes(size=n))+
  geom_errorbar(aes(ymin=iqr_lower_PCN, ymax=iqr_upper_PCN))+
  geom_errorbarh(aes(xmin=iqr_lower_size, xmax=iqr_upper_size))+
  scale_x_log10()+
  scale_y_log10()+
  stat_smooth(method="lm", se=FALSE, colour="grey")+
  ggtitle("(b) Ramiro-Martínez dataset (PTUs)")+
  theme(legend.position = c(0.8, 0.8))
# Try a segmented regression on this data:
ramiro.martinez.PTU.summary$log10_size = log10(ramiro.martinez.PTU.summary$median_size)
ramiro.martinez.PTU.summary$log10_PCN = log10(ramiro.martinez.PTU.summary$median_PCN)

ramiro.martinez.PTU.lm = lm(log10_PCN ~ log10_size, data=ramiro.martinez.PTU.summary)
ramiro.martinez.PTU.segmented.lm = segmented(
  ramiro.martinez.PTU.lm,
  seg.Z = ~log10_size,
  psi = list(log10_size = 4.25))
AIC(ramiro.martinez.PTU.lm, ramiro.martinez.PTU.segmented.lm)
summary(ramiro.martinez.PTU.segmented.lm)
# Add segmented model predictions to the dataset so we can plot
ramiro.martinez.PTU.summary$segmented.lm.pre = predict(ramiro.martinez.PTU.segmented.lm, data=ramiro.martinez.PTU.summary$log10_size)

#pdf('figures/ramiro-martinez-PTUs.pdf', width=6, height=4)
p.ramiro.martinez.PTUs.fit = p.ramiro.martinez.PTUs + 
 # geom_vline(xintercept = 10**ramiro.martinez.PTU.segmented.lm$psi[2], linetype='dashed', size=2, colour='grey')+
  geom_line(data=ramiro.martinez.PTU.summary, aes(10**log10_size, 10**segmented.lm.pre), colour="red", size=1)+
  theme_bw()+
  theme(legend.position=c(0.8,0.8))+
  theme(panel.grid=element_blank())+
  xlab("Size")+
  ylab("Copy number")

# Add: values of k for each, AIC comparison

# Maddamsetti - individual plasmids
maddamsetti.lm = lm(log10_PCN ~ log10_Size, data=maddamsetti)
maddamsetti.segmented.lm = segmented(
  maddamsetti.lm,
  seg.Z = ~log10_Size,
  psi = list(log10_Size = 4.25))
AIC(mdadamsetti.segmented.lm, maddamsetti.lm)
p.maddamsetti = ggplot(maddamsetti, aes(Size, PCN))+
  geom_point(size=0.2)+
  scale_x_log10()+
  scale_y_log10()+
  stat_smooth(method="lm", se=FALSE, colour="grey")+
  ggtitle("(a) Maddamsetti dataset")+
  xlab("Size")+
  ylab("Copy number")
# Add segmented model predictions to the dataset so we can plot
maddamsetti$segmented.lm.pre = predict(maddamsetti.segmented.lm, data=maddamsetti$log10_Size)

p.maddamsetti.fit = p.maddamsetti + 
  #geom_vline(xintercept = 10**maddamsetti.segmented.lm$psi[2], linetype='dashed', size=2, colour='grey')+
  geom_line(data=maddamsetti, aes(10**log10_Size, 10**segmented.lm.pre), colour="red", size=1)+
  theme_bw()+
  theme(panel.grid=element_blank())

# Limits
x.limits = c(1000, max(maddamsetti$Size))
y.limits = c(min(min(ramiro.martinez$PCN), min(maddamsetti$PCN)),
             max(max(ramiro.martinez$PCN), max(maddamsetti$PCN)))
pdf("figures/Fig3-regression.pdf", width=8, height=4)
cowplot::plot_grid( p.maddamsetti.fit+
                      scale_x_log10(limits=x.limits, labels = scales::label_log())+
                      scale_y_log10(limits=y.limits, labels = scales::label_log()), 
                    p.ramiro.martinez.PTUs.fit+
                      scale_x_log10(limits=x.limits, labels = scales::label_log())+
                      scale_y_log10(limits=y.limits, labels = scales::label_log()),
                    rel_widths=c(1,1))
dev.off()

# Breakpoints
# Maddamsetti
10**maddamsetti.segmented.lm$psi[2]
# Ramiro-Martinez
10**ramiro.martinez.PTU.segmented.lm$psi[2]
