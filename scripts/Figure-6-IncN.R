source('load-data.R')



betalactamase.plasmids = amrfinder.amr$Contig.id[which(amrfinder.amr$Class=="BETA-LACTAM")]
maddamsetti$betalactamase =  ifelse(maddamsetti$SeqID %in% betalactamase.plasmids, 
                                    "yes", "no")

IncN.plasmids = maddamsetti[which(maddamsetti$rep.types=="IncN"),]
IncN.plasmids$amr = ifelse(IncN.plasmids$amr=="yes", "With AMR genes", "Without AMR genes")


maddamsetti$novick = ifelse(maddamsetti$PCN>6, "Multicopy\n(>6x)",
                            ifelse(maddamsetti$PCN>2 & maddamsetti$PCN<6, "Oligocopy\n(2-6x)", "Unit copy\n(<2x)"))

p.novick = ggplot(maddamsetti, aes(novick, Size))+
  scale_y_log10(labels=scales::label_log())+
  geom_boxplot()+
  xlab("")+
  theme_bw()
ggsave(p.novick, file="figures/Novick-categories.pdf", width=6, height=4)


pdf("figures/Fig-6-IncN.pdf", width=8, height=4)
p.IncN = ggplot(IncN.plasmids, aes(Size, PCN, group=amr))+
  scale_y_log10()+
  stat_smooth(method="lm", colour="grey")+
  geom_point()+
  facet_wrap(~ amr, nrow=1) +
  scale_x_log10(breaks = c(4e4, 5e4, 7e4,1e5),
                labels = scales::parse_format()(c("4 %*% 10^4", "5 %*% 10^4", "7 %*% 10^4", "1 %*% 10^5")))+
  ylab("Copy number")+
  theme_bw()
p.IncN
dev.off()


# Consider plasmids <100kb
AMR.comparison = maddamsetti %>%   filter(Size>20000 & Size<100000) %>%
  mutate(size_bin = floor(Size / 1000) * 5000) %>%
  group_by(size_bin, Genus, betalactamase) %>%
  summarise(median_PCN = median(PCN, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = betalactamase, values_from = median_PCN) %>%
  mutate(PCN_diff = yes - no) %>%
  filter(!is.na(PCN_diff)) 

AMR.comparison.key.genera= AMR.comparison %>% filter(Genus %in% head(names(sort(table(AMR.comparison$Genus), decreasing = TRUE)), 3))
ggplot(AMR.comparison.key.genera, aes(size_bin, PCN_diff))+
  geom_bar(stat="identity")+
  facet_wrap(~Genus)+
  ylim(c(-5,5))+
  theme_bw()

table((AMR.comparison %>%filter(Genus=="Staphylococcus"))[,"PCN_diff"]>0)
table((AMR.comparison %>%filter(Genus=="Escherichia"))[,"PCN_diff"]>0)
table((AMR.comparison %>%filter(Genus=="Klebsiella"))[,"PCN_diff"]>0)

AMR.comparison %>%
  summarise(
    n_bins = n(),
    n_AMR_higher = sum(PCN_diff > 0),
    n_nonAMR_higher = sum(PCN_diff < 0),
    n_equal = sum(PCN_diff == 0),
    median_delta = median(PCN_diff),
    mean_delta = mean(PCN_diff)
  )



maddamsetti %>%   filter(Size>20000 & Size<100000 & Genus=="Escherichia") %>%
  mutate(size_bin = floor(Size / 1000) * 1000) %>%
  group_by(size_bin, amr) %>%
  summarise(median_PCN = median(PCN, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = amr, values_from = median_PCN) %>%
  mutate(PCN_diff = yes - no) %>%
  filter(!is.na(PCN_diff)) %>%
  summarise(
    n_bins = n(),
    n_AMR_higher = sum(PCN_diff > 0),
    n_nonAMR_higher = sum(PCN_diff < 0),
    n_equal = sum(PCN_diff == 0),
    median_delta = median(PCN_diff),
    mean_delta = mean(PCN_diff)
  )
