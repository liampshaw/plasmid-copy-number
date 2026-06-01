maddamsetti = read.csv('data/Maddamsetti-Fig1BC-Source-Data.csv', header=T)
ramiro.martinez = read.csv('data/ramiro-supplementary-dataset-1-per-plasmid.csv', header=T)
plsdb = read.csv('data/plsdb-typing.csv', header=T)
rownames(plsdb) = plsdb$NUCCORE_ACC

amrfinder = read.csv('data/amrfinder-results.tsv', header=T, sep='\t')
amrfinder.amr = amrfinder[which(amrfinder$Type=="AMR"),]
# Group by contig ID
amr.summary = data.frame(amrfinder.amr %>% group_by(Contig.id) %>%
                           summarise(n=length(Contig.id)))
rownames(amr.summary) = amr.summary$Contig.id
maddamsetti$amr = ifelse(maddamsetti$SeqID %in% rownames(amr.summary), "yes", "no")


# So we can use consistent names
maddamsetti$PCN = maddamsetti$InitialCopyNumberEstimate
maddamsetti$log10_PCN = log10(maddamsetti$PCN)
maddamsetti$Size = maddamsetti$replicon_length 
maddamsetti$log10_Size = log10(maddamsetti$Size)

ramiro.martinez$SeqID = ramiro.martinez$Contig

# Drop Maddamsetti plasmids that aren't in PLSDB
# Filter out missing (~5%)
maddamsetti = maddamsetti[which(maddamsetti$SeqID %in% rownames(plsdb)),]
maddamsetti$rep.types = plsdb[maddamsetti$SeqID, "rep_type.s."]
maddamsetti$relaxase.types = plsdb[maddamsetti$SeqID, "relaxase_type.s."]
maddamsetti$mpf.types = plsdb[maddamsetti$SeqID, "mpf_type"]