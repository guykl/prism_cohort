##nmds code 
library("vegan", lib.loc="~/R/win-library/3.5")
setwd("C:\\Users\\azolty\\Desktop\\RWD\\cross_omics")


OTU=read.("METAGENOMICS_TAXA_DATA.xlsx", header=TRUE)

OTUtrans=t(OTU)
tiff("rarefaction_no_labes.tiff", height = 15, width = 20, units = 'cm', compression = "lzw", res = 500)
rarecurve(OTUtrans, step=100, xlab="Number of sequences", ylab="OTUs", label=F, font=4)
dev.off()

OTUmetadata=read.table("metadata2_NO73.txt", header=TRUE, row.names=1)

transformedOTU=decostand(OTUtrans, "hell")
disDNA <- vegdist(transformedOTU, method="bray")
nMDSDNA <- monoMDS(disDNA, model = "global")

nMDSscoresDNA<-scores(nMDSDNA)
with(OTUmetadata, levels(shapeSoilHabitat))

shapeR=c("Control"=8, "Brachypodium-soil"=15, "Brachypodium-rhizosphere"=17, "Brachypodium-root"=19, "Arabidopsis-rhizosphere"=24, "Arabidopsis-soil"=22, "Arabidopsis-root"=21)
colR=c("T1"="#66c2a5", "T2"="#fc8d62", "T3"="#8da0cb", "T4"="#e78ac3", "T5"="#a6d854", "T6"="#ffd92f", "T7"="#e5c494")


ordination_score<-as.data.frame(nMDSscoresDNA)
mergedscores= merge(ordination_score, OTUmetadata, by=0)
head(mergedscores)

tiff("nmdsALL_final.tiff", height = 15, width = 22, units = 'cm', compression = "lzw", res = 800)
ggplot(mergedscores, aes(x=MDS1, y=MDS2, color=Treatment, shape=shapeSoilHabitat, font=8)) + geom_point(size=3)+ theme_bw()+ scale_shape_manual(values=shapeR)+ scale_color_manual(values=colR)+ labs(color="Treatment", shape="Plant and niche")
dev.off()