# Gut microbiome structure and metabolic activity in inflammatory bowel disease
# Fransoza 2019, Nature Microbiology

library(readxl)
# load dataset
metab_assoc <- read_excel('W:/Reh-Jer/Innovative Technologies/Cross-omics/crossomics/Franzosa_IBD/Processed_data/METABOLITES_ASSOCIATION_cut.xlsx')
raw_data <- read_excel('W:/Reh-Jer/Innovative Technologies/Cross-omics/crossomics/Franzosa_IBD/Processed_data/METABOLOMICS_SAMPLES.xlsx')

# extract metadata
metadata <- raw_data[1:7,]
names <- metadata$`# Feature / Sample`
metadata <- metadata[,-1]
rownames(metadata) <- names
metadata <- t(metadata)
metadata <- data.frame(metadata[1:155,])

# join the two
library(dplyr)
join1 <- metab_assoc[,1:2]
join2 <- raw_data[-(1:7),]

# rename columns
colnames(join1) = c('Sample', 'Metabolite')
colnames(join2)[1] = c('Sample')

# merge the two datasets
df_raw <- data.frame(inner_join(join1, join2, by = "Sample"))

# manipulate df
df <- df_raw[,-(1:2)]
rownames(df) <- df_raw$Sample
df <- data.frame(t(df))
df <- df[1:155,] # extract PRISM cohort only (n=155)

# PCA
df <- data.frame(sapply(df, as.numeric)) # columns to numeric
pca <- prcomp(df, center = TRUE, scale = FALSE, retx = TRUE)
autoplot(pca, data = metadata, colour = 'Diagnosis') + theme_classic()

# PCoA
library(vegan)
dists <- vegdist(df, method = "bray")
pcoa <- data.frame(cmdscale(dists))
pcoa_metadata <- cbind(pcoa,metadata)
ggplot(pcoa_metadata, aes(x=X1, y=X2, color=Diagnosis, shape=Diagnosis)) + geom_point() + theme_classic()
