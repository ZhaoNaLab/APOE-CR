library(tidyverse)
library(openxlsx)
library(GUniFrac)
library(ggpubr)

raw_metadata <- as.data.frame(read_tsv("/research/bsi/archive/microbiome/s216158.apoE_mouse_microbiome//processing/metadata.txt"))
metadata <- raw_metadata
rownames(metadata) <- gsub("-", ".", raw_metadata$SampleID)
metadata$SampleID <- rownames(metadata)


species <- read.table("/research/labs/microbiome/chia/m141127/humann_nf_test/bracken_fixed_again_mpa/combined.species.ct.txt", sep="\t", header=TRUE, comment.char="")
colnames(species) <- gsub(".k2.025_bracken.species.report", "", colnames(species))
colnames(species) <- gsub("Qiao.", "", colnames(species))
colnames(species) <- gsub("_.*", "", colnames(species))
rownames(species) <- species$X.Classification
species <- subset(species, select = -c(X.Classification))

genus <- read.table("/research/labs/microbiome/chia/m141127/humann_nf_test/bracken_fixed_again_mpa/combined.genus.ct.txt", sep="\t", header=TRUE, comment.char="")
colnames(genus) <- gsub(".k2.025_bracken.species.report", "", colnames(genus))
colnames(genus) <- gsub("Qiao.", "", colnames(genus))
colnames(genus) <- gsub("_.*", "", colnames(genus))
rownames(genus) <- genus$X.Classification
genus <- subset(genus, select = -c(X.Classification))




#Batch2 only, ND vs CR
metadata_b2 <- metadata[metadata$Batch=="Batch2" & metadata$Group %in% c("Normal", "CR"),]
species_b2 <- species[,rownames(metadata[metadata$Batch=="Batch2" & metadata$Group %in% c("Normal", "CR"),])]
genus_b2 <- genus[,rownames(metadata[metadata$Batch=="Batch2" & metadata$Group %in% c("Normal", "CR"),])]

metadata_b2_e2 <- metadata_b2[metadata_b2$genotype=="E2",]
species_b2_e2 <- species_b2[,rownames(metadata_b2[metadata_b2$genotype=="E2",])]
species_b2_e2 <- species_b2_e2[rowSums(species_b2_e2) > 0,]

genus_b2_e2 <- genus_b2[,rownames(metadata_b2[metadata_b2$genotype=="E2",])]
genus_b2_e2 <- genus_b2_e2[rowSums(genus_b2_e2) > 0,]


metadata_b2_e3 <- metadata_b2[metadata_b2$genotype=="E3",]
species_b2_e3 <- species_b2[,rownames(metadata_b2[metadata_b2$genotype=="E3",])]
species_b2_e3 <- species_b2_e3[rowSums(species_b2_e3) > 0,]

genus_b2_e3 <- genus_b2[,rownames(metadata_b2[metadata_b2$genotype=="E3",])]
genus_b2_e3 <- genus_b2_e3[rowSums(genus_b2_e3) > 0,]



metadata_b2_e4 <- metadata_b2[metadata_b2$genotype=="E4",]
species_b2_e4 <- species_b2[,rownames(metadata_b2[metadata_b2$genotype=="E4",])]
species_b2_e4 <- species_b2_e4[rowSums(species_b2_e4) > 0,]

genus_b2_e4 <- genus_b2[,rownames(metadata_b2[metadata_b2$genotype=="E4",])]
genus_b2_e4 <- genus_b2_e4[rowSums(genus_b2_e4) > 0,]


## Batch2, ND vs CR, E2 genotype

set.seed(123)
ZicoSeq.e2 <- ZicoSeq(meta.dat = metadata_b2_e2,
                       feature.dat = as.matrix(species_b2_e2),
                       grp.name = 'Group', adj.name = NULL, feature.dat.type = "count",
                       # Filter to remove rare taxa
                       prev.filter = 0, mean.abund.filter = 0,
                       max.abund.filter = 0.002, min.prop = 0,
                       # Winsorization to replace outliers
                       is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                       # Posterior sampling
                       is.post.sample = TRUE, post.sample.no = 25,
                       # Use the square-root transformation
                       link.func = list(function (x) x^0.5), stats.combine.func = max,
                       # Permutation-based multiple testing correction
                       perm.no = 999,  strata = NULL,
                       # Reference-based multiple stage normalization
                       ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                       # Family-wise error rate control
                       is.fwer = TRUE, verbose = TRUE, return.feature.dat = T)


df = data.frame(t(ZicoSeq.e2$feature.dat), check.names = FALSE)
prev = data.frame(df > 0, check.names = FALSE)
df$var = ZicoSeq.e2$meta.dat[[ZicoSeq.e2$grp.name]]
prev$var = ZicoSeq.e2$meta.dat[[ZicoSeq.e2$grp.name]]

group_means = df %>%
  mutate_if(is.numeric, ~ . / rowSums(select(df, where(is.numeric)))) %>% 
  group_by(var) %>%
  summarise_all(mean) %>%
  pivot_longer(cols=-var) %>%
  pivot_wider(names_from = var, values_from = value)

colnames(group_means) = c("Taxa", paste0(colnames(group_means)[-1], "_mean"))


group_prev = prev %>%
  group_by(var) %>%
  summarise_all(sum) %>%
  pivot_longer(cols=-var) %>%
  pivot_wider(names_from = var, values_from = value)

colnames(group_prev) = c("Taxa", paste0(colnames(group_prev)[-1], "_preval"))


zico_tab <- merge(group_means, group_prev,by="Taxa", all=TRUE)
zico_tab$R2 <- ZicoSeq.e2$R2[zico_tab$Taxa,]
zico_tab$padj <- ZicoSeq.e2$p.adj.fdr[zico_tab$Taxa]

zico_tab <- zico_tab[c("Taxa", "padj", "R2", colnames(zico_tab)[!colnames(zico_tab) %in% c("Taxa", "padj", "R2")])  ]

zico_tab_e2.sp <- tibble(zico_tab) %>% mutate(across(where(is.numeric), signif, 10))



ZicoSeq.e2 <- ZicoSeq(meta.dat = metadata_b2_e2,
                      feature.dat = as.matrix(genus_b2_e2),
                      grp.name = 'Group', adj.name = NULL, feature.dat.type = "count",
                      # Filter to remove rare taxa
                      prev.filter = 0, mean.abund.filter = 0,
                      max.abund.filter = 0.002, min.prop = 0,
                      # Winsorization to replace outliers
                      is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                      # Posterior sampling
                      is.post.sample = TRUE, post.sample.no = 25,
                      # Use the square-root transformation
                      link.func = list(function (x) x^0.5), stats.combine.func = max,
                      # Permutation-based multiple testing correction
                      perm.no = 999,  strata = NULL,
                      # Reference-based multiple stage normalization
                      ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                      # Family-wise error rate control
                      is.fwer = TRUE, verbose = TRUE, return.feature.dat = T)


df = data.frame(t(ZicoSeq.e2$feature.dat), check.names = FALSE)
prev = data.frame(df > 0, check.names = FALSE)
df$var = ZicoSeq.e2$meta.dat[[ZicoSeq.e2$grp.name]]
prev$var = ZicoSeq.e2$meta.dat[[ZicoSeq.e2$grp.name]]

group_means = df %>%
  mutate_if(is.numeric, ~ . / rowSums(select(df, where(is.numeric)))) %>% 
  group_by(var) %>%
  summarise_all(mean) %>%
  pivot_longer(cols=-var) %>%
  pivot_wider(names_from = var, values_from = value)

colnames(group_means) = c("Taxa", paste0(colnames(group_means)[-1], "_mean"))


group_prev = prev %>%
  group_by(var) %>%
  summarise_all(sum) %>%
  pivot_longer(cols=-var) %>%
  pivot_wider(names_from = var, values_from = value)

colnames(group_prev) = c("Taxa", paste0(colnames(group_prev)[-1], "_preval"))


zico_tab <- merge(group_means, group_prev,by="Taxa", all=TRUE)
zico_tab$R2 <- ZicoSeq.e2$R2[zico_tab$Taxa,]
zico_tab$padj <- ZicoSeq.e2$p.adj.fdr[zico_tab$Taxa]

zico_tab <- zico_tab[c("Taxa", "padj", "R2", colnames(zico_tab)[!colnames(zico_tab) %in% c("Taxa", "padj", "R2")])  ]

zico_tab_e2.ge <- tibble(zico_tab) %>% mutate(across(where(is.numeric), signif, 10))




## Batch2, ND vs CR, E3 genotype


set.seed(123)
ZicoSeq.e3 <- ZicoSeq(meta.dat = metadata_b2_e3,
                       feature.dat = as.matrix(species_b2_e3),
                       grp.name = 'Group', adj.name = NULL, feature.dat.type = "count",
                       # Filter to remove rare taxa
                       prev.filter = 0, mean.abund.filter = 0,
                       max.abund.filter = 0.002, min.prop = 0,
                       # Winsorization to replace outliers
                       is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                       # Posterior sampling
                       is.post.sample = TRUE, post.sample.no = 25,
                       # Use the square-root transformation
                       link.func = list(function (x) x^0.5), stats.combine.func = max,
                       # Permutation-based multiple testing correction
                       perm.no = 999,  strata = NULL,
                       # Reference-based multiple stage normalization
                       ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                       # Family-wise error rate control
                       is.fwer = TRUE, verbose = TRUE, return.feature.dat = T)

df = data.frame(t(ZicoSeq.e3$feature.dat), check.names = FALSE)
prev = data.frame(df > 0, check.names = FALSE)
df$var = ZicoSeq.e3$meta.dat[[ZicoSeq.e3$grp.name]]
prev$var = ZicoSeq.e3$meta.dat[[ZicoSeq.e3$grp.name]]

group_means = df %>%
  mutate_if(is.numeric, ~ . / rowSums(select(df, where(is.numeric)))) %>% 
  group_by(var) %>%
  summarise_all(mean) %>%
  pivot_longer(cols=-var) %>%
  pivot_wider(names_from = var, values_from = value)

colnames(group_means) = c("Taxa", paste0(colnames(group_means)[-1], "_mean"))


group_prev = prev %>%
  group_by(var) %>%
  summarise_all(sum) %>%
  pivot_longer(cols=-var) %>%
  pivot_wider(names_from = var, values_from = value)

colnames(group_prev) = c("Taxa", paste0(colnames(group_prev)[-1], "_preval"))


zico_tab <- merge(group_means, group_prev,by="Taxa", all=TRUE)
zico_tab$R2 <- ZicoSeq.e3$R2[zico_tab$Taxa,]
zico_tab$padj <- ZicoSeq.e3$p.adj.fdr[zico_tab$Taxa]

zico_tab <- zico_tab[c("Taxa", "padj", "R2", colnames(zico_tab)[!colnames(zico_tab) %in% c("Taxa", "padj", "R2")])  ]

zico_tab_e3.sp <- tibble(zico_tab) %>% mutate(across(where(is.numeric), signif, 10))


set.seed(123)
ZicoSeq.e3 <- ZicoSeq(meta.dat = metadata_b2_e3,
                      feature.dat = as.matrix(genus_b2_e3),
                      grp.name = 'Group', adj.name = NULL, feature.dat.type = "count",
                      # Filter to remove rare taxa
                      prev.filter = 0, mean.abund.filter = 0,
                      max.abund.filter = 0.002, min.prop = 0,
                      # Winsorization to replace outliers
                      is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                      # Posterior sampling
                      is.post.sample = TRUE, post.sample.no = 25,
                      # Use the square-root transformation
                      link.func = list(function (x) x^0.5), stats.combine.func = max,
                      # Permutation-based multiple testing correction
                      perm.no = 999,  strata = NULL,
                      # Reference-based multiple stage normalization
                      ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                      # Family-wise error rate control
                      is.fwer = TRUE, verbose = TRUE, return.feature.dat = T)

df = data.frame(t(ZicoSeq.e3$feature.dat), check.names = FALSE)
prev = data.frame(df > 0, check.names = FALSE)
df$var = ZicoSeq.e3$meta.dat[[ZicoSeq.e3$grp.name]]
prev$var = ZicoSeq.e3$meta.dat[[ZicoSeq.e3$grp.name]]

group_means = df %>%
  mutate_if(is.numeric, ~ . / rowSums(select(df, where(is.numeric)))) %>% 
  group_by(var) %>%
  summarise_all(mean) %>%
  pivot_longer(cols=-var) %>%
  pivot_wider(names_from = var, values_from = value)

colnames(group_means) = c("Taxa", paste0(colnames(group_means)[-1], "_mean"))


group_prev = prev %>%
  group_by(var) %>%
  summarise_all(sum) %>%
  pivot_longer(cols=-var) %>%
  pivot_wider(names_from = var, values_from = value)

colnames(group_prev) = c("Taxa", paste0(colnames(group_prev)[-1], "_preval"))


zico_tab <- merge(group_means, group_prev,by="Taxa", all=TRUE)
zico_tab$R2 <- ZicoSeq.e3$R2[zico_tab$Taxa,]
zico_tab$padj <- ZicoSeq.e3$p.adj.fdr[zico_tab$Taxa]

zico_tab <- zico_tab[c("Taxa", "padj", "R2", colnames(zico_tab)[!colnames(zico_tab) %in% c("Taxa", "padj", "R2")])  ]

zico_tab_e3.ge <- tibble(zico_tab) %>% mutate(across(where(is.numeric), signif, 10))



## Batch2, ND vs CR, E4 genotype

set.seed(123)
ZicoSeq.e4 <- ZicoSeq(meta.dat = metadata_b2_e4,
                      feature.dat = as.matrix(species_b2_e4),
                      grp.name = 'Group', adj.name = NULL, feature.dat.type = "count",
                      # Filter to remove rare taxa
                      prev.filter = 0, mean.abund.filter = 0,
                      max.abund.filter = 0.002, min.prop = 0,
                      # Winsorization to replace outliers
                      is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                      # Posterior sampling
                      is.post.sample = TRUE, post.sample.no = 25,
                      # Use the square-root transformation
                      link.func = list(function (x) x^0.5), stats.combine.func = max,
                      # Permutation-based multiple testing correction
                      perm.no = 999,  strata = NULL,
                      # Reference-based multiple stage normalization
                      ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                      # Family-wise error rate control
                      is.fwer = TRUE, verbose = TRUE, return.feature.dat = T)

df = data.frame(t(ZicoSeq.e4$feature.dat), check.names = FALSE)
prev = data.frame(df > 0, check.names = FALSE)
df$var = ZicoSeq.e4$meta.dat[[ZicoSeq.e4$grp.name]]
prev$var = ZicoSeq.e4$meta.dat[[ZicoSeq.e4$grp.name]]

group_means = df %>%
  mutate_if(is.numeric, ~ . / rowSums(select(df, where(is.numeric)))) %>% 
  group_by(var) %>%
  summarise_all(mean) %>%
  pivot_longer(cols=-var) %>%
  pivot_wider(names_from = var, values_from = value)

colnames(group_means) = c("Taxa", paste0(colnames(group_means)[-1], "_mean"))


group_prev = prev %>%
  group_by(var) %>%
  summarise_all(sum) %>%
  pivot_longer(cols=-var) %>%
  pivot_wider(names_from = var, values_from = value)

colnames(group_prev) = c("Taxa", paste0(colnames(group_prev)[-1], "_preval"))


zico_tab <- merge(group_means, group_prev,by="Taxa", all=TRUE)
zico_tab$R2 <- ZicoSeq.e4$R2[zico_tab$Taxa,]
zico_tab$padj <- ZicoSeq.e4$p.adj.fdr[zico_tab$Taxa]

zico_tab <- zico_tab[c("Taxa", "padj", "R2", colnames(zico_tab)[!colnames(zico_tab) %in% c("Taxa", "padj", "R2")])  ]

zico_tab_e4.sp <- tibble(zico_tab) %>% mutate(across(where(is.numeric), signif, 10))

set.seed(123)
ZicoSeq.e4 <- ZicoSeq(meta.dat = metadata_b2_e4,
                      feature.dat = as.matrix(genus_b2_e4),
                      grp.name = 'Group', adj.name = NULL, feature.dat.type = "count",
                      # Filter to remove rare taxa
                      prev.filter = 0, mean.abund.filter = 0,
                      max.abund.filter = 0.002, min.prop = 0,
                      # Winsorization to replace outliers
                      is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                      # Posterior sampling
                      is.post.sample = TRUE, post.sample.no = 25,
                      # Use the square-root transformation
                      link.func = list(function (x) x^0.5), stats.combine.func = max,
                      # Permutation-based multiple testing correction
                      perm.no = 999,  strata = NULL,
                      # Reference-based multiple stage normalization
                      ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                      # Family-wise error rate control
                      is.fwer = TRUE, verbose = TRUE, return.feature.dat = T)

df = data.frame(t(ZicoSeq.e4$feature.dat), check.names = FALSE)
prev = data.frame(df > 0, check.names = FALSE)
df$var = ZicoSeq.e4$meta.dat[[ZicoSeq.e4$grp.name]]
prev$var = ZicoSeq.e4$meta.dat[[ZicoSeq.e4$grp.name]]

group_means = df %>%
  mutate_if(is.numeric, ~ . / rowSums(select(df, where(is.numeric)))) %>% 
  group_by(var) %>%
  summarise_all(mean) %>%
  pivot_longer(cols=-var) %>%
  pivot_wider(names_from = var, values_from = value)

colnames(group_means) = c("Taxa", paste0(colnames(group_means)[-1], "_mean"))


group_prev = prev %>%
  group_by(var) %>%
  summarise_all(sum) %>%
  pivot_longer(cols=-var) %>%
  pivot_wider(names_from = var, values_from = value)

colnames(group_prev) = c("Taxa", paste0(colnames(group_prev)[-1], "_preval"))


zico_tab <- merge(group_means, group_prev,by="Taxa", all=TRUE)
zico_tab$R2 <- ZicoSeq.e4$R2[zico_tab$Taxa,]
zico_tab$padj <- ZicoSeq.e4$p.adj.fdr[zico_tab$Taxa]

zico_tab <- zico_tab[c("Taxa", "padj", "R2", colnames(zico_tab)[!colnames(zico_tab) %in% c("Taxa", "padj", "R2")])  ]

zico_tab_e4.ge <- tibble(zico_tab) %>% mutate(across(where(is.numeric), signif, 10))

## Adjust for genotype instead of splitting by genotype


set.seed(123)

species_b2 <- species_b2[rowSums(species_b2) > 0,]
ZicoSeq.adj <- ZicoSeq(meta.dat = metadata_b2,
                      feature.dat = as.matrix(species_b2),
                      grp.name = 'Group', adj.name = "genotype", feature.dat.type = "count",
                      # Filter to remove rare taxa
                      prev.filter = 0, mean.abund.filter = 0,
                      max.abund.filter = 0.002, min.prop = 0,
                      # Winsorization to replace outliers
                      is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                      # Posterior sampling
                      is.post.sample = TRUE, post.sample.no = 25,
                      # Use the square-root transformation
                      link.func = list(function (x) x^0.5), stats.combine.func = max,
                      # Permutation-based multiple testing correction
                      perm.no = 999,  strata = NULL,
                      # Reference-based multiple stage normalization
                      ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                      # Family-wise error rate control
                      is.fwer = TRUE, verbose = TRUE, return.feature.dat = T)

df = data.frame(t(ZicoSeq.adj$feature.dat), check.names = FALSE)
prev = data.frame(df > 0, check.names = FALSE)
df$var = ZicoSeq.adj$meta.dat[[ZicoSeq.adj$grp.name]]
prev$var = ZicoSeq.adj$meta.dat[[ZicoSeq.adj$grp.name]]

group_means = df %>%
  mutate_if(is.numeric, ~ . / rowSums(select(df, where(is.numeric)))) %>% 
  group_by(var) %>%
  summarise_all(mean) %>%
  pivot_longer(cols=-var) %>%
  pivot_wider(names_from = var, values_from = value)

colnames(group_means) = c("Taxa", paste0(colnames(group_means)[-1], "_mean"))


group_prev = prev %>%
  group_by(var) %>%
  summarise_all(sum) %>%
  pivot_longer(cols=-var) %>%
  pivot_wider(names_from = var, values_from = value)

colnames(group_prev) = c("Taxa", paste0(colnames(group_prev)[-1], "_preval"))


zico_tab <- merge(group_means, group_prev,by="Taxa", all=TRUE)
zico_tab$R2 <- ZicoSeq.adj$R2[zico_tab$Taxa,]
zico_tab$padj <- ZicoSeq.adj$p.adj.fdr[zico_tab$Taxa]

zico_tab <- zico_tab[c("Taxa", "padj", "R2", colnames(zico_tab)[!colnames(zico_tab) %in% c("Taxa", "padj", "R2")])  ]

zico_tab_adj.sp <- tibble(zico_tab) %>% mutate(across(where(is.numeric), signif, 10))


set.seed(123)
genus_b2 <- genus_b2[rowSums(genus_b2) > 0,]

ZicoSeq.adj <- ZicoSeq(meta.dat = metadata_b2,
                       feature.dat = as.matrix(genus_b2),
                       grp.name = 'Group', adj.name = "genotype", feature.dat.type = "count",
                       # Filter to remove rare taxa
                       prev.filter = 0, mean.abund.filter = 0,
                       max.abund.filter = 0.002, min.prop = 0,
                       # Winsorization to replace outliers
                       is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                       # Posterior sampling
                       is.post.sample = TRUE, post.sample.no = 25,
                       # Use the square-root transformation
                       link.func = list(function (x) x^0.5), stats.combine.func = max,
                       # Permutation-based multiple testing correction
                       perm.no = 999,  strata = NULL,
                       # Reference-based multiple stage normalization
                       ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                       # Family-wise error rate control
                       is.fwer = TRUE, verbose = TRUE, return.feature.dat = T)

df = data.frame(t(ZicoSeq.adj$feature.dat), check.names = FALSE)
prev = data.frame(df > 0, check.names = FALSE)
df$var = ZicoSeq.adj$meta.dat[[ZicoSeq.adj$grp.name]]
prev$var = ZicoSeq.adj$meta.dat[[ZicoSeq.adj$grp.name]]

group_means = df %>%
  mutate_if(is.numeric, ~ . / rowSums(select(df, where(is.numeric)))) %>% 
  group_by(var) %>%
  summarise_all(mean) %>%
  pivot_longer(cols=-var) %>%
  pivot_wider(names_from = var, values_from = value)

colnames(group_means) = c("Taxa", paste0(colnames(group_means)[-1], "_mean"))


group_prev = prev %>%
  group_by(var) %>%
  summarise_all(sum) %>%
  pivot_longer(cols=-var) %>%
  pivot_wider(names_from = var, values_from = value)

colnames(group_prev) = c("Taxa", paste0(colnames(group_prev)[-1], "_preval"))


zico_tab <- merge(group_means, group_prev,by="Taxa", all=TRUE)
zico_tab$R2 <- ZicoSeq.adj$R2[zico_tab$Taxa,]
zico_tab$padj <- ZicoSeq.adj$p.adj.fdr[zico_tab$Taxa]

zico_tab <- zico_tab[c("Taxa", "padj", "R2", colnames(zico_tab)[!colnames(zico_tab) %in% c("Taxa", "padj", "R2")])  ]

zico_tab_adj.ge <- tibble(zico_tab) %>% mutate(across(where(is.numeric), signif, 10))



list_of_datasets <- list("E2, species, ND vs CR" = zico_tab_e2.sp, "E2, genus, ND vs CR" = zico_tab_e2.ge, 
                         "E3, species, ND vs CR" = zico_tab_e3.sp, "E3, genus, ND vs CR" = zico_tab_e3.ge,
                         "E4, species, ND vs CR" = zico_tab_e4.sp, "E4, genus, ND vs CR" = zico_tab_e4.ge,
                         "GT adj, species, ND vs CR" = zico_tab_adj.sp, "GT adj, genus, ND vs CR" = zico_tab_adj.ge)

write.xlsx(list_of_datasets, file = "zicoseq_batch2_nd_vs_cr.xlsx")



## PCoA plots
ordinate_plot <- function(dist, meta, voi){
  obj <- cmdscale(as.dist(dist), k=2, eig=T)
  pve <- round(obj$eig[1:2]/sum(abs(obj$eig))*100, 1)
  y <- cbind(obj$points[, 1], obj$points[, 2])
  xlab <- paste0('PC1(', pve[1], '%)')
  ylab <- paste0('PC2(', pve[2], '%)')
  
  colnames(y) <- c("PC1", "PC2")
  y <- as.data.frame(y)
  y$type <- factor(meta[rownames(y), voi])
  centroids <- aggregate(cbind(PC1,PC2)~type,data=y,mean)
  y <- merge(y, centroids, by="type", suffixes=c("",".centroid"))
  mtest <- melt(y, id=c('PC1', 'PC2', 'PC1.centroid', 'PC2.centroid', measure=c('type')))
  obj <- ggplot(mtest, aes(x=PC1, y=PC2, color=type)) +
    geom_point(size=3) +
    geom_point(data=centroids,aes(x=PC1,y=PC2,color=type)) +
    geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=type)) +
    stat_ellipse() +
    theme_bw() + xlab(xlab) + ylab(ylab)
  return(obj)
}

set.seed(123)
rarefy_out <- Rarefy(t(species_b2))
species_b2.rff <- t(rarefy_out$otu.tab.rff)

b2_e2.rff <- species_b2.rff[,rownames(metadata_b2[metadata_b2$genotype=="E2",])]
b2_e2.rff <- b2_e2.rff[rowSums(b2_e2.rff) > 0,]

b2_e3.rff <- species_b2.rff[,rownames(metadata_b2[metadata_b2$genotype=="E3",])]
b2_e3.rff <- b2_e3.rff[rowSums(b2_e3.rff) > 0,]

b2_e4.rff <- species_b2.rff[,rownames(metadata_b2[metadata_b2$genotype=="E4",])]
b2_e4.rff <- b2_e4.rff[rowSums(b2_e4.rff) > 0,]


BC.dist <- vegdist(t(species_b2.rff), method="bray")

b2_e2.rff.dist <- vegdist(t(b2_e2.rff), method="bray")
b2_e3.rff.dist <- vegdist(t(b2_e3.rff), method="bray")
b2_e4.rff.dist <- vegdist(t(b2_e4.rff), method="bray")

all_gt.plt <- ordinate_plot(BC.dist, metadata_b2, "genotype")
all_diet.plt <- ordinate_plot(BC.dist, metadata_b2, "Group")
all_gender.plt <- ordinate_plot(BC.dist, metadata_b2, "gender")

b2_e2.plt <- ordinate_plot(b2_e2.rff.dist, metadata_b2_e2, "Group")
b2_e3.plt <- ordinate_plot(b2_e3.rff.dist, metadata_b2_e3, "Group")
b2_e4.plt <- ordinate_plot(b2_e4.rff.dist, metadata_b2_e4, "Group")

ggpar(all_gt.plt, palette="jco", legend.title = "Genotype") + theme_pubr() + 
  ggpar(all_diet.plt, palette="jco", legend.title = "Diet") + theme_pubr() +
  ggpar(all_gender.plt, palette="jco", legend.title = "gender") + theme_pubr() 

ggpar(b2_e2.plt, palette="jco", title = "E2", legend.title = "Diet") + theme_pubr(legend="left") + 
  ggpar(b2_e3.plt, palette="jco", title = "E3", legend.title = "Diet") + theme_pubr(legend="none") +
  ggpar(b2_e4.plt, palette="jco", title = "E4", legend.title = "Diet") + theme_pubr(legend="none")


metadata.tmp <- metadata_b2
metadata.tmp$Group2 <- paste0(metadata_b2$genotype, "-", metadata_b2$Group)

all_gt_diet.plt <- ordinate_plot(BC.dist, metadata.tmp, "Group2")
ggpar(all_gt_diet.plt, palette="jco", legend.title = "Group") + theme_pubr()




