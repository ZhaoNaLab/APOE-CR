source('Analysis/util.r')
library(arsenal)
library(rrcov)
source('/research/labs/moleneurosci/bug/projects/Plotting/single-cell/heatmap.R')

col.data = colData(x) %>% as_tibble

# distribution of zeros
df = assay(x)
table(rowSums(df == 0))
df[rowSums(df == 0) > 0,]

p = col.data %>%
  select(lipidomics.ID, Mouse.weight : CFC.cued, genotype.diet, sex) %>%
  mutate(across(Mouse.weight : CFC.cued, ~as.numeric(.x))) %>%
  pivot_longer( Mouse.weight : CFC.cued, names_to = 'Characteristic') %>%
  mutate(
    value = ifelse(is.na(as.numeric(value)), NA, as.numeric(value))
  ) %>%
  ggplot(aes(x=genotype.diet, y=value)) +
  geom_jitter(aes(fill=sex), shape=21, size=2, width=0.1) + 
  facet_wrap(~Characteristic, scales='free_y') +
  theme(text = element_text(size = 10))
cowplot::save_plot(
  plot = p,
  filename = file.path("Results", 'eda', 'eda - boxplots.png'),
  base_height = 10,
  dpi = 400
)

h = col.data %>%
  select(Mouse.weight : CFC.cued) %>%
  select(where(is.double)) %>%
  drop_na() %>%
  cor %>%
  as.matrix %>%
  Heatmap
cowplot::save_plot(
  plot = draw_heatmap(h),
  filename = file.path('Results', 'eda', 'eda - correlation heatmap.png'),
  base_height=8)

by_covariate = col.data %>%
  select(Mouse.weight : CFC.cued, sex, Group, lipidomics.ID, genotype.diet) %>%
  select(where(is.double), sex, Group, lipidomics.ID, genotype.diet) %>%
  drop_na() %>%
  pivot_longer(Mouse.weight : OFA.distance, names_to = 'covariate') %>%
  group_by(covariate) %>%
  nest()

covariate_model <- function(df) {
  lm(value ~ sex + Group, data=df)
}

regressions = by_covariate %>%
  group_by(covariate) %>%
  mutate(
    model = purrr::map(data, ~lm(value ~ sex, data=.x)),
    glanced = purrr::map(model, broom::glance),
    augmented = purrr::map(model, broom::augment, se_fit=TRUE, interval='confidence'),
    tidied = purrr::map(model, broom::tidy)
  )

p = regressions %>%
  select(augmented, data) %>% unnest() %>%
  ggplot(aes(x=genotype.diet, y=.resid)) +
  geom_jitter(aes(fill=sex), shape=21, size=2, width=0.1) + 
  facet_wrap(~covariate, scales='free_y') +
  theme(text = element_text(size = 10))
cowplot::save_plot(
  plot = p,
  filename = file.path("Results", 'eda', 'eda - boxplots sex adjusted.png'),
  base_height = 10,
  dpi = 400
)

h = regressions %>%
  select(data, augmented) %>% unnest() %>%
  select(.std.resid, lipidomics.ID) %>%
  pivot_wider(names_from=covariate, values_from=.std.resid, id_cols=lipidomics.ID) %>%
  select(-lipidomics.ID) %>%
  cor %>%
  as.matrix %>%
  Heatmap()
cowplot::save_plot(
  plot = draw_heatmap(h),
  filename = file.path('Results', 'eda', 'eda - correlation heatmap sex adjusted.png'),
  base_height=8)

x = readRDS("Data/speciesSE.rds")
row.data = rowData(x) %>% as.data.frame
keep = !rowSums(assay(x) == 0)
summary(keep)
fit = lm(t(log(assay(x)[keep,])) ~ sex + Group, data=colData(x))
lipid.abbrs = rowData(x)[keep,'Lipid.Abbr']
corr = cor(fit$residuals)
corr = cor(t(assay(x, 'lognorm')[keep,]))

keep = row.data[rownames(corr),]$Lipid.Abbr %in% c('CL', 'CAR')

dendro = hclust(dist(corr[keep,keep]))
order = cutree(dendro, h=4)

tmp = arsenal::tableby(order ~ Lipid.Abbr, data = data.frame(Lipid.Abbr = lipid.abbrs[keep]))
summary(tmp, text=TRUE)

h = Heatmap(
  corr[keep,keep], 
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  row_split = order, #lipid.abbrs[keep],
  column_split = order, #lipid.abbrs[keep],
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_column_slices = TRUE,
  row_title_rot = 0,
  heatmap_legend_param = list(
    title = 'Normalized gene expression',
    direction = 'vertical',
    legend_height = unit(3, 'cm'),
    at = c(-1, 0, 1)
  ),
  rect_gp = gpar(col='white', lwd=1)
)

#### Plot PCA Pairs ####
x = readRDS("Data/serum_speciesSE.rds")
plot_pca_pairs(x, 'prop.of.class')
plot_pca_pairs(x, 'prop.of.class', log=TRUE)
plot_pca_pairs(x, 'lognorm')

#x = x[,colnames(x) != 'BU-LF-L22']
x = x[,colnames(x) != 'Q-LF-SLP2']
plot_pca_pairs(x, 'prop.of.class', 'pairs-wo-2')
plot_pca_pairs(x, 'prop.of.class', 'pairs-wo-2', log=TRUE)
plot_pca_pairs(x, 'lognorm', 'pairs-wo-2')

x = readRDS('Data/serum_classSE.rds')
plot_pca_pairs(x, 'lipid.class')
plot_pca_pairs(x, 'lipid.class.prop')
plot_pca_pairs(x, 'lipid.class', log=TRUE)
plot_pca_pairs(x, 'lipid.class.prop', log=TRUE)

#x = x[,colnames(x) != 'BU-LF-L22']
x = x[,colnames(x) != 'Q-LF-SLP2']
plot_pca_pairs(x, 'lipid.class', 'pairs-wo-2')
plot_pca_pairs(x, 'lipid.class.prop', 'pairs-wo-2')
plot_pca_pairs(x, 'lipid.class', 'pairs-wo-2', log=TRUE)
plot_pca_pairs(x, 'lipid.class.prop', 'pairs-wo-2', log=TRUE)

##### PCA bi-plots #####
x = readRDS("Data/serum_speciesSE.rds")
plot_pca_biplot(x, 'prop.of.class')
plot_pca_biplot(x, 'prop.of.class', log=TRUE)
plot_pca_biplot(x, 'lognorm')

#x = x[,colnames(x) != 'BU-LF-L22']
x = x[,colnames(x) != 'Q-LF-SLP2']
plot_pca_biplot(x, 'prop.of.class', 'biplots-wo-2')
plot_pca_biplot(x, 'prop.of.class', 'biplots-wo-2', log=TRUE)
plot_pca_biplot(x, 'lognorm', 'biplots-wo-2')

x = readRDS('Data/serum_classSE.rds')
plot_pca_biplot(x, 'lipid.class')
plot_pca_biplot(x, 'lipid.class.prop')
plot_pca_biplot(x, 'lipid.class', log=TRUE)
plot_pca_biplot(x, 'lipid.class.prop', log=TRUE)

#x = x[,colnames(x) != 'BU-LF-L22']
#plot_pca_biplot(x, 'lipid.class', 'biplots-wo-22')
#plot_pca_biplot(x, 'lipid.class.prop', 'biplots-wo-22')
#plot_pca_biplot(x, 'lipid.class', 'biplots-wo-22', log=TRUE)
#plot_pca_biplot(x, 'lipid.class.prop', 'biplots-wo-22', log=TRUE)

x = readRDS("Data/classSE.rds")
pca = PcaClassic(t(log(assay(x, 'lipid.class'))), k=10, scale = TRUE)
summary(pca)
plot(pca)

bind_rows(pca$loadings %>%
  as.data.frame %>%
  arrange(desc(PC5)) %>% 
  head(5),
pca$loadings %>%
  as.data.frame %>%
  arrange(desc(PC5)) %>% 
  tail(5)
)

pca$loadings %>%
  as.data.frame %>%
  arrange(desc(PC4 * PC5))

biplot(pca, c(4,5))

x = readRDS("Data/speciesSE.rds")
keep = ! rowSums(assay(x) == 0)
pca = PcaClassic(t(log(assay(x, 'prop.of.class')[keep,])), k=10, scale = TRUE)
#pca = PcaClassic(t(assay(x, 'lipid.class.prop')), k=10)
summary(pca)
plot(pca)

pca$loadings['PE D18:1-22:6',]

bind_rows(pca$loadings %>%
  as.data.frame %>%
  arrange(desc(PC4)) %>% 
  head(5),
pca$loadings %>%
  as.data.frame %>%
  arrange(desc(PC4)) %>% 
  tail(5)
)

ggpairs(
  t(assay(x, 'lognorm'))[,
      rownames(x) %in% names(order)[order %in% c(1,2)] &
        (rowData(x)$Lipid.Abbr == 'CL' | rownames(x) == 'CAR 18:0')
    ] %>% as.data.frame
)

fit = lm(t(sqrt(assay(x))) ~ sex + Group + sex:Group + genotype + sex:genotype, data=colData(x))
assay(x,'resids') = t(resid(fit))
plot_species(x,'prop.of.class') %>%
  filter(Lipid.Abbr == 'ST') %>%
  ggplot(aes(x = interaction(genotype, Group), y = value, color=sex)) + 
  geom_quasirandom(width=0.2) + 
  facet_wrap(~Lipid.Species, scales = 'free_y')

dat = t(assay(x, 'lipid'))[,rownames(x) %in% names(order)[order %in% c(1,2)] &
                            (rowData(x)$Lipid.Abbr == 'CL' | rownames(x) == 'CAR 18:0')]
aggregate(dat ~ genotype + Group, data = colData(x), FUN="mean")

dat = t(assay(x, 'lipid'))[, (rowData(x)$Lipid.Abbr == 'CL')]
aggregate(dat ~ genotype + Group, data = colData(x), FUN="mean")

plot_species_pairs(c('CL', 'CAR'), 'Group', 
                   keep = rownames(x) %in% names(order)[order %in% c(1,2)] # & 
                     )#(rowData(x)$Lipid.Abbr == 'CL' | rownames(x) == 'CAR 18:0'))

t(log(assay(x)[]))

assay(x, 'lognorm')[rowData(x)$Lipid.Abbr %in% c('CBS'),] %>%
  as.data.frame() %>% rownames_to_column('Lipid.Species') %>%
  pivot_longer(-Lipid.Species, names_to='SampleID') %>%
  left_join(colData(x) %>% as.data.frame %>% rownames_to_column('SampleID')) %>%
  left_join(rowData(x) %>% as.data.frame %>% rownames_to_column('Lipid.Species')) %>%
  ggplot(aes(
    x = interaction(Group, genotype),
    y = value,
    fill = sex)
  ) + 
  geom_jitter(
    shape=21,
    size=1,
    width=0.1
  ) + 
  facet_wrap(~Lipid.Species, scales='free_y') + 
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))

library(rrcov)
x = readRDS("Data/serum_classSE.rds")
pca = PcaClassic(x = t(assay(x, 'lipid.class')))
plot(pca)
cbind(col.data, pca$scores) %>%
  ggplot(aes(PC5, PC6, fill=sex)) + 
  geom_point(size=3, shape=21) + 
  #geom_text(aes(label = lipidomics.ID), nudge_x = 0.001, nudge_y = 0.001) + 
  #scale_fill_viridis_c(option='inferno') + 
  xlab('PC1 (59.3%)') + 
  ylab('PC2 (27.6%)') 
  #scale_fill_brewer(palette='Paired') #+ 
  #scale_fill_manual(values = genotype.diet.colors) + 

pca$loadings[,'PC3'] # CE

p = assay(x, 'lipid.class.prop') %>%
  as.data.frame %>%
  rownames_to_column('lipid.class') %>%
  pivot_longer(-lipid.class, names_to='id') %>%
  mutate(id = str_remove(id, 'BU-LF-L')) %>%
  ggplot(aes(x=lipid.class, y=value)) + 
  geom_jitter(width=0.1) + 
  geom_text(aes(label = id)) + 
  facet_wrap(~lipid.class, scales='free')
cowplot::save_plot(
  plot = p,
  filename = 'Results/eda/eda - class distributions.png',
  base_height=8
)


# Test correlation of lipids with covariates
x = readRDS("Data/serum_speciesSE.rds")
col_data = col.data %>%
  select(Mouse.weight : Pancreat) %>%
  mutate(Stomach.weight = as.numeric(Stomach.weight)) %>%
  mutate(Stomach.weight = ifelse(is.na(Stomach.weight), median(Stomach.weight, na.rm=TRUE), Stomach.weight))
df = assay(x, 'prop.of.class')
fit = lm(t(df) ~ sex + genotype + diet, data=col.data)
Heatmap(cor(as.matrix(col_data[,c('Mouse.weight', 
                          'Brain.weight',
                          'Stomach.weight',
                          'spleen.weight',
                          'Pancreat'
                          )]), 
    resid(fit)))

cor.test(as.matrix(assay(x)), col.data$Mouse.weight)


##### Visualize the distribution of lipid classes #####
x = readRDS("Data/classSE.rds")
col.data = colData(x) %>% as.data.frame %>% set_levels()
df = assay(x, 'lipid.class.prop') %>% as.matrix
df = df[rownames(df) != 'TC',]

t(df) %>% 
  as.data.frame %>%
  bind_cols(col.data) %>%
  pivot_longer(BMP : ST, names_to='Lipid.Abbr') %>%
  filter(lipidomics.ID != '22') %>%
  filter(! Lipid.Abbr %in% c('FC', 'PE', 'PC')) %>%
  mutate(Lipid.Abbr = fct_reorder(as.factor(Lipid.Abbr), .x=value, .fun=sum)) %>%
  ggplot(
    aes(x = interaction(lipidomics.ID, Group, genotype),
        y = value,
        fill = Lipid.Abbr)) + 
  geom_col(position='fill') + 
  theme(axis.text.x = element_blank()) +
  coord_flip() + 
  scale_fill_manual(values = as.vector(Lipid.Abbr.colors[!Lipid.Abbr.colors %in% c('TC')]))

df = assay(x)
df = df[rownames(df) != 'TC',] # TC is already FC + cholesterol ester
df = df[,colnames(df) != 'BU-LF-22']
dr_data = DirichletReg::DR_data(t(df), base=12) # Remove PC -- highly abundant phospholipid of the cell membrane that I don't think should change much between experimental groups.

fit0 = DirichletReg::DirichReg(dr_data ~ sex + Group + Brain.weight | sex + Group + Brain.weight, data = col.data, model = 'alternative')
fit1 = DirichletReg::DirichReg(dr_data ~ sex + Group | sex + Group, data = col.data, model = 'alternative')
fit2 = DirichletReg::DirichReg(dr_data ~ sex + Group + genotype | sex + Group + genotype, data = col.data, model = 'alternative')
fit3 = DirichletReg::DirichReg(dr_data ~ sex + Group + genotype + Group:genotype | sex + Group + genotype + Group:genotype, data = col.data, model = 'alternative')
fit4 = DirichletReg::DirichReg(dr_data ~ sex * Group * genotype | sex * Group * genotype, data = col.data, model = 'alternative')
fit5 = DirichletReg::DirichReg(dr_data ~ sex * Group * genotype | sex + Group * genotype, data = col.data, model = 'alternative')
fit6 = DirichletReg::DirichReg(dr_data ~ sex * Group * genotype | sex + Group + genotype, data = col.data, model = 'alternative')
anova(fit1,fit2,fit3,fit4,fit5,fit6) #  rule out fit1-2
anova(fit3, fit4, fit5, fit6) # use fit3 or fit4
AIC(fit3)
AIC(fit4)
AIC(fit5)
AIC(fit6)

s = summary(fit3)
coefs = as.matrix(s$coef.mat)
term = rownames(coefs)
coef = str_remove(names(s$coefficients), ':.*$')
coefs = bind_cols(
  coef = coef,
  term = term,
  coefs,
)
coefs <- coefs %>%
  filter(coef != '(phi)') %>%
  group_by(term) %>%
  mutate(adj.p.val = p.adjust(`Pr(>|z|)`, method='bonferroni'))

coefs %>%
  filter(term == 'GroupCR:genotypeE2') %>%
  arrange(`Pr(>|z|)`)

coefs %>%
  filter(term == 'genotypeE2') %>%
  arrange(`Pr(>|z|)`)

#coefs %>%
#  filter(coef == 'CE')

dr_data %>% as.matrix %>%
  bind_cols(as.data.frame(colData(x))) %>%
  pivot_longer(BMP : ST, names_to='Lipid.Abbr') %>%
  filter(Lipid.Abbr == 'CE' & genotype %in% c('E2', 'E3', 'E4')) %>%
  ggplot(
    aes(x = genotype,
        y = value,
        fill = sex)
  ) + 
  geom_quasirandom(width=0.2, shape=21) + 
  #scale_fill_manual(values = Lipid.Abbr.colors) + 
  facet_wrap(~Group, scales='free_y') + 
  theme(axis.text.x = element_text(angle = 90, hjust=0.5, vjust=0.5)) + 
  scale_y_continuous(trans='log10')

fit4.fitted = fitted(fit4) %>% as.matrix %>% as.data.frame
rownames(fit4.fitted) = colnames(df)
fit4.fitted %>%
  rownames_to_column('SampleID') %>%
  pivot_longer(-SampleID, names_to = 'Lipid.Abbr', values_to='fitted') %>%
  left_join(dr_data %>% as.matrix %>%
              bind_cols(as.data.frame(colData(x))) %>%
              rownames_to_column('SampleID') %>%
              pivot_longer(BMP : ST, names_to='Lipid.Abbr', values_to='value'),
            by = c('Lipid.Abbr', "SampleID")) %>%
  filter(Lipid.Abbr == 'CE') %>%
  ggplot() + 
  geom_text(
    mapping = aes(
      x = interaction(Group, genotype),
      y = value,
      label = SampleID,
      #fill = sex
    ),
    #width=0.2, 
    #shape=21
  ) + 
  geom_point(
    data = . %>% distinct(SampleID, Lipid.Abbr, fitted, .keep_all=TRUE),
    mapping = aes(x = interaction(Group, genotype),
                  y = fitted,
                  fill = sex),
    shape=25
  ) + 
  #scale_fill_manual(values = Lipid.Abbr.colors) + 
  facet_wrap(~sex, scales='free_y') + 
  theme(axis.text.x = element_text(angle = 90, hjust=0.5, vjust=0.5)) + 
  scale_y_continuous(trans='log10')

ci = confint(fit4)$ci[[1]]
ci$CE
 
fit = lm(t(log(assay(x)[keep,])) ~ sex + Group + genotype + Group:genotype, data=colData(x))
#design = model.matrix(t(log(assay(x))) ~ sex + Group + genotype + Group:genotype + Group:sex + genotype:sex, data = colData(x)) 
#fitted = design[,c(1:3,9)] %*% coef(fit)[c(1:3,9),]
assay(x, 'resids') = t(fit$residuals)
assay(x, 'fitted') = t(fitted(fit))

p = plot_class(x = readRDS("Data/serum_classSE.rds"), assay='lipid.class') %>%
  ggplot(
    aes(x=interaction(Group, genotype, sep = ' '), 
        y=log(value), 
        fill=interaction(Group, genotype, sep=' '),
        shape = sex)
    ) + 
  geom_quasirandom(width=0.1, size=2) + 
  facet_wrap(~Lipid.Abbr, scales='free_y') + 
  scale_shape_manual(values=c(21,25)) + 
  scale_fill_manual(values = as.vector(genotype.diet.colors)) + 
  guides(fill = guide_legend(title = 'Diet Genotype',
                             override.aes = list(shape = 21))) + 
  ylab('log(sum of class)')
cowplot::save_plot(
  plot = p,
  filename = file.path('Results', 'eda', 'Class', 'quasi-log', 'classes.png'),
  base_height=10,
  base_asp=2
)

plot_class(x = readRDS("Data/serum_classSE.rds"), assay = 'lipid.class') %>%
  mutate(Lipid.Abbr = as.character(Lipid.Abbr)) %>%
  group_by(Lipid.Abbr) %>%
  group_map(.f=~{
    p = .x %>%
      ggplot(aes(x=genotype, y=log(value), fill=genotype)) + 
      geom_quasirandom(width=0.1, shape=21) + 
      facet_wrap(~Group) + 
      labs(title = paste(.y, 'log(sum class)'))
    cowplot::save_plot(plot=p,
                       filename=file.path('Results', 'eda', 'Class', 'quasi-log', paste0(.y, '.png')), 
                       base_height=5)
  })

plot_species(readRDS("Data/serum_speciesSE.rds"), assay='lognorm') %>%
  mutate(Lipid.Abbr = as.character(Lipid.Abbr)) %>%
  group_by(Lipid.Abbr) %>%
  group_walk(.f=~{
    p = .x %>%
      ggplot(aes(x=interaction(Group, genotype, sep=' '), y=value, fill=interaction(genotype, Group, sep=' '), shape=sex)) + 
      geom_quasirandom(width=0.1, size=3) + 
      scale_shape_manual(values=c(21,25)) + 
      guides(fill = guide_legend(title='Genotype Group', override.aes = list(shape = c(21)))) + 
      scale_fill_manual(values=genotype.diet.colors) + 
      facet_wrap(~Lipid.Species, scales='free_y') + 
      theme(axis.text.x = element_text(angle=90, hjust=0.5)) + 
      labs(title = paste('Lipid species - lognorm', .y))
    cowplot::save_plot(plot=p,
                       filename=paste0('Results/eda/Species/quasi-lognorm/', .y, '.png'), 
                       base_height=10)
  })

plot_species(readRDS("Data/serum_speciesSE.rds"), assay='prop.of.class') %>%
  #filter(Lipid.Abbr == 'OX') %>%
  #filter(SampleID != 'BU-LF-L22') %>%
  mutate(Lipid.Abbr = as.character(Lipid.Abbr)) %>%
  group_by(Lipid.Abbr) %>%
  group_walk(.f=~{
    p = .x %>%
      ggplot(aes(x=interaction(Group, genotype, sep=' '), y=value, fill=interaction(genotype, Group, sep=' '), shape=sex)) + 
      geom_quasirandom(width=0.1, size=3) + 
      scale_shape_manual(values=c(21,25)) + 
      guides(fill = guide_legend(title='Genotype Group', override.aes = list(shape = c(21)))) + 
      scale_fill_manual(values=genotype.diet.colors) + 
      facet_wrap(~Lipid.Species, scales='free_y') + 
      theme(axis.text.x = element_text(angle=90, hjust=0.5)) + 
      labs(title = paste('Lipid species - prop.of.class', .y))
    cowplot::save_plot(plot=p,
                       filename=paste0('Results/eda/Species/quasi-prop.of.class/', .y, '.png'), 
                       base_height=10)
  })

p = plot_species(x=readRDS("Data/serum_speciesSE.rds"), assay = 'lipid') %>%
  mutate(genotype = relevel(genotype, 'E2')) %>%
  ggplot(aes(x = interaction(lipidomics.ID, sex),
             y = value, 
             fill = fct_reorder(Lipid.Abbr, value, sum))) + 
  geom_col(position='stack') + 
  coord_flip() + 
  facet_wrap(~interaction(genotype, Group, sep=' '), scales='free_y', ncol = 1, switch = TRUE) + 
  scale_fill_manual(values=Lipid.Abbr.colors) + 
  guides(fill = guide_legend(title='Lipid Class')) + 
  ylab('Mouse')
cowplot::save_plot(
  plot = p,
  filename = file.path('Results', 'eda', 'Class', 'Bars', 'stack classes.png'),
  base_height=8,
  base_asp=1.6
)

p = plot_species(x=readRDS("Data/serum_speciesSE.rds"), assay = 'lipid') %>%
  mutate(genotype = relevel(genotype, 'E2')) %>%
  ggplot(aes(x = interaction(lipidomics.ID, sex),
             y = value, 
             fill = fct_reorder(Lipid.Abbr, value, sum))) + 
  geom_col(position='fill') + 
  coord_flip() + 
  facet_wrap(~interaction(genotype, Group, sep=' '), scales='free_y', ncol = 1, switch = TRUE) + 
  scale_fill_manual(values=Lipid.Abbr.colors) + 
  guides(fill = guide_legend(title='Lipid Class')) + 
  ylab('Mouse')
cowplot::save_plot(
  plot = p,
  filename = file.path('Results', 'eda', 'Class', 'Bars', 'fill classes.png'),
  base_height=8,
  base_asp=1.6
)

# Dendrograms
pdf(file = 'Results/eda/Species/dendro/prop.of.class.pdf', width = 60, height=20)
dendro = plot_dendro(x = readRDS("Data/serum_speciesSE.rds"), assay = 'prop.of.class')
dev.off()

pdf(file = 'Results/eda/Species/dendro/lognorm.pdf', width = 60, height=20)
dendro = plot_dendro(x = readRDS("Data/serum_speciesSE.rds"), assay = 'lognorm')
dev.off()

pdf(file = 'Results/eda/Class/dendro/log lipid.class dendro.pdf', width = 20, height=10)
dendro = plot_dendro(x = readRDS("Data/serum_classSE.rds"), assay = 'lipid.class', log=TRUE)
dev.off()

#order = cutree(dendro, h=5.5)
#tmp = arsenal::tableby(Lipid.Abbr ~ factor(order),
#                 data = bind_cols(
#                   row.data[names(order),],
#                   order = order
#                 ))
#summary(tmp, text=TRUE)

# Heatmaps
# Species
h = plot_heatmap(readRDS("Data/serum_speciesSE.rds"), 'lognorm')
cowplot::save_plot(plot = draw_heatmap(h), 
                   filename = 'Results/eda/Species/correlations-lognorm/all-species.png',
                   base_height=40,
                   base_asp=1.2)
h = plot_heatmap(readRDS("Data/serum_speciesSE.rds"), 'prop.of.class')
cowplot::save_plot(plot = draw_heatmap(h), 
                   filename = 'Results/eda/Species/correlations-prop.of.class/all-species.png',
                   base_height=40,
                   base_asp=1.2)

for (abbr in unique(rowData(x)$Lipid.Abbr)) {
  h = plot_heatmap(x, 'lognorm', 
                   Lipid.Abbr = abbr, 
                   order = TRUE,
                   h=3)
  cowplot::save_plot(plot = draw_heatmap(h),
                     filename = paste0('Results/eda/Species/correlations-lognorm/', abbr, '.png'),
                     base_height=8,
                     base_asp=1.2)
}
for (abbr in unique(rowData(x)$Lipid.Abbr)) {
  h = plot_heatmap(x, 'prop.of.class', 
                   Lipid.Abbr = abbr, 
                   order = TRUE,
                   h=3)
  cowplot::save_plot(plot = draw_heatmap(h),
                     filename = paste0('Results/eda/Species/correlations-prop.of.class/', abbr, '.png'),
                     base_height=8,
                     base_asp=1.2)
}
# Class
h = plot_heatmap(x = readRDS("Data/serum_classSE.rds"), 'lipid.class', order = TRUE, h=2, show_names = TRUE)
cowplot::save_plot(plot = draw_heatmap(h),
                   filename = 'Results/eda/Class/correlations/lipid.class correlations.png',
                   base_height=5,
                   base_asp=1.2)
# Class
h = plot_heatmap(x = readRDS("Data/serum_classSE.rds"), 'lipid.class.prop', order = TRUE, h=2.25, show_names = TRUE)
cowplot::save_plot(plot = draw_heatmap(h),
                   filename = 'Results/eda/Class/correlations/lipid.class.prop correlations.png',
                   base_height=5,
                   base_asp=1.2)

# Visualize which lipids appear to change the most. 
# This seems to be CE lipids. 
plot_class(x[,x$`lipidomics-ID` != '22'], assay = 'lipid.class') %>%
  #pivot_wider(names_from='Lipid.Abbr', values_from='value') %>%
  ggplot() + 
  geom_boxplot(aes(x = Group, y = value, fill=Group)) + 
  geom_quasirandom(
    #data = plot_class(x[,x$`lipidomics-ID` != '22'], assay = 'lipid.class'),
    aes(x = Group, y = value, fill=Group), width=0.1, shape=21
  ) + 
  facet_wrap(~Lipid.Abbr, scales='free')
  
  #geom_density(aes(x = log(value))) +  
  
  #geom_point(aes(x = Brain.weight, y=value, fill=Group), shape=21) + 
  ##geom_path(aes(x = Brain.weight, y=log(value), group=genotype), linetype='dashed') + 
  #geom_smooth(aes(x = Brain.weight, y=value, color=Group), method='lm', alpha=0.1, linetype='dashed') + 
  #facet_wrap(~Lipid.Abbr, scales='free')
  

  ggplot(aes(x = , y = `Lyso CL`)) + 
  geom_point(aes(fill=Group), shape=21) + 
  geom_smooth(aes(color=Group), method='lm', linetype='dashed')
  