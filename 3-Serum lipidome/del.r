source('Analysis/util.r')
source("Analysis/Normalizations.r")
source('../../projects/Plotting/single-cell/heatmap.R')

x = readRDS("Data/serum_speciesSE.rds")
lipid.assay.longer = assay(x, "prop.of.class") %>% as.data.frame %>% rownames_to_column("Lipid.Species") %>%
  pivot_longer(-Lipid.Species, names_to='ID') %>%
  left_join(rowData(x) %>% as.data.frame %>% rownames_to_column('Lipid.Species')) %>%
  left_join(colData(x) %>% as.data.frame %>% rownames_to_column('ID')) %>%
  filter(!Lipid.Abbr %in% c('TC', 'FC')) %>%
  mutate(diet = Group)

#### ####
#### Visualization of lipid - lipid correlations ####
dat = assay(x, 'lipid')[!rownames(x) %in% c('TC', 'FC'),]
corr = cor(t(dat), method='kendall')
column_split = rowData(x)$Lipid.Abbr[!rownames(x) %in% c("TC", "FC")]
h = Heatmap(as.matrix(corr),
        row_split = column_split, 
        column_split=column_split, 
        row_labels = rep('', nrow(corr)), 
        column_labels=rep('', ncol(corr)))
cowplot::save_plot(plot=draw_heatmap(h), filename='Plots/lipid correlation.png', base_height=20)

#### Visualization of lipid class - lipid class correlations ####
lipid.assay.longer %>%
  dplyr::select(value, diet, genotype, Lipid.Abbr, ID) %>%
  pivot_wider(values_from=value, names_from=Lipid.Abbr) %>%
  ggpairs(columns=Lipid.Abbr.levels[!Lipid.Abbr.levels %in% c("TC", "FC")], mapping=aes(color=diet))
p
cowplot::save_plot(plot=p, filename='Plots/Lipid class ggpairs.png', base_height=20, base_asp=1.3)

#### View sums of values, summing over class ####
lipid.assay.sum = assay(xb, 'lipid.class') %>% 
  t() %>% as_tibble(rownames = "ID") %>% 
  pivot_longer(-ID, names_to='Lipid.Abbr') %>%
  left_join(colData(xb) %>% as_tibble(rownames='ID'), by='ID') %>%
  mutate(group = paste(genotype, Group, sex) %>% as.factor) %>%
  mutate(ID.group = factor(ID.group, ID.group.levels),
         Lipid.Abbr = factor(Lipid.Abbr, Lipid.Abbr.levels))

p = lipid.assay.sum %>%
  ggplot() + 
  geom_col(aes(x=ID.group, y=value, fill=Lipid.Abbr, group=Lipid.Abbr), color='black', position = 'stack') + 
  scale_fill_manual(values=Lipid.Abbr.colors[levels(lipid.assay.sum$Lipid.Abbr)]) + 
  theme(axis.text.x = element_text(angle=45, hjust=1)) #+ 
cowplot::save_plot(plot=p, filename='Plots/Lipids summed by class, by ID-group - col.png', base_height=8)

p = lipid.assay.sum %>%
  ggplot() + 
  geom_col(aes(x=group, y=value, fill=Lipid.Abbr, group=Lipid.Abbr), color='black', position = 'stack') +
  scale_fill_manual(values=Lipid.Abbr.colors[levels(lipid.assay.sum$Lipid.Abbr)]) + 
  theme(axis.text.x = element_text(angle=45, hjust=1)) #+ 
cowplot::save_plot(plot=p, filename='Plots/Lipids summed by class, by group - col.png', base_height=8)

p = lipid.assay.sum %>%
  ggplot() + 
  geom_smooth(aes(x=group, y=value, color=Lipid.Abbr, group=Lipid.Abbr), alpha=0.5, fill='gray70', linetype='dashed', size=0.5, show.legend=FALSE) + 
  geom_quasirandom(aes(x=group, y=value, fill=Lipid.Abbr), width=0.1, size=2, shape=21) + 
  scale_color_manual(values=Lipid.Abbr.colors) + 
  scale_fill_manual(values=Lipid.Abbr.colors[levels(lipid.assay.sum$Lipid.Abbr)]) + 
  theme(axis.text.x = element_text(angle=45, hjust=1)) + 
  scale_y_log10()
cowplot::save_plot(plot=p, filename='Plots/Lipid summed by class, by group - quasirandom.png', base_height=8)

p = lipid.assay.sum %>%
  filter(!Lipid.Abbr %in% c('TC', 'FC', 'PC', 'PE')) %>%
  #filter(!Lipid.Abbr %in% c('PS', 'CBS', 'CL', 'PI', 'ST', 'SM', 'CE')) %>%
  #filter(!Lipid.Abbr %in% c('CER', 'PA', 'LPE', 'LPC', 'PG')) %>%
  ggplot() + 
  geom_col(aes(x=ID.group, y=value, fill=Lipid.Abbr), color='gray20') + 
  scale_fill_manual(values=Lipid.Abbr.colors[levels(lipid.assay.sum$Lipid.Abbr)]) + 
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x='Sample', y='Lipid summed by class (nmol / pg protein)')
cowplot::save_plot(plot=p, filename='Plots/Lipid class compositions.png', base_height=10)

##### Linear models ####
library(modelr)
#### --- Model Lipid Class ---- ####
x = readRDS("Data/serum_classSE.rds")

x = x[, colnames(x) != 'Q-LF-SLP2']
all_glanced = tibble()
#all_glanced = readxl::read_excel('Results/glance all normalizations.xlsx', sheet = 'class')

regressions = run_lm(x, 'lipid.class', 'Lipid.Class', formula = 'sex + genotype + diet + genotype:diet', trans=identity)
all_glanced <- all_glanced %>%
  bind_rows(regressions %>% select(assay, group, glanced, formula, trans) %>% unnest(glanced))

regressions = run_lm(x, 'lipid.class', 'Lipid.Class', formula = 'sex + genotype + diet + genotype:diet', trans=log)
all_glanced <- all_glanced %>%
  bind_rows(regressions %>% select(assay, group, glanced, formula, trans) %>% unnest(glanced))

regressions = run_lm(x, 'lipid.class', 'Lipid.Class', formula = 'sex + genotype + diet + genotype:diet + sex:diet + genotype:sex + genotype:diet:sex', trans=identity)
all_glanced <- all_glanced %>%
  bind_rows(regressions %>% select(assay, group, glanced, formula, trans) %>% unnest(glanced))

regressions = run_lm(x, 'lipid.class', 'Lipid.Class', formula = 'sex + genotype + diet + genotype:diet + sex:diet + genotype:sex + genotype:diet:sex', trans=log)
all_glanced <- all_glanced %>%
  bind_rows(regressions %>% select(assay, group, glanced, formula, trans) %>% unnest(glanced))

regressions = run_lm(x, 'lipid.class.prop', 'Lipid.Class', formula = 'sex + genotype + diet + genotype:diet', trans=identity)
all_glanced <- all_glanced %>%
  bind_rows(regressions %>% select(assay, group, glanced, formula, trans) %>% unnest(glanced))

regressions = run_lm(x, 'lipid.class.prop', 'Lipid.Class', formula = 'sex + genotype + diet + genotype:diet', trans=log)
all_glanced <- all_glanced %>%
  bind_rows(regressions %>% select(assay, group, glanced, formula, trans) %>% unnest(glanced))

regressions = run_lm(x, 'lipid.class.prop', 'Lipid.Class', formula = 'sex + genotype + diet + genotype:diet + sex:diet + genotype:sex + genotype:diet:sex', trans=identity)
all_glanced <- all_glanced %>%
  bind_rows(regressions %>% select(assay, group, glanced, formula, trans) %>% unnest(glanced))

regressions = run_lm(x, 'lipid.class.prop', 'Lipid.Class', formula = 'sex + genotype + diet + genotype:diet + sex:diet + genotype:sex + genotype:diet:sex', trans=log)
all_glanced <- all_glanced %>%
  bind_rows(regressions %>% select(assay, group, glanced, formula, trans) %>% unnest(glanced))

regressions = run_lm(x, 'lipid.class', 'Lipid.Class', formula = 'sex + genotype + diet + genotype:diet + sex:diet + genotype:sex + genotype:diet:sex', trans=log)
all_glanced <- all_glanced %>%
  bind_rows(regressions %>% select(assay, group, glanced, formula, trans) %>% unnest(glanced))

regressions = run_lm(x, 'lipid.class.prop', 'Lipid.Class', formula = 'sex + genotype + diet + genotype:diet + sex:diet + genotype:sex + genotype:diet:sex', trans=log)
all_glanced <- all_glanced %>%
  bind_rows(regressions %>% select(assay, group, glanced, formula, trans) %>% unnest(glanced))

#### --- Visualize glanced --- ####

p = all_glanced %>%
  #filter((assay == 'lipid.class.prop' & trans == 'identity')) %>% # |
  #         (assay == 'lipid.class' & trans == 'log' & !str_detect(formula, 'genotype:sex'))) %>%
  ggplot(aes(x=interaction(assay, formula, trans), y=AIC)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_quasirandom(width=0.2) +
  facet_wrap(~assay+trans, scales='free', nrow=1) + 
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + 
  ggtitle('Lipid Class Models: AIC')
cowplot::save_plot(plot=p, filename='Results/model-selection-wo2/Class - AIC.png', base_height=10, base_asp=2)

p = all_glanced %>%
  #filter((assay == 'lipid.class.prop' & trans == 'identity')) %>% # |
  #         (assay == 'lipid.class' & trans == 'log' & !str_detect(formula, 'genotype:sex'))) %>%
  ggplot(aes(x=interaction(assay, formula, trans), y=BIC)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_quasirandom(width=0.2) +
  facet_wrap(~assay + trans, scales='free', nrow=1) + 
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + 
  ggtitle('Lipid Class Models: BIC')
cowplot::save_plot(plot=p, filename='Results/model-selection-wo2/Class - BIC.png', base_height=10, base_asp=2)

p = all_glanced %>%
  #filter((assay == 'lipid.class.prop' & trans == 'identity')) %>% # |
  #         (assay == 'lipid.class' & trans == 'log' & !str_detect(formula, 'genotype:sex'))) %>%
  ggplot(aes(x=interaction(assay, formula, trans), y=r.squared)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_quasirandom(width=0.2) +
  facet_wrap(~trans, scales='free_x') + 
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + 
  ggtitle('Lipid Class Models: R squared')
cowplot::save_plot(plot=p, filename='Results/model-selection-wo2/Class - r.squared.png', base_height=10, base_asp=2)

p = all_glanced %>%
  #filter((assay == 'lipid.class.prop' & trans == 'identity')) %>% # |
  #         (assay == 'lipid.class' & trans == 'log' & !str_detect(formula, 'genotype:sex'))) %>%
  ggplot(aes(x=interaction(assay, formula, trans), y=r.squared)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_quasirandom(width=0.2) +
  facet_wrap(~trans, scales='free_x') + 
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + 
  ggtitle('Lipid Class Models: Adjusted R squared')
cowplot::save_plot(plot=p, filename='Results/model-selection-wo2/Class - adj.r.squared.png', base_height=10, base_asp=2)

model_glance[['class']] = all_glanced

#### Visualize augmented #####
augmented = regressions %>% select(data, augmented, trans) %>% unnest()

# Fitted vs value
p = fit_vs_value(augmented, group='Lipid.Abbr')
cowplot::save_plot(
  plot = p,
  filename = file.path('Results-wo2', model_name(regressions), 'fit vs value.png'),
  base_height = 10, dpi=400
)

# Residuals vs fitted
p = stdresid_vs_fitted(augmented, group='Lipid.Class')
cowplot::save_plot(
  plot = p,
  filename = file.path('Results-wo2', model_name(regressions), 'stdresid-vs-fitted.png'),
  base_height=10
)

# Plot confidence intervals
regressions %>%
  unnest(augmented, data) %>%
  plot_ci(group="Lipid.Class")

# qq plot
p = regressions %>%
  unnest(augmented, data) %>%
  ggplot(aes(sample=.resid)) + 
  stat_qq() + 
  stat_qq_line() + 
  facet_wrap(~Lipid.Class, scales='free')
cowplot::save_plot(
  plot = p,
  filename = file.path('Results-wo2', model_name(regressions), 'qq.png'),
  base_height=10
)

# ANOVA
anova_list = run_anova(regressions, 'Lipid.Class')

# Post-Hoc linear model
lm_table = run_posthoc(regression, 'Lipid.Class')

#### ####
#### ---- Model Lipid Species ---- ####
x = readRDS("Data/serum_speciesSE.rds")
x = x[, colnames(x) != 'Q-LF-SLP2']

all_glanced = tibble()
#all_glanced = readxl::read_excel('Results/glance all normalizations.xlsx', sheet='species')
regressions = run_lm(x, 'lipid', 'Lipid.Species', formula = 'sex + genotype + diet + genotype:diet', trans=identity)
all_glanced <- all_glanced %>%
  bind_rows(regressions %>% select(assay, group, glanced, formula, trans) %>% unnest(glanced))

regressions = run_lm(x, 'lipid', 'Lipid.Species', formula = 'sex + genotype + diet + genotype:diet', trans=log)
all_glanced <- all_glanced %>%
  bind_rows(regressions %>% select(assay, group, glanced, formula, trans) %>% unnest(glanced))

regressions = run_lm(x, 'lipid', 'Lipid.Species', formula = 'sex + genotype + diet + genotype:diet + sex:diet + genotype:sex + genotype:diet:sex', trans=identity)
all_glanced <- all_glanced %>%
  bind_rows(regressions %>% select(assay, group, glanced, formula, trans) %>% unnest(glanced))

regressions = run_lm(x, 'lipid', 'Lipid.Species', formula = 'sex + genotype + diet + genotype:diet + sex:diet + genotype:sex + genotype:diet:sex', trans=log)
all_glanced <- all_glanced %>%
  bind_rows(regressions %>% select(assay, group, glanced, formula, trans) %>% unnest(glanced))

regressions = run_lm(x, 'prop.of.class', 'Lipid.Species', formula = 'sex + genotype + diet + genotype:diet', trans=identity)
all_glanced <- all_glanced %>%
  bind_rows(regressions %>% select(assay, group, glanced, formula, trans) %>% unnest(glanced))

regressions = run_lm(x, 'prop.of.class', 'Lipid.Species', formula = 'sex + genotype + diet + genotype:diet', trans=log)
all_glanced <- all_glanced %>%
  bind_rows(regressions %>% select(assay, group, glanced, formula, trans) %>% unnest(glanced))

regressions = run_lm(x, 'prop.of.class', 'Lipid.Species', formula = 'sex + genotype + diet + genotype:diet + sex:diet + genotype:sex + genotype:diet:sex', trans=identity)
all_glanced <- all_glanced %>%
  bind_rows(regressions %>% select(assay, group, glanced, formula, trans) %>% unnest(glanced))

regressions = run_lm(x, 'prop.of.class', 'Lipid.Species', formula = 'sex + genotype + diet + genotype:diet + sex:diet + genotype:sex + genotype:diet:sex', trans=log)
all_glanced <- all_glanced %>%
  bind_rows(regressions %>% select(assay, group, glanced, formula, trans) %>% unnest(glanced))

#### ---- Visualize glanced --- ####
p = all_glanced %>%
  #filter((assay == 'lipid.class.prop' & trans == 'identity')) %>% # |
  #         (assay == 'lipid.class' & trans == 'log' & !str_detect(formula, 'genotype:sex'))) %>%
  ggplot(aes(x=interaction(assay, formula, trans), y=AIC)) + 
  geom_boxplot() + 
  geom_quasirandom(width=0.2) +
  facet_wrap(~assay + trans, scales='free', nrow=1) + 
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + 
  ggtitle('Lipid Species Models: AIC')
cowplot::save_plot(plot=p, filename='Results/model-selection-wo2/Species - AIC.png', base_height=10, base_asp=2)

p = all_glanced %>%
  #filter((assay == 'lipid.class.prop' & trans == 'identity')) %>% # |
  #         (assay == 'lipid.class' & trans == 'log' & !str_detect(formula, 'genotype:sex'))) %>%
  ggplot(aes(x=interaction(assay, formula, trans), y=BIC)) + 
  geom_boxplot() + 
  facet_wrap(~assay+trans, scales='free', nrow=1) + 
  geom_quasirandom(width=0.2) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + 
  ggtitle('Lipid Species Models: BIC')
cowplot::save_plot(plot=p, filename='Results/model-selection-wo2/Species - BIC.png', base_height=10, base_asp=2)

p = all_glanced %>%
  #filter((assay == 'lipid.class.prop' & trans == 'identity')) %>% # |
  #         (assay == 'lipid.class' & trans == 'log' & !str_detect(formula, 'genotype:sex'))) %>%
  ggplot(aes(x=interaction(assay, formula, trans), y=r.squared)) + 
  geom_boxplot() + 
  facet_wrap(~trans, scales='free_x') + 
  geom_quasirandom(width=0.2) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + 
  ggtitle('Lipid Species Models: R squared')
cowplot::save_plot(plot=p, filename='Results/model-selection-wo2/Species - r.squared.png', base_height=10, base_asp=2)

p = all_glanced %>%
  #filter((assay == 'lipid.class.prop' & trans == 'identity')) %>% # |
  #         (assay == 'lipid.class' & trans == 'log' & !str_detect(formula, 'genotype:sex'))) %>%
  ggplot(aes(x=interaction(assay, formula, trans), y=adj.r.squared)) + 
  geom_boxplot() + 
  facet_wrap(~trans, scales='free_x') + 
  geom_quasirandom(width=0.2) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + 
  ggtitle('Lipid Species Models: Adjusted R squared')
cowplot::save_plot(plot=p, filename='Results/model-selection-wo2/Species - adj.r.squared.png', base_height=10, base_asp=2)

model_glance[['species']] = all_glanced
writexl::write_xlsx(model_glance, 'Results/glance all normalizations-wo2.xlsx')

#### --- Visualize augmented --- ####
augmented = regressions %>% select(data, augmented, trans) %>% unnest()
# Fitted vs value
p.list = augmented %>%
  group_by(Lipid.Abbr) %>%
  group_map(.f=~{
    fit_vs_value(.x, group='Lipid.Species') + 
    labs(title = .y)
  })
walk(p.list, ~{
  cowplot::save_plot(
    plot = .x,
    filename = file.path('Results', model_name(regressions), 'fitted-vs-value', paste0(.x$labels$title[[1]], '.png')),
    base_height=10
  )
})

# Residuals vs fitted
p.list = augmented %>%
  group_by(Lipid.Abbr) %>%
  group_map(~{
    stdresid_vs_fitted(.x, group='Lipid.Species') + 
    labs(title = .y)
  })
walk(p.list, ~{
  cowplot::save_plot(
    plot = .x,
    filename = file.path('Results', model_name(regressions), 'stdresid-vs-fitted', paste0(.x$labels$title[[1]], '.png')),
    base_height=10
  )
})

# Confidence intervals
regressions %>%
  unnest(augmented, data) %>%
  plot_ci(group='Lipid.Abbr')

# QQ plots
p.list = regressions %>%
  unnest(augmented, data) %>%
  group_by(Lipid.Abbr) %>%
  group_map(~{
    p = .x %>%      
      ggplot(aes(sample=.resid)) + 
      stat_qq() + 
      stat_qq_line() + 
      facet_wrap(~Lipid.Species, scales='free') + 
      labs(title = .y) 
    p
  })
walk(p.list, ~{
  cowplot::save_plot(
    plot = .x,
    filename = file.path('Results', model_name(regressions), 'qq', paste0(.x$labels$title[[1]], '.png')),
    base_height=10
  )
})

# ANOVA
anova_table = run_anova(regressions, 'Lipid.Species')

anova_table$`genotype x diet` %>%
  left_join(row.data) %>%
  ggplot(aes(x = fct_reorder(Lipid.Abbr, `F value`), y = `F value`)) + 
  geom_quasirandom(width=0.2)

anova_table$`genotype` %>%
  left_join(row.data) %>%
  ggplot(aes(x = `Pr(>F)`)) + 
  geom_histogram() + 
  facet_wrap(~Lipid.Abbr)

# post-hoc table
lm_table = run_posthoc(regression, 'Lipid.Species')

lm_table$genotypeE2 %>% 
  left_join(row.data) %>% 
  ggplot(aes(x = Estimate, y = -log10(`Pr(>|t|)`), fill=Lipid.Abbr)) +
  geom_point(shape=21,size=2) +
  scale_fill_manual(values = Lipid.Abbr.colors) +
  facet_wrap(~Lipid.Abbr)

lm_table$dietCR %>%
  left_join(
  anova_table$diet %>% select(Lipid.Species, p_val_adj)
)

lm_table$genotypeE2 %>%
  filter(Lipid.Species %in% (anova_table$genotype %>% filter(p_val_adj < 0.05) %>% select(Lipid.Species) %>% unlist))

lm_table$genotypeE4 %>%
  filter(Lipid.Species %in% (anova_table$genotype %>% filter(p_val_adj < 0.05) %>% select(Lipid.Species) %>% unlist))

lm_table$`genotypeE2 x dietCR` %>%
  filter(Lipid.Species %in% (anova_table$`genotype x diet` %>% head(n=10) %>% select(Lipid.Species) %>% unlist))

lm_table$`genotypeE4 x dietCR` %>%
  filter(Lipid.Species %in% (anova_table$`genotype x diet` %>% head(n=10) %>% select(Lipid.Species) %>% unlist))

# enrichment test
library(clusterProfiler)
TERM2GENE = rowData(x) %>% as.data.frame %>% rownames_to_column('gene') %>% rename(term=Lipid.Abbr) %>% select(term, gene)
TERM2GENE %<>% mutate(term = ifelse(gene %in% c('TC', 'FC'), 'CH', term),
                     term = ifelse(term == 'CE', 'CH', term))
TERM2NAME = rowData(x) %>% as.data.frame %>% rename(term=Lipid.Abbr,name=Lipid.Class) %>% select(term, name)
TERM2NAME %<>% mutate(term = ifelse(term %in% c('TC', 'FC'), 'CH', term),
                      term = ifelse(term == 'CE', 'CH', term))
anova_ora = map(anova_table, ~{
  #gene = .x %>% filter(`F value` > 20) %>% select(Lipid.Species) %>% unlist %>% unname
  gene = .x %>% filter(p_val_adj < 0.05) %>% select(Lipid.Species) %>% unlist %>% unname
  #gene = .x %>% head(20) %>% select(Lipid.Species) %>% unlist %>% unname
  enricher(gene, TERM2GENE=TERM2GENE, TERM2NAME=TERM2NAME, minGSSize = 1)
})
names(anova_ora) = names(anova_table)
map(anova_ora, ~if(!is.null(.x)) { .x@result })
anova_ora_tables = map(anova_ora, ~if(!is.null(.x)) {as.data.frame(.x@result)})
anova_ora_tables = anova_ora_tables[unlist(map(anova_ora_tables, ~!is.null(.x)))]
writexl::write_xlsx(x = anova_ora_tables,
                    path = file.path("Results", model_name(regressions), 'anova_ora.xlsx'))
