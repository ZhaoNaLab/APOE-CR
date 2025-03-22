library(SummarizedExperiment)
library(magrittr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(batchtools)
library(forcats)
library(ggpubr)

library(ComplexHeatmap)
library(patchwork)
library(ggplot2)
library(ggbeeswarm)
library(RColorBrewer)
library(GGally)
library(ggforce)

if (!'x' %in% ls()) {
  x = readRDS("Data/ox_initialSE.rds")
}
x$ID.group = paste(colData(x)$genotype, colData(x)$Group, str_sub(colData(x)$sex, end=1), str_remove(colnames(x), '[A-Z-]*')) %>% as.factor
ID.group.levels = levels(x$ID.group)
Lipid.Abbr.levels = unique(rowData(x)$Lipid.Abbr)
Lipid.Class.levels = unique(rowData(x)$Lipid.Class)
Lipid.Abbr.colors = c(brewer.pal(5,'RdPu'),
                      brewer.pal(5, 'PuBu'),
                      brewer.pal(5, 'BuGn'),
                      brewer.pal(5, 'OrRd'))
names(Lipid.Abbr.colors) = Lipid.Abbr.levels
x$genotype.diet = paste(x$genotype, x$Group) %>% as.factor
genotype.diet.colors = c(brewer.pal(n=9,'Greens')[c(2,7)],
                         brewer.pal(n=9,'Blues')[c(2,7)],
                         brewer.pal(n=9,"Reds")[c(2,7)])
names(genotype.diet.colors) = levels(x$genotype.diet)

set_levels <- function(object) {
  object %>%
    mutate(
      diet = Group,
      genotype = factor(genotype, c('E3', 'E2', 'E4')),
      diet = factor(diet, c('AL', 'CR')),
      sex = factor(sex, c('Female', 'Male')),
      genotype.diet = interaction(genotype, diet, sep=' ')  
    )
}

# x : classSE.rds or speciesSE.rds
# assay: assays
# group: Lipid.Species or Lipid.Class
run_lm <- function(x, assay, group, formula, trans) {
  source('Analysis/util.r')
  library(rlang)
  library(purrr)
  library(broom)
  
  ##### ------ Setup data ------ #####
  lipid.assay.longer <- assay(x, assay) %>% 
    as.data.frame %>% rownames_to_column(group) %>%
    pivot_longer(-sym(group), names_to='ID') %>%
    left_join(rowData(x) %>% as.data.frame %>% rownames_to_column(group)) %>%
    left_join(colData(x) %>% as.data.frame %>% rownames_to_column('ID')) %>% 
    set_levels()
  
  #### --- Apply a filter --- ####
  lipid.assay.longer <- lipid.assay.longer %>%
    group_by(!!sym(group)) %>%
    mutate(num_zeros = sum(value == 0)) %>%
    filter(num_zeros == 0) %>%
    mutate(num_zeros = NULL)
  
  if (group == 'Lipid.Species' & assay == 'prop.of.class') {
    by_lipid = lipid.assay.longer %>% 
      filter(!Lipid.Abbr %in% c('TC', 'FC')) %>%
      group_by(!!sym(group)) %>%
      nest()
  } else {
    by_lipid = lipid.assay.longer %>% 
      group_by(!!sym(group)) %>%
      nest()
  }
  ##### ------ Do Modeling ----- #####
  library(modelr)
  by_lipid_model = function(df, formula, trans) {
    lm(as.formula(paste0('value ~ ', formula)), data=df %>% mutate(value = trans(value)))
  }
  ##### ----- return result ----- #####
  regressions = by_lipid %>%
    mutate(
      assay=assay,
      group=group,
      formula = formula,
      trans = as_name(enquo(trans)),
      model = purrr::map(data, by_lipid_model, formula=formula, trans=trans),
      glanced = purrr::map(model, broom::glance),
      augmented = purrr::map(model, broom::augment, se_fit=TRUE, interval='confidence'),
      tidied = purrr::map(model, broom::tidy)
    )
  return(regressions) 
}

# Visualize the output of run_lm
run_visualization <- function(regressions) {
  group = regressions[1,'group'][[1]]
  ##### ------ Plots ----- #####
  p = regressions %>%
    unnest(glanced) %>%
    arrange(r.squared) %>%
    unnest(data) %>%
    dplyr::select(-c(lipidomics.ID:CFC.cued)) %>%
    distinct(r.squared, .keep_all=TRUE) %>%
    ggplot(aes(fct_reorder(Lipid.Abbr, r.squared), r.squared)) + 
    geom_quasirandom(width=0.2) + 
    theme(axis.text.x=element_text(angle=45, hjust=1))
  cowplot::save_plot(plot=p, filename=file.path('Results', model_name(regressions), 'r squared by class.png'), base_height=6)
  
  ###### ------ Plot bad fits  ------ #######
  bad_fit <- regressions %>% 
    unnest(glanced) %>%
    filter(r.squared < 0.1)
  if (nrow(bad_fit) > 0) {
    
    p = bad_fit %>% unnest(data, augmented) %>% 
      ggplot(aes(.fitted, .resid)) + 
      geom_point() + 
      facet_wrap(group, scales='free') + 
      geom_hline(yintercept = 0, linetype = 'dashed')
    
    p = bad_fit %>% unnest(augmented, data) %>%
      ggplot() + 
      geom_boxplot(aes(x=genotype, value, fill=Group)) + 
      facet_wrap(group, scales='free')
    cowplot::save_plot(plot=p,
                       filename=file.path('Results', model_name(regressions), 'bad fits boxplots.png'),
                       base_asp=1.6,
                       base_height=12)
      
    p = bad_fit %>% unnest(augmented) %>%
      ggplot(aes(sample=.std.resid)) + 
      stat_qq() + geom_abline() + 
      facet_wrap(group, scales='free')
    cowplot::save_plot(plot=p,
                     filename=file.path('Results', model_name(regressions), 'bad fits qq.png'),
                     base_asp=1.6,
                     base_height=12)
  }
 
  ##### ---------- Coef plots ----------- ##### 
  library(coefplot)
  lipid_models = regressions$model
  coef_plots = lapply(lipid_models, coefplot)
  coef_plots = purrr::map2(coef_plots, by_lipid[[group]], ~{.x$labels$title = .y; .x})
  if (group == 'Lipid.Species' & assay == 'prop.of.class') {
    coef_plots.list = tibble(
      p=coef_plots,
      as_tibble(rowData(x)) %>% 
        filter(!Lipid.Abbr %in% c('FC', 'TC'))
      ) %>%
        group_by(Lipid.Abbr) %>%
        group_split
  } else {
    coef_plots.list = tibble(
      p=coef_plots,
      as_tibble(rowData(x))
      ) %>%
        group_by(Lipid.Abbr) %>%
        group_split
  }
  names(coef_plots.list) = purrr::map(coef_plots.list, ~{.x$Lipid.Abbr[[1]]}) %>% unlist
  
  purrr::walk2(coef_plots.list, names(coef_plots.list), .f = ~{
    cowplot::save_plot(plot=Reduce(`+`, .x$p),
                       filename=file.path('Results', model_name(regressions), 'coefplots', paste0(.y, '.png')),
                       base_asp=1.6,
                       base_height=12)
    })
}

# Plot fitted value vs values -- to see where models might be breaking down.
fit_vs_value <- function(object, group) {
  object %>%
  mutate(
    value = case_when(
      trans == 'log' ~ log(value),
      TRUE ~ value
    )
  ) %>%
  mutate(id = str_extract(ID, '[0-9]+')) %>%
  ggplot(aes(value, .fitted)) + 
  geom_point(aes(fill=genotype.diet, shape=sex), size=3) + 
  scale_shape_manual(values=c(21,25)) + 
  scale_fill_manual(values=genotype.diet.colors) + 
  guides(fill = guide_legend(override.aes = list(shape = c(21,25)))) + 
  geom_text(aes(label=id)) + 
  geom_abline(intercept = 0, slope=1, linetype='dashed') + 
  facet_wrap(group, scales='free') + 
  theme(text=element_text(size=10))
}

# Plot standardized residuals vs fitted values
stdresid_vs_fitted <- function(object, group) {
  object %>%
    mutate(
      value = case_when(
        trans == 'log' ~ log(value),
        TRUE ~ value
      )
    ) %>%
    mutate(id = str_extract(ID, '[0-9]+')) %>%
    ggplot(aes(.fitted, .std.resid)) + 
    geom_point(aes(fill=genotype.diet, shape=sex), size=3) + 
    scale_shape_manual(values=c(21,25)) + 
    scale_fill_manual(values=genotype.diet.colors) + 
    guides(fill = guide_legend(override.aes = list(shape = c(21,25)))) + 
    geom_text(aes(label=id)) + 
    geom_hline(yintercept = 0, linetype='dashed') + 
    facet_wrap(group, scales='free') + 
    theme(text=element_text(size=10))
}

# Output anova summary statistics for given regression model
run_anova <- function(regressions, group) {
  anova_summaries = purrr::map(regressions$model, ~summary(aov(.x)))
  names(anova_summaries) = regressions[[group]]
  anova_table = purrr::map2(anova_summaries, names(anova_summaries), ~{
    .x[[1]] <- .x[[1]] %>% as.data.frame %>%
      rownames_to_column('coef') %>%
      mutate(coef = coef %>% str_trim %>% str_squish)
    .x[[1]][[group]] = .y
    .x[[1]]
  }) %>% bind_rows
  anova_list = anova_table %>%
    filter(coef != 'Residuals') %>%
    group_by(coef) %>%
    mutate(p_val_adj = p.adjust(`Pr(>F)`)) %>%
    arrange(`Pr(>F)`) %>%
    group_split(.keep = TRUE)
  names(anova_list) = lapply(anova_list, function(x) str_replace_all(unique(x$coef), ':', ' x '))
  writexl::write_xlsx(anova_list, file.path('Results', model_name(regressions), 'anova_table.xlsx'))
  return(anova_list)
}

# Test contrasts
run_posthoc <- function(regression, group) {
  lm_summaries = purrr::map(regressions$model, ~summary(lm(.x)))
  names(lm_summaries) = regressions[[group]]
  lm_table = purrr::map2(lm_summaries, names(lm_summaries), ~{
    .x$coefficients <- .x$coefficients %>% as.data.frame %>%
      rownames_to_column('coef') %>%
      mutate(coef = coef %>% str_trim %>% str_squish)
    .x$coefficients[[group]] = .y
    .x$coefficients
  }) %>% bind_rows
  lm_list = lm_table %>%
    filter(coef != '(Intercept)') %>%
    rowwise() %>% 
    mutate(
      p_val_posthoc = case_when(
        str_detect(coef, 'genotype') ~ p.adjust(`Pr(>|t|)`, n=2, method='bonferroni'),
        TRUE ~ `Pr(>|t|)`
      )
    ) %>%
    arrange(`Pr(>|t|)`) %>%
    group_by(coef) %>%
    group_split(.keep = TRUE)
  names(lm_list) = lapply(lm_list, function(x) str_replace_all(unique(x$coef), ':', ' x '))
  writexl::write_xlsx(lm_list, file.path('Results', model_name(regressions), 'lm_table.xlsx'))
  return(lm_list)
}

# Give output of run_lm a name to save in Results
model_name = function(regressions) {
  model = regressions[1,'model'][[1]]
  model_class = class(model[[1]])
  
  formula = regressions[1,'formula'][[1]]
  if (formula == "sex + genotype + diet + genotype:diet + sex:diet + genotype:sex + genotype:diet:sex") {
    formula = 'sex*genotype*diet'
  } else if (formula == "sex + genotype + diet + genotype:diet") {
    formula = "sex+genotype*diet"
  }
  
  trans = regressions[1,'trans'][[1]]
  assay = regressions[1,'assay'][[1]]
  group = regressions[1,'group'][[1]]
  model_name = file.path(
    paste0(model_class, '-', trans, '(x)~', formula),
    assay,
#    group,
    ''
  )
  return(model_name)
}

# Plot confidence intervals
plot_ci <- function(object, group) {
  object %>%
    group_by(!!sym(group)) %>%
    group_walk(.f = ~{
      p = .x %>%
        mutate(
          value = case_when(
            trans == 'log' ~ log(value),
            TRUE ~ value
          )
        ) %>%
        ggplot() + 
        geom_jitter(aes(x=genotype.diet,
                        y=value,
                        group = interaction(genotype.diet, sex),
                        shape=sex,
                        fill=genotype.diet),
                    size=3,
                    position = position_jitterdodge(jitter.width=0.5, dodge.width=1)) + 
        geom_errorbar(aes(x = genotype.diet,
                          y = .fitted,
                          group = interaction(genotype.diet, sex),
                          ymin = .lower,
                          ymax = .upper),
                      position = position_dodge(width=1),
                      width = 0.5
                      ) + 
        scale_shape_manual(values=c(21,25)) + 
        scale_fill_manual(values=genotype.diet.colors) + 
        guides(fill = guide_legend(override.aes = list(shape = c(21,25)))) + 
        theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
      
      if (group == 'Lipid.Abbr') {
        p = p + 
          facet_wrap(~Lipid.Species, scales='free_y')
      }
      
      cowplot::save_plot(
        plot = p,
        filename = file.path('Results', model_name(regressions), 'ConfidenceIntervals', paste(.y, 'confidence intervals.png')),
        base_height=10,
        dpi = 600
      )
    }
  )
}

# Source plotting scripts.
source('Analysis/plot.r')