source("/research/labs/moleneurosci/bug/projects/Plotting/single-cell/heatmap.R")
plot_class <- function(x = readRDS("Data/classSE.rds"), assay = 'lipid.class') {
  library(rlang)
  col.data = colData(x) %>% as.data.frame %>% set_levels()
  df = assay(x, i=assay) %>% as.matrix
  df = df[rownames(df) != 'TC',]
  
  beg = rownames(x)[1]
  end = rownames(x)[nrow(x)]
  
  t(df) %>% 
    as.data.frame %>%
    bind_cols(col.data) %>%
    pivot_longer(!!sym(beg) : !!sym(end), names_to='Lipid.Abbr') %>%
    mutate(Lipid.Abbr = fct_reorder(as.factor(Lipid.Abbr), .x=value, .fun=sum))
}

plot_species <- function(x = readRDS("Data/speciesSE.rds"), assay = 'lognorm') {
  col.data = colData(x) %>% as.data.frame %>% set_levels()
  df = assay(x, i=assay) %>% as.matrix
  df = df[rownames(df) != 'TC',]
  
  t(df) %>% 
    as.data.frame %>%
    rownames_to_column('SampleID') %>%
    pivot_longer(-'SampleID', names_to="Lipid.Species") %>%
    left_join(rowData(x) %>% as.data.frame %>% rownames_to_column('Lipid.Species')) %>%
    left_join(col.data %>% rownames_to_column('SampleID'), by = 'SampleID') %>%
    mutate(Lipid.Abbr = fct_reorder(as.factor(Lipid.Abbr), .x=value, .fun=sum))
} 
plot_pairs <- function(keep, group, x = readRDS("Data/speciesSE.rds"), assay = 'lognorm') {
  col.data = colData(x) %>% as.data.frame %>% set_levels()
  row.data = rowData(x) %>% as.data.frame %>% rownames_to_column('Lipid.Species')
  df = assay(x, i=assay) %>% as.matrix
  df = df[keep,]
  df = t(df) %>% as.data.frame
  df[,'group'] = colData(x)[,group]
  GGally::ggpairs(data = df, columns = colnames(df)[-ncol(df)], mapping = aes(color = group))
}
plot_species_pairs <- function( abbr, group, x = readRDS("Data/speciesSE.rds"), keep = rep(T, nrow(x)), assay = 'lognorm'
) {
  dat = bind_cols(
    t(assay(x, i=assay))[, (rowData(x)$Lipid.Abbr %in% c(abbr)) & keep],
    group = colData(x)[,group]
  )
  GGally::ggpairs(
    dat,
    columns = colnames(dat)[-ncol(dat)],
    mapping = aes(color = group)
  )
}
plot_pca_pairs <- function(x, assay, outdir = 'pairs', log=FALSE) {
  type = ifelse(nrow(x) == length(unique(rowData(x)$Lipid.Abbr)),
                'Class', 'Species')
  keep = !rowSums(assay(x) == 0)
  
  if (log==TRUE) {
    pca = PcaClassic(t(log(assay(x, i=assay)[keep,])), k=10, scale=TRUE)
  } else {
    pca = PcaClassic(t(assay(x, i=assay)[keep,]), k=10, scale=TRUE)
  }
  
  p = ggpairs(
    data = bind_cols(pca$scores[x$Group == 'CR',], as.data.frame(colData(x)[x$Group == 'CR',])), 
    columns = paste0('PC',1:10), 
    mapping = aes(color = genotype)
  )
  cowplot::save_plot(plot=p,
                     filename = file.path('Results/eda', type, outdir, paste0(ifelse(log, 'log ', ''), assay, ' CR pairs.png')),
                     base_height=20)
  p = ggpairs(
    data = bind_cols(pca$scores[x$Group == 'AL',], as.data.frame(colData(x)[x$Group == 'AL',])), 
    columns = paste0('PC',1:10), 
    mapping = aes(color = genotype)
  )
  cowplot::save_plot(plot=p, 
                     filename = file.path('Results/eda', type, outdir, paste0(ifelse(log, 'log ', ''), assay, ' AL pairs.png')),
                     base_height=20)
   
}
plot_pca_biplot <- function(x, assay, outdir = 'biplots', log=FALSE, lipid.abbr = unique(rowData(x)$Lipid.Abbr)) {
  type = ifelse(nrow(x) == length(unique(rowData(x)$Lipid.Abbr)),
                'Class', 'Species')
  keep = !rowSums(assay(x) == 0)
  x <- x[rowData(x)$Lipid.Abbr %in% lipid.abbr, ]
  
  if (log==TRUE) {
    pca = PcaClassic(t(log(assay(x, i=assay)[keep,])), k=10, scale=TRUE)
  } else {
    pca = PcaClassic(t(assay(x, i=assay)[keep,]), k=10, scale=TRUE)
  }

  p.list = purrr::map2(paste0('PC', 1:7), paste0('PC', 2:8), ~{
    bind_cols(pca$scores, as.data.frame(colData(x))) %>%
    ggplot(aes(!!sym(.x), !!sym(.y), fill=genotype)) + 
    facet_wrap(~Group) + 
    geom_point(shape=21, size=2, show.legend = .x %in% c('PC3', 'PC6')) + 
    geom_mark_ellipse(alpha=0.3, show.legend = .x %in% c('PC3', 'PC6'))
  }) 
  cowplot::save_plot(
    plot=wrap_plots(p.list, ncol=3),
    filename = file.path('Results/eda', type, outdir, paste0(ifelse(log, 'log ', ''), assay, ' biplots.png')),
    base_height=10,
    base_asp=2
  )
}

plot_heatmap <- function(
  x = readRDS("Data/classSE.rds"),
  assay = 'lipid.class',
  Lipid.Abbr = unique(rowData(x)$Lipid.Abbr),
  order=FALSE,
  h=7,
  show_names = FALSE
) {
  row.data = rowData(x) %>% as.data.frame
  keep = !rowSums(assay(x) == 0)
  x = x[keep,]
  if (assay == 'prop.of.class') {
    keep = !rownames(x) %in% c('TC', 'FC')
    x = x[keep,]
    row.data = row.data[keep,]
  }
  lipid.abbrs = rowData(x)$Lipid.Abbr
  corr = cor(t(assay(x, i=assay)), method = 'spearman')
  keep = row.data[rownames(corr), 'Lipid.Abbr'] %in% Lipid.Abbr
  dendro = hclust(dist(corr[keep,keep], ))
  if (order == TRUE) {
    order = cutree(dendro, h=h)
    row_split = column_split = order
    if (show_names == TRUE) {
      show_row_names = show_column_names = TRUE
    } else {
      show_row_names = show_column_names = FALSE
    }
  } else {
    row_split = column_split = lipid.abbrs[keep]
    show_row_names = show_column_names = FALSE
  }
  
  h = 
    make_heatmap(
      corr[keep,keep], 
      col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
      row_split = row_split,
      column_split = column_split,
      width = sum(keep) * unit(3, 'mm'),
      height = sum(keep) * unit(3, 'mm'),
      show_row_names = show_row_names,
      show_column_names = show_column_names,
      cluster_column_slices = TRUE,
      row_title_rot = 0,
      show_row_dend = TRUE,
      show_column_dend = TRUE,
      heatmap_legend_param = list(
        title = 'Normalized gene expression',
        direction = 'vertical',
        legend_height = unit(3, 'cm'),
        at = c(-1, 0, 1)
      ),
      rect_gp = gpar(col='white', lwd=0)
    )
  h
}

plot_dendro <- function(
  x = readRDS("Data/classSE.rds"),
  assay = 'lipid.class',
  log=FALSE,
  Lipid.Abbr = unique(rowData(x)$Lipid.Abbr),
  h=7,
  show_names = FALSE
) {
  row.data = rowData(x) %>% as.data.frame
  keep = !rowSums(assay(x) == 0)
  x = x[keep,]
  if (str_detect(assay, 'prop')) {
    keep = !rownames(x) %in% c('TC', 'FC')
    x = x[keep,]
    row.data = row.data[keep,]
  }
  lipid.abbrs = rowData(x)$Lipid.Abbr
  if (log) {
    corr = cor(t(log(assay(x, i=assay))), method = 'spearman')
  } else {
    corr = cor(t(assay(x, i=assay)), method = 'spearman')
  }
  keep = row.data[rownames(corr), 'Lipid.Abbr'] %in% Lipid.Abbr
  dendro = hclust(dist(corr[keep,keep], ))
  plot(dendro)
  dendro
}
