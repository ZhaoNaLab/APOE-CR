source('Analysis/util.r')
source('../../projects/Plotting/single-cell/heatmap.R')

#### Normalizations #####
# as a first pass, perform log normalization.
setup_normalization <- function(
  initial='Data/initialSE.rds', 
  species='Data/speciesSE.rds', 
  class='Data/classSE.rds'
) {
  x = readRDS(initial)
  lipid.assay = assay(x)
  lognorm.assay = log(assay(x))
  lipid.assay.longer = lipid.assay %>% as.data.frame %>% rownames_to_column("Lipid.Species") %>%
    pivot_longer(-Lipid.Species, names_to='ID') %>%
    left_join(rowData(x) %>% as.data.frame %>% rownames_to_column('Lipid.Species')) %>%
    left_join(colData(x) %>% as.data.frame %>% rownames_to_column('ID')) %>%
    mutate( value = ifelse(str_detect(Lipid.Unit, 'nmol'), value*1e3, value))
  prop.of.class.assay = lipid.assay.longer %>%  # A proportion of each class
    group_by(ID, Lipid.Abbr) %>%
    mutate(value = value / sum(value)) %>%
    pivot_wider(values_from=value, names_from=Lipid.Species, id_cols = ID) %>% 
    column_to_rownames('ID') %>% t
  prop.of.total.assay = lipid.assay.longer %>%  # A proportion of total 
    group_by(ID) %>%
    mutate(value = value / sum(value)) %>%
    pivot_wider(values_from=value, names_from=Lipid.Species, id_cols = ID) %>% 
    column_to_rownames('ID') %>% t
  lipid.class.assay = lipid.assay.longer %>% # a sum of each class, in pmol
    group_by(ID, Lipid.Abbr) %>%
    summarise(value = sum(value)) %>%
    pivot_wider(values_from=value, names_from=Lipid.Abbr, id_cols = ID) %>% 
    column_to_rownames('ID') %>% t
  lipid.class.prop.assay = lipid.assay.longer %>% # the proportion of each class
    group_by(ID, Lipid.Abbr) %>%
    summarise(value = sum(value)) %>%
    group_by(ID) %>%
    summarise(value = value / sum(value), Lipid.Abbr) %>%
    pivot_wider(values_from=value, names_from=Lipid.Abbr, id_cols = ID) %>% 
    column_to_rownames('ID') %>% t
    
  assays(x) = list(lipid=lipid.assay,
                   lognorm=lognorm.assay,
                   prop.of.class=prop.of.class.assay,
                   prop.of.total=prop.of.total.assay)
  saveRDS(x, species)
  
  xb = SummarizedExperiment(assays=list(lipid.class = lipid.class.assay[,colnames(x)],
                                        lipid.class.prop = lipid.class.prop.assay[, colnames(x)]),
                            colData = colData(x))
  rowData(xb)$Lipid.Abbr = rownames(xb)
  saveRDS(xb, class)
}
setup_normalization()

assay_split = list(E2 = assay(x)[,x$genotype == 'E2'],
                   E3 = assay(x)[,x$genotype == 'E3'],
                   E4 = assay(x)[,x$genotype == 'E4'])

# 12/8/2022
# prepare rownames
assay_split_new = purrr::map(assay_split, ~{
  rownames(.x) = str_replace_all(rownames(.x), 'Lyso ', 'L')
  rownames(.x) = str_replace_all(rownames(.x), 'CER', 'Cer')
  rownames(.x)[rownames(.x) == 'FC'] = 'ST 27:1;O'
  .x <- .x[rownames(.x) != "TC",]
  #tmp = str_split_fixed(rownames(.x), ' ', 2)
  #rownames(.x) = paste0(tmp[,1], '(', tmp[,2], ')')
  keep = str_detect(rownames(.x), 'CE CE')
  rownames(.x)[keep] = paste('CE', str_trim(str_remove(rownames(.x)[keep], 'CE CE')))
  rownames(.x) = str_replace(rownames(.x), 'BMP', 'LBPA')
  #rownames(.x) = str_remove(str_replace(rownames(.x), '\\(CE', '('), ' ')
  .x
})
rownames(assay_split_new[[1]])
dat = t(assay_split_new$E2) %>% 
  as.data.frame %>%
  rownames_to_column('Sample')
dat$Group = colData(x)[dat$Sample,]$Group
write.csv(dat, 'Data/E2.csv', row.names = T)
write.csv(assay_split_new$E3, 'Data/E3.csv', row.names = F)
write.csv(assay_split_new$E4, 'Data/E4.csv', row.names = T)

#### Normalizations -- includying oxylsterols #####
# as a first pass, perform log normalization.
setup_normalization(
  'Data/ox_initialSE.rds',
  'Data/ox_speciesSE.rds',
  'Data/ox_classSE.rds'
)

#### Normalization -- for serum lipidomics #####
setup_normalization(
  'Data/serum_initialSE.rds',
  'Data/serum_speciesSE.rds',
  'Data/serum_classSE.rds'
)