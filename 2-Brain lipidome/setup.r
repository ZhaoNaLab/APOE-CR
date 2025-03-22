source('Analysis/util.r')
result = readxl::read_excel('Data/Result_48 Cortical Powder_Dr.Bu_12.31.2021-raw.xlsx', range='A4:AY333')
result %<>%
  rowwise() %>%
  mutate(
    group = ifelse(...3 == 'Mass', ...1, NA)
  ) %>%
  ungroup() %>%
  relocate(group, .after='...3')
for (i in 1:nrow(result)) {
  if (!is.na(result[i,'group'])) {
    group = result[i,'group']
  }
  result[i,'group'] = group
}
result %<>% mutate(group = str_squish(group))
result %<>% mutate(
  Lipid.Unit = str_extract(group, '\\((p|n)mol\\/mg protein\\)'),
  Lipid.Class = str_remove(group, Lipid.Unit) %>% str_remove('\\(\\)') %>% str_trim,
  Lipid.Class = ifelse(Lipid.Class == 'Cholesterol Ester(CE)', 'Cholesterol Ester (CE)', Lipid.Class),
  Lipid.Abbr = str_extract(Lipid.Class, ' (\\([A-Z].*\\))') %>% str_trim %>% str_remove('\\(') %>% str_remove('\\)'),
  Lipid.Class = str_remove(Lipid.Class, Lipid.Abbr) %>% str_remove('\\(\\)') %>% str_trim,
  Lipid.Abbr = case_when(
    Lipid.Class == 'Free Cholesterol' ~ 'FC',
    Lipid.Class == 'Total Cholesterol' ~ 'TC',
    TRUE ~ Lipid.Abbr
  ),
  group = NULL
) %>% relocate(Lipid.Unit, Lipid.Class, Lipid.Abbr, .after=`...3`)
colnames(result)[1:3] = c('Lipid.Species', 'Lipid.m/z', 'Lipid.Mass')
colnames(result)[7:ncol(result)] = result[1,7:ncol(result)]
result %<>% mutate(Lipid.Species = ifelse(is.na(Lipid.Species), Lipid.Abbr, paste0(Lipid.Abbr, ' ', Lipid.Species)))
data <- result %>% filter(!is.na(`Lipid.m/z`) & `Lipid.m/z` != 'm/z')
data %<>% mutate(across(c(`Lipid.m/z`, Lipid.Mass, starts_with('BU')), ~as.numeric(.x))) %>% as.data.frame %>%
  column_to_rownames('Lipid.Species')

data_behavior = readxl::read_excel('Data/lipidomics-mouse-behavior.xlsx') %>% as.data.frame
rownames(data_behavior) = paste0("BU-LF-L", data_behavior$`lipidomics-ID`)

col_lipid = data %>% dplyr::select(-starts_with('Lipid'))
row_lipid = data %>% dplyr::select(starts_with('Lipid'))

x = SummarizedExperiment(assays = list(lipid = col_lipid),
                         rowData = row_lipid,
                         colData = data_behavior)
saveRDS(x, 'Data/initialSE.rds')

# Add oxysterol information
x = readRDS("Data/initialSE.rds")
df = assay(x)
row_lipid = rowData(x) %>% as.data.frame
col_lipid = colData(x) %>% as.data.frame

ox = readxl::read_excel("Data/Result_48 Cortical Powder with OS_Dr.Bu.xlsx", sheet = 'Oxysterols', skip=4)
row_ox = ox %>% 
  select(`Oxysterols (nmol/mg protein)` : `Molecular Weight`) %>%
  rename(
    Lipid.m.z = `Molecular Weight`,
    Lipid.Mass = `Mass m/z`,
    ) %>%
  mutate(
    Lipid.Class = 'Oxysterol',
    Lipid.Abbr = 'OX',
    Lipid.Unit = '(nmol/mg protein)'
  ) %>%
  column_to_rownames('Oxysterols (nmol/mg protein)')
df_ox = ox[,-c(1:3)] %>% as.data.frame
rownames(df_ox) = rownames(row_ox)
colnames(df_ox) = colnames(df)
df_new = bind_rows(
  df,
  df_ox
)
row_new = bind_rows(
  row_lipid,
  row_ox
)
x_new = SummarizedExperiment(
  assays = list(lipid = df_new),
  rowData = row_new,
  colData = col_lipid
)
saveRDS(x_new, 'Data/ox_initialSE.rds')

##### Next add plasma lipid assay #####
result = readxl::read_excel('Data/Results 48 serum Dr Wenhai Qiao Mayo.xlsx', range='A4:AY336')
result = result[-c(1:3),]
result %<>%
  rowwise() %>%
  mutate(
    group = ifelse(...2 == 'MASS' | ...2 == 'NL', ...1, NA)
  ) %>%
  ungroup() %>%
  relocate(group, .after='...2')
for (i in 1:nrow(result)) {
  if (!is.na(result[i,'group'])) {
    group = result[i,'group']
  }
  result[i,'group'] = group
}
result %<>% mutate(group = str_squish(group))
raw_data = result %>% 
  rename(Lipid.Species = `...1`, Lipid.Mass = `...2`, `Lipid.m/z` = `Index`) %>%
  #select(Lipid.Species : Lipid.MW) %>% 
  mutate(
    Lipid.Unit = str_extract(group, '\\([a-z/ ].*\\)') %>% str_remove('\\(') %>% str_remove('\\)'),
    Lipid.Abbr = str_extract(group, '\\([A-Z]*\\)') %>% str_remove('\\(') %>% str_remove('\\)'),
    Lipid.Class = str_remove(group, '\\(.*') %>% str_squish,
    Lipid.Species = paste(Lipid.Abbr, Lipid.Species)
  ) %>%
  filter(Lipid.Mass != 'MASS' & Lipid.Mass != 'NL' & !is.na(Lipid.Species) & Lipid.Mass != 'Sum')
row_data <- raw_data %>%
  select(Lipid.Species, Lipid.Class, Lipid.Abbr, Lipid.Unit, Lipid.Mass, `Lipid.m/z`) %>%
  column_to_rownames('Lipid.Species')
lipid_data <- raw_data %>% select(`#1` : `#48`) %>% 
  mutate(across(.fns=~as.numeric(.x))) %>%
  as.data.frame
col_data <- readxl::read_excel('Data/serum lipidomics-mouse-info.xlsx') %>% as.data.frame
rownames(col_data) = paste0('Q-LF-SLP', col_data$`lipidomics-ID`)
rownames(lipid_data) = rownames(row_data)
colnames(lipid_data) = rownames(col_data)

x = SummarizedExperiment(
  assays = list(lipid = lipid_data),
  rowData = row_data,
  colData = col_data
)
saveRDS(x, 'Data/serum_initialSE.rds')
