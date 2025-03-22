%% lifestyle_micorbiome_5_cral_heatmaps.m WQ 05312023
%% DABs hierarchical clustering heatmap  
clear;
close all hidden;
clc;

%% mat_path 
mat_path = 'D:\Microbiome\Ana-mat';
cd(mat_path);
load APOE-microbiome-abundance.mat;
load microbiome_B2_CRAL_cmp.mat;
load microbiome_B2_cmp_log2FC.mat;
save_fig_path = 'D:\Microbiome\Ana-fig';

%% find the dabs for CR vs. AL
dab_CRAL_all = cmp_B2_CRAL.GT_adj_species_dab;
dab_CRAL_E2E3E4 = unique([cmp_B2_CRAL.E2_species_dab;cmp_B2_CRAL.E3_species_dab;cmp_B2_CRAL.E4_species_dab]); 
dab_CRAL_species = unique(dab_CRAL_E2E3E4);

%% unbias hierarchical clustering for combined dabs
abd_B2_species = abd.species(:,120:end);

sample_info_B2 = abd.B2_microbiome_mouse_info;
sample_names_B2 = abd.sample_name(120:end); 
sample_info_B2_CRAL_ids= [find(contains(sample_info_B2(:,2),'AL'));...
    find(startsWith(sample_info_B2(:,2),'C'))]; 
sample_info_B2_CRAL = sample_info_B2(sample_info_B2_CRAL_ids,:);
sample_names_B2_CRAL = cellstr(sample_names_B2(sample_info_B2_CRAL_ids));
abd_B2_species_CRAL = abd_B2_species(:,sample_info_B2_CRAL_ids); 
abd_species_name = abd.species_name;

[~, ~, idb5] = intersect(dab_CRAL_species,abd_species_name,'stable');

dab_abd = abd_B2_species_CRAL(idb5,:);
row_name = dab_CRAL_species; 

collumn_name1 = sample_info_B2_CRAL(:,1);% genotype
collumn_name2 = sample_info_B2_CRAL(:,2);% diet
collumn_name3 = cellfun(@(x) x(1), sample_info_B2_CRAL(:,3),'un',0);% sex
collumn_name4 = cellfun(@(x) x(4:end),sample_names_B2_CRAL,'UniformOutput',0);% sample names
collumn_name = strcat(collumn_name4,collumn_name1,'-',collumn_name2,'-',collumn_name3);

CGobj  = clustergram(dab_abd,'Standardize','row','colormap',redbluecmap);  %% 
set(CGobj,'RowLabels',row_name,'ColumnLabels',collumn_name); 
y_label = flip(CGobj.RowLabels);
y_label_neo = cell(length(y_label),1);

for i = 1:length(y_label)
    temp = split(y_label(i),'|');
    y_label_neo(i) = temp(end);
end
x_label = CGobj.ColumnLabels; 

%% unbias hierarchical clustering for dab : cr vs al 
dab_CRAL = cmp_B2_CRAL.GT_adj_species_dab;
abd_B2_species = abd.species(:,120:end);

sample_info_B2 = abd.B2_microbiome_mouse_info;
sample_info_B2_CRAL_ids= [find(contains(sample_info_B2(:,2),'AL'));...
    find(startsWith(sample_info_B2(:,2),'C'))]; 
sample_info_B2_CRAL = sample_info_B2(sample_info_B2_CRAL_ids,:);
abd_B2_species_CRAL = abd_B2_species(:,sample_info_B2_CRAL_ids); 
abd_species_name = abd.species_name;

[identified11, ida11, idb11] = intersect(dab_CRAL,abd_species_name,'stable');
dab_abd = abd_B2_species_CRAL(idb11,:);

row_name = dab_CRAL; 

collumn_name1 = sample_info_B2_CRAL(:,1);% genotype
collumn_name2 = sample_info_B2_CRAL(:,2);% diet
collumn_name3 = sample_info_B2_CRAL(:,3);% sex
collumn_name = strcat(collumn_name1,'-',collumn_name2,'-',collumn_name3);

CGobj  = clustergram(dab_abd,'Standardize','row','colormap',redbluecmap);  %% 
set(CGobj,'RowLabels',row_name   ,'ColumnLabels',collumn_name); 
y_label = flip(CGobj.RowLabels);
x_label = CGobj.ColumnLabels; 
y_label_neo = cell(length(y_label),1);

for i = 1:length(y_label)
    temp = split(y_label(i),'|');
    y_label_neo(i) = temp(end);
end
x_label = CGobj.ColumnLabels; 

%% unbias hierarchical clustering for dab : cr vs al --E2, E3, E4
sample_info_B2_CRAL_E2_ids= find(contains(sample_info_B2_CRAL(:,1),'E2'));
sample_info_B2_CRAL_E3_ids= find(contains(sample_info_B2_CRAL(:,1),'E3'));
sample_info_B2_CRAL_E4_ids= find(contains(sample_info_B2_CRAL(:,1),'E4'));

sample_info_B2_CRAL_E2 = sample_info_B2_CRAL(sample_info_B2_CRAL_E2_ids,:);
sample_info_B2_CRAL_E3 = sample_info_B2_CRAL(sample_info_B2_CRAL_E3_ids,:);
sample_info_B2_CRAL_E4 = sample_info_B2_CRAL(sample_info_B2_CRAL_E4_ids,:);

[identified6, ida6, idb6] = intersect(cmp_B2_CRAL.E2_species_dab,abd_species_name,'stable');
dab_abd_E2 = abd_B2_species_CRAL(idb6,sample_info_B2_CRAL_E2_ids);
[identified7, ida7, idb7] = intersect(cmp_B2_CRAL.E3_species_dab,abd_species_name,'stable');
dab_abd_E3 = abd_B2_species_CRAL(idb7,sample_info_B2_CRAL_E3_ids);
[identified8, ida8, idb8] = intersect(cmp_B2_CRAL.E4_species_dab,abd_species_name,'stable');
dab_abd_E4 = abd_B2_species_CRAL(idb8,sample_info_B2_CRAL_E4_ids);
clear row_name;

row_name = cmp_B2_CRAL.E2_species_dab;
collumn_name1 = sample_info_B2_CRAL_E2(:,1);% genotype
collumn_name2 = sample_info_B2_CRAL_E2(:,2);% diet
collumn_name3 = sample_info_B2_CRAL_E2(:,3);% sex
collumn_name = strcat(collumn_name1,'-',collumn_name2,'-',collumn_name3);
CGobj_E2  = clustergram(dab_abd_E2,'Standardize','row','colormap',redbluecmap);  %% 
set(CGobj_E2,'RowLabels',row_name,'ColumnLabels',collumn_name); 
y_label_E2 = flip(CGobj.RowLabels);
x_label_E2 = CGobj.ColumnLabels; 
y_label_neo_E2 = cell(length(y_label),1);

for i = 1:length(y_label)
    temp = split(y_label(i),'|');
    y_label_neo_E2(i) = temp(end);
end
x_label = CGobj.ColumnLabels; 

row_name = cmp_B2_CRAL.E3_species_dab;
collumn_name1 = sample_info_B2_CRAL_E3(:,1);% genotype
collumn_name2 = sample_info_B2_CRAL_E3(:,2);% diet
collumn_name3 = sample_info_B2_CRAL_E3(:,3);% sex
collumn_name = strcat(collumn_name1,'-',collumn_name2,'-',collumn_name3);

CGobj_E3  = clustergram(dab_abd_E3,'Standardize','row','colormap',redbluecmap);  %% 
set(CGobj_E3,'RowLabels',row_name,'ColumnLabels',collumn_name); 
y_label_E3 = flip(CGobj.RowLabels);
x_label_E3 = CGobj.ColumnLabels; 
y_label_neo_E3 = cell(length(y_label),1);

for i = 1:length(y_label)
    temp = split(y_label(i),'|');
    y_label_neo_E3(i) = temp(end);
end
x_label = CGobj.ColumnLabels; 

row_name = cmp_B2_CRAL.E4_species_dab;
collumn_name1 = sample_info_B2_CRAL_E4(:,1);% genotype
collumn_name2 = sample_info_B2_CRAL_E4(:,2);% diet
collumn_name3 = sample_info_B2_CRAL_E4(:,3);% sex
collumn_name = strcat(collumn_name1,'-',collumn_name2,'-',collumn_name3);
CGobj_E4  = clustergram(dab_abd_E4,'Standardize','row','colormap',redbluecmap);  %% 
set(CGobj_E4,'RowLabels',row_name,'ColumnLabels',collumn_name); 
y_label_E4 = flip(CGobj.RowLabels);
x_label_E4 = CGobj.ColumnLabels; 

y_label_neo_E4 = cell(length(y_label),1);

for i = 1:length(y_label)
    temp = split(y_label(i),'|');
    y_label_neo_E4(i) = temp(end);
end
x_label = CGobj.ColumnLabels; 

%% hierarchical clustering organized sample sequence: used in the paper. organzied based on
%% Line 20 paragraph
cd('D:\Microbiome\Data analysis-Stephen\2nd round');
[a,~,~] = xlsread('hierarchical sample seq.xlsx','AL-CR');

sample_info_B2 = abd.B2_microbiome_mouse_info;
sample_IDs_B2 = abd.B2_microbiome_mouse_ID;
[identified9, ida9, idb9] = intersect(a(:,1),sample_IDs_B2,'stable');

sample_info_B2_CRAL = sample_info_B2(idb9,:);
sample_names_B2_CRAL = cellstr(sample_names_B2(idb9));
abd_species_name = abd.species_name;

abd_B2_species = abd.species(:,120:end);
abd_B2_species_CRAL = abd_B2_species(:,idb9); 
[~, ~, idb2] = intersect(dab_CRAL_species,abd_species_name,'stable');
dab_abd = abd_B2_species_CRAL(idb2,:);

row_name = dab_CRAL_species; 

collumn_name1 = sample_info_B2_CRAL(:,1);% genotype
collumn_name2 = sample_info_B2_CRAL(:,2);% diet
collumn_name3 = cellfun(@(x) x(1), sample_info_B2_CRAL(:,3),'un',0);% sex
collumn_name4 = cellfun(@(x) x(4:end),sample_names_B2_CRAL,'UniformOutput',0);% sample names
collumn_name = strcat(collumn_name4,collumn_name1,'-',collumn_name2,'-',collumn_name3);

CGobj  = clustergram(dab_abd,'Standardize','row','cluster','column','colormap',redbluecmap);  %% 
set(CGobj,'RowLabels',row_name,'ColumnLabels',collumn_name); 
y_label = flip(CGobj.RowLabels);
y_label_neo = cell(length(y_label),1);
y_label_phyla = cell(length(y_label),1);

for i = 1:length(y_label)
    temp = split(y_label(i),'|');
    y_label_neo(i) = temp(end);
    y_label_phyla(i) = temp(2);
end
x_label = CGobj.ColumnLabels; 

%% save y label  
cmp_B2_CRAL.dab_cral_ylabel_noall = y_label;
cmp_B2_CRAL.dab_cral_ylabel_noall_neo = y_label_neo;

cd(mat_path);
save('microbiome_B2_CRAL_cmp.mat',"cmp_B2_CRAL");

%% adj p for all and apoe2 3 4
CRAL_all_bac = cmp_B2_CRAL.GT_adj_species_list;
[idt1,~,idb1] = intersect(y_label,CRAL_all_bac,'stable');
adjp(:,1) = cmp_B2_CRAL.GT_adj_species_adjp(idb1);

CRAL_E2_bac = cmp_B2_CRAL.E2_species_list;
[idt2,ida2,idb2] = intersect(y_label,CRAL_E2_bac,'stable');
adjp(ida2,2) = cmp_B2_CRAL.E2_species_adjp(idb2);

CRAL_E3_bac = cmp_B2_CRAL.E3_species_list;
[idt3,ida3,idb3] = intersect(y_label,CRAL_E3_bac,'stable');
adjp(ida3,3) = cmp_B2_CRAL.E3_species_adjp(idb3);

CRAL_E4_bac = cmp_B2_CRAL.E4_species_list;
[idt4,ida4,idb4] = intersect(y_label,CRAL_E4_bac,'stable');
adjp(ida4,4) = cmp_B2_CRAL.E4_species_adjp(idb4);

%%
CRAL_bac = cmp_B2_log2FC.bac_list;
CRAL_all_log2FC = cmp_B2_log2FC.CR_AL_log2FC;
[~,~,idb] = intersect(y_label,CRAL_bac,'stable');
log2FC(:,1) = CRAL_all_log2FC(idb);

CRAL_bac_APOE = cmp_B2_log2FC.bac_list_APOE;
CRAL_E2_log2FC = cmp_B2_log2FC.E2_CR_AL_log2FC;
[~,~,idb] = intersect(y_label,CRAL_bac_APOE,'stable');
log2FC(:,2) = CRAL_E2_log2FC(idb);

CRAL_E3_log2FC = cmp_B2_log2FC.E3_CR_AL_log2FC;
[~,~,idb] = intersect(y_label,CRAL_bac_APOE,'stable');
log2FC(:,3) = CRAL_E3_log2FC(idb);

CRAL_E4_log2FC = cmp_B2_log2FC.E4_CR_AL_log2FC;
[~,~,idb] = intersect(y_label,CRAL_bac_APOE,'stable');
log2FC(:,4) = CRAL_E4_log2FC(idb);

%% fold change heatmap of the combined dabs
lfc_map = log2FC(:,2:4);
figure ('position',[0.0010    0.0418    1.5336    0.7408]*1000);
subplot(121);
xvalues = {'E2' 'E3' 'E4'};
yvalues = y_label_neo;
heatmap(xvalues,yvalues,lfc_map,'MissingDataColor',[0.5 0.5 0.5]); %% ,
colormap(gca,bluewhitered(256)); 
colorbar;
ylabel('LFC');
tit = 'CR-AL-LFC';
title(tit);

subplot(122);
xvalues = {'E2' 'E3' 'E4'};
yvalues = y_label_neo;
adjp_map = -log10(adjp(:,2:4));
adjp_map(isinf(adjp_map)) = NaN;
adjp_map(isnan(adjp_map)) = NaN;

adjp_map(adjp(:,2:4)>0.05)=0;
id1 = intersect(find(adjp(:,2:4)<0.05),find(adjp_map<2));
adjp_map(id1)=1;%% *
id2 = intersect(find(adjp_map>2),find(adjp_map<3));
adjp_map(id2)=2;%% **
id3 = intersect(find(adjp_map>3),find(adjp_map<4));
adjp_map(id3)=3;%% ***
id4 = find(adjp_map>4);
adjp_map(id4)=4;%% ****
heatmap(xvalues,yvalues,adjp_map,'MissingDataColor',[0.5 0.5 0.5]); %% ,
colormap(gca,bluewhitered(256)); 
colorbar;
ylabel('adjp');
tit = 'CR-AL-adjp';
title(tit);
fig_fn = 'CRAL-logfc-adjp-updated_e234_noall.emf';
cd(save_fig_path);
saveas(gcf,fig_fn);


