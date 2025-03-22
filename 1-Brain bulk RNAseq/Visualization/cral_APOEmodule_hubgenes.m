%% lifestyle_RNAseq_5_APOEmodule_hubgenes.m   %% 
%% generate module related APOE related for IPA; hub gene network plot
%% module of interests: salmon, brown, green. 
clear;
close all;
clc;

%% 
cd('D:\Cortical RNAseq\new_analysis_flowcell_coverage_extractiondate_adjusted\Ana_mat\WGCNA');
load ('WGCNA_AL_CR.mat');
APOE_sigmodule_CRAL.module_oi = {'lightcyan' 'pink' 'yellow' 'lightyellow' 'black' 'tan'};
n_module = length(APOE_sigmodule_CRAL.module_oi); 

for i = 1:n_module
module_name = APOE_sigmodule_CRAL.module_oi(i);
id = find(startsWith(WGCNA_AL_CR.gene_module_identity,module_name));
APOE_sigmodule_CRAL.gene_info{i} = WGCNA_AL_CR.gene_IDs(id,:);
mmp_id = find(contains(WGCNA_AL_CR.column_terms,strcat('.',module_name)));
APOE_sigmodule_CRAL.mmp {i}= WGCNA_AL_CR.gene_mm_p(id,mmp_id-10);

%%%% hub genes
n_hub = 10;
[sort_p,sort_id] = sort(APOE_sigmodule_CRAL.mmp {i}(:,2),'ascend');
if i==3||i==6
APOE_sigmodule_CRAL.hub_gene{i} = APOE_sigmodule_CRAL.gene_info{i}(sort_id(1:n_hub),:);
elseif i==1
APOE_sigmodule_CRAL.hub_gene{i} = APOE_sigmodule_CRAL.gene_info{i}(sort_id([1:8,10:11]),:);
elseif i==2
APOE_sigmodule_CRAL.hub_gene{i} = APOE_sigmodule_CRAL.gene_info{i}(sort_id([4:6,11,13,15,18:20,23]),:);
elseif i==4
APOE_sigmodule_CRAL.hub_gene{i} = APOE_sigmodule_CRAL.gene_info{i}(sort_id([1,5,8:15]),:);
elseif i==5
APOE_sigmodule_CRAL.hub_gene{i} = APOE_sigmodule_CRAL.gene_info{i}(sort_id([1,2,5,6,8:13]),:);
end
%%%%
cd('D:\Cortical RNAseq\new_analysis_flowcell_coverage_extractiondate_adjusted\Ana_mat');
load ('gene_expression_adjusted.mat');
gene_list2 = expression_adjusted.gene_info(:,2);
CR_AL_expr = expression_adjusted.AL_CR_data;
[identified,ida,idb] = intersect(APOE_sigmodule_CRAL.hub_gene{i}(:,2),gene_list2,'stable');
APOE_sigmodule_CRAL.hub_gene_expr{i} = CR_AL_expr(idb,:);

%%%%
save_fig_path = 'D:\Cortical RNAseq\new_analysis_flowcell_coverage_extractiondate_adjusted\Ana-fig\CRND-wgcna';
APOE_sigmodule_CRAL.hubds{i} = module_hub_gene_network(APOE_sigmodule_CRAL.hub_gene{i}(:,2),APOE_sigmodule_CRAL.hub_gene_expr{i},...
    char(module_name),save_fig_path);
end

%%%%
cd('D:\Cortical RNAseq\new_analysis_flowcell_coverage_extractiondate_adjusted\Ana_mat')
mat_fn = 'APOE_sigmodule_CRAL.mat';
save(mat_fn,"APOE_sigmodule_CRAL");

%% hub gene FC and p
cd('D:\Cortical RNAseq\new_analysis_flowcell_coverage_extractiondate_adjusted\Ana_mat')
load('CR_AL_DEG.mat');
gene_list3 = cellstr(NDCRANOVA.GeneName);









