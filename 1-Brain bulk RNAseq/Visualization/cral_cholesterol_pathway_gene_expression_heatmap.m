%% APOE_lifestyle_cholesterol_pathway_gene_expression_heatmap.m 
close all;
clear all;
clc;

%%  gene list on the figure 18 genes
chlesterol_genes = {'Hmgcs1' 'Hmgcr' 'Mvk' 'Mvd' 'Pmvk' 'Fdps' 'Fdft1' 'Sqle' 'Lss' 'Cyp51' 'Tm7sf2' 'Msmo1' 'Nsdhl' 'Hsd17b7' 'Dhcr24' 'Ebp' 'Sc5d' 'Dhcr7'};

%% gene expression fold changes and p values
cd("D:\Cortical RNAseq\new_analysis_flowcell_coverage_extractiondate_adjusted\Ana_mat");
load('CR_AL_DEG.mat');
gene_list = cellstr(NDCRANOVA.GeneName);

[identified,ida,idb] = intersect(chlesterol_genes,gene_list,'stable');

all_FDR_FC = [NDCRANOVA.stepuppvalueCRVsND(idb) NDCRANOVA.FoldChangeCRVsND(idb)];
E2_FDR_FC = [NDCRANOVA.stepuppvalueE2CRVsE2ND(idb) NDCRANOVA.FoldChangeE2CRVsE2ND(idb)];
E3_FDR_FC = [NDCRANOVA.stepuppvalueE3CRVsE3ND(idb) NDCRANOVA.FoldChangeE3CRVsE3ND(idb)];
E4_FDR_FC = [NDCRANOVA.stepuppvalueE4CRVsE4ND(idb) NDCRANOVA.FoldChangeE4CRVsE4ND(idb)];

%% 
fcq = [E2_FDR_FC E3_FDR_FC E4_FDR_FC];
z = fcq(:,[2 4 6]);
q = -log10(fcq(:,[1 3 5]));
zthr = 0;
pthr = 1.3;
xlabels = {'E2' 'E3' 'E4'};
fig_fn = 'Fig_cholesterol_syn_genes_log2fcfdr_dotplot.emf';
colorbar_bin=0.5;
ylabels = chlesterol_genes';

savefigpath = 'D:\Cortical RNAseq\new_analysis_flowcell_coverage_extractiondate_adjusted\Ana-fig';
z (q<1.3) = nan; q(q<1.3)=nan;
zpfcp_dotplot_neo(z,q,zthr,pthr,xlabels,ylabels,fig_fn,colorbar_bin);
cd(savefigpath);
saveas(gcf,fig_fn);



%% expression heatmap
cd("D:\Cortical RNAseq\new_analysis_flowcell_coverage_extractiondate_adjusted\Ana_mat");
load("gene_expression_adjusted.mat");

gene_exp_list = expression_adjusted.gene_info(:,2);

[identified,ida,idb] = intersect(chlesterol_genes,gene_exp_list,'stable');

exp = expression_adjusted.AL_CR_data(idb,:);
E2IDs = find(strcmp(expression_adjusted.AL_CR_sample_info(:,3),'E2'));
E3IDs = find(strcmp(expression_adjusted.AL_CR_sample_info(:,3),'E3'));
E4IDs = find(strcmp(expression_adjusted.AL_CR_sample_info(:,3),'E4'));

ALIDs = find(strcmp(expression_adjusted.AL_CR_sample_info(:,4),'AL'));
CRIDs = find(strcmp(expression_adjusted.AL_CR_sample_info(:,4),'CR'));

MIDs = find(strcmp(expression_adjusted.AL_CR_sample_info(:,5),'Male'));
FIDs = find(strcmp(expression_adjusted.AL_CR_sample_info(:,5),'Female'));

neoids = [
    mintersect(E2IDs, ALIDs,MIDs); mintersect(E2IDs, ALIDs,FIDs);...
    mintersect(E2IDs, CRIDs,MIDs); mintersect(E2IDs, CRIDs,FIDs);...
    mintersect(E3IDs, ALIDs,MIDs); mintersect(E3IDs, ALIDs,FIDs);...
    mintersect(E3IDs, CRIDs,MIDs); mintersect(E3IDs, CRIDs,FIDs);...
    mintersect(E4IDs, ALIDs,MIDs); mintersect(E4IDs, ALIDs,FIDs);...
    mintersect(E4IDs, CRIDs,MIDs); mintersect(E4IDs, CRIDs,FIDs)];

exp_neo = exp(:,neoids);
expnormz = zscore(exp_neo');

figure('Position',[488.0000   41.8000  560.0000  740.8000]);
heatmap(neoids,ylabels,expnormz');
colormap(linspecer);
colormap(flip(brewermap([],"RdGY")));
fig_fn = 'Fig_cholesterol_syn_genes_expnormz.emf';
cd(savefigpath);
saveas(gcf,fig_fn);



















