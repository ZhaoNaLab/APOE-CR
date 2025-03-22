%% APOE_lifestyle_pathway_dotplot.m  WQ03122024
close all;
clear;
clc;

%% 
save_fig_path = 'D:\Cortical RNAseq\RNAseq-Validation\prediction-summary';
cd(save_fig_path);
[a0,b0,c0] = xlsread('CRAL-all-e234-pathway-zscore.xls','final');

%% p=rawp
z = a0(:,1:4);
p = a0(:,5:8);
zthr = 2;
pthr = 1.3;
xlabels = {'All' 'E2' 'E3' 'E4'};
ylabels = b0(2:end,1);
fig_fn = 'CR-AL-all-e234-pathway-dotplot.emf';
colorbar_bin = 1;
z (p<1.3) = nan;p(p<1.3)=nan;

%%
zpfcp_dotplot_neo(z,p,zthr,pthr,xlabels,ylabels,fig_fn,colorbar_bin);
saveas(gcf,fig_fn);

%% 



