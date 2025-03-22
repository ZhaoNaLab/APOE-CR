%% lifestyle_RNAseq_4_cral_ipa_pathway.m   %% 
%% generate module related CRAL gene expression fold change/p values for IPA; hub gene network plot
%% module of interests: salmon, brown, green. 
clear;
close all;
clc;

%% green brown salmon module IPA pathway figures
datapath = 'D:\Cortical RNAseq\new_analysis_flowcell_coverage_extractiondate_adjusted\WGCNA\ND_CR\IPA';
save_fig_path = 'D:\Cortical RNAseq\new_analysis_flowcell_coverage_extractiondate_adjusted\Ana-fig\CRND-wgcna';
cd(datapath);
brown_fn = 'Brown-only genes-pathway.xls';
salmon_fn = 'Salmon-only genes-pathway.xls';
green_fn = 'Green-only genes-pathway.xls';
fn_names = {brown_fn salmon_fn green_fn};
n_fn = length(fn_names);
p_thr = 2;
for i = 1:n_fn
temp = char(fn_names(i));
cd(datapath);
[a,b,~] = xlsread(temp);
p = a(:,1);
p_sig_id = find(p>p_thr);
p_sig = p(p_sig_id);
overlap_ratio = a(p_sig_id,2); 
n_pathway = length(p_sig);
pathway_name = b(p_sig_id+2,1);

figure('position',[0.0010    0.0410    1.5360    0.7488]*1000);
for j = 1:n_pathway
    h = bar(j,p_sig(j),0.4);
    set(h,'facecolor',[1,1,1]);
    hold on;
    text(j,p_sig(j)+0.2,num2cell(overlap_ratio(j)*100));   
    xtickangle(gca,90);
end
plot([0 n_pathway+0.5],[p_thr p_thr],':k');

set(gca,'xlim',[0 n_pathway+0.5],'xtick',1:n_pathway,'xticklabel',pathway_name);
ylabel('-log(p-value)');
fig_fn = sprintf('%s-IPA-p-%s.emf',temp(1:end-4),string(p_thr));
title(fig_fn(1:end-4));
cd(save_fig_path);
saveas(gcf,fig_fn);
end




%% 























