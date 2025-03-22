%% lifestyle_RNAseq_6_cpur_heatmaps.m   %% 
%%% CRAL_IPA comparison
clear;
close all;
clc;

%%
cd('D:\Cortical RNAseq\new_analysis_flowcell_coverage_extractiondate_adjusted\DEG');
[a1,b1,~] = xlsread('CR-AL-pathway-cmp.xls');
[a2,b2,~] = xlsread('CR-AL-ur-cmp.xls');
[a3,b3,~] = xlsread('CR-AL-pathway-cmp-rawp-zscore=2.xls');

%% pathway heatmap
figure ('position',[0.0010    0.0410    1.5360    0.7488]*1000);
subplot(121)
% xvalues = {'all' 'E2','E3','E4'};
% yvalues = b1(3:end,1);
% h = heatmap(xvalues,yvalues,a1,'Colormap',redbluecmap); %% 
% h.Title = 'CR-AL-APOE-pathway-cmp-FDR1.3';
% h.XLabel = 'CR vs AL';
% h.YLabel = 'IPA pathway';
xvalues = {'all' 'E2','E3','E4'};
yvalues = b3(3:end,1);
% h = heatmap(xvalues,yvalues,a3,'Colormap',redbluecmap); %% 
h = heatmap(xvalues,yvalues,a3,'MissingDataColor',[0.5 0.5 0.5]);
h2=colormap(gca,bluewhitered(256)); 
h3=colormap(gca);
colorbar;

h.Title = 'CR-AL-APOE-pathway-cmp-rawp-z=2';
h.XLabel = 'CR vs AL';
h.YLabel = 'IPA pathway';

subplot(122)
xvalues = {'all' 'E2','E3','E4'};
yvalues = b2(3:end,1);
yvalues = uplow_cap(yvalues);

% h = heatmap(xvalues,yvalues,a2,'Colormap',redbluecmap); %% 
h = heatmap(xvalues,yvalues,a2,'MissingDataColor',[0.5 0.5 0.5]);
colormap(gca,bluewhitered(256)); 
colorbar;
h.Title = 'CR-AL-APOE-ur-cmp-z2';
h.XLabel = 'CR vs AL';
h.YLabel = 'Upstream regulator';

fn = 'CR-AL-APOE-pathway-ur-cmp-redblue.emf';
cd('D:\Cortical RNAseq\new_analysis_flowcell_coverage_extractiondate_adjusted\Ana-fig');
saveas(gcf,fn);

%%
% cd('D:\Cortical RNAseq\new_analysis_flowcell_coverage_extractiondate_adjusted\DEG');
% [a3,b4,~] = xlsread('CR-AL-match analysis.xls');
% 
% row_name = b4(2,3:7);
% collumn_name1 = b4(3:end,1);% genotype
% collumn_name2 = b4(3:end,2);% diet
% collumn_name = strcat(collumn_name1,'-',collumn_name2);
% CGobj  = clustergram(a3','Standardize','column','colormap',redbluecmap);  %% 
% set(CGobj,'RowLabels', row_name   ,'ColumnLabels',collumn_name); 
% y_label = flip(CGobj.RowLabels);
% x_label = CGobj.ColumnLabels; 


