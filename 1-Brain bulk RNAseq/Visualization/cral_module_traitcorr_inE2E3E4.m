%% lifestyle_RNAseq_3_E3black_module_inE2E4.m
%% purpose: calculate the MEs of the black module genes in the E2 and E4 cohorts.
clear;
close all;
clc;

% %% E3_module_gene file 
% data_path = 'D:\Cortical RNAseq\new_analysis_flowcell_coverage_extractiondate_adjusted\WGCNA\ND_CR\E3';
% gene_info_fn = 'geneInfo.csv';
% [a1,b1,c1] = xlsread(gene_info_fn);
% 
% %%  
% E3_WGCNA_AL_CR.gene_info = c1(2:end,2:8);
% E3_WGCNA_AL_CR.gene_info_term = c1(1,2:8);
% E3_WGCNA_AL_CR.gene_module_identity = c1(2:end,8); 
% E3_WGCNA_AL_CR.gene_mm_p = a1(:,8:end);
% E3_WGCNA_AL_CR.gene_gs_p = a1(:,9:10);
% E3_WGCNA_AL_CR.gene_IDs = b1(2:end,2:3);
% 
% cd('D:\Cortical RNAseq\new_analysis_flowcell_coverage_extractiondate_adjusted\Ana_mat\WGCNA');
% mat_fn = 'E3_WGCNA_AL_CR.mat';
% save(mat_fn,"E3_WGCNA_AL_CR");
% 
% %% 
% module_oi = 'black';
% id_oi = find(contains(E3_WGCNA_AL_CR.gene_module_identity,'black')); 
% gene_oi = E3_WGCNA_AL_CR.gene_IDs(id_oi,2);

%% module_gene file 
data_path = 'D:\Cortical RNAseq\new_analysis_flowcell_coverage_extractiondate_adjusted\WGCNA\ND_CR\allsamples';
cd(data_path);
gene_info_fn = 'geneInfo.csv';
[a1,b1,c1] = xlsread(gene_info_fn);

%%  
WGCNA_AL_CR.gene_IDs = b1(2:end,2:3);
WGCNA_AL_CR.gene_info_term = c1(1,2:3);
WGCNA_AL_CR.gene_module_identity = c1(2:end,8); 
WGCNA_AL_CR.gene_chr = a1(2:end,4); 
WGCNA_AL_CR.column_terms = b1(1,:);
WGCNA_AL_CR.gene_mm_p = a1(:,8:end);
WGCNA_AL_CR.gene_gs_p = a1(:,6:7);
WGCNA_AL_CR.modules = unique(WGCNA_AL_CR.gene_module_identity);
WGCNA_AL_CR.n_modules = length(WGCNA_AL_CR.modules );
WGCNA_AL_CR.module_figure = {'blue' 'yellow' 'salmon' 'tan' 'black' 'grey60' 'brown'...
'turquoise' 'greenyellow' 'lightgreen' 'green' 'royalblue' 'lightcyan'...
'lightyellow' 'magenta' 'midnightblue' 'purple' 'red' 'cyan' 'pink' 'grey'};

cd('D:\Cortical RNAseq\new_analysis_flowcell_coverage_extractiondate_adjusted\Ana_mat\WGCNA');
mat_fn = 'WGCNA_AL_CR.mat';
save(mat_fn,"WGCNA_AL_CR");

%%   
cd('D:\Cortical RNAseq\new_analysis_flowcell_coverage_extractiondate_adjusted\Ana_mat');
load('gene_expression_adjusted.mat'); 

gene_list = expression_adjusted.gene_info (:,2);
AL_CR_exprs = expression_adjusted.AL_CR_data;
E2_AL_CR_ID = find(contains(expression_adjusted.AL_CR_sample_info(:,3),'E2'));
E3_AL_CR_ID = find(contains(expression_adjusted.AL_CR_sample_info(:,3),'E3'));
E4_AL_CR_ID = find(contains(expression_adjusted.AL_CR_sample_info(:,3),'E4'));
E2_AL_CR_meta = expression_adjusted.AL_CR_sample_info(E2_AL_CR_ID,:);
E3_AL_CR_meta = expression_adjusted.AL_CR_sample_info(E3_AL_CR_ID,:);
E4_AL_CR_meta = expression_adjusted.AL_CR_sample_info(E4_AL_CR_ID,:);
AL_CR_seq_ID = expression_adjusted.AL_CR_sample_info(:,1);

all_binary = zeros(size(AL_CR_exprs,2),1);
CR_ids= find(contains(expression_adjusted.AL_CR_sample_info(:,4),'CR'));
all_binary(CR_ids)=1; 

E2_binary = zeros(length(E2_AL_CR_ID),1);
[~,ida2,idb2] = intersect(E2_AL_CR_ID,find(contains(expression_adjusted.AL_CR_sample_info(:,4),'CR')),'stable');
E2_binary(ida2)=1;

E3_binary = zeros(length(E3_AL_CR_ID),1);
[~,ida3,idb3] = intersect(E3_AL_CR_ID,find(contains(expression_adjusted.AL_CR_sample_info(:,4),'CR')),'stable');
E3_binary(ida3)=1;

E4_binary = zeros(length(E4_AL_CR_ID),1);
[~,ida4,idb4] = intersect(E4_AL_CR_ID,find(contains(expression_adjusted.AL_CR_sample_info(:,4),'CR')),'stable');
E4_binary(ida4)=1;

%% all MEs
cd('D:\Cortical RNAseq\new_analysis_flowcell_coverage_extractiondate_adjusted\WGCNA\ND_CR\allsamples');
A = importdata('MEs.txt');
MEs = A.data;
sample_ID = A.textdata(2:end);
temp1 = A.textdata(1);
temp2 = split(temp1, ' ');
modules_YR = regexp(temp2,'[^""]*','match','once');
modules_YR = regexp(modules_YR,'[^ME]*','match','once');

[identified1, ida5,idb5] = intersect(AL_CR_seq_ID,sample_ID,'stable');
[identified2, ida6,idb6] = intersect(WGCNA_AL_CR.module_figure,modules_YR,'stable');
all_MEs = MEs(idb5,idb6)';

%% generate r p matrix: all samples, E2, E3, E4
r_matrix = zeros(WGCNA_AL_CR.n_modules,4); p_matrix = zeros(WGCNA_AL_CR.n_modules,4); 
E2_MEs = zeros(WGCNA_AL_CR.n_modules,length(E2_AL_CR_ID));
E3_MEs = zeros(WGCNA_AL_CR.n_modules,length(E3_AL_CR_ID));
E4_MEs = zeros(WGCNA_AL_CR.n_modules,length(E4_AL_CR_ID));
MEsv2 = zeros(WGCNA_AL_CR.n_modules,size(AL_CR_exprs,2));
for  i = 1:WGCNA_AL_CR.n_modules  
module_oi = WGCNA_AL_CR.modules (i);  %% 'green' is the lipid synthesis module   
id_oi = find(contains(WGCNA_AL_CR.gene_module_identity,module_oi)); 
gene_oi = WGCNA_AL_CR.gene_IDs(id_oi,2);

[identified,ida,idb] = intersect(gene_oi,gene_list,'stable'); %% note, some were not identified, try EN id instead.

E2_AL_CR_oi = AL_CR_exprs(idb,E2_AL_CR_ID);
E3_AL_CR_oi = AL_CR_exprs(idb,E3_AL_CR_ID);
E4_AL_CR_oi = AL_CR_exprs(idb,E4_AL_CR_ID);
AL_CR_oi = AL_CR_exprs (idb,:);

%%%% calculate ME 
[~,scores,~] = pca (E2_AL_CR_oi');
E2_MEs (i,:) = scores(:,1);
[~,scores,~] = pca (E3_AL_CR_oi');
E3_MEs (i,:)= scores(:,1);
[~,scores,~] = pca (E4_AL_CR_oi');
E4_MEs (i,:)= scores(:,1);
[~,scores,~] = pca (AL_CR_oi');
MEsv2 (i,:)= scores(:,1);

%%%% ME correlation with the meta
[r,p] = corrcoef(all_binary,MEsv2(i,:));p_matrix (i,1) = p(1,2);r_matrix (i,1)= r(1,2);
[r,p] = corrcoef(E2_binary,E2_MEs(i,:));p_matrix (i,2) = p(1,2);r_matrix (i,2)= r(1,2);
[r,p] = corrcoef(E3_binary,E3_MEs(i,:));p_matrix (i,3) = p(1,2);r_matrix (i,3)= r(1,2);
[r,p] = corrcoef(E4_binary,E4_MEs(i,:));p_matrix (i,4) = p(1,2);r_matrix (i,4)= r(1,2);
end

%% 
[identified3, ida7,idb7] = intersect(WGCNA_AL_CR.module_figure,WGCNA_AL_CR.modules,'stable');
r_matrix_neo = r_matrix(idb7,:);
p_matrix_neo = p_matrix(idb7,:);

for  i = 1:WGCNA_AL_CR.n_modules 
[r,p] = corrcoef(all_binary,all_MEs(i,:));
p_matrix_neo (i,1) = p(1,2);
r_matrix_neo (i,1)= r(1,2);
end 

fdr_matrix = zeros(size(p_matrix_neo));
hom_bof_matrix = zeros(size(p_matrix_neo));

for i = 1:size(p_matrix_neo,2)
fdr_matrix(:,i) = multicmp(p_matrix_neo(:,i),'fdr',0.05);
hom_bof_matrix(:,i) = multicmp(p_matrix_neo(:,i),'down',0.05);
end

%%%% heatmap clustergram
figure ('position',[0.0010    0.0410    1.5360    0.7488]*1000);
subplot(121)
xvalues = {'all' 'E2','E3','E4'};
yvalues = WGCNA_AL_CR.module_figure;
h = heatmap(xvalues,yvalues,r_matrix_neo,'Colormap',redbluecmap);

h.Title = 'Module Trait corr';
h.XLabel = 'CR vs AL';
h.YLabel = 'moduleME';
title('CR-AL-module-trait-r');

% figure;
subplot(122)
% h = heatmap(xvalues,yvalues,p_matrix_neo,'Colormap',redbluecmap);
p_map = p_matrix_neo;
p_map(p_map>0.05)=0;
id1 = intersect(find(p_map<0.05),find(p_map>0.01));
p_map(id1)=1;%% *
id2 = intersect(find(p_map<0.01),find(p_map>0.001));
p_map(id2)=2;%% **
id3 = intersect(find(p_map<0.001),find(p_map>0.0001));
p_map(id3)=3;%% ***
id4 = intersect(find(p_map<0.0001),find(p_map>0));
p_map(id4)=4;%% ****
h = heatmap(xvalues,yvalues,p_map);
colormap(gca,bluewhitered(256)); 
colorbar;

h.Title = 'Module Trait corr';
h.XLabel = 'CR vs AL';
h.YLabel = 'moduleME';
title('CR-AL-module-trait-p');
% fig_fn = 'CR_AL_module_trait_rp
% .emf'; 
fig_fn = 'CR_AL_module_trait_rp_RBWheatmap.emf'; 

cd('D:\Cortical RNAseq\new_analysis_flowcell_coverage_extractiondate_adjusted\Ana-fig\CRND-wgcna');
saveas(gcf,fig_fn);

%%
figure ('position',[0.0010    0.0410    1.5360    0.7488]*1000);

subplot(121)
p_map = fdr_matrix;
p_map(p_map>0.05)=0;
id1 = intersect(find(p_map<0.05),find(p_map>0.01));
p_map(id1)=1; %% *
id2 = intersect(find(p_map<0.01),find(p_map>0.001));
p_map(id2)=2; %% **
id3 = intersect(find(p_map<0.001),find(p_map>0.0001));
p_map(id3)=3; %% ***
id4 = intersect(find(p_map<0.0001),find(p_map>0));
p_map(id4)=4; %% ****
h = heatmap(xvalues,yvalues,p_map,'Colormap',redbluecmap);
h.Title = 'Module Trait corr';
h.XLabel = 'CR vs AL';
h.YLabel = 'moduleME';
title('CR-AL-module-trait-fdr');

subplot(122)
p_map = hom_bof_matrix;
p_map(p_map>0.05)=0;
id1 = intersect(find(p_map<0.05),find(p_map>0.01));
p_map(id1)=1; %% *
id2 = intersect(find(p_map<0.01),find(p_map>0.001));
p_map(id2)=2; %% **
id3 = intersect(find(p_map<0.001),find(p_map>0.0001));
p_map(id3)=3; %% ***
id4 = intersect(find(p_map<0.0001),find(p_map>0));
p_map(id4)=4; %% ****

h = heatmap(xvalues,yvalues,p_map,'Colormap',redbluecmap);
h.Title = 'Module Trait corr';
h.XLabel = 'CR vs AL';
h.YLabel = 'moduleME';
title('CR-AL-module-trait-hom-bof');

fig_fn = 'CR_AL_module_trait_fdr_hom_bof.emf'; 
cd('D:\Cortical RNAseq\new_analysis_flowcell_coverage_extractiondate_adjusted\Ana-fig\CRND-wgcna');
saveas(gcf,fig_fn);

%% module behavioral correlation analysis
%% module behavioral correlation analysis
cd('D:\Cortical RNAseq\new_analysis_flowcell_coverage_extractiondate_adjusted\Ana_mat');
fn = 'Lifestyle-RNA-seq-sample-Behavior.xlsx';
[a1,b1,~] = xlsread(fn);
seqid = a1(:,1);
mouseid = a1(:,2);
wgcnaidstemp = cellfun(@(x) x(7:end),AL_CR_seq_ID,'UniformOutput',false);
wgcnaids = cellfun(@(x) str2double(x),wgcnaidstemp);
[idt,ida,idb] = intersect(wgcnaids,seqid,'stable');
wgcnamouseids = mouseid(idb);

cd('D:\Manuscript\3 Raw data\Health');
fn = '0-Behavior-OFA-cral.xlsx';
[a2,b2,c2] = xlsread(fn,'sctratio-removeoutlier');
ofamouseid = a2(:,1);
ofasctratio = a2(:,14);

fn = '0-Behavior-CFC-cral.xlsx';
[a3,b3,c3] = xlsread(fn,'cued-no outliers');
cuedmouseid = a3(2:end,1);
cuedfreezeperc = a3(2:end,end);

[idt1,ida1,idb1] = intersect(wgcnamouseids,ofamouseid,'stable');
seqofa = ofasctratio(idb1,:);

[idt2,ida2,idb2] = intersect(wgcnamouseids,cuedmouseid,'stable');
seqcued = cuedfreezeperc(idb2,:);

n_modules = WGCNA_AL_CR.n_modules;

%% 
rpofa = zeros(n_modules,2);
for  i = 1:n_modules
[x,y] = corrcoef(seqofa,all_MEs(i,ida1),"Rows","complete");
rpofa(i,:) = [x(1,2) y(2,1)];
end
clear x y;

%%
rpcued = zeros(n_modules,2);
for  i = 1:n_modules
[x,y] = corrcoef(seqcued,all_MEs(i,ida2),"Rows","complete");
rpcued(i,:) = [x(1,2) y(2,1)];
end

%% 
neoids = [E2_AL_CR_ID;E3_AL_CR_ID;E4_AL_CR_ID];
neoall_MEs = [E2_MEs,E3_MEs,E4_MEs];
all_MEsv2 (:,neoids) = neoall_MEs;

%%
rpofav2 = zeros(n_modules,2);
for  i = 1:n_modules
[x,y] = corrcoef(seqofa,all_MEsv2(i,ida1),"Rows","complete");
rpofav2(i,:) = [x(1,2) y(2,1)];
end
clear x y;
rpcuedv2 = zeros(n_modules,2);

for  i = 1:n_modules
[x,y] = corrcoef(seqcued,all_MEsv2(i,ida2),"Rows","complete");
rpcuedv2(i,:) = [x(1,2) y(2,1)];
end

%% correlation in E2 E3 E4
%% 
rpofae2 = zeros(n_modules,2);
rpofae3 = zeros(n_modules,2);
rpofae4 = zeros(n_modules,2);
for  i = 1:n_modules
[x,y] = corrcoef(alcrofa(E2_AL_CR_ID),all_MEs(i,E2_AL_CR_ID));
rpofae2(i,:) = [x(1,2) y(2,1)];
[x,y] = corrcoef(alcrofa(E3_AL_CR_ID),all_MEs(i,E3_AL_CR_ID));
rpofae3(i,:) = [x(1,2) y(2,1)];
[x,y] = corrcoef(alcrofa(E4_AL_CR_ID),all_MEs(i,E4_AL_CR_ID));
rpofae4(i,:) = [x(1,2) y(2,1)];
end
clear x y;

%%
rpcuede2 = zeros(n_modules,2);
rpcuede3 = zeros(n_modules,2);
rpcuede4 = zeros(n_modules,2);
for  i = 1:n_modules
[x,y] = corrcoef(alcrcued(E2_AL_CR_ID),all_MEs(i,E2_AL_CR_ID),"Rows","complete");
rpcuede2(i,:) = [x(1,2) y(2,1)];
[x,y] = corrcoef(alcrcued(E3_AL_CR_ID),all_MEs(i,E3_AL_CR_ID),"Rows","complete");
rpcuede3(i,:) = [x(1,2) y(2,1)];
[x,y] = corrcoef(alcrcued(E4_AL_CR_ID),all_MEs(i,E4_AL_CR_ID),"Rows","complete");
rpcuede4(i,:) = [x(1,2) y(2,1)];
end
clear x y;


%% correlation in E2 E3 E4
%% 
rpofae2 = zeros(n_modules,2);
rpofae3 = zeros(n_modules,2);
rpofae4 = zeros(n_modules,2);
for  i = 1:n_modules
[x,y] = corrcoef(alcrofa(E2_AL_CR_ID),E2_MEs(i,:));
rpofae2(i,:) = [x(1,2) y(2,1)];
[x,y] = corrcoef(alcrofa(E3_AL_CR_ID),E3_MEs(i,:));
rpofae3(i,:) = [x(1,2) y(2,1)];
[x,y] = corrcoef(alcrofa(E4_AL_CR_ID),E4_MEs(i,:));
rpofae4(i,:) = [x(1,2) y(2,1)];
end
clear x y;

%%
rpcuede2 = zeros(n_modules,2);
rpcuede3 = zeros(n_modules,2);
rpcuede4 = zeros(n_modules,2);
for  i = 1:n_modules
[x,y] = corrcoef(alcrcued(E2_AL_CR_ID),E2_MEs(i,:),"Rows","complete");
rpcuede2(i,:) = [x(1,2) y(2,1)];
[x,y] = corrcoef(alcrcued(E3_AL_CR_ID),E3_MEs(i,:),"Rows","complete");
rpcuede3(i,:) = [x(1,2) y(2,1)];
[x,y] = corrcoef(alcrcued(E4_AL_CR_ID),E4_MEs(i,:),"Rows","complete");
rpcuede4(i,:) = [x(1,2) y(2,1)];
end
clear x y;

%% 
r_matrix_neo = [rpofa(:,1) rpcued(:,1)];
p_matrix_neo = [rpofa(:,2) rpcued(:,2)];

fdr_matrix = zeros(size(p_matrix_neo));
hom_bof_matrix = zeros(size(p_matrix_neo));

for i = 1:size(p_matrix_neo,2)
fdr_matrix(:,i) = multicmp(p_matrix_neo(:,i),'fdr',0.05);
hom_bof_matrix(:,i) = multicmp(p_matrix_neo(:,i),'down',0.05);
end

%%%% heatmap clustergram
figure ('position',[0.0010    0.0410    1.5360    0.7488]*1000);
subplot(121)
xvalues = {'ofa' 'cued'};
yvalues = WGCNA_AL_CR.module_figure;
h = heatmap(xvalues,yvalues,r_matrix_neo,'Colormap',redbluecmap);

h.Title = 'Module Trait corr';
h.XLabel = 'CR vs AL';
h.YLabel = 'moduleME';
title('CR-AL-module-ofacued-r');

% figure;
subplot(122)
% h = heatmap(xvalues,yvalues,p_matrix_neo,'Colormap',redbluecmap);
p_map = p_matrix_neo;
p_map(p_map>0.05)=0;
id1 = intersect(find(p_map<0.05),find(p_map>0.01));
p_map(id1)=1;%% *
id2 = intersect(find(p_map<0.01),find(p_map>0.001));
p_map(id2)=2;%% **
id3 = intersect(find(p_map<0.001),find(p_map>0.0001));
p_map(id3)=3;%% ***
id4 = intersect(find(p_map<0.0001),find(p_map>0));
p_map(id4)=4;%% ****
h = heatmap(xvalues,yvalues,p_map);
colormap(gca,bluewhitered(256)); 
colorbar;

h.Title = 'Module ofacued corr';
h.XLabel = 'CR vs AL';
h.YLabel = 'moduleME';
title('CR-AL-module-trait-p');
fig_fn = 'CR_AL_module_ofacued_rp_RBWheatmap.emf'; 

cd('D:\Cortical RNAseq\new_analysis_flowcell_coverage_extractiondate_adjusted\Ana-fig\CRND-wgcna');
saveas(gcf,fig_fn);




















