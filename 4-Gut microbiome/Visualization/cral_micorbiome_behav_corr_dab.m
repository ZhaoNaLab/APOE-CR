%% lifestyle_micorbiome_6_cral_behav_corr_updated.m WQ 05312023
%% correlation between abundance,diversity and behaviroal readout. @ species level

clear;
close all hidden;
clc;

%% load mat
save_mat_path = 'D:\Microbiome\Ana-mat';
save_fig_path = 'D:\Microbiome\Ana-fig';
cd(save_mat_path);
load APOE-microbiome-abundance.mat;
load microbiome_B2_CRAL_cmp.mat;

%% mouse info and abundance for B2
B2_microbiome_mouse_info = abd.B2_microbiome_mouse_info;
B2_microbiome_mouse_id = abd.B2_microbiome_mouse_ID;
B2_abd = abd.species(:,120:end);
species_name = abd.species_name;
dab_CRAL_species = cmp_B2_CRAL.dab_cral_ylabel_noall; 

% B2_abd = abd.genus(:,120:end);
% species_name = abd.genus_name;
% dab_CRAL_species = cmp_B2_CRAL.dab_cral_ylabel_genus; 

[idt0,ida0,idb0] = intersect(dab_CRAL_species,species_name,'stable');

%%  
B2_microbiome_CRAL_id= [find(contains(B2_microbiome_mouse_info(:,2),'AL'));...
    find(startsWith(B2_microbiome_mouse_info(:,2),'C'))]; 
B2_CRAL_info = B2_microbiome_mouse_info(B2_microbiome_CRAL_id,:); 
B2_abd_CRAL = B2_abd(idb0,B2_microbiome_CRAL_id); 

%% 
B2_abd_CRAL_neo = B2_abd_CRAL;
B2_CRAL_species_name = dab_CRAL_species;
n_species_CRAL = length(B2_CRAL_species_name);

%E2
E2_IDs = find(contains(B2_CRAL_info (:,1),'E2'));
B2_CRAL_info_E2 = B2_CRAL_info(E2_IDs,:);
B2_abd_CRAL_neo_E2 = B2_abd_CRAL_neo(:,E2_IDs);

%E3
E3_IDs = find(contains(B2_CRAL_info (:,1),'E3'));
B2_CRAL_info_E3 = B2_CRAL_info(E3_IDs,:);
B2_abd_CRAL_neo_E3 = B2_abd_CRAL_neo(:,E3_IDs);

%E4
E4_IDs = find(contains(B2_CRAL_info (:,1),'E4'));
B2_CRAL_info_E4 = B2_CRAL_info(E4_IDs,:);
B2_abd_CRAL_neo_E4 = B2_abd_CRAL_neo(:,E4_IDs);

%% behaviroral readout
behaviral_readout_path = 'D:\Microbiome\Mouse info\Behavioral readouts- Lipid panels';
cd(behaviral_readout_path);
fn = 'B1-B2-microbiome-behavior-readouts-Lipid-panel.xlsx';
[a1,b1,~] = xlsread(fn,'B2-OFA');
[a2,b2,~] = xlsread(fn,'B2-SA');
[a3,b3,~] = xlsread(fn,'B2-CFC-Contextual');
[a4,b4,~] = xlsread(fn,'B2-CFC-Cued');
[a5,b5,~] = xlsread(fn,'B2-Serum-Lipid Panel');

%%  
OFA_eartag = b1(2:end,5); 
SA_eartag = b2(2:end,5); 
context_eartag = b3(2:end,5); 
cued_eartag = b4(2:end,5); 
lipid_eartag = b5(2:end,5); 

%% 
[idt6,ida6,idb6] = intersect(OFA_eartag,B2_CRAL_info_E2(:,4),'stable');
[idt7,ida7,idb7] = intersect(OFA_eartag,B2_CRAL_info_E3(:,4),'stable');
[idt8,ida8,idb8] = intersect(OFA_eartag,B2_CRAL_info_E4(:,4),'stable');

[idt9,ida9,idb9] = intersect(SA_eartag,B2_CRAL_info_E2(:,4),'stable');
[idt10,ida10,idb10] = intersect(SA_eartag,B2_CRAL_info_E3(:,4),'stable');
[idt11,ida11,idb11] = intersect(SA_eartag,B2_CRAL_info_E4(:,4),'stable');

[idt12,ida12,idb12] = intersect(context_eartag,B2_CRAL_info_E2(:,4),'stable');
[idt13,ida13,idb13] = intersect(context_eartag,B2_CRAL_info_E3(:,4),'stable');
[idt14,ida14,idb14] = intersect(context_eartag,B2_CRAL_info_E4(:,4),'stable');

[idt15,ida15,idb15] = intersect(cued_eartag,B2_CRAL_info_E2(:,4),'stable');
[idt16,ida16,idb16] = intersect(cued_eartag,B2_CRAL_info_E3(:,4),'stable');
[idt17,ida17,idb17] = intersect(cued_eartag,B2_CRAL_info_E4(:,4),'stable');

OFA_B2_CRAL_E2 = a1(ida6,7);OFA_B2_CRAL_E3 = a1(ida7,7);OFA_B2_CRAL_E4 = a1(ida8,7);
OFA_mouse_info_E2 = b1(ida6+1,2:5);OFA_mouse_info_E3 = b1(ida7+1,2:5);OFA_mouse_info_E4 = b1(ida8+1,2:5);
OFA_abd_E2 = B2_abd_CRAL_neo_E2(:,idb6);OFA_abd_E3 = B2_abd_CRAL_neo_E3(:,idb7);
OFA_abd_E4 = B2_abd_CRAL_neo_E4(:,idb8); 

SA_B2_CRAL_E2 = a1(ida9,7);SA_B2_CRAL_E3 = a1(ida10,7);SA_B2_CRAL_E4 = a1(ida11,7);
SA_mouse_info_E2 = b1(ida9+1,2:5);SA_mouse_info_E3 = b1(ida10+1,2:5);SA_mouse_info_E4 = b1(ida11+1,2:5);
SA_abd_E2 = B2_abd_CRAL_neo_E2(:,idb9);SA_abd_E3 = B2_abd_CRAL_neo_E3(:,idb10);
SA_abd_E4 = B2_abd_CRAL_neo_E4(:,idb11); 

context_B2_CRAL_E2 = a1(ida12,7);context_B2_CRAL_E3 = a1(ida13,7);context_B2_CRAL_E4 = a1(ida14,7);
context_mouse_info_E2 = b1(ida12+1,2:5);context_mouse_info_E3 = b1(ida13+1,2:5);context_mouse_info_E4 = b1(ida14+1,2:5);
context_abd_E2 = B2_abd_CRAL_neo_E2(:,idb12);context_abd_E3 = B2_abd_CRAL_neo_E3(:,idb13);
context_abd_E4 = B2_abd_CRAL_neo_E4(:,idb14); 

cued_B2_CRAL_E2 = a1(ida12,7);cued_B2_CRAL_E3 = a1(ida13,7);cued_B2_CRAL_E4 = a1(ida14,7);
cued_mouse_info_E2 = b1(ida12+1,2:5);cued_mouse_info_E3 = b1(ida13+1,2:5);cued_mouse_info_E4 = b1(ida14+1,2:5);
cued_abd_E2 = B2_abd_CRAL_neo_E2(:,idb12);cued_abd_E3 = B2_abd_CRAL_neo_E3(:,idb13);
cued_abd_E4 = B2_abd_CRAL_neo_E4(:,idb14); 

%% 
[idt1,ida1,idb1] = intersect(OFA_eartag,B2_CRAL_info(:,4),'stable');
[idt2,ida2,idb2] = intersect(SA_eartag,B2_CRAL_info(:,4),'stable');
[idt3,ida3,idb3] = intersect(context_eartag,B2_CRAL_info(:,4),'stable');
[idt4,ida4,idb4] = intersect(cued_eartag,B2_CRAL_info(:,4),'stable');
[idt5,ida5,idb5] = intersect(lipid_eartag,B2_CRAL_info(:,4),'stable');

%%
OFA_B2_CRAL = a1(ida1,7);
OFA_mouse_info = b1(ida1+1,2:4);
OFA_abd = B2_abd_CRAL_neo(:,idb1);

SA_B2_CRAL = a2(ida2,7);
SA_mouse_info = b2(ida2+1,2:4);
SA_abd = B2_abd_CRAL_neo(:,idb2);

Context_B2_CRAL =  a3(ida3,8);
Context_mouse_info = b3(ida3+1,2:4);
Context_abd = B2_abd_CRAL_neo(:,idb3);

Cued_B2_CRAL = a4(ida4,8);
Cued_mouse_info = b4(ida4+1,2:4);
Cued_abd = B2_abd_CRAL_neo(:,idb4);

lipid_B2_CRAL = a5(ida5,7:10);
lipid_mouse_info = b5(ida5+1,2:4);
lipid_abd = B2_abd_CRAL_neo(:,idb5);

lipid_terms = b5(1,7:10);

%% correlation
p_nocorrect = zeros(n_species_CRAL,8); %% ofa,sa,cfc-co,cfc-cu,cholesterol, tg, hdl,ldl
fdr_nocorrect = zeros(n_species_CRAL,8); %% ofa,sa,cfc-co,cfc-cu,cholesterol, tg, hdl,ldl
r_nocorrect = zeros(n_species_CRAL,8); %% ofa,sa,cfc-co,cfc-cu,cholesterol, tg, hdl,ldl
% [r_nocorrect(:,1),p_nocorrect(:,1)] = corr(OFA_abd',OFA_B2_CRAL,'type','Spearman');
% [r_nocorrect(:,2),p_nocorrect(:,2)] = corr(SA_abd',SA_B2_CRAL,'type','Spearman');
% [r_nocorrect(:,3),p_nocorrect(:,3)] = corr(Context_abd',Context_B2_CRAL,'type','Spearman');
% [r_nocorrect(:,4),p_nocorrect(:,4)] = corr(Cued_abd',Cued_B2_CRAL,'type','Spearman');
% [r_nocorrect(:,5:8),p_nocorrect(:,5:8)] = corr(lipid_abd',lipid_B2_CRAL,'type','Spearman');

[r_nocorrect(:,1),p_nocorrect(:,1)] = corr(OFA_abd',OFA_B2_CRAL,'type','Pearson');
[r_nocorrect(:,2),p_nocorrect(:,2)] = corr(SA_abd',SA_B2_CRAL,'type','Pearson');
[r_nocorrect(:,3),p_nocorrect(:,3)] = corr(Context_abd',Context_B2_CRAL,'type','Pearson');
[r_nocorrect(:,4),p_nocorrect(:,4)] = corr(Cued_abd',Cued_B2_CRAL,'type','Pearson');
[r_nocorrect(:,5:8),p_nocorrect(:,5:8)] = corr(lipid_abd',lipid_B2_CRAL,'type','Pearson');


for i = 1:8
fdr_nocorrect(:,i) = multicmp(p_nocorrect(:,i),'down',0.05);
end

r_map = [r_nocorrect(:,1) r_nocorrect(:,4)];
figure ('position',[0.0010    0.0418    1.5336    0.7408]*1000);
subplot(121);
xvalues = {'OFA' 'Cued'};
% yvalues = cmp_B2_CRAL.dab_cral_ylabel_genus_neo;
yvalues = cmp_B2_CRAL.dab_cral_ylabel_noall_neo;

heatmap(xvalues,yvalues,r_map,'MissingDataColor',[0.5 0.5 0.5]); %% ,
colormap(gca,bluewhitered(256)); 
colorbar;
% ylabel('correlation coefficient (spearman)');
ylabel('correlation coefficient (pearson)');

tit = 'r';
title(tit);

subplot(122);
xvalues = {'OFA' 'Cued'};
% p_map = [fdr_nocorrect(:,1) fdr_nocorrect(:,4)];
p_map = [p_nocorrect(:,1) p_nocorrect(:,4)];
p_map(p_map>0.05)=0;
id1 = intersect(find(p_map<0.05),find(p_map>0.01));
p_map(id1)=1;%% *
id2 = intersect(find(p_map<0.01),find(p_map>0.001));
p_map(id2)=2;%% **
id3 = intersect(find(p_map<0.001),find(p_map>0.0001));
p_map(id3)=3;%% ***
id4 = intersect(find(p_map<0.0001),find(p_map>0));
p_map(id4)=4;%% ****

heatmap(xvalues,yvalues,p_map,'MissingDataColor',[0.5 0.5 0.5]); %% ,
colormap(gca,bluewhitered(256));
colorbar;
ylabel('p');
tit = 'p';
title(tit);
% fig_fn = 'CRAL-abd-behav-corr-spearman_bof_nocorrection_dabnoall.emf';
% fig_fn = 'CRAL-genus-abd-behav-corr-spearman_fdr_nocorrection_dabnoall.emf';
% fig_fn = 'CRAL-abd-behav-corr-pearson_bof_nocorrection_dabnoall.emf';
fig_fn = 'CRAL-abd-behav-corr-pearson_rawp_nocorrection_dabnoall.emf';

suptitle(fig_fn(1:end-4));

cd(save_fig_path);
saveas(gcf,fig_fn);

%%  adjust apoe
p_correct = zeros(n_species_CRAL,2); %% ofa,cfc-cu,
fdr_correct = zeros(n_species_CRAL,2); %% ofa,cfc-cu,
r_correct = zeros(n_species_CRAL,2); %% ofa,cfc-cu,

%% correction of behaviors and abundance
ds = table(OFA_mouse_info(:,1),OFA_mouse_info(:,2),...
    OFA_mouse_info(:,3),OFA_B2_CRAL,'VariableNames',[abd.microbiome_mouse_info_terms(1:3) 'OFA']);

% mdl = fitlm(ds,'OFA ~ 1+apoe');
mdl = fitlm(ds,'OFA ~ 1+apoe*diet*sex');
temp2 = mdl.Residuals(:,1);
OFA_B2_CRAL_resid= table2array(temp2);

OFA_abd_resid = zeros(size(OFA_abd));
for i = 1:n_species_CRAL
temp1 = OFA_abd(i,:)';
ds = table(OFA_mouse_info(:,1),OFA_mouse_info(:,2),...
    OFA_mouse_info(:,3),temp1, 'VariableNames',[abd.microbiome_mouse_info_terms(1:3) 'abundance']);
% mdl = fitlm(ds,'abundance ~ 1+apoe');
mdl = fitlm(ds,'abundance ~ 1+apoe*diet*sex');

temp2 = mdl.Residuals(:,1);
OFA_abd_resid (i,:) = table2array(temp2);
end

ds = table(Cued_mouse_info(:,1),Cued_mouse_info(:,2),...
    Cued_mouse_info(:,3),Cued_B2_CRAL,'VariableNames',[abd.microbiome_mouse_info_terms(1:3) 'Cued']);
% mdl = fitlm(ds,'Cued ~ 1+apoe');
mdl = fitlm(ds,'Cued ~ 1+apoe*diet*sex');

temp2 = mdl.Residuals(:,1);
Cued_B2_CRAL_resid= table2array(temp2);

Cued_abd_resid = zeros(size(Cued_abd));
for i = 1:n_species_CRAL
temp1 = Cued_abd(i,:)';
ds = table(Cued_mouse_info(:,1),Cued_mouse_info(:,2),...
    Cued_mouse_info(:,3),temp1, 'VariableNames',[abd.microbiome_mouse_info_terms(1:3) 'abundance']);
% mdl = fitlm(ds,'abundance ~ 1+apoe');
mdl = fitlm(ds,'abundance ~ 1+apoe*diet*sex');

temp2 = mdl.Residuals(:,1);
Cued_abd_resid (i,:) = table2array(temp2);
end

[r_correct(:,1),p_correct(:,1)] = corr(OFA_abd_resid',OFA_B2_CRAL_resid,"type","Spearman");
[r_correct(:,2),p_correct(:,2)] = corr(Cued_abd_resid',Cued_B2_CRAL_resid,"type","Spearman");

for i = 1:2
fdr_correct(:,i) = multicmp(p_correct(:,i),'fdr',0.05);
end

%% heatmap for r and p
r_map = r_correct(:,1:2);
figure ('position',[0.0010    0.0418    1.5336    0.7408]*1000);
subplot(121);
xvalues = {'OFA' 'Cued'};

heatmap(xvalues,yvalues,r_map,'MissingDataColor',[0.5 0.5 0.5]); %% ,
colormap(gca,bluewhitered(256)); 
colorbar;
ylabel('correlation coefficient (spearman)');

tit = 'r';
title(tit);

subplot(122);
xvalues = {'OFA' 'Cued'};
p_map = p_correct(:,1:2);
% p_map = fdr_correct(:,1:2);

p_map(p_map>0.05)=0;
id1 = intersect(find(p_map<0.05),find(p_map>0.01));
p_map(id1)=1;%% *
id2 = intersect(find(p_map<0.01),find(p_map>0.001));
p_map(id2)=2;%% **
id3 = intersect(find(p_map<0.001),find(p_map>0.0001));
p_map(id3)=3;%% ***
id4 = intersect(find(p_map<0.0001),find(p_map>0));
p_map(id4)=4;%% ****

heatmap(xvalues,yvalues,p_map,'MissingDataColor',[0.5 0.5 0.5]); %%
colormap(gca,bluewhitered(256)); 
colorbar;
ylabel('p');
tit = 'p';
title(tit);
% % fig_fn = 'CRAL-abd-behav-corr-Spearman-apoe-residue-nobh-dabnoall.emf';
fig_fn = 'CRAL-abd-behav-corr-Spearman-apoe+diet+sex-residue-nobh-dabnoall.emf';



% % fig_fn = 'CRAL-abd-behav-corr-Spearman-apoe-residue-bh-dabnoall.emf';
% % fig_fn = 'CRAL-abd-behav-corr-Spearman-residue-nobh-genus.emf';
% % fig_fn = 'CRAL-abd-behav-corr-Spearman-residue-bh-genus.emf';
% 
suptitle(fig_fn(1:end-4));

cd(save_fig_path);
saveas(gcf,fig_fn);

% 



