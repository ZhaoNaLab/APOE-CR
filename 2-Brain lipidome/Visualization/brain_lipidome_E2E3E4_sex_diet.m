%% brain_lipidome_6_modeling_stratified_APOE_Species_sex_diet.m 09262023
%% brain_lipidome.mat pair wised comparison adjust sex. 
%%%% 
close all hidden;
clear;
clc;

%% raw  data
cd('D:\Cotical lipidome\Qiao-analysis\ana_mat');
load('brain_lipidome.mat');
  
%% meta data
mouseinfo = brain_lipid.mouse_info;
n_sample = brain_lipid.n_sample;
n_feature = size(mouseinfo,1);

lipid_class = brain_lipid.class_abbr;
class_abd = brain_lipid.class_level;
class_abd_log = log(class_abd);
n_class = length(lipid_class);

class_species = brain_lipid.lipis_class_abbr;
lipid_species = brain_lipid.lipid_list;
species_abd = brain_lipid.lipid_level+10^(-6);

[xid, yid]= ind2sub(size(species_abd),find(species_abd==0));
filtered_IDs = setdiff(1:brain_lipid.n_lipid,unique(xid));
species_abd_log = log(species_abd(filtered_IDs,:));
n_species = length(filtered_IDs);
filtered_lipid_species = lipid_species(filtered_IDs);
filtered_lipid_class_species = class_species(filtered_IDs);

brain_lipid.E2_IDs = find(contains(brain_lipid.mouse_info(:,1),'E2')); 
brain_lipid.E3_IDs = find(contains(brain_lipid.mouse_info(:,1),'E3')); 
brain_lipid.E4_IDs = find(contains(brain_lipid.mouse_info(:,1),'E4')); 
brain_lipid.mouse_info_E2 = brain_lipid.mouse_info(brain_lipid.E2_IDs,:);
brain_lipid.mouse_info_E3 = brain_lipid.mouse_info(brain_lipid.E3_IDs,:);
brain_lipid.mouse_info_E4 = brain_lipid.mouse_info(brain_lipid.E4_IDs,:);
brain_lipid.class_abd_log_E2 = class_abd_log(:,brain_lipid.E2_IDs);
brain_lipid.class_abd_log_E3 = class_abd_log(:,brain_lipid.E3_IDs);
brain_lipid.class_abd_log_E4 = class_abd_log(:,brain_lipid.E4_IDs);
brain_lipid.species_abd_log_E2 = species_abd_log(:,brain_lipid.E2_IDs);
brain_lipid.species_abd_log_E3 = species_abd_log(:,brain_lipid.E3_IDs);
brain_lipid.species_abd_log_E4 = species_abd_log(:,brain_lipid.E4_IDs);

%% F, anova, ml tables E2
F_table = zeros((n_species),3);
CRAL_lm = zeros((n_species)*2,4);
CRAL_anova = zeros((n_species)*2,5);

for i = 1:n_species
temp = species_abd_log(i,brain_lipid.E2_IDs);
ds = table(brain_lipid.mouse_info_E2(:,1),brain_lipid.mouse_info_E2(:,2),brain_lipid.mouse_info_E2(:,3),temp',...
    'VariableNames',{'apoe' 'diet' 'sex' 'abd'});
ds.diet = ordinal(brain_lipid.mouse_info_E2(:,2));
ds.diet = reordercats(ds.diet,{'AL' 'CR'});
ds.sex = ordinal(brain_lipid.mouse_info_E2(:,3)); %% Male vs Female
ds.sex = reordercats(ds.sex,{'Female' 'Male'});

mdl = fitlm(ds,'abd ~ 1+sex+diet');
terms = mdl.CoefficientNames;
lm_term_oi = terms(2:end);
n_lm_terms_oi = length(lm_term_oi);
coefficients_of_interest = mdl.Coefficients(2:end,:);
CRAL_lm ((i-1)*n_lm_terms_oi+1:i*n_lm_terms_oi,:) = table2array(coefficients_of_interest);

tbl = anova(mdl);
anova_terms = tbl.Properties.RowNames;
anova_terms_oi = anova_terms(1:end-1);
n_anova_terms_oi = length(anova_terms_oi);
CRAL_anova((i-1)*n_anova_terms_oi+1:i*n_anova_terms_oi,:) = table2array(tbl(1:end-1,:));
F_table(i,:) = tbl.F;
end     

%% % plot mean F
mean_F = mean(F_table);
terms = tbl.Properties.RowNames;
figure;
b=bar(mean_F);
set(gca,'XTickLabel',terms);
hold on;
plot([0 n_anova_terms_oi+2],[1.25 1.25],'-k');%% Yingxue uses 1.25 as a threshold (Don't know why?)
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom');
title('E2');

save_fig_path= 'D:\Cotical lipidome\Qiao-analysis\ana_fig';
fig_fn = 'CR_AL_brain_species_Mean_F_sex+diet-E2.emf';
cd(save_fig_path);
saveas(gcf,fig_fn);

%% save anova and lm table
species_list_temp = repmat(filtered_lipid_species(1:end),n_lm_terms_oi,1);
lipid_species_neo = species_list_temp(:);
coef = repmat(lm_term_oi',n_species,1);

fdr = zeros(length(CRAL_lm),1);
hom_bof_p = zeros(length(CRAL_lm),1);

for k = 1:n_lm_terms_oi
temp  = CRAL_lm (k:n_lm_terms_oi:end,end);
fdr(k:n_lm_terms_oi:end) = multicmp(temp,'fdr',0.05);
hom_bof_p(k:n_lm_terms_oi:end) = multicmp(temp,'down',0.05);
end

lm_tbl = table(coef,lipid_species_neo,CRAL_lm(:,1),CRAL_lm(:,2),CRAL_lm(:,3),CRAL_lm(:,4),fdr,hom_bof_p,...
    'VariableNames',{'terms','lipid_species','Estimate','SE','tStat','pValue','BH_p','bof_p'});
cd('D:\Cotical lipidome\Qiao-analysis\ana_mat');
save('CRAL-brain-lipid-species-lm-E2-sex-diet.mat', "lm_tbl");
fn = 'CRAL-brain-lipid-species-lm-E2-sex-diet.xlsx';
writetable (lm_tbl,fn);

species_list_temp = repmat(filtered_lipid_species(1:end),n_anova_terms_oi,1);
lipid_species_neo = species_list_temp(:);
coef = repmat(anova_terms_oi,n_species,1);

fdr = zeros(length(CRAL_anova),1);
hom_bof_p = zeros(length(CRAL_anova),1);

for k = 1:n_anova_terms_oi
temp  = CRAL_anova (k:n_anova_terms_oi:end,end);
fdr(k:n_anova_terms_oi:end) = multicmp(temp,'fdr',0.05);
hom_bof_p(k:n_anova_terms_oi:end) = multicmp(temp,'down',0.05);
end

anova_tbl = table(coef,lipid_species_neo,CRAL_anova(:,1),CRAL_anova(:,2),CRAL_anova(:,3),CRAL_anova(:,4),...
    CRAL_anova (:,5),fdr,hom_bof_p,'VariableNames',{'terms','lipid_species','SumSq','DF','MeanSq','F','pValues','BH_p','bof_p'});
cd('D:\Cotical lipidome\Qiao-analysis\ana_mat');

save('CRAL-brain-lipid-species-anova-E2-sex-diet.mat', "anova_tbl");
fn = 'CRAL-brain-lipid-species-anova-E2-sex-diet.xlsx';
writetable (anova_tbl,fn);

%% F, anova, ml tables E3
F_table = zeros((n_species),3);
CRAL_lm = zeros((n_species)*2,4);
CRAL_anova = zeros((n_species)*2,5);

for i = 1:n_species
temp = species_abd_log(i,brain_lipid.E3_IDs);
ds = table(brain_lipid.mouse_info_E3(:,1),brain_lipid.mouse_info_E3(:,2),brain_lipid.mouse_info_E3(:,3),temp',...
    'VariableNames',{'apoe' 'diet' 'sex' 'abd'});
ds.diet = ordinal(brain_lipid.mouse_info_E3(:,2));
ds.diet = reordercats(ds.diet,{'AL' 'CR'});
ds.sex = ordinal(brain_lipid.mouse_info_E3(:,3)); %% Male vs Female
ds.sex = reordercats(ds.sex,{'Female' 'Male'});

mdl = fitlm(ds,'abd ~ 1+sex+diet');
terms = mdl.CoefficientNames;
lm_term_oi = terms(2:end);
n_lm_terms_oi = length(lm_term_oi);
coefficients_of_interest = mdl.Coefficients(2:end,:);
CRAL_lm ((i-1)*n_lm_terms_oi+1:i*n_lm_terms_oi,:) = table2array(coefficients_of_interest);

tbl = anova(mdl);
anova_terms = tbl.Properties.RowNames;
anova_terms_oi = anova_terms(1:end-1);
n_anova_terms_oi = length(anova_terms_oi);
CRAL_anova((i-1)*n_anova_terms_oi+1:i*n_anova_terms_oi,:) = table2array(tbl(1:end-1,:));
F_table(i,:) = tbl.F;
end     

%% % plot mean F
mean_F = mean(F_table);
terms = tbl.Properties.RowNames;
figure;
b=bar(mean_F);
set(gca,'XTickLabel',terms);
hold on;
plot([0 n_anova_terms_oi+2],[1.25 1.25],'-k');%% Yingxue uses 1.25 as a threshold (Don't know why?)
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom');
title('E3');

save_fig_path= 'D:\Cotical lipidome\Qiao-analysis\ana_fig';
fig_fn = 'CR_AL_brain_species_Mean_F_E3-sex-diet.emf';
cd(save_fig_path);
saveas(gcf,fig_fn);

%% save anova and lm table
species_list_temp = repmat(filtered_lipid_species(1:end),n_lm_terms_oi,1);
lipid_species_neo = species_list_temp(:);
coef = repmat(lm_term_oi',n_species,1);

fdr = zeros(length(CRAL_lm),1);
hom_bof_p = zeros(length(CRAL_lm),1);

for k = 1:n_lm_terms_oi
temp  = CRAL_lm (k:n_lm_terms_oi:end,end);
fdr(k:n_lm_terms_oi:end) = multicmp(temp,'fdr',0.05);
hom_bof_p(k:n_lm_terms_oi:end) = multicmp(temp,'down',0.05);
end

lm_tbl = table(coef,lipid_species_neo,CRAL_lm(:,1),CRAL_lm(:,2),CRAL_lm(:,3),CRAL_lm(:,4),fdr,hom_bof_p,...
    'VariableNames',{'terms','lipid_species','Estimate','SE','tStat','pValue','BH_p','bof_p'});
cd('D:\Cotical lipidome\Qiao-analysis\ana_mat');
save('CRAL-brain-lipid-species-lm-E3-sex-diet.mat', "lm_tbl");
fn = 'CRAL-brain-lipid-species-lm-E3-sex-diet.xlsx';
writetable (lm_tbl,fn);

species_list_temp = repmat(filtered_lipid_species(1:end),n_anova_terms_oi,1);
lipid_species_neo = species_list_temp(:);
coef = repmat(anova_terms_oi,n_species,1);

fdr = zeros(length(CRAL_anova),1);
hom_bof_p = zeros(length(CRAL_anova),1);

for k = 1:n_anova_terms_oi
temp  = CRAL_anova (k:n_anova_terms_oi:end,end);
fdr(k:n_anova_terms_oi:end) = multicmp(temp,'fdr',0.05);
hom_bof_p(k:n_anova_terms_oi:end) = multicmp(temp,'down',0.05);
end

anova_tbl = table(coef,lipid_species_neo,CRAL_anova(:,1),CRAL_anova(:,2),CRAL_anova(:,3),CRAL_anova(:,4),...
    CRAL_anova (:,5),fdr,hom_bof_p,'VariableNames',{'terms','lipid_species','SumSq','DF','MeanSq','F','pValues','BH_p','bof_p'});
cd('D:\Cotical lipidome\Qiao-analysis\ana_mat');

save('CRAL-brain-lipid-species-anova-E3-sex-diet.mat', "anova_tbl");
fn = 'CRAL-brain-lipid-species-anova-E3-sex-diet.xlsx';
writetable (anova_tbl,fn);

%% F, anova, ml tables E4
F_table = zeros((n_species),3);
CRAL_lm = zeros((n_species)*2,4);
CRAL_anova = zeros((n_species)*2,5);

for i = 1:n_species
temp = species_abd_log(i,brain_lipid.E4_IDs);
ds = table(brain_lipid.mouse_info_E4(:,1),brain_lipid.mouse_info_E4(:,2),brain_lipid.mouse_info_E4(:,3),temp',...
    'VariableNames',{'apoe' 'diet' 'sex' 'abd'});
ds.diet = ordinal(brain_lipid.mouse_info_E4(:,2));
ds.diet = reordercats(ds.diet,{'AL' 'CR'});
ds.sex = ordinal(brain_lipid.mouse_info_E4(:,3)); %% Male vs Female
ds.sex = reordercats(ds.sex,{'Female' 'Male'});

mdl = fitlm(ds,'abd ~ 1+sex+diet');
terms = mdl.CoefficientNames;
lm_term_oi = terms(2:end);
n_lm_terms_oi = length(lm_term_oi);
coefficients_of_interest = mdl.Coefficients(2:end,:);
CRAL_lm ((i-1)*n_lm_terms_oi+1:i*n_lm_terms_oi,:) = table2array(coefficients_of_interest);

tbl = anova(mdl);
anova_terms = tbl.Properties.RowNames;
anova_terms_oi = anova_terms(1:end-1);
n_anova_terms_oi = length(anova_terms_oi);
CRAL_anova((i-1)*n_anova_terms_oi+1:i*n_anova_terms_oi,:) = table2array(tbl(1:end-1,:));
F_table(i,:) = tbl.F;
end     

%% % plot mean F
mean_F = mean(F_table);
terms = tbl.Properties.RowNames;
figure;
b=bar(mean_F);
set(gca,'XTickLabel',terms);
hold on;
plot([0 n_anova_terms_oi+2],[1.25 1.25],'-k');%% Yingxue uses 1.25 as a threshold (Don't know why?)
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom');
title('E4');
save_fig_path= 'D:\Cotical lipidome\Qiao-analysis\ana_fig';
fig_fn = 'CR_AL_brain_species_Mean_F_E4-sex-diet.emf';
cd(save_fig_path);
saveas(gcf,fig_fn);

%% save anova and lm table
species_list_temp = repmat(filtered_lipid_species(1:end),n_lm_terms_oi,1);
lipid_species_neo = species_list_temp(:);
coef = repmat(lm_term_oi',n_species,1);

fdr = zeros(length(CRAL_lm),1);
hom_bof_p = zeros(length(CRAL_lm),1);

for k = 1:n_lm_terms_oi
temp  = CRAL_lm (k:n_lm_terms_oi:end,end);
fdr(k:n_lm_terms_oi:end) = multicmp(temp,'fdr',0.05);
hom_bof_p(k:n_lm_terms_oi:end) = multicmp(temp,'down',0.05);
end

lm_tbl = table(coef,lipid_species_neo,CRAL_lm(:,1),CRAL_lm(:,2),CRAL_lm(:,3),CRAL_lm(:,4),fdr,hom_bof_p,...
    'VariableNames',{'terms','lipid_species','Estimate','SE','tStat','pValue','BH_p','bof_p'});
cd('D:\Cotical lipidome\Qiao-analysis\ana_mat');
save('CRAL-brain-lipid-species-lm-E4-sex-diet.mat', "lm_tbl");
fn = 'CRAL-brain-lipid-species-lm-E4-sex-diet.xlsx';
writetable (lm_tbl,fn);

species_list_temp = repmat(filtered_lipid_species(1:end),n_anova_terms_oi,1);
lipid_species_neo = species_list_temp(:);
coef = repmat(anova_terms_oi,n_species,1);

fdr = zeros(length(CRAL_anova),1);
hom_bof_p = zeros(length(CRAL_anova),1);

for k = 1:n_anova_terms_oi
temp  = CRAL_anova (k:n_anova_terms_oi:end,end);
fdr(k:n_anova_terms_oi:end) = multicmp(temp,'fdr',0.05);
hom_bof_p(k:n_anova_terms_oi:end) = multicmp(temp,'down',0.05);
end

anova_tbl = table(coef,lipid_species_neo,CRAL_anova(:,1),CRAL_anova(:,2),CRAL_anova(:,3),CRAL_anova(:,4),...
    CRAL_anova (:,5),fdr,hom_bof_p,'VariableNames',{'terms','lipid_species','SumSq','DF','MeanSq','F','pValues','BH_p','bof_p'});
cd('D:\Cotical lipidome\Qiao-analysis\ana_mat');

save('CRAL-brain-lipid-species-anova-E4-sex-diet.mat', "anova_tbl");
fn = 'CRAL-brain-lipid-species-anova-E4-sex-diet.xlsx';
writetable (anova_tbl,fn);
%% 
cd('D:\Cotical lipidome\Clark-analysis-summary\Results\lm_table_summary');
[a,b,c] = xlsread('lm_table_all.xlsx');

lmterms = b(1,:);
lm_tbl = table(b(2:end,1),a(:,1),a(:,2),a(:,3),a(:,4),b(2:end,6),a(:,6),'VariableNames',lmterms);
cd('D:\Cotical lipidome\Qiao-analysis\ana_mat');

save('CRAL-brain-lipid-species-lm-betaforall.mat', "lm_tbl");
fn = 'CRAL-brain-lipid-species-lm-betaforall.xlsx';
writetable (lm_tbl,fn);


