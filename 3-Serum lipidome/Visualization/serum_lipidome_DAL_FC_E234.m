%% serum_lipid_7_FC_rawadb_species_anoava_allE234.m
%%%%% FC heatmap 
%%%%
close all hidden;
clear;
clc;

%%
mat_path = 'D:\serum lipidome\ana_mat';
cd(mat_path);
load('CRAL-serum-lipid-species-anova.mat');
save_fig_path = 'D:\serum lipidome\ana Fig\modify0';
p_thr = 0.05;

%% CR vs AL all, 86 DAL
terms_anova = unique(anova_tbl.terms);
term_oi = terms_anova(3); %% diet
bof_p = anova_tbl.bof_p(strcmp(anova_tbl.terms,term_oi));
species = anova_tbl.lipid_species(strcmp(anova_tbl.terms,term_oi));
n_species = length(species);

%%%% significant species
sig_IDs = find(bof_p<p_thr);
species_sig = species(sig_IDs);
species_sig_bofp = bof_p(sig_IDs);
n_sig = length(sig_IDs);

%% E2 58 DAL E3 4 DAL E4 8 DAL
cd(mat_path);
load CRAL-serum-lipid-species-anova-E2-sex-diet.mat;
terms_anovae2 = unique(anova_tbl.terms);
term_oi = terms_anovae2(1); %% diet
bof_pe2 = anova_tbl.bof_p(strcmp(anova_tbl.terms,term_oi));
speciese2 = anova_tbl.lipid_species(strcmp(anova_tbl.terms,term_oi));

%%%% significant species
sig_IDse2 = find(bof_pe2<p_thr);
species_sige2 = speciese2(sig_IDse2);
species_sig_bofpe2 = bof_pe2(sig_IDse2);
n_sige2 = length(sig_IDse2);

load CRAL-serum-lipid-species-anova-E3-sex-diet.mat;
terms_anovae3 = unique(anova_tbl.terms);
term_oi = terms_anovae3(1); %% diet
bof_pe3 = anova_tbl.bof_p(strcmp(anova_tbl.terms,term_oi));
speciese3 = anova_tbl.lipid_species(strcmp(anova_tbl.terms,term_oi));

%%%% significant species
sig_IDse3 = find(bof_pe3<p_thr);
species_sige3 = speciese3(sig_IDse3);
species_sig_bofpe3 = bof_pe3(sig_IDse3);
n_sige3 = length(sig_IDse3);

load CRAL-serum-lipid-species-anova-E4-sex-diet.mat;
terms_anovae4 = unique(anova_tbl.terms);
term_oi = terms_anovae4(1); %% diet
bof_pe4 = anova_tbl.bof_p(strcmp(anova_tbl.terms,term_oi));
speciese4 = anova_tbl.lipid_species(strcmp(anova_tbl.terms,term_oi));

%%%% significant species
sig_IDse4 = find(bof_pe4<p_thr);
species_sige4 = speciese4(sig_IDse4);
species_sig_bofpe4 = bof_pe4(sig_IDse4);
n_sige4 = length(sig_IDse4);
%%% do the test below before other codes
%%%% speceis, speciese2, speciese3, speciese4 share the same list.
% [idt,ida,idb] = intersect(speciese2,speciese3,'stable');
% [idt,ida,idb] = intersect(speciese2,speciese4,'stable');
% find(ida==idb);
% [idt,ida,idb] = intersect(species,speciese3,'stable');

%%% 
custom_DAL_list = unique([species_sige2;species_sige3;species_sige4]);
bofp234 = [bof_pe2,bof_pe3,bof_pe4];

%% expression hierarichical heatmap
cd(mat_path);
load serum_lipidome.mat;
lipidabd = serum_lipid.lipid_level;
mouseinfo = serum_lipid.mouse_info_w2;
mouse_id = [1,3:48];
lipid = serum_lipid.lipid_list;

[~, ~, idb] = intersect(custom_DAL_list,lipid,'stable');
dal_abd = lipidabd(idb,:);
row_name = custom_DAL_list; 

collumn_name1 = mouseinfo(1,:);% genotype
collumn_name2 = mouseinfo(2,:);% diet
collumn_name3 = cellfun(@(x) x(1), mouseinfo(3,:),'un',0);% sex
collumn_name4 = num2str(mouse_id');% sample id
collumn_name = strcat(collumn_name4,collumn_name1','-',collumn_name2','-',collumn_name3');

CGobj  = clustergram(dal_abd,'Standardize','row','colormap',redbluecmap);  %% 
set(CGobj,'RowLabels',row_name,'ColumnLabels',collumn_name); 
x_label = CGobj.ColumnLabels; 
[idt0,ida0,idb0] = intersect(x_label,collumn_name,'stable');

new_xlabel = x_label([find(contains(x_label,'E2-AL-F')),find(contains(x_label,'E2-AL-M')),...
    find(contains(x_label,'E2-CR-F')),find(contains(x_label,'E2-CR-M')),...
    find(contains(x_label,'E3-AL-F')),find(contains(x_label,'E3-AL-M')),...
    find(contains(x_label,'E3-CR-F')),find(contains(x_label,'E3-CR-M')),...
    find(contains(x_label,'E4-AL-F')),find(contains(x_label,'E4-AL-M')),...
    find(contains(x_label,'E4-CR-F')),find(contains(x_label,'E4-CR-M'))]);

[~,~,idb] = intersect(new_xlabel,collumn_name,'stable');
dal_abdre = dal_abd(:,idb);

CGobj  = clustergram(dal_abdre,'Standardize','row','colormap',redbluecmap,'cluster','column');  %%

set(CGobj,'RowLabels',row_name,'ColumnLabels',new_xlabel); 
y_label = flip(CGobj.RowLabels);
x_label = CGobj.ColumnLabels; 

%% E234 heatmap FC
cd(mat_path);
load serum_lipidome.mat;
lipidlist = serum_lipid.lipid_list;
fc = serum_lipid.fc;
fce2 = serum_lipid.fce2;
fce3 = serum_lipid.fce3;
fce4 = serum_lipid.fce4;

fccmb = [fce2,fce3,fce4];

%% fc map; bof p map
[idt0,ida0,idb0] = intersect(y_label,lipidlist,'stable');
fc_sigs = fccmb(idb0,:);
[idt00,ida00,idb00] = intersect(y_label,speciese2,'stable');
bofp_sigs234 = bofp234(idb00,:);
[idt000,ida000,idb000] = intersect(y_label,species,'stable');
bofp_sigsall = bof_p(idb000,:);
bofp_sigs = bofp_sigs234;

%% 
figure ('position',[0.0010    0.0410    1.5360    0.7488]*1000);
minfc = min(fc_sigs(:));maxfc = max(fc_sigs(:));
subplot(121);

xvalues = {'E2' 'E3' 'E4'};
yvalues = y_label;
heatmap(xvalues,yvalues,fc_sigs,'MissingDataColor',[0.5 0.5 0.5]); %% ,
colormap(gca,bluewhitered(256)); 
colorbar;
tit = 'APOE234-DAL-FC';
title(tit);
ylabel('FC');
cmap = colormap;

subplot(122);
p_map = -log10(bofp_sigs);

p_map(p_map<-log10(p_thr))=0;
id1 = intersect(find(p_map>-log10(p_thr)),find(p_map<2));
p_map(id1)=1;%% *

id2 = intersect(find(p_map>2),find(p_map<3));
p_map(id2)=2;%% **
id3 = intersect(find(p_map>3),find(p_map<4));
p_map(id3)=3;%% ***
id4 = find(p_map>4);
p_map(id4)=4;%% ****
heatmap(xvalues,yvalues,p_map,'MissingDataColor',[0.5 0.5 0.5]); %% ,
colormap(gca,bluewhitered(256));
colorbar;
tit = 'APOE234-DAL-bofp';
title(tit);
ylabel('bofp');

fig_fn = 'APOE234-DAL-anova-based-rawabd-species-FC-bofp-heatmap.emf';
cd(save_fig_path);
saveas(gcf,fig_fn);




