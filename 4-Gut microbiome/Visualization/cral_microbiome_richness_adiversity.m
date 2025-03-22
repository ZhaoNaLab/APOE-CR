%% lifestyle_micorbiome_2_richness_adiversity_B1.m WQ 05242023
%% alpha diversity calculation at different levels in different batches. 
%% observed alpha diversity or richness index (R), 
%% shannon index (H), evenness (EH),inversed simpson index (D) will be calculated. 
%% there are six levels. 

clear;
close all;
clc;

%% 
save_fig_path = 'D:\Microbiome\Ana-fig';
mat_path = 'D:\Microbiome\Ana-mat';
cd(mat_path);
load APOE-microbiome-abundance.mat;
load APOE-microbiome-count.mat;

%%
sample_info = [abd.B1_microbiome_mouse_info;abd.B2_microbiome_mouse_info];
sample_name = abd.sample_name';
levels = {'phylum' 'class' 'order' 'family' 'genus' 'species'};

%% Richness R = number of types; shannon index (H) = -sum(pi*lnpi); Evennese (EH) = H / ln(S); InvSimpson index (D) = 1/sum((pi2))
% if a species was not detected in all samples, it will be excluded. 
% checking: the number of types and abd values were the same with Stephen's
% numbers in B1.B2 has less types at different levels.
n_sample = length(abd.phylum); 
R = zeros(n_sample,6); %% phylum, class, order, family, genus, species.
H = zeros(n_sample,6); %% phylum, class, order, family, genus, species.
EH = zeros(n_sample,6); %% phylum, class, order, family, genus, species.
D = zeros(n_sample,6); %% phylum, class, order, family, genus, species.

for i = 1:n_sample
    temp2 = abd.phylum(:,i);
    temp3 = cnt.phylum_count(:,i);
R (i,1) = length(find(temp3));
H (i,1) = -sum(temp2(temp2~=0).*log(temp2(temp2~=0)));
EH(i,1) = H (i,1)/log(sum(temp3));
D (i,1) = 1/sum(temp2.^2); 

    temp2 = abd.class(:,i);
    temp3 = cnt.class_count(:,i);
R (i,2) = length(find(temp3));
H (i,2) = -sum(temp2(temp2~=0).*log(temp2(temp2~=0)));
EH(i,2) = H (i,1)/log(sum(temp3));
D (i,2) = 1/sum(temp2.^2); 

    temp2 = abd.family(:,i);
    temp3 = cnt.family_count(:,i);
R (i,3) = length(find(temp3));
H (i,3) = -sum(temp2(temp2~=0).*log(temp2(temp2~=0)));
EH(i,3) = H (i,3)/log(sum(temp3));
D (i,3) = 1/sum(temp2.^2); 

    temp2 = abd.genus(:,i);
    temp3 = cnt.genus_count(:,i);
R (i,4) = length(find(temp3));
H (i,4) = -sum(temp2(temp2~=0).*log(temp2(temp2~=0)));
EH(i,4) = H (i,4)/log(sum(temp3));
D (i,4) = 1/sum(temp2.^2); 

    temp2 = abd.order(:,i);
    temp3 = cnt.order_count(:,i);
R (i,5) = length(find(temp3));
H (i,5) = -sum(temp2(temp2~=0).*log(temp2(temp2~=0)));
EH(i,5) = H (i,5)/log(sum(temp3));
D (i,5) = 1/sum(temp2.^2); 

    temp2 = abd.species(:,i);
    temp3 = cnt.species_count(:,i);
R (i,6) = length(find(temp3));
H (i,6) = -sum(temp2(temp2~=0).*log(temp2(temp2~=0)));
EH(i,6) = H (i,6)/log(sum(temp3));
D (i,6) = 1/sum(temp2.^2); 
end 

%% save to abundance mat
abd.richness_R = R;
abd.shannon_diversity_H = H;
abd.evenness_EH = EH;
abd.invSimpson_D = D;

cd('D:\Microbiome\Ana-mat');
mat_fn = 'APOE-microbiome-anundance.mat';
save(mat_fn,"abd");


%% Figures at different levels B1
temp1 = abd.B1_microbiome_mouse_info;
E2_AL_M_id = intersect(intersect(find(contains(temp1(:,1),'E2')),find(contains(temp1(:,2),'AL'))),find(contains(temp1(:,3),'Male')));
E2_AL_F_id = intersect(intersect(find(contains(temp1(:,1),'E2')),find(contains(temp1(:,2),'AL'))),find(contains(temp1(:,3),'Female')));
R_E2_AL_M = R (E2_AL_M_id,:); R_E2_AL_F = R (E2_AL_F_id,:); R_E2_AL = [R_E2_AL_M; R_E2_AL_F];
H_E2_AL_M = H (E2_AL_M_id,:); H_E2_AL_F = H (E2_AL_F_id,:); H_E2_AL = [H_E2_AL_M; H_E2_AL_F];
EH_E2_AL_M = EH (E2_AL_M_id,:); EH_E2_AL_F = EH (E2_AL_F_id,:); EH_E2_AL = [EH_E2_AL_M; EH_E2_AL_F];
D_E2_AL_M = D (E2_AL_M_id,:); D_E2_AL_F = D (E2_AL_F_id,:); D_E2_AL = [D_E2_AL_M; D_E2_AL_F];
[R_E2_AL_ave, R_E2_AL_err] = mean_se([R_E2_AL_M; R_E2_AL_F]);
[H_E2_AL_ave, H_E2_AL_err] = mean_se([H_E2_AL_M; H_E2_AL_F]);
[EH_E2_AL_ave, EH_E2_AL_err] = mean_se([EH_E2_AL_M; EH_E2_AL_F]);
[D_E2_AL_ave, D_E2_AL_err] = mean_se([D_E2_AL_M; D_E2_AL_F]);
E2_CR_M_id = intersect(intersect(find(contains(temp1(:,1),'E2')),find(startsWith(temp1(:,2),'CR'))),find(contains(temp1(:,3),'Male')));
E2_CR_F_id = intersect(intersect(find(contains(temp1(:,1),'E2')),find(startsWith(temp1(:,2),'CR'))),find(contains(temp1(:,3),'Female')));
R_E2_CR_M = R (E2_CR_M_id,:); R_E2_CR_F = R (E2_CR_F_id,:); R_E2_CR = [R_E2_CR_M; R_E2_CR_F];
[R_E2_CR_ave, R_E2_CR_err] = mean_se([R_E2_CR_M; R_E2_CR_F]); 
H_E2_CR_M = H (E2_CR_M_id,:); H_E2_CR_F = H (E2_CR_F_id,:); H_E2_CR = [H_E2_CR_M; H_E2_CR_F];
[H_E2_CR_ave, H_E2_CR_err] = mean_se([H_E2_CR_M; H_E2_CR_F]); 
EH_E2_CR_M = EH (E2_CR_M_id,:); EH_E2_CR_F = EH (E2_CR_F_id,:); EH_E2_CR = [EH_E2_CR_M; EH_E2_CR_F];
[EH_E2_CR_ave, EH_E2_CR_err] = mean_se([EH_E2_CR_M; EH_E2_CR_F]); 
D_E2_CR_M = D (E2_CR_M_id,:); D_E2_CR_F = D (E2_CR_F_id,:); D_E2_CR = [D_E2_CR_M; D_E2_CR_F];
[D_E2_CR_ave, D_E2_CR_err] = mean_se([D_E2_CR_M; D_E2_CR_F]); 
E2_HFD_M_id = intersect(intersect(find(contains(temp1(:,1),'E2')),find(contains(temp1(:,2),'H'))),find(contains(temp1(:,3),'Male')));
E2_HFD_F_id = intersect(intersect(find(contains(temp1(:,1),'E2')),find(contains(temp1(:,2),'H'))),find(contains(temp1(:,3),'Female')));
R_E2_HFD_M = R (E2_HFD_M_id,:); R_E2_HFD_F = R (E2_HFD_F_id,:); R_E2_HFD = [R_E2_HFD_M; R_E2_HFD_F];
[R_E2_HFD_ave, R_E2_HFD_err] = mean_se([R_E2_HFD_M; R_E2_HFD_F]);
H_E2_HFD_M = H (E2_HFD_M_id,:); H_E2_HFD_F = H (E2_HFD_F_id,:); H_E2_HFD = [H_E2_HFD_M; H_E2_HFD_F];
[H_E2_HFD_ave, H_E2_HFD_err] = mean_se([H_E2_HFD_M; H_E2_HFD_F]);
EH_E2_HFD_M = EH (E2_HFD_M_id,:); EH_E2_HFD_F = EH (E2_HFD_F_id,:); EH_E2_HFD = [EH_E2_HFD_M; EH_E2_HFD_F];
[EH_E2_HFD_ave, EH_E2_HFD_err] = mean_se([EH_E2_HFD_M; EH_E2_HFD_F]);
D_E2_HFD_M = D (E2_HFD_M_id,:); D_E2_HFD_F = D (E2_HFD_F_id,:); D_E2_HFD = [D_E2_HFD_M; D_E2_HFD_F];
[D_E2_HFD_ave, D_E2_HFD_err] = mean_se([D_E2_HFD_M; D_E2_HFD_F]);

E3_AL_M_id = intersect(intersect(find(contains(temp1(:,1),'E3')),find(contains(temp1(:,2),'AL'))),find(contains(temp1(:,3),'Male')));
E3_AL_F_id = intersect(intersect(find(contains(temp1(:,1),'E3')),find(contains(temp1(:,2),'AL'))),find(contains(temp1(:,3),'Female')));
R_E3_AL_M = R (E3_AL_M_id,:); R_E3_AL_F = R (E3_AL_F_id,:); R_E3_AL = [R_E3_AL_M; R_E3_AL_F];
H_E3_AL_M = H (E3_AL_M_id,:); H_E3_AL_F = H (E3_AL_F_id,:); H_E3_AL = [H_E3_AL_M; H_E3_AL_F];
EH_E3_AL_M = EH (E3_AL_M_id,:); EH_E3_AL_F = EH (E3_AL_F_id,:); EH_E3_AL = [EH_E3_AL_M; EH_E3_AL_F];
D_E3_AL_M = D (E3_AL_M_id,:); D_E3_AL_F = D (E3_AL_F_id,:); D_E3_AL = [D_E3_AL_M; D_E3_AL_F];
[R_E3_AL_ave, R_E3_AL_err] = mean_se([R_E3_AL_M; R_E3_AL_F]);
[H_E3_AL_ave, H_E3_AL_err] = mean_se([H_E3_AL_M; H_E3_AL_F]);
[EH_E3_AL_ave, EH_E3_AL_err] = mean_se([EH_E3_AL_M; EH_E3_AL_F]);
[D_E3_AL_ave, D_E3_AL_err] = mean_se([D_E3_AL_M; D_E3_AL_F]);
E3_CR_M_id = intersect(intersect(find(contains(temp1(:,1),'E3')),find(startsWith(temp1(:,2),'CR'))),find(contains(temp1(:,3),'Male')));
E3_CR_F_id = intersect(intersect(find(contains(temp1(:,1),'E3')),find(startsWith(temp1(:,2),'CR'))),find(contains(temp1(:,3),'Female')));
R_E3_CR_M = R (E3_CR_M_id,:); R_E3_CR_F = R (E3_CR_F_id,:); R_E3_CR = [R_E3_CR_M; R_E3_CR_F];
[R_E3_CR_ave, R_E3_CR_err] = mean_se([R_E3_CR_M; R_E3_CR_F]); 
H_E3_CR_M = H (E3_CR_M_id,:); H_E3_CR_F = H (E3_CR_F_id,:); H_E3_CR = [H_E3_CR_M; H_E3_CR_F];
[H_E3_CR_ave, H_E3_CR_err] = mean_se([H_E3_CR_M; H_E3_CR_F]); 
EH_E3_CR_M = EH (E3_CR_M_id,:); EH_E3_CR_F = EH (E3_CR_F_id,:); EH_E3_CR = [EH_E3_CR_M; EH_E3_CR_F];
[EH_E3_CR_ave, EH_E3_CR_err] = mean_se([EH_E3_CR_M; EH_E3_CR_F]); 
D_E3_CR_M = D (E3_CR_M_id,:); D_E3_CR_F = D (E3_CR_F_id,:); D_E3_CR = [D_E3_CR_M; D_E3_CR_F];
[D_E3_CR_ave, D_E3_CR_err] = mean_se([D_E3_CR_M; D_E3_CR_F]); 
E3_HFD_M_id = intersect(intersect(find(contains(temp1(:,1),'E3')),find(contains(temp1(:,2),'H'))),find(contains(temp1(:,3),'Male')));
E3_HFD_F_id = intersect(intersect(find(contains(temp1(:,1),'E3')),find(contains(temp1(:,2),'H'))),find(contains(temp1(:,3),'Female')));
R_E3_HFD_M = R (E3_HFD_M_id,:); R_E3_HFD_F = R (E3_HFD_F_id,:); R_E3_HFD = [R_E3_HFD_M; R_E3_HFD_F];
[R_E3_HFD_ave, R_E3_HFD_err] = mean_se([R_E3_HFD_M; R_E3_HFD_F]);
H_E3_HFD_M = H (E3_HFD_M_id,:); H_E3_HFD_F = H (E3_HFD_F_id,:); H_E3_HFD = [H_E3_HFD_M; H_E3_HFD_F];
[H_E3_HFD_ave, H_E3_HFD_err] = mean_se([H_E3_HFD_M; H_E3_HFD_F]);
EH_E3_HFD_M = EH (E3_HFD_M_id,:); EH_E3_HFD_F = EH (E3_HFD_F_id,:); EH_E3_HFD = [EH_E3_HFD_M; EH_E3_HFD_F];
[EH_E3_HFD_ave, EH_E3_HFD_err] = mean_se([EH_E3_HFD_M; EH_E3_HFD_F]);
D_E3_HFD_M = D (E3_HFD_M_id,:); D_E3_HFD_F = D (E3_HFD_F_id,:); D_E3_HFD = [D_E3_HFD_M; D_E3_HFD_F];
[D_E3_HFD_ave, D_E3_HFD_err] = mean_se([D_E3_HFD_M; D_E3_HFD_F]);

E4_AL_M_id = intersect(intersect(find(contains(temp1(:,1),'E4')),find(contains(temp1(:,2),'AL'))),find(contains(temp1(:,3),'Male')));
E4_AL_F_id = intersect(intersect(find(contains(temp1(:,1),'E4')),find(contains(temp1(:,2),'AL'))),find(contains(temp1(:,3),'Female')));
R_E4_AL_M = R (E4_AL_M_id,:); R_E4_AL_F = R (E4_AL_F_id,:); R_E4_AL = [R_E4_AL_M; R_E4_AL_F];
H_E4_AL_M = H (E4_AL_M_id,:); H_E4_AL_F = H (E4_AL_F_id,:); H_E4_AL = [H_E4_AL_M; H_E4_AL_F];
EH_E4_AL_M = EH (E4_AL_M_id,:); EH_E4_AL_F = EH (E4_AL_F_id,:); EH_E4_AL = [EH_E4_AL_M; EH_E4_AL_F];
D_E4_AL_M = D (E4_AL_M_id,:); D_E4_AL_F = D (E4_AL_F_id,:); D_E4_AL = [D_E4_AL_M; D_E4_AL_F];
[R_E4_AL_ave, R_E4_AL_err] = mean_se([R_E4_AL_M; R_E4_AL_F]);
[H_E4_AL_ave, H_E4_AL_err] = mean_se([H_E4_AL_M; H_E4_AL_F]);
[EH_E4_AL_ave, EH_E4_AL_err] = mean_se([EH_E4_AL_M; EH_E4_AL_F]);
[D_E4_AL_ave, D_E4_AL_err] = mean_se([D_E4_AL_M; D_E4_AL_F]);
E4_CR_M_id = intersect(intersect(find(contains(temp1(:,1),'E4')),find(startsWith(temp1(:,2),'CR'))),find(contains(temp1(:,3),'Male')));
E4_CR_F_id = intersect(intersect(find(contains(temp1(:,1),'E4')),find(startsWith(temp1(:,2),'CR'))),find(contains(temp1(:,3),'Female')));
R_E4_CR_M = R (E4_CR_M_id,:); R_E4_CR_F = R (E4_CR_F_id,:); R_E4_CR = [R_E4_CR_M; R_E4_CR_F];
[R_E4_CR_ave, R_E4_CR_err] = mean_se([R_E4_CR_M; R_E4_CR_F]); 
H_E4_CR_M = H (E4_CR_M_id,:); H_E4_CR_F = H (E4_CR_F_id,:); H_E4_CR = [H_E4_CR_M; H_E4_CR_F];
[H_E4_CR_ave, H_E4_CR_err] = mean_se([H_E4_CR_M; H_E4_CR_F]); 
EH_E4_CR_M = EH (E4_CR_M_id,:); EH_E4_CR_F = EH (E4_CR_F_id,:); EH_E4_CR = [EH_E4_CR_M; EH_E4_CR_F];
[EH_E4_CR_ave, EH_E4_CR_err] = mean_se([EH_E4_CR_M; EH_E4_CR_F]); 
D_E4_CR_M = D (E4_CR_M_id,:); D_E4_CR_F = D (E4_CR_F_id,:); D_E4_CR = [D_E4_CR_M; D_E4_CR_F];
[D_E4_CR_ave, D_E4_CR_err] = mean_se([D_E4_CR_M; D_E4_CR_F]); 
E4_HFD_M_id = intersect(intersect(find(contains(temp1(:,1),'E4')),find(contains(temp1(:,2),'H'))),find(contains(temp1(:,3),'Male')));
E4_HFD_F_id = intersect(intersect(find(contains(temp1(:,1),'E4')),find(contains(temp1(:,2),'H'))),find(contains(temp1(:,3),'Female')));
R_E4_HFD_M = R (E4_HFD_M_id,:); R_E4_HFD_F = R (E4_HFD_F_id,:); R_E4_HFD = [R_E4_HFD_M; R_E4_HFD_F];
[R_E4_HFD_ave, R_E4_HFD_err] = mean_se([R_E4_HFD_M; R_E4_HFD_F]);
H_E4_HFD_M = H (E4_HFD_M_id,:); H_E4_HFD_F = H (E4_HFD_F_id,:); H_E4_HFD = [H_E4_HFD_M; H_E4_HFD_F];
[H_E4_HFD_ave, H_E4_HFD_err] = mean_se([H_E4_HFD_M; H_E4_HFD_F]);
EH_E4_HFD_M = EH (E4_HFD_M_id,:); EH_E4_HFD_F = EH (E4_HFD_F_id,:); EH_E4_HFD = [EH_E4_HFD_M; EH_E4_HFD_F];
[EH_E4_HFD_ave, EH_E4_HFD_err] = mean_se([EH_E4_HFD_M; EH_E4_HFD_F]);
D_E4_HFD_M = D (E4_HFD_M_id,:); D_E4_HFD_F = D (E4_HFD_F_id,:); D_E4_HFD = [D_E4_HFD_M; D_E4_HFD_F];
[D_E4_HFD_ave, D_E4_HFD_err] = mean_se([D_E4_HFD_M; D_E4_HFD_F]);
%%
figure ('position',[1 41 1536 748.8]);
for i = 1:6
subplot(2,3,i);
h=bar([1 2 3],[R_E2_AL_ave(i),R_E2_CR_ave(i),R_E2_HFD_ave(i);R_E3_AL_ave(i),R_E3_CR_ave(i),R_E3_HFD_ave(i);...
    R_E4_AL_ave(i),R_E4_CR_ave(i),R_E4_HFD_ave(i);]);
hold on;
h(1).FaceColor='none';
h(2).FaceColor='g';
h(3).FaceColor='b';
errorbar(h(1).XEndPoints,[R_E2_AL_ave(i),R_E3_AL_ave(i),R_E4_AL_ave(i)],...
    [R_E2_AL_err(i),R_E3_AL_err(i),R_E4_AL_err(i)],'LineStyle','none','Color','k','LineWidth',1);
errorbar(h(2).XEndPoints,[R_E2_CR_ave(i),R_E3_CR_ave(i),R_E4_CR_ave(i)],...
    [R_E2_CR_err(i),R_E3_CR_err(i),R_E4_CR_err(i)],'LineStyle','none','Color','k','LineWidth',1);
errorbar(h(3).XEndPoints,[R_E2_HFD_ave(i),R_E3_HFD_ave(i),R_E4_HFD_ave(i)],...
    [R_E2_HFD_err(i),R_E3_HFD_err(i),R_E4_HFD_err(i)],'LineStyle','none','Color','k','LineWidth',1);

% jitter. Starting in R2020b you can do this easily with scatter:
scatter(repmat(h(1).XEndPoints(1), length(R_E2_AL_M(:,i)),1),R_E2_AL_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(1).XEndPoints(1), length(R_E2_AL_F(:,i)),1),R_E2_AL_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(2).XEndPoints(1), length(R_E2_CR_M(:,i)),1),R_E2_CR_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(2).XEndPoints(1), length(R_E2_CR_F(:,i)),1),R_E2_CR_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(3).XEndPoints(1), length(R_E2_HFD_M(:,i)),1),R_E2_HFD_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(3).XEndPoints(1), length(R_E2_HFD_F(:,i)),1),R_E2_HFD_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);

scatter(repmat(h(1).XEndPoints(2), length(R_E3_AL_M(:,i)),1),R_E3_AL_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(1).XEndPoints(2), length(R_E3_AL_F(:,i)),1),R_E3_AL_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(2).XEndPoints(2), length(R_E3_CR_M(:,i)),1),R_E3_CR_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(2).XEndPoints(2), length(R_E3_CR_F(:,i)),1),R_E3_CR_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(3).XEndPoints(2), length(R_E3_HFD_M(:,i)),1),R_E3_HFD_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(3).XEndPoints(2), length(R_E3_HFD_F(:,i)),1),R_E3_HFD_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);

scatter(repmat(h(1).XEndPoints(3), length(R_E4_AL_M(:,i)),1),R_E4_AL_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(1).XEndPoints(3), length(R_E4_AL_F(:,i)),1),R_E4_AL_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(2).XEndPoints(3), length(R_E4_CR_M(:,i)),1),R_E4_CR_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(2).XEndPoints(3), length(R_E4_CR_F(:,i)),1),R_E4_CR_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(3).XEndPoints(3), length(R_E4_HFD_M(:,i)),1),R_E4_HFD_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(3).XEndPoints(3), length(R_E4_HFD_F(:,i)),1),R_E4_HFD_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);

xticklabels(["E2" "E3" "E4"]);
tit = sprintf('Richness-%s',char(levels(i)));
ylabel (tit);
set(gca,'TickDir','out','box','off');
end

fig_fn = 'B1-Richness-oberved alpha-6 levels.emf';
suptitle(fig_fn(1:end-4));
cd('D:\Microbiome\Ana-fig\'); 
saveas(gcf,fig_fn);

%%
figure ('position',[1 41 1536 748.8]);
for i = 1:6
subplot(2,3,i);
h=bar([1 2 3],[H_E2_AL_ave(i),H_E2_CR_ave(i),H_E2_HFD_ave(i);H_E3_AL_ave(i),H_E3_CR_ave(i),H_E3_HFD_ave(i);...
    H_E4_AL_ave(i),H_E4_CR_ave(i),H_E4_HFD_ave(i);]);
hold on;
h(1).FaceColor='none';
h(2).FaceColor='g';
h(3).FaceColor='b';
errorbar(h(1).XEndPoints,[H_E2_AL_ave(i),H_E3_AL_ave(i),H_E4_AL_ave(i)],...
    [H_E2_AL_err(i),H_E3_AL_err(i),H_E4_AL_err(i)],'LineStyle','none','Color','k','LineWidth',1);
errorbar(h(2).XEndPoints,[H_E2_CR_ave(i),H_E3_CR_ave(i),H_E4_CR_ave(i)],...
    [H_E2_CR_err(i),H_E3_CR_err(i),H_E4_CR_err(i)],'LineStyle','none','Color','k','LineWidth',1);
errorbar(h(3).XEndPoints,[H_E2_HFD_ave(i),H_E3_HFD_ave(i),H_E4_HFD_ave(i)],...
    [H_E2_HFD_err(i),H_E3_HFD_err(i),H_E4_HFD_err(i)],'LineStyle','none','Color','k','LineWidth',1);

% jitter. Starting in R2020b you can do this easily with scatter:
scatter(repmat(h(1).XEndPoints(1), length(H_E2_AL_M(:,i)),1),H_E2_AL_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(1).XEndPoints(1), length(H_E2_AL_F(:,i)),1),H_E2_AL_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(2).XEndPoints(1), length(H_E2_CR_M(:,i)),1),H_E2_CR_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(2).XEndPoints(1), length(H_E2_CR_F(:,i)),1),H_E2_CR_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(3).XEndPoints(1), length(H_E2_HFD_M(:,i)),1),H_E2_HFD_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(3).XEndPoints(1), length(H_E2_HFD_F(:,i)),1),H_E2_HFD_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);

scatter(repmat(h(1).XEndPoints(2), length(H_E3_AL_M(:,i)),1),H_E3_AL_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(1).XEndPoints(2), length(H_E3_AL_F(:,i)),1),H_E3_AL_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(2).XEndPoints(2), length(H_E3_CR_M(:,i)),1),H_E3_CR_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(2).XEndPoints(2), length(H_E3_CR_F(:,i)),1),H_E3_CR_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(3).XEndPoints(2), length(H_E3_HFD_M(:,i)),1),H_E3_HFD_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(3).XEndPoints(2), length(H_E3_HFD_F(:,i)),1),H_E3_HFD_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);

scatter(repmat(h(1).XEndPoints(3), length(H_E4_AL_M(:,i)),1),H_E4_AL_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(1).XEndPoints(3), length(H_E4_AL_F(:,i)),1),H_E4_AL_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(2).XEndPoints(3), length(H_E4_CR_M(:,i)),1),H_E4_CR_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(2).XEndPoints(3), length(H_E4_CR_F(:,i)),1),H_E4_CR_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(3).XEndPoints(3), length(H_E4_HFD_M(:,i)),1),H_E4_HFD_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(3).XEndPoints(3), length(H_E4_HFD_F(:,i)),1),H_E4_HFD_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);

xticklabels(["E2" "E3" "E4"]);
tit = sprintf('Shannon-%s',char(levels(i)));
ylabel (tit);
set(gca,'TickDir','out','box','off');
end

fig_fn = 'B1-Shannon-oberved alpha-6 levels.emf';
suptitle(fig_fn(1:end-4));
cd('D:\Microbiome\Ana-fig\'); 
saveas(gcf,fig_fn);

%%
figure ('position',[1 41 1536 748.8]);
for i = 1:6
subplot(2,3,i);
h=bar([1 2 3],[EH_E2_AL_ave(i),EH_E2_CR_ave(i),EH_E2_HFD_ave(i);EH_E3_AL_ave(i),EH_E3_CR_ave(i),EH_E3_HFD_ave(i);...
    EH_E4_AL_ave(i),EH_E4_CR_ave(i),EH_E4_HFD_ave(i);]);
hold on;
h(1).FaceColor='none';
h(2).FaceColor='g';
h(3).FaceColor='b';
errorbar(h(1).XEndPoints,[EH_E2_AL_ave(i),EH_E3_AL_ave(i),EH_E4_AL_ave(i)],...
    [EH_E2_AL_err(i),EH_E3_AL_err(i),EH_E4_AL_err(i)],'LineStyle','none','Color','k','LineWidth',1);
errorbar(h(2).XEndPoints,[EH_E2_CR_ave(i),EH_E3_CR_ave(i),EH_E4_CR_ave(i)],...
    [EH_E2_CR_err(i),EH_E3_CR_err(i),EH_E4_CR_err(i)],'LineStyle','none','Color','k','LineWidth',1);
errorbar(h(3).XEndPoints,[EH_E2_HFD_ave(i),EH_E3_HFD_ave(i),EH_E4_HFD_ave(i)],...
    [EH_E2_HFD_err(i),EH_E3_HFD_err(i),EH_E4_HFD_err(i)],'LineStyle','none','Color','k','LineWidth',1);

% jitter. Starting in R2020b you can do this easily with scatter:
scatter(repmat(h(1).XEndPoints(1), length(EH_E2_AL_M(:,i)),1),EH_E2_AL_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(1).XEndPoints(1), length(EH_E2_AL_F(:,i)),1),EH_E2_AL_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(2).XEndPoints(1), length(EH_E2_CR_M(:,i)),1),EH_E2_CR_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(2).XEndPoints(1), length(EH_E2_CR_F(:,i)),1),EH_E2_CR_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(3).XEndPoints(1), length(EH_E2_HFD_M(:,i)),1),EH_E2_HFD_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(3).XEndPoints(1), length(EH_E2_HFD_F(:,i)),1),EH_E2_HFD_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);

scatter(repmat(h(1).XEndPoints(2), length(EH_E3_AL_M(:,i)),1),EH_E3_AL_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(1).XEndPoints(2), length(EH_E3_AL_F(:,i)),1),EH_E3_AL_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(2).XEndPoints(2), length(EH_E3_CR_M(:,i)),1),EH_E3_CR_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(2).XEndPoints(2), length(EH_E3_CR_F(:,i)),1),EH_E3_CR_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(3).XEndPoints(2), length(EH_E3_HFD_M(:,i)),1),EH_E3_HFD_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(3).XEndPoints(2), length(EH_E3_HFD_F(:,i)),1),EH_E3_HFD_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);

scatter(repmat(h(1).XEndPoints(3), length(EH_E4_AL_M(:,i)),1),EH_E4_AL_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(1).XEndPoints(3), length(EH_E4_AL_F(:,i)),1),EH_E4_AL_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(2).XEndPoints(3), length(EH_E4_CR_M(:,i)),1),EH_E4_CR_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(2).XEndPoints(3), length(EH_E4_CR_F(:,i)),1),EH_E4_CR_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(3).XEndPoints(3), length(EH_E4_HFD_M(:,i)),1),EH_E4_HFD_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(3).XEndPoints(3), length(EH_E4_HFD_F(:,i)),1),EH_E4_HFD_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);

xticklabels(["E2" "E3" "E4"]);
tit = sprintf('evenness-%s',char(levels(i)));
ylabel (tit);
set(gca,'TickDir','out','box','off');
end

fig_fn = 'B1-evenness-6 levels.emf';
suptitle(fig_fn(1:end-4));
cd('D:\Microbiome\Ana-fig\'); 
saveas(gcf,fig_fn);

%%
figure ('position',[1 41 1536 748.8]);
for i = 1:6
subplot(2,3,i);
h=bar([1 2 3],[D_E2_AL_ave(i),D_E2_CR_ave(i),D_E2_HFD_ave(i);D_E3_AL_ave(i),D_E3_CR_ave(i),D_E3_HFD_ave(i);...
    D_E4_AL_ave(i),D_E4_CR_ave(i),D_E4_HFD_ave(i);]);
hold on;
h(1).FaceColor='none';
h(2).FaceColor='g';
h(3).FaceColor='b';
errorbar(h(1).XEndPoints,[D_E2_AL_ave(i),D_E3_AL_ave(i),D_E4_AL_ave(i)],...
    [D_E2_AL_err(i),D_E3_AL_err(i),D_E4_AL_err(i)],'LineStyle','none','Color','k','LineWidth',1);
errorbar(h(2).XEndPoints,[D_E2_CR_ave(i),D_E3_CR_ave(i),D_E4_CR_ave(i)],...
    [D_E2_CR_err(i),D_E3_CR_err(i),D_E4_CR_err(i)],'LineStyle','none','Color','k','LineWidth',1);
errorbar(h(3).XEndPoints,[D_E2_HFD_ave(i),D_E3_HFD_ave(i),D_E4_HFD_ave(i)],...
    [D_E2_HFD_err(i),D_E3_HFD_err(i),D_E4_HFD_err(i)],'LineStyle','none','Color','k','LineWidth',1);

% jitter. Starting in R2020b you can do this easily with scatter:
scatter(repmat(h(1).XEndPoints(1), length(D_E2_AL_M(:,i)),1),D_E2_AL_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(1).XEndPoints(1), length(D_E2_AL_F(:,i)),1),D_E2_AL_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(2).XEndPoints(1), length(D_E2_CR_M(:,i)),1),D_E2_CR_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(2).XEndPoints(1), length(D_E2_CR_F(:,i)),1),D_E2_CR_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(3).XEndPoints(1), length(D_E2_HFD_M(:,i)),1),D_E2_HFD_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(3).XEndPoints(1), length(D_E2_HFD_F(:,i)),1),D_E2_HFD_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);

scatter(repmat(h(1).XEndPoints(2), length(D_E3_AL_M(:,i)),1),D_E3_AL_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(1).XEndPoints(2), length(D_E3_AL_F(:,i)),1),D_E3_AL_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(2).XEndPoints(2), length(D_E3_CR_M(:,i)),1),D_E3_CR_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(2).XEndPoints(2), length(D_E3_CR_F(:,i)),1),D_E3_CR_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(3).XEndPoints(2), length(D_E3_HFD_M(:,i)),1),D_E3_HFD_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(3).XEndPoints(2), length(D_E3_HFD_F(:,i)),1),D_E3_HFD_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);

scatter(repmat(h(1).XEndPoints(3), length(D_E4_AL_M(:,i)),1),D_E4_AL_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(1).XEndPoints(3), length(D_E4_AL_F(:,i)),1),D_E4_AL_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(2).XEndPoints(3), length(D_E4_CR_M(:,i)),1),D_E4_CR_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(2).XEndPoints(3), length(D_E4_CR_F(:,i)),1),D_E4_CR_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(3).XEndPoints(3), length(D_E4_HFD_M(:,i)),1),D_E4_HFD_M(:,i),15,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);
scatter(repmat(h(3).XEndPoints(3), length(D_E4_HFD_F(:,i)),1),D_E4_HFD_F(:,i),15,'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'LineWidth',1,'XJitter','randn','XJitterWidth',.1);

xticklabels(["E2" "E3" "E4"]);
tit = sprintf('InvSimpson-%s',char(levels(i)));
ylabel (tit);
set(gca,'TickDir','out','box','off');
end

fig_fn = 'B1-InvSimpson-6 levels.emf';
suptitle(fig_fn(1:end-4));
cd('D:\Microbiome\Ana-fig\'); 
saveas(gcf,fig_fn);

%% alternative one
for i = 1:6
ave = [D_E2_AL_ave(i),D_E2_CR_ave(i),D_E2_HFD_ave(i);D_E3_AL_ave(i),D_E3_CR_ave(i),D_E3_HFD_ave(i);...
    D_E4_AL_ave(i),D_E4_CR_ave(i),D_E4_HFD_ave(i);];
err = [D_E2_AL_err(i),D_E2_CR_err(i),D_E2_HFD_err(i);D_E3_AL_err(i),D_E3_CR_err(i),D_E3_HFD_err(i);...
    D_E4_AL_err(i),D_E4_CR_err(i),D_E4_HFD_err(i)];
A_num = [length(D_E2_AL_M(:,i)),length(D_E2_CR_M(:,i)),length(D_E2_HFD_M(:,i));...
    length(D_E3_AL_M(:,i)),length(D_E3_CR_M(:,i)),length(D_E3_HFD_M(:,i));...
    length(D_E4_AL_M(:,i)),length(D_E4_CR_M(:,i)),length(D_E4_HFD_M(:,i))];
B_num = [length(D_E2_AL_F(:,i)),length(D_E2_CR_F(:,i)),length(D_E2_HFD_F(:,i));...
    length(D_E3_AL_F(:,i)),length(D_E3_CR_F(:,i)),length(D_E3_HFD_F(:,i));...
    length(D_E4_AL_F(:,i)),length(D_E4_CR_F(:,i)),length(D_E4_HFD_F(:,i))];
raw_data.rd11 = D_E2_AL(:,i);raw_data.rd12 = D_E2_CR(:,i);raw_data.rd13 = D_E2_HFD(:,i);
raw_data.rd21 = D_E3_AL(:,i);raw_data.rd22 = D_E3_CR(:,i);raw_data.rd23 = D_E3_HFD(:,i);
raw_data.rd31 = D_E4_AL(:,i);raw_data.rd32 = D_E4_CR(:,i);raw_data.rd33 = D_E4_HFD(:,i);

xticklelabel_list = ["E2" "E3" "E4"];
y_label = sprintf('InvSimpson-%s',char(levels(i)));
fig_fn = strcat(y_label,'.emf'); 
save_fig_path ='D:\'; 
output = bar_witherrorbar_datapoint(ave,err,A_num,B_num,raw_data,xticklelabel_list,y_label,fig_fn,save_fig_path);
end
















