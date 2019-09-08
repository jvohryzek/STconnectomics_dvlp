% This program computes Figure2 of 'Dynamic spatio-temporal patterns of brain connectivity reorganize across development' Vohryzek et al. in prepation
% Alessandra Griffa 2016, edited by Jakub Vohryzek 2019
load('/Users/jakub/Matlab/September_2018/SD_TD_v2/SD_TD_gpo_gpy_1000_Globalwise_trial4_corr_v3_PERMUTATION_entprob.mat')

%load('/Volumes/JakubExtHD/localadmin/matlab/September_2018/STconnectome_dvlp_matlab/workSpace1.mat')
opt.FontSize = 18;
fontSize = 22;
i = 1 % global measure
% % SD_gpAtp2 = SD_prob_gpAtp2;
% % SD_gpBtp2 = SD_prob_gpBtp2;
% % 
% % SD_gpA = SD_prob_gpA;
% % SD_gpB = SD_prob_gpB;

TD_alltp = (TD_gpAtp2-TD_gpBtp2);
SD_alltp = (SD_gpAtp2-SD_gpBtp2);

if (TD_gpA(i)-TD_gpB(i)) <= 0
    p_TD(1,i) = (sum((TD_gpA(i)-TD_gpB(i)) > (TD_gp1(4,i,:) - TD_gp2(4,i,:)))/nb_permutations); % one-tailed
elseif (TD_gpA(i)-TD_gpB(i)) > 0
    p_TD(1,i) = (sum((TD_gpA(i)-TD_gpB(i)) < (TD_gp1(4,i,:) - TD_gp2(4,i,:)))/nb_permutations); % one-tailed
end
if (SD_gpA(i)-SD_gpB(i)) <= 0
    p_SD(1,i) = (sum((SD_gpA(i)-SD_gpB(i)) > (SD_gp1(4,i,:) - SD_gp2(4,i,:)))/nb_permutations); % one-tailed
elseif (SD_gpA(i)-SD_gpB(i)) > 0
    p_SD(1,i) = (sum((SD_gpA(i)-SD_gpB(i)) < (SD_gp1(4,i,:) - SD_gp2(4,i,:)))/nb_permutations); % one-tailed
end
%% PLOTTING
H1 = figure;
h1 = histfit(squeeze(TD_gp1(4,1,:)) - squeeze(TD_gp2(4,1,:)),30); hold on; plot([TD_alltp TD_alltp],[0 100],'r--', 'LineWidth', 3);%hold on;plot([0.05 0.05],[0 800], 'LineWidth', 2);
h1(1).FaceColor = [.8 .8 .8];h1(1).EdgeColor = 'none';
h1(2).Color = [.1 .1 .1];
title(sprintf(' p=%.3f',p_TD(1)),'FontSize',fontSize);
xlabel('Diff. in STD btw the two groups','FontSize',fontSize)
set(H1,'Units','Inches','Position',[0, 0, 15, 15],'PaperUnits','Inches','PaperSize',[15, 15])
H1 = fancy_figure(H1, opt);
saveas(H1,strcat(SaveFile,'Figure2_STD'),'eps')
saveas(H1,strcat(SaveFile,'Figure2_STD'),'jpg')

H2 = figure
h2 = histfit(squeeze(SD_gp1(4,1,:)) - squeeze(SD_gp2(4,1,:)),30); hold on; plot([SD_alltp SD_alltp],[0 100],'r--', 'LineWidth', 3);%hold on; plot([0.05 0.05],[0 800], 'LineWidth', 2);
h2(1).FaceColor = [.8 .8 .8]; h2(1).EdgeColor = 'none';
h2(2).Color = [.1 .1 .1];
title(sprintf(' p=%.3f',p_SD(1)),'FontSize',fontSize);
xlabel('Diff. in SD btw the two groups','FontSize',fontSize)
set(H2,'Units','Inches','Position',[0, 0, 15, 15],'PaperUnits','Inches','PaperSize',[15, 15])
H2 = fancy_figure(H2, opt);

saveas(H2,strcat(SaveFile,'Figure2_SD'),'eps')
saveas(H2,strcat(SaveFile,'Figure2_SD'),'jpg')
