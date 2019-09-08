% This program computes Figure4 part 2 of 'Dynamic spatio-temporal patterns of brain connectivity reorganize across development' Vohryzek et al. in prepation
% Alessandra Griffa 2016, edited by Jakub Vohryzek 2019

load('/Users/jakub/Matlab/September_2018/SD_TD_v2/SD_TD_gpo_gpy_1000_nodewise_LONGIT_corr_trial42019.mat')
SaveFile = '/Volumes/JakubExtHD/localadmin/matlab/September_2018/STconnectome_dvlp_matlab/Figures/';

scale = 4;
path_correspondance = sprintf('/Users/jakub/Matlab/September_2018/SD_TD_v2/correspondance_yeo_scale%d.txt', scale);
correspondance = load(path_correspondance, '-ascii'); % linking functional systems to Lausanne 2008
correspondance = correspondance(correspondance > 0); % excluding the subcortical regions indexed by 0 
nFS = max(correspondance);
labelFS = {'VIS','SM','DA','VA','LIM','FP','DM'};
markerFS = {'^','v','o','p','h','s','d'};
%% initialising
trial = 4;
nb_struct = 7;
numROIs = 448;
numPerm = 1000;
typeTest = 'two-sided';
p_SD = zeros(trial,nb_struct); p_TD = zeros(trial,nb_struct);
fontSize = 34;
hColor = [120 18 134; 70 30 180; 0 118 14; 196 58 250; 220 248 164; 230 148 34; 205 62 78]./256;

%% Correlation plots 23/02/2019
opt.FontSize = 28;

HcorrTD = figure,
for i=1:7
    plot(TD_gpA(correspondance == i),TD_gpB(correspondance == i),markerFS{i},'Color',hColor(i,:),'MarkerFaceColor',hColor(i,:),'MarkerSize',8)
    hold on
end
hold on
plot(0.68:0.01:0.88,0.68:0.01:0.88,'k--','LineWidth',4)

xlabel('STD Adults')
ylabel('STD Children')
title(strcat('r = ',num2str(round(corr(TD_gpA',TD_gpB','Type','Spearman'),2))))
legend(labelFS,'Location','southeast','FontSize',32)
HcorrTD = fancy_figure(HcorrTD, opt);
pos = get(HcorrTD,'Position');
set(HcorrTD,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

HcorrSDprob = figure, 
for i=1:7
    plot(SD_prob_gpA(correspondance == i),SD_prob_gpB(correspondance == i),markerFS{i},'Color',hColor(i,:),'MarkerFaceColor',hColor(i,:),'MarkerSize',8)
    hold on
end
hold on
plot(1.3:0.1:2.4,1.3:0.1:2.4,'k--','LineWidth',4)
xlabel('SD Adults')
ylabel('SD Children')
title(strcat('r = ',num2str(round(corr(SD_prob_gpA',SD_prob_gpB','Type','Spearman'),2))))
legend(labelFS,'Location','southeast','FontSize',32)
HcorrSDprob = fancy_figure(HcorrSDprob, opt);
pos = get(HcorrSDprob,'Position');
set(HcorrSDprob,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

saveas(HcorrSDprob,strcat(SaveFile,'Figure4_corrSD'),'svg')
saveas(HcorrSDprob,strcat(SaveFile,'Figure4_corrSD'),'jpg')
saveas(HcorrTD,strcat(SaveFile,'Figure4_corrSTD'),'svg')
saveas(HcorrTD,strcat(SaveFile,'Figure4_corrSTD'),'jpg')
