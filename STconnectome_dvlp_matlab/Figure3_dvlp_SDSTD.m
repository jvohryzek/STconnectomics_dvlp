%%%%%%%%
% This program computes Figure3 of 'Dynamic spatio-temporal patterns of brain connectivity reorganize across development' Vohryzek et al. in prepation
% Alessandra Griffa 2016, edited by Jakub Vohryzek 2019
% saving to Figuresv2.1
%%%%%%%%
clear all;
load('/Volumes/JakubExtHD/localadmin/matlab/September_2018/STconnectome_dvlp_matlab/workSpace7.mat')
SaveFile = '/Volumes/JakubExtHD/localadmin/matlab/September_2018/STconnectome_dvlp_matlab/Figures/';


%% INITIALISE

nb_struct = 7;
nb_permutations = 1000;
typeTest = 'one-sided';
p_SD = zeros(1,nb_struct); p_TD = zeros(1,nb_struct);
hColor = [120 18 134; 70 30 180; 0 118 14; 196 58 250; 220 248 164; 230 148 34; 205 62 78]./256;
labelFS = {'VIS','SM','DA','VA','LIM','FP','DM'};

entProb = 1;

if entProb
    
SD_gpA=SD_prob_gpA;
SD_gpB=SD_prob_gpB;

SD_gpAtp2 = SD_prob_gpAtp2;
SD_gpBtp2 = SD_prob_gpBtp2;

SDent_gpAtp2 = SDentprob_gpAtp2;
SDent_gpBtp2 = SDentprob_gpBtp2;
end
%% PLOTTING
opt.FontSize = 28;
nameCondition = ' SD '; nameGroup = ' Diff (A-CH)';
[H1] = histplotdifference(SD_gpAtp2,SD_gpBtp2,SD_gp1,SD_gp2,nb_permutations,typeTest,hColor,labelFS,nameCondition,nameGroup);
H1 = fancy_figure(H1, opt);
pos = get(H1,'Position');
set(H1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

nameCondition = ' STD '; nameGroup = ' Diff (A-CH)';
[H2] = histplotdifference(TD_gpAtp2,TD_gpBtp2,TD_gp1,TD_gp2,nb_permutations,typeTest,hColor,labelFS,nameCondition,nameGroup);
H2 = fancy_figure(H2, opt);
pos = get(H2,'Position');
set(H2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

saveas(H1,strcat(SaveFile,'Figure3_SD'),'svg')
saveas(H1,strcat(SaveFile,'Figure3_SD'),'eps')
saveas(H2,strcat(SaveFile,'Figure3_STD'),'svg')
saveas(H2,strcat(SaveFile,'Figure3_STD'),'eps')
