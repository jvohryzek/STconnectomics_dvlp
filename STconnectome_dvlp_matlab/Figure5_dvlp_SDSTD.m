% This program computes Figure5 of 'Dynamic spatio-temporal patterns of brain connectivity reorganize across development' Vohryzek et al. in prepation
% Alessandra Griffa 2016, edited by Jakub Vohryzek 2019

path_CCfilt = '/Volumes/JakubExtHD/localadmin/matlab/September_2018/STconnectome_dvlp_matlab/INPUT/connected_components_filtered_scrub/';

opt.FontSize = 28
subj_struct = dir(fullfile(path_CCfilt,'*3T*'));

count = 0;
for i = 1:1:length(subj_struct)
    count = count + 1;
    subj{count} = fullfile(path_CCfilt,subj_struct(i).name);
    code{count} = subj_struct(i).name(1:end-16);
end
numSubj = length(code);

%% Load data
load('/Volumes/JakubExtHD/localadmin/matlab/September_2018/STconnectome_dvlp_matlab/edges_measures/INPUT/shortFibers.mat')
load('/Volumes/JakubExtHD/localadmin/matlab/September_2018/STconnectome_dvlp_matlab/edges_measures/INPUT/longFibers.mat')
load('/Volumes/JakubExtHD/localadmin/matlab/September_2018/STconnectome_dvlp_matlab/edges_measures/INPUT/intersystemFibers.mat')
age = [];
for j = 1:1:numSubj % SUBJECT LOOP
        if exist(strcat(path_CCfilt, code{j},'_CCfiltscrub.mat'),'file') == 2
            load(strcat(path_CCfilt, code{j},'_CCfiltscrub.mat'));
            display(['SUBJ: ', code{j}])
            age = [age, str2double(CCfilt{1,1}.subject_age)];
        end
end

% Longitudinal data
for i = 1:length(code)
    Longitud_data(i,1) = str2num(code{i}(1:4));
end
[C,IA,IC] = unique(Longitud_data);

for i=1:length(C)
    
    tmp_fl = long_fibers(IC==i);
    tmp_sl = short_fibers(IC==i);
    tmp_isl = InterSystem_fibers(IC==i);
    tmp_age = age(IC==i);

    if length(age(IC==i)) >= 2
        Longitud_long_fibers(i) = mean(tmp_fl);
        Longitud_short_fibers(i) = mean(tmp_sl);
        Longitud_InterSystem_fibers(i) = mean(tmp_isl);
        Longitud_age(i) = mean(tmp_age); 
    else

        Longitud_long_fibers(i) = tmp_fl;
        Longitud_short_fibers(i) = tmp_sl;
        Longitud_InterSystem_fibers(i) = tmp_isl;
        Longitud_age(i) = tmp_age;

    end
    clear tmp_age tmp_fl tmp_sl tmp_isl tmp_FD tmp_DVARS

end


%% Plotting
SaveFile = '/Volumes/JakubExtHD/localadmin/matlab/September_2018/STconnectome_dvlp_matlab/Figures/';
fontSize = 18;
h1 = figure;

[R_mlsf, p_mlsf] = Plot_Fit_ConfidenceInterval(Longitud_age,Longitud_short_fibers,1,0.05)
display(strcat(sprintf('Longitudinal Corrected Short Fibers r = %.3f:',R_mlsf(1,2)),sprintf(' p = %.3f', p_mlsf(1,2))))

xlabel('Age (years)','FontSize',fontSize);ylabel('Short Fibers','FontSize',fontSize);
title(strcat(sprintf('r = %.3f',R_mlsf(1,2)),sprintf(' (p = %.3f)', p_mlsf(1,2))),'FontSize',fontSize)
ylim([0.25 0.5]);xlim([5 34]);
h1 = fancy_figure(h1, opt);
pos = get(h1,'Position');
set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

saveas(h1,strcat(SaveFile,'Figure5_ShortFibers'),'svg')
saveas(h1,strcat(SaveFile,'Figure5_ShortFibers'),'pdf')

h2 = figure;
display(strcat(sprintf('Longitudinal Corrected Long Fibers r = %.3f:',R_mlsf(1,2)),sprintf(' p = %.3f', p_mlsf(1,2))))

[R_mllf, p_mllf] = Plot_Fit_ConfidenceInterval(Longitud_age,Longitud_long_fibers,1,0.05)
xlabel('Age (years)','FontSize',fontSize);ylabel('Long Fibers','FontSize',fontSize);
title(strcat(sprintf('r = %.3f',R_mllf(1,2)),sprintf(' (p = %.3f)', p_mllf(1,2))),'FontSize',fontSize)
ylim([0.2 0.45]);xlim([5 34]);
h2 = fancy_figure(h2, opt);
% set(h2,'Units','Inches','Position',[0, 0, 30, 10],'PaperUnits','Inches','PaperSize',[30, 10])
%set(h2,'Units','Inches');
pos = get(h2,'Position');
set(h2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

saveas(h2,strcat(SaveFile,'Figure5_LongFibers'),'svg')
saveas(h2,strcat(SaveFile,'Figure5_LongFibers'),'pdf')

h3 = figure;
display(strcat(sprintf('Longitudinal Corrected Inter-System Fibers r = %.3f:',R_mlsf(1,2)),sprintf(' p = %.3f', p_mlsf(1,2))))

[R_mlis, p_mlis, H2] = Plot_Fit_ConfidenceInterval(Longitud_age,Longitud_InterSystem_fibers,1,0.05)
xlabel('Age (years)','FontSize',fontSize);ylabel('Inter-System Fibers','FontSize',fontSize);
title(strcat(sprintf('r = %.3f',R_mlis(1,2)),sprintf(' (p = %.3f)', p_mlis(1,2))),'FontSize',fontSize)
ylim([0.2 0.5]);xlim([5 34]);
h3 = fancy_figure(h3, opt);
pos = get(h3,'Position');
set(h3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

saveas(h3,strcat(SaveFile,'Figure5_IntersystemFibers'),'svg')
saveas(h3,strcat(SaveFile,'Figure5_IntersystemFibers'),'pdf')
