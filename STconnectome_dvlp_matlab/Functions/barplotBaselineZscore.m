%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HISTPLOTDIFFERENCE
%
% This functions plots baseline and zscore for each measure and each group
% summary_________________________________________________
% Authors: Jakub Vohryzek
% Hagmann Group
% CHUV-UNIL
% July 2018
% Version $1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H] = barplotBaselineZscore(Group,permGroup,errorGroup,saveFile,fileName,hColor,labelFS,nameCondition,nameGroup,savePlot)

trial = 1;
labelFSorder = cell(1,7);
barFontSize = 24;
H = figure
opt.FontSize = 18;
% baseline
[evalBaseline, ixBaseline] = sort(Group,'descend');
xBaseline = 1:7;

% zscore
for i=1:7
    tmp = zscore([squeeze(permGroup(trial,i,:));Group(i)]);
    zscoreGroup(i) = tmp(end);   
end
[evalZscore, ixZscore] = sort(zscoreGroup,'descend');
xZscore = 1:7;

for b = 1:7
    subplot(1,2,1)
	% Plot one single bar as a separate bar series.
	handlebar(b) = bar(xBaseline(b), evalBaseline(b), 'BarWidth', 0.9);
    hold on
    errorbar(xBaseline(b),evalBaseline(b),std(errorGroup{ixBaseline(b)}),'Color',hColor(ixBaseline(b),:),'LineWidth',3)

	% Apply the color to this bar series.
	set(handlebar(b), 'FaceColor', hColor(ixBaseline(b),:),'EdgeColor',[1 1 1]);
	% Place text atop the bar
	% barTopper = sprintf('y(%d) = %.3f', x(b), eval(b));
	% text(x(b)-0.2, eval(b)+3, char(labelFS{ix}), 'FontSize', barFontSize,'Color',hColor(ix(b),:));
	% set(gca, 'XTickLabels', labelFS{b});
    
    title(strcat('Baseline:', nameCondition,nameGroup))
    hold on;
    labelFSorder{b} = labelFS{ixBaseline(b)};
    
    switch nameCondition
        case ' Spatio-Temporal Diversity ' % for TD values between 0 and 1
    set(gca, 'ylim',[floor(min(10*evalBaseline))/10, inf],'XTick',1:7,'XTickLabels', labelFSorder);
        case ' System Diversity '
    set(gca, 'ylim',[floor(min(evalBaseline)), inf],'XTick',1:7,'XTickLabels', labelFSorder);
    end
    subplot(1,2,2)
	% Plot one single bar as a separate bar series.
	handlebarZscore(b) = bar(xZscore(b), evalZscore(b), 'BarWidth', 0.9);
	% Apply the color to this bar series.
	set(handlebarZscore(b), 'FaceColor', hColor(ixZscore(b),:),'EdgeColor',[1 1 1]);
	% Place text atop the bar
	%barTopper = sprintf('y(%d) = %.3f', x(b), eval(b));
	%text(xZscore(b)-0.2, evalZscore(b)+3, char(labelFS{ixZscore}), 'FontSize', barFontSize,'Color',hColor(ixZscore(b),:));
	%set(gca, 'XTickLabels', labelFS{b});
    hold on;
    labelFSorderZscore{b} = labelFS{ixZscore(b)}; % reordering labels
    title(strcat('Zscore Baseline:', nameCondition,nameGroup))
    
    if min(evalZscore) > 0
    set(gca, 'ylim',[floor(min(evalZscore)), inf],'XTick',1:7,'XTickLabels', labelFSorderZscore,'FontSize',32);
    else
    set(gca, 'ylim',[-inf ceil(max(evalZscore))],'XTick',1:7,'XTickLabels', labelFSorderZscore,'FontSize',32);
    end
        
    view([90 90])

end

%set(gca, 'XTickMode', 'Auto');
%view([90 90])
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
H = fancy_figure(H, opt);
if savePlot
    saveas(H,strcat(saveFile,fileName),'svg')
    saveas(H,strcat(saveFile,fileName),'jpg') 
end
end

