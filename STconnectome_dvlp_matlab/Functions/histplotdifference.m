%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HISTPLOTDIFFERENCE
%
% This functions qunatifies the zscore differences from the permuted
% distributions for each FS and plots a bar graph with the zscore difference
% summary_________________________________________________
% Authors: Jakub Vohryzek
% Hagmann Group
% CHUV-UNIL
% July 2018
% Version $1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H] = histplotdifference(groupA,groupB,permGroupA,permGroupB,numPerm,typeTest,hColor,labelFS,nameCondition,nameGroup)

H = figure,
opt.FontSize = 18;
nb_struct = length(groupA);
fontSize = 22;
labelFS = {'VIS','SM','DA','VA','LIM','FP','DM'};
for i=1:nb_struct
    switch typeTest
        
        case 'one-sided'
            if (groupA(i)-groupB(i)) <= 0
                groupP(1,i) = (sum((groupA(i)-groupB(i)) > (permGroupA(1,i,:) - permGroupB(1,i,:)))/numPerm); % one-tailed
            elseif (groupA(i)-groupB(i)) > 0
                groupP(1,i) = (sum((groupA(i)-groupB(i)) < (permGroupA(1,i,:) - permGroupB(1,i,:)))/numPerm); % one-tailed
            end
            
        case 'two-sided'
            groupP(1,i) = (sum(abs(groupA(i)-groupB(i)) > abs(permGroupA(1,i,:) - permGroupB(1,i,:)))/numPerm); % two-tailed   
    end
    ax = subplot(2,4,i); h2 = histfit(squeeze(permGroupA(1,i,:) - permGroupB(1,i,:)),25);
    
    hold on;
    plot([groupA(i)-groupB(i),groupA(i)-groupB(i)],[0 max(h2(1).YData)],'r--', 'LineWidth', 2);%hold on; plot([0.05 0.05],[0 800], 'LineWidth', 2);
    ax.LineWidth = 3;
    h2(1).EdgeColor = 'none';
    h2(1).FaceColor = hColor(i,:);
    alpha(.8)
    h2(2).Color = [.1 .1 .1];
    title(sprintf(' %s: p=%.3f',labelFS{i},groupP(1,i)),'FontSize',fontSize);
    set(ax, 'FontSize', 24)
    tmp = zscore([squeeze(permGroupA(1,i,:) - permGroupB(1,i,:));(groupA(i)-groupB(i))]);
    diffZscore(i) = tmp(end);
end
subplot(2,4,8)
labelFSorder = cell(1,7);
[eval, ix] = sort(diffZscore,'descend');
x = 1:7;
for b = 1 : 7
	% Plot one single bar as a separate bar series.
	handlebar(b) = bar(x(b), eval(b), 'BarWidth', 0.9);
	% Apply the color to this bar series.
	set(handlebar(b), 'FaceColor', hColor(ix(b),:),'EdgeColor','none');
	set(gca, 'XTickLabels', labelFS{b});
    hold on;
    labelFSorder{b} = labelFS{ix(b)};
end
set(gca, 'XTickMode', 'Auto');
set(gca, 'XTick',[1:7],'XTickLabels', labelFSorder,'FontSize',32);
view([90 90])

set(H,'Units','Inches','Position',[0, 0, 30, 10],'PaperUnits','Inches','PaperSize',[30, 10])
%title(strcat('Difference:', nameCondition,nameGroup))
Htitle=get(gca,'Title');
Htitle.String = strcat(nameCondition,nameGroup);

H = fancy_figure(H, opt);

end

