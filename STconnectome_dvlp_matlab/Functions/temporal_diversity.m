%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPATIO-TEMPORAL DIVERSITY
%
% This function caclulates Spatio-Temporal Diversity.
%
% Input Parameters:
%   list_FMs       : list of all the FM to compute the system diversity from
%   ts_CCs         : Index of the CCs affiliation to Nodes/FS/Global.
%   nROIs          : number of regions
%__________________________________________________
% Authors: Jakub Vohryzek and Emeline Mullier
% Hagmann Group
% CHUV-UNIL
% July 2018
% Version $1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [TD,TDpairwise] = temporal_diversity(listFMs, tsCCs)

    
    [nCCs, nROIs] = size(tsCCs);
    TD = zeros(1,nROIs);
    TDpairwise = {};
    for i=1:nROIs
        TD(i) = mean(pdist(listFMs(find(tsCCs(:,i)==1),:),'cosine'));
        TDpairwise{i} = pdist(listFMs(find(tsCCs(:,i)==1),:),'cosine');
    end
end
%% Calculating FM from CC and normalising by l2-norm
%     for i=1:nCCs
%         ID = list_CCs(i).anatomical_id;
%         for j=1:length(ID)
%             FM(i,ID(j)) = FM(i,ID(j)) + 1;
%             %FM(i,ID(j)) =  1;
%         end     
%     end
%     
%     FM = FM./norm(FM,2);