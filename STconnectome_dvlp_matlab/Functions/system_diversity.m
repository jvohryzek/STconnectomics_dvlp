%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYSTEM DIVERSITY
%
% This function caclulates System Diversity.
%
% Input Parameters:
%   list_CCs       : list of all the CC to compute the system diversity from
%   ts_CCs         : Index of the CCs affiliation to Nodes/FS/Global.
%   scale          : Parcelation
%__________________________________________________
% Authors: Jakub Vohryzek and Emeline Mullier
% Hagmann Group
% CHUV-UNIL
% July 2018
% Version $1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [SD_jaccard, SDent, SD_prob,SDentprob, Ent_jaccard, Ent_prob, jaccard,prob] = system_diversity(list_CCs, ts_CCs, scale)
    % yeo functional correspondance
    path_correspondance = sprintf('/Volumes/JakubExtHD/localadmin/matlab/September_2018/SD_TD_v2/correspondance_yeo_scale%d.txt', scale);
    correspondance = load(path_correspondance, '-ascii'); % linking functional systems to Lausanne 2008
    correspondance = correspondance(correspondance > 0); % excluding the subcortical regions indexed by 0 
    
    nCCs = size(ts_CCs, 1); % number of CCs
   
    nROIs = size(ts_CCs, 2); % number of ROIS
    
    SD_jaccard = zeros(1,nROIs);
    Ent_jaccard  = zeros(1,nCCs);
    jaccard = zeros(nCCs,7);
    
    SD_prob = zeros(1,nROIs);
    Ent_prob  = zeros(1,nCCs);
    prob = zeros(nCCs,7);

    for i = 1:nCCs
        % STEP 1:
        %%% JACCARD - taking into consideration the functional system size
        %%% A&B/A|B
        clear unique_node 
        unique_node  = unique(list_CCs(i).anatomical_id);
        list_nodes = zeros(size(correspondance));
        list_nodes(unique_node) = 1;
        
        for ind = 1:7
            jaccard(i, ind) = (sum((list_nodes) & (correspondance == ind)))./(sum((list_nodes) | (correspondance == ind)));
            %%% PROBABILITY - number of nodes in the functional system
            %%% divided by number of nodes in the CCs
            prob(i,ind) = length(find(correspondance(unique(list_CCs(i).anatomical_id)) == ind))./length(unique(list_CCs(i).anatomical_id));
        end
        
       % STEP 2: 
       e_ind = jaccard(i,:).*log2(jaccard(i,:));
       e_ind(isnan(e_ind)) = 0;
       Ent_jaccard(i) = -sum(e_ind);
       
       e_ind2 = prob(i,:).*log2(prob(i,:));
       e_ind2(isnan(e_ind2)) = 0;
       Ent_prob(i) = -sum(e_ind2);

    end  
    % STEP 3:
    SDent = {};
    SDentprob = {};
    for i=1:nROIs
        SD_prob(i) = mean(Ent_prob(find(ts_CCs(:,i))));
        SD_jaccard(i) = mean(Ent_jaccard(find(ts_CCs(:,i))));
        SDent{i} = Ent_jaccard(find(ts_CCs(:,i)));
        SDentprob{i} = Ent_prob(find(ts_CCs(:,i)));
    end

end