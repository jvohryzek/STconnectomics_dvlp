%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ts_CCs_per_ROI
%
% This function caclulates ts_CCs_per_ROI.
%
% Input Parameters:
%   list_CCs       : list of all the CC to compute the system diversity from
%   nROIs          : number of regions
%__________________________________________________
% Authors: Jakub Vohryzek and Emeline Mullier
% Hagmann Group
% CHUV-UNIL
% July 2018
% Version $1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mat_CCs_per_ROIs = ts_CCs_per_ROI(list_CCs, nROIs)

    nCCs = length(list_CCs);
    mat_CCs_per_ROIs = zeros(nCCs, nROIs);

    for i=1:nCCs
        mat_CCs_per_ROIs(i,unique(list_CCs(i).anatomical_id)) = 1;        
    end
        
end
