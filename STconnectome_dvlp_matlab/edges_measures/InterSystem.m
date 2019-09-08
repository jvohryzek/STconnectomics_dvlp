%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of InterSystem fibers
%
% Analysing CCfilt with their functional network affiliation
%__________________________________________________
% Author: Jakub Vohryzek
% Hagmann Group
% CHUV-UNIL
% September 2018
% Version $1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;

% load necessary data
path_CCfilt = '/Volumes/JakubExtHD/localadmin/matlab/September_2018/STconnectome_dvlp_matlab/INPUT/connected_components_filtered_scrub/';
%% Initialising

yeo_matrix_all =[];
subj_young = 0;subj_old = 0;
age = [];yeo_matrix_individual = [];
%% Computing the subject names

subj_struct = dir(fullfile(path_CCfilt,'*3T*'));
count = 0;
for i = 1:1:length(subj_struct)
    count = count + 1;
    subj{count} = fullfile(path_CCfilt,subj_struct(i).name);
    code{count} = subj_struct(i).name(1:end-16);
end
%% Cortical indices
load('/Volumes/JakubExtHD/localadmin/matlab/September_2018/STconnectome_dvlp_matlab/edges_measures/index_CORTICAL_Laus2008_scale4.mat')
ix = ix;

%% Computing Yeo matrix
edge_intra = [];
edge_inter = [];
for j = 1:1:length(code)
    if exist(strcat(path_CCfilt, code{j},'_CCfiltscrub.mat'),'file') == 2
        load(strcat(path_CCfilt, code{j},'_CCfiltscrub.mat'));
        clear matrix
        yeo_matrix = zeros(7,7,length(CCfilt));
        for i = 1:1:length(CCfilt)
            % each node associated with gien functional network
            yeo = CCfilt{i}.functional_network;

            % this for loop extracts an array for nodes with given func-network
            clear f_network_index
            f_network_index = str2num(yeo(:,end));

            % here we loop through the seven functional systems to get unique nodes
            for n = 1:1:7
                unique_node_yeo = CCfilt{i}.nodes(f_network_index == n);
                % we loop all the edges and associate the different f-network
                % to them if 2 edge belongs to one network if 1 to two different

                for t= 1:1:length(CCfilt{i}.edges)
                                        
                    intseectionEdge2Yeo = intersect(CCfilt{i}.edges(t,:),unique_node_yeo);
                    if length(intseectionEdge2Yeo) == 2
                        % adding condition to exclude connections between the same
                        % region in different tps
                       if CCfilt{i}.anatomical_id(CCfilt{i}.nodes == intseectionEdge2Yeo(:,1)) ~= CCfilt{i}.anatomical_id(CCfilt{i}.nodes == intseectionEdge2Yeo(:,2))
                          yeo_matrix(n,n,i) = yeo_matrix(n,n,i) + 1;
                       end
                    elseif length(intseectionEdge2Yeo) == 1
                        for m =  1:1:7
                            unique_node_yeo_2 = CCfilt{i}.nodes(f_network_index == m);
                            y = intersect(CCfilt{i}.edges(t,:),unique_node_yeo_2);
                            if length(y) == 1
                                if intseectionEdge2Yeo ~= y
                                   yeo_matrix(n,m,i) = yeo_matrix(n,m,i) + 1;
                                end
                           end
                        end
                    end
                end
            end
        end

        yeo_matrix_all = cat(3,yeo_matrix_all,yeo_matrix);
        yeo_matrix_individual = cat(3,yeo_matrix_individual,sum(yeo_matrix,3))
    end
end

%% Saving
for mn = 1:87
Inter_diag(mn) = (sum(sum(triu(yeo_matrix_individual(:,:,mn),1))))./(sum(sum(triu(yeo_matrix_individual(:,:,mn),0))))
end
InterSystem_fibers = Inter_diag';
save('/Volumes/JakubExtHD/localadmin/matlab/September_2018/STconnectome_dvlp_matlab/edges_measures/INPUT/intersystemFibers.mat','InterSystem_fibers')
