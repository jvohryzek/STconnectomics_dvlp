%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of Short and Long fibers
%
%__________________________________________________
% Author: Jakub Vohryzek
% Hagmann Group
% CHUV-UNIL
% September 2018
% Version $1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input

path_CCfilt = '/Volumes/JakubExtHD/localadmin/matlab/September_2018/STconnectome_dvlp_matlab/INPUT/connected_components_filtered_scrub/'
%% Initialising
num_edge_cohort_fl_short = [];
num_edge_cohort_fl_long = [];
num_edge_cohort = [];
ICV_orig = [];
    
%% Cortical indices
load('index_CORTICAL_Laus2008_scale4.mat')
ix = ix;

%% Loading SC
load('SCgroup50_cortical.mat')
SC = double(SC);

%% Loading length SC
load('lengthSC.mat')

fl_concate = L(ix,ix,:);
fl_concate(fl_concate == 0) = NaN;

fl_average = mean(fl_concate,3,'omitnan');
fl_average(isnan(fl_average)) = 0;
fl_average = fl_average.*SC; % (recurrence_th50_71);
%% ICV
load('ICV_aseg.mat')
templateICV = mean(str2num(ICV));
%% Extracting subject names

subj_struct = dir(fullfile(path_CCfilt,'*3T*'));
count = 0;
for i = 1:1:length(subj_struct)
    count = count + 1;
    subj{count} = fullfile(path_CCfilt,subj_struct(i).name);
    code{count} = subj_struct(i).name(1:end-16);
end
%% Analysis

num_subj = length(code);    

for j = 1:1:num_subj % Subject loop
        if exist(strcat(path_CCfilt, code{j},'_CCfiltscrub.mat'),'file') == 2
            load(strcat(path_CCfilt, code{j},'_CCfiltscrub.mat'));
            display(['SUBJ: ', code{j}])
            fl_threshold_short = ((double(CCfilt{1, 1}.ICV)/templateICV)^(1/3)).*20; %%% threshold(th)%%%;
            fl_threshold_long = ((double(CCfilt{1, 1}.ICV)/templateICV)^(1/3)).*42;  %%% threshold(th)%%%;

            edge_fl_short = []; edge_fl_long = [];
            edge_subj = [];
            for i = 1:1:length(CCfilt) % CC LOOP
                display(['Subj: ', code{j}, ' , CC: ', num2str(i), ', # of Edges: ',num2str(length(CCfilt{i}.edges)) ])

                %% Indexing edges
                edge_CC = [];
                for t = 1:1:length(CCfilt{i}.edges) % EDGE LOOP
                    index_nn = CCfilt{i}.edges(t,1);
                    index_mm = CCfilt{i}.edges(t,2);
                    edge_anat_ID_n = CCfilt{i}.anatomical_id((CCfilt{i}.nodes == index_nn));
                    edge_anat_ID_m = CCfilt{i}.anatomical_id((CCfilt{i}.nodes == index_mm));
                    edge_fl = fl_average(edge_anat_ID_n,edge_anat_ID_m);
                    edge_CC = [edge_CC,edge_fl]; % all edges for a CC
                    %% Separating long/short edges
                    if edge_fl < fl_threshold_short
                        edge_fl_short = [edge_fl_short,edge_fl];
                    elseif edge_fl > fl_threshold_long
                        edge_fl_long = [edge_fl_long,edge_fl];
                    end        
                    
                end
                %% Percentage of long/short fibers per subject
                edge_subj = [edge_subj, edge_CC] ; % all edges for a subj
            end
            % ICV extraction
            ICV_orig = [ICV_orig, CCfilt{1,1}.ICV];
            num_edge_cohort = [num_edge_cohort,size(nonzeros(edge_subj),1)];
            num_edge_cohort_fl_short = [num_edge_cohort_fl_short, size(nonzeros(edge_fl_short),1)];
            num_edge_cohort_fl_long  = [num_edge_cohort_fl_long, size(nonzeros(edge_fl_long),1)];

        end      
end

long_fibers_norm = num_edge_cohort_fl_long./num_edge_cohort;
short_fibers_norm = num_edge_cohort_fl_short./num_edge_cohort;
%% Saving

long_fibers = (num_edge_cohort_fl_long./num_edge_cohort)';
short_fibers = (num_edge_cohort_fl_short./num_edge_cohort)';
save('/Volumes/JakubExtHD/localadmin/matlab/September_2018/STconnectome_dvlp_matlab/edges_measures/INPUT/shortFibers.mat','short_fibers')
save('/Volumes/JakubExtHD/localadmin/matlab/September_2018/STconnectome_dvlp_matlab/edges_measures/INPUT/longFibers.mat','long_fibers')
