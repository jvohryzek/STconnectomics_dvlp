%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of SD and STD for Functional Systems
%
%__________________________________________________
% Author: Jakub Vohryzek
% Hagmann Group
% CHUV-UNIL
% September 2018
% Version $1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Input
% load necessary data
path_FMfilt = '/Volumes/JakubExtHD/localadmin/matlab/September_2018/STconnectome_dvlp_matlab/INPUT/feature_matrix_filtered_scrub/';
path_CCfilt = '/Volumes/JakubExtHD/localadmin/matlab/September_2018/STconnectome_dvlp_matlab/INPUT/connected_components_filtered_scrub/';

% Loading correspondance
scale = 4;
path_correspondance = sprintf('correspondance_yeo_scale%d.txt', scale);
correspondance = load(path_correspondance, '-ascii'); % linking functional systems to Lausanne 2008
correspondance = correspondance(correspondance > 0); % excluding the subcortical regions indexed by 0 

% loading the longitudinal index
load('longitudinalDataIndex.mat')

%% Initialisation
listSubs =[]; listSubs_gpo =[]; listSubs_gpy =[];
% CCs for SD
listCCs =[]; listCCs_gpo =[]; listCCs_gpy =[];
% FMs for TD
listFMs =[]; listFMs_gpo =[]; listFMs_gpy =[];
age_upper_threshold = 20; age_lower_threshold = 14;
age = [];nROIs = 448;
% Longitudinal
listCCs_gpytp2  = []; listFMs_gpytp2 = []; listSubs_gpytp2 = []; 
listCCs_gpotp2  = []; listFMs_gpotp2 = [];listSubs_gpotp2 = [];
age_otp2 = []; age_ytp2 = [];

%% Analysis

subj_struct = dir(fullfile(path_CCfilt,'*3T*'));
count = 0;
for ii = 1:length(subj_struct)
    count = count + 1;
    subj{count} = fullfile(path_CCfilt,subj_struct(ii).name);
    code{count} = subj_struct(ii).name(1:end-16);
    load(strcat(path_CCfilt, code{ii},'_CCfiltscrub.mat'));
    load(strcat(path_FMfilt, code{ii},'_FMscrub.mat'));
    age = [age, str2double(CCfilt{1,1}.subject_age)];

    display(strcat('Subj','_',code{ii}))
    
    lastTP_sbj = [];
    lastTP_sbj = lastTP(ii)*ones(1,length(CCfilt)); % longit index per CC
    for j=1:length(CCfilt)
        
        listCCs  = [listCCs CCfilt{j}];
        listFMs  = [listFMs; FMsubj(j,:)];
        listSubs = [listSubs ii]; 
        
        if str2double(CCfilt{1,1}.subject_age) > age_upper_threshold
            listCCs_gpo = [listCCs_gpo CCfilt{j}];
            listFMs_gpo = [listFMs_gpo; FMsubj(j,:)];
            listSubs_gpo = [listSubs_gpo; ii]; 
            
            if lastTP_sbj(j) == 1
                listCCs_gpotp2  = [listCCs_gpotp2 CCfilt{j}];
                listFMs_gpotp2 = [listFMs_gpotp2; FMsubj(j,:)];
                listSubs_gpotp2 = [listSubs_gpotp2; ii]; 
                age_otp2 = [age_otp2 str2double(CCfilt{1,1}.subject_age)];
            end
        elseif str2double(CCfilt{1,1}.subject_age) < age_lower_threshold
            listCCs_gpy  = [listCCs_gpy CCfilt{j}];
            listFMs_gpy = [listFMs_gpy; FMsubj(j,:)];
            listSubs_gpy = [listSubs_gpy; ii]; 
            if lastTP_sbj(j) == 1
                listCCs_gpytp2  = [listCCs_gpytp2 CCfilt{j}];
                listFMs_gpytp2 = [listFMs_gpytp2; FMsubj(j,:)];
                listSubs_gpytp2 = [listSubs_gpytp2; ii]; 
                age_ytp2 = [age_ytp2 str2double(CCfilt{1,1}.subject_age)];
            end
        end
    end  
end
%%
%%%%%%%%%%%%%%% TRANSLATING LISTS TO A COMMON CODE %%%%%%%%%%%%%%%%%%%

listCCs_gpA = listCCs_gpo;listFMs_gpA = listFMs_gpo; 
listCCs_gpB = listCCs_gpy;listFMs_gpB = listFMs_gpy;

% longitudinal 
listCCs_gpAtp2 = listCCs_gpotp2; listFMs_gpAtp2 = listFMs_gpotp2;
listCCs_gpBtp2 = listCCs_gpytp2; listFMs_gpBtp2 = listFMs_gpytp2;

                    %%%%%%%%%%%%%%% COMMON CODE %%%%%%%%%%%%%%%%%%%
%% Computation System and Temporal Diversity
ts_CC_gpA = ts_CCs_per_ROI(listCCs_gpA, nROIs);
ts_CC_gpB = ts_CCs_per_ROI(listCCs_gpB, nROIs);

% longitudinal
ts_CC_gpAtp2 = ts_CCs_per_ROI(listCCs_gpAtp2, nROIs);
ts_CC_gpBtp2 = ts_CCs_per_ROI(listCCs_gpBtp2, nROIs);

nFS = 7;
nb_struct = nFS; %% nROIs, nFS, 1

ts_gpA_struc = zeros(size(ts_CC_gpA,1),nb_struct); ts_gpB_struc = zeros(size(ts_CC_gpB,1),nb_struct);
% longitudinal
ts_gpAtp2_struc = zeros(size(ts_CC_gpAtp2,1),nb_struct); ts_gpBtp2_struc = zeros(size(ts_CC_gpBtp2,1),nb_struct);

%%%%% FOR AT LEAST repCC% OF FS IN THE CCS %%%%%%
percFS = 0.2;
for i=1:nb_struct
    repCC = percFS*sum(correspondance == i); % xx percent of a given FS => accounting for FS size
    ts_gpA_struc(:,i)  = double(sum(ts_CC_gpA(:, correspondance == i),2)>repCC);
    ts_gpB_struc(:,i)  = double(sum(ts_CC_gpB(:, correspondance == i),2)>repCC);

    % longitudinal
    ts_gpAtp2_struc(:,i)  = double(sum(ts_CC_gpAtp2(:, correspondance == i),2)>repCC);
    ts_gpBtp2_struc(:,i)  = double(sum(ts_CC_gpBtp2(:, correspondance == i),2)>repCC);
end
 
[SD_gpA,SDent_gpA,SD_prob_gpA,SDentprob_gpA] = system_diversity(listCCs_gpA, ts_gpA_struc, scale);
[SD_gpB,SDent_gpB,SD_prob_gpB,SDentprob_gpB] = system_diversity(listCCs_gpB, ts_gpB_struc, scale);

[TD_gpA,TDpairwise_gpA] = temporal_diversity(listFMs_gpA, ts_gpA_struc);
[TD_gpB,TDpairwise_gpB] = temporal_diversity(listFMs_gpB, ts_gpB_struc);

% longitudinal

[SD_gpAtp2,SDent_gpAtp2,SD_prob_gpAtp2,SDentprob_gpAtp2] = system_diversity(listCCs_gpAtp2, ts_gpAtp2_struc, scale);
[SD_gpBtp2,SDent_gpBtp2,SD_prob_gpBtp2,SDentprob_gpBtp2] = system_diversity(listCCs_gpBtp2, ts_gpBtp2_struc, scale);
[TD_gpAtp2,TDpairwise_gpAtp2] = temporal_diversity(listFMs_gpAtp2, ts_gpAtp2_struc);
[TD_gpBtp2,TDpairwise_gpBtp2] = temporal_diversity(listFMs_gpBtp2, ts_gpBtp2_struc);

%% Permutation of the CCs
nb_permutations = 1000;
SD_gp1 = zeros(1,nb_struct,nb_permutations); SD_gp2 = zeros(1,nb_struct,nb_permutations);
TD_gp1 = zeros(1,nb_struct,nb_permutations); TD_gp2 = zeros(1,nb_struct,nb_permutations);
p_SD = zeros(1,nb_struct); p_TD = zeros(1,nb_struct);

% longitudinal
CCs_all_sub = [];FMs_all_sub=[];list_node_all=[];
    
CCs_all_sub = [listCCs_gpAtp2 listCCs_gpBtp2]; % concatenated CCs for old and young
FMs_all_sub = [listFMs_gpAtp2; listFMs_gpBtp2]; % concatenated FMs for old and young
list_node_all = [ts_CC_gpAtp2; ts_CC_gpBtp2]; % index of all the regions participation in CCs or FMs

ts_struct = zeros(size(list_node_all,1),nb_struct);

for i=1:nb_struct
    repCCperm = percFS*sum(correspondance == i); % xx percent of a given FS => accounting for FS size
    ts_struct(:,i)  = double(sum(list_node_all(:, correspondance==i),2)>repCCperm);
end

for i=1:nb_permutations
    [x_CCs, y_struct] = ind2sub(size(ts_struct),find(ts_struct));
    i
    parfor j=1:max(y_struct) % number of unique global
        tmp1 = randperm(length(find(y_struct == j)));
        temp2 = x_CCs(find(y_struct == j));
        perm_grp = temp2(tmp1);
        perm_grp1 = perm_grp(1:round(length(perm_grp)/2));
        perm_grp2 = perm_grp((round(length(perm_grp)/2)+1:end));
        CCs_gp1  = CCs_all_sub(perm_grp1);
        CCs_gp2  = CCs_all_sub(perm_grp2);
        FMs_gp1  = FMs_all_sub(perm_grp1,:);
        FMs_gp2  = FMs_all_sub(perm_grp2,:);
        init1 = ones(length(perm_grp1),1);
        init2 = ones(length(perm_grp2),1);
        [SD_gp1(1,j,i),~,SD_prob_gp1(1,j,i)] = system_diversity(CCs_gp1, init1, scale);
        [SD_gp2(1,j,i),~,SD_prob_gp2(1,j,i)] = system_diversity(CCs_gp2, init2, scale);
        TD_gp1(1,j,i) = temporal_diversity(FMs_gp1, init1);
        TD_gp2(1,j,i) = temporal_diversity(FMs_gp2, init2);
    end
end

save(strcat('workSpace',num2str(nFS)))