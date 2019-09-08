%% This program performs scrubbing by removing the CC containing frames 2 discard and recomputing the FM
% Alessandra Griffa 2016, edited by Jakub Vohryzek 2019
% CCfilt = filtered by min & max nb of components --> here min tp = 2,
% min # ROIs = 6 vs max tp = 196 (all) and max ROIs = 448 (all)
% CC{x}.width : #of tp (1tp = 2.4s)
% CC{x}.anatomical_id : ROIs (which ROIs are involved)
% CC{x}.nodes : unique ID for each noed over all tp (at the moment for all sbj)
% CC{x}.subject : sbj number
% CC{x}.PH : sbj ID
% CC{x}.tp: subject space (1-196)
% CC{x}.heigth: max # ROIs?
% FMsubj : 448 ROI x nb of CCfilt for this sbj (information compressed in
% time --> sum & normalization tp that each node is active)
%% paths
project_dir1 = '/Volumes/JakubExtHD/localadmin/PycharmProjects/Lausanne/STconnectome_dvlp_matlab/OUTPUT/cc_filtered/';
project_dir2 = '/Volumes/JakubExtHD/localadmin/PycharmProjects/Lausanne/STconnectome_dvlp_matlab/OUTPUT/feature_matrix_filtered/';
project_dir3 = '/Volumes/JakubExtHD/RAD/ELLIEZ_DEVELOPMENT/DATA/';

%% load sbj processing order

subj_struct = dir(fullfile(project_dir3,'*3T*'));
count = 0;
for i = 1:1:length(subj_struct)
    count = count + 1;
    subj{count} = fullfile(project_dir3,subj_struct(i).name);
    code{count} = subj_struct(i).name;
end

%% parameters
nROI = 448;
count =0;
%% load CC
for s = 1:length(code)
    if exist([ project_dir3 code{s} '/T1/CMP/fMRI/preprocessing_Jakub_v4/tp_after_scrubbing.mat'],'file') == 2 && exist(strcat(project_dir1,code{s},'_CCfilt.mat'),'file') == 2;
    % load frames to remove (scrubbing)
    load([ project_dir3 code{s} '/T1/CMP/fMRI/preprocessing_Jakub_v4/tp_after_scrubbing.mat']); % for the excluded sbj already
    % load CCfilt
    load(strcat(project_dir1,code{s},'_CCfilt.mat'));
    % load FMfilt:
    load(strcat(project_dir2,code{s},'_FMfilt.mat'));
    count = count +1
    display(code{s})
    %% HOMEMADE SCRUBBING: remove CC which contain one of the frames where an
    %% artifact is suspected and save it
    savedir1 =  [project_dir1 'connected_components_filtered_scrub/'];
    if ~exist(savedir1)
        mkdir(savedir1)
    end
    savedir2 =  [project_dir2 'feature_matrix_filtered_scrub/'];
    if ~exist(savedir2)
        mkdir(savedir2)
    end
    savenameCC = [savedir1 code{s} '_CCfiltscrub.mat'];
    savenameFM = [savedir2 code{s} '_FMscrub.mat'];
    
    FMsubj=FMfilt;
    f2d = nonzeros(index)';
    if ~isempty(f2d) % there are frames to discard for this sbj
        CCfilt_noscrub = CCfilt;
        for f = 1: length(f2d)
            for cc = 1: length(CCfilt_noscrub)
                if find(unique(CCfilt_noscrub{cc}.tp) == f2d(f)) % if tp with artifacts are included in this component
                    CCfilt{cc} = []; % remove this CC
                    fprintf('CC %i removed for sbj %i because it contains frame %i\n', cc, s, f2d(f));
                end
            end
        end
        finalCCnr =find(~cellfun(@isempty, CCfilt));
        CCfilt = CCfilt(finalCCnr);
        save(savenameCC, 'CCfilt', 'CCfilt_noscrub');
        
        % compute new FM
        FM = zeros(size(CCfilt,1),nROI);
        for cc = 1: length(CCfilt)
            active_nodes = unique(CCfilt{cc}.anatomical_id);
            for a = 1: length(active_nodes)
                FM(cc,active_nodes(a)) = length(find(CCfilt{cc}.anatomical_id == active_nodes(a)));
            end
            % normalize by norm along each row
            FMnorm(cc,:) = FM(cc,:)./norm(FM(cc,:));
        end
        FMsubj=[];
        FMsubj = FMnorm;
        clear FMnorm
        save(savenameFM, 'FMsubj')
    else
        save(savenameCC, 'CCfilt');
        save(savenameFM, 'FMsubj');
    end
    
    end
    
end