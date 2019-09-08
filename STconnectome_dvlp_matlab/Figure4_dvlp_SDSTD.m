% This program computes Figure4 part 1 of 'Dynamic spatio-temporal patterns of brain connectivity reorganize across development' Vohryzek et al. in prepation
% Alessandra Griffa 2016, edited by Jakub Vohryzek 2019
% SD can be run as SD (Jaccard ME) or SD_prob (probability Alessandra)
%% INPUT
% load('/Users/jakub/Matlab/September_2018/SD_TD_v2/SD_TD_gpo_gpy_10000_nodewise_LONGIT_corr_trial4v2.mat')
load('/Users/jakub/Matlab/September_2018/SD_TD_v2/SD_TD_gpo_gpy_1000_nodewise_LONGIT_corr_trial42019.mat')
numPerm = 1000;
numROIs = 448;
typeTest = 'one-sided';

%% Calculating the zscore colormap cut-off
% adults
% SD =SD_prob_gpAtp2;
% STD =TD_gpAtp2;
% children
SD =SD_prob_gpBtp2;
STD =TD_gpBtp2; 

Surfr = Read_Surface('/Volumes/JakubExtHD/localadmin/matlab/Collaboration_Elliez/plot_kit/needed_data/fsaverage/surf/rh.inflated');
Surfl = Read_Surface('/Volumes/JakubExtHD/localadmin/matlab/Collaboration_Elliez/plot_kit/needed_data/fsaverage/surf/lh.inflated');

[txtr,ctabr,colorsr] = read_cfiles('/Volumes/JakubExtHD/localadmin/matlab/Collaboration_Elliez/plot_kit/needed_data/fsaverage/label/rh.lausanne2008.scale4.annot');
[txtl,ctabl,colorsl] = read_cfiles('/Volumes/JakubExtHD/localadmin/matlab/Collaboration_Elliez/plot_kit/needed_data/fsaverage/label/lh.lausanne2008.scale4.annot');
[txty_rh,ctaby_rh,colorsy_rh] = read_cfiles('/Volumes/JakubExtHD/localadmin/matlab/Collaboration_Elliez/plot_kit/needed_data/fsaverage/label/rh.Yeo2011_7Networks_N1000.annot');
[txty_lh,ctaby_lh,colorsy_lh] = read_cfiles('/Volumes/JakubExtHD/localadmin/matlab/Collaboration_Elliez/plot_kit/needed_data/fsaverage/label/lh.Yeo2011_7Networks_N1000.annot');

%% YEO
Surfr.SurfData.FaceVertexCData = colorsy_rh;
Surfl.SurfData.FaceVertexCData = colorsy_lh;


Surfl.SurfData.vertices(:,1,1) = Surfl.SurfData.vertices(:,1,1)-65;
Surf = Compound_Surf([Surfr; Surfl]);
Surf_yeo{1,1} = Surf;
% Surface_Viewer(Surf_yeo);
%% TD
Surfr.SurfData.FaceVertexCData = colorsr;
Surfl.SurfData.FaceVertexCData = colorsl;

%%% Contours indexes
[Trip] = Vert_Neibp(double(Surfl.SurfData.faces),size(Surfl.SurfData.vertices,1),size(Surfl.SurfData.faces,1));
Temp = sum(Trip); Trip(:,Temp==0) = []; temp = Trip(:,3:end); indz = find(temp == 0); temp(indz) = 1; labid = txty_lh;
temp1 = labid(temp); temp1(indz) =  max(temp1(:))+1; NewMat = temp1-repmat(min(temp1')',[1 size(temp1,2)]); NewMat(indz) = 0;
a = logical(sum(logical(NewMat)')'); ind_lh = find(a);
[Trip] = Vert_Neibp(double(Surfr.SurfData.faces),size(Surfr.SurfData.vertices,1),size(Surfr.SurfData.faces,1));
Temp = sum(Trip); Trip(:,Temp==0) = []; temp = Trip(:,3:end); indz = find(temp == 0); temp(indz) = 1; labid = txty_rh;
temp1 = labid(temp); temp1(indz) =  max(temp1(:))+1; NewMat = temp1-repmat(min(temp1')',[1 size(temp1,2)]); NewMat(indz) = 0;
a = logical(sum(logical(NewMat)')'); ind_rh = find(a); 

ind_ctxr = [];
Is = zeros(length(txtr),1);
for i = 1:size(ctabr.table,1)-1
    ind = find(txtr == ctabr.table(i+1,5));
    Is(ind) = STD(i); 
    ind_ctxr = [ind_ctxr; ind];   
end
ind_subctxr = find(~ismember([1:length(txtr)], ind_ctxr));
%% Artificially adding maxlim and minlim to the colormap to unify the baseline
limHighTD = max([TD_gpAtp2,TD_gpBtp2]);
limLowTD = min([TD_gpAtp2,TD_gpBtp2]); 
%%

% creates rescaled colors => necessary later for creating the colorbar
Surfr.Is = Is;%[limLowTD;Is;limHighTD];
Colors = Surf_ColorAdj(Surfr,'winter');

Surfr.SurfData.FaceVertexCData = Colors; 
% changes some of the vertices of the contours to white and subcortical to
% black
Surfr.SurfData.FaceVertexCData(ind_rh,:) = ones(size(Surfr.SurfData.FaceVertexCData(ind_rh,:))); 
Surfr.SurfData.FaceVertexCData(ind_subctxr,:) = zeros(size(Surfr.SurfData.FaceVertexCData(ind_subctxr,:))); 

% mapping measure values to the appropriate vertices
ind_ctxl = [];
Is = zeros(length(txtl),1);
for i = 1:size(ctabl.table,1)-1
    ind = find(txtl == ctabl.table(i+1,5));
    Is(ind) = STD(i+size(ctabr.table,1)-1);
    ind_ctxl = [ind_ctxl; ind];   
end
ind_subctxl = find(~ismember([1:length(txtl)], ind_ctxl));
% creates rescaled colors => necessary later for creating the colorbar
Surfl.Is =Is;% [limLowTD;Is;limHighTD]; 
Colors = Surf_ColorAdj(Surfl,'winter');

Surfl.SurfData.FaceVertexCData = Colors; 
% changes some of the vertices of the contours to white and subcortical to
% black
Surfl.SurfData.FaceVertexCData(ind_lh,:) = ones(size(Surfl.SurfData.FaceVertexCData(ind_lh,:))); 
Surfl.SurfData.FaceVertexCData(ind_subctxl,:) = zeros(size(Surfl.SurfData.FaceVertexCData(ind_subctxl,:))); 

Surf = Compound_Surf([Surfr; Surfl]);

Surf_total{1,1} = Surfl; 
Surf_total{2,1} = Surfl; 
Surf_total{3,1} = Surf; 
Surf_total{4,1} = Surfr; 
Surf_total{5,1} = Surfr; 

%% SD

[Trip] = Vert_Neibp(double(Surfl.SurfData.faces),size(Surfl.SurfData.vertices,1),size(Surfl.SurfData.faces,1));
Temp = sum(Trip); Trip(:,Temp==0) = []; temp = Trip(:,3:end); indz = find(temp == 0); temp(indz) = 1; labid = txty_lh;
temp1 = labid(temp); temp1(indz) =  max(temp1(:))+1; NewMat = temp1-repmat(min(temp1')',[1 size(temp1,2)]); NewMat(indz) = 0;
a = logical(sum(logical(NewMat)')'); ind_lh = find(a);
[Trip] = Vert_Neibp(double(Surfr.SurfData.faces),size(Surfr.SurfData.vertices,1),size(Surfr.SurfData.faces,1));
Temp = sum(Trip); Trip(:,Temp==0) = []; temp = Trip(:,3:end); indz = find(temp == 0); temp(indz) = 1; labid = txty_rh;
temp1 = labid(temp); temp1(indz) =  max(temp1(:))+1; NewMat = temp1-repmat(min(temp1')',[1 size(temp1,2)]); NewMat(indz) = 0;
a = logical(sum(logical(NewMat)')'); ind_rh = find(a); 

ind_ctxr = [];
Is = zeros(length(txtr),1);
for i = 1:size(ctabr.table,1)-1
    ind = find(txtr == ctabr.table(i+1,5));
    Is(ind) = SD(i);  
    ind_ctxr = [ind_ctxr; ind];   
end
ind_subctxr = find(~ismember([1:length(txtr)], ind_ctxr));
%% Artificially adding maxlim and minlim to the colormap to unify the baseline
limHighSD = max([SD_gpAtp2,SD_gpBtp2]);
limLowSD = min([SD_gpAtp2,SD_gpBtp2]); 
%%
% creates rescaled colors => necessary later for creating the colorbar
Surfr.Is = Is; % [limLowSD;Is;limHighSD];
Colors = Surf_ColorAdj(Surfr,'winter');

Surfr.SurfData.FaceVertexCData = Colors; 
% changes some of the vertices of the contours to white and subcortical to
% black
Surfr.SurfData.FaceVertexCData(ind_rh,:) = ones(size(Surfr.SurfData.FaceVertexCData(ind_rh,:))); 
Surfr.SurfData.FaceVertexCData(ind_subctxr,:) = zeros(size(Surfr.SurfData.FaceVertexCData(ind_subctxr,:))); 

ind_ctxl = [];
Is = zeros(length(txtl),1);
for i = 1:size(ctabl.table,1)-1
    ind = find(txtl == ctabl.table(i+1,5));
    Is(ind) = SD(i+size(ctabr.table,1)-1);  
    ind_ctxl = [ind_ctxl; ind];   
end
ind_subctxl = find(~ismember([1:length(txtl)], ind_ctxl));
% creates rescaled colors => necessary later for creating the colorbar
Surfl.Is = Is;%[limLowSD;Is;limHighSD]; 
Colors = Surf_ColorAdj(Surfl,'winter');
Surfl.SurfData.FaceVertexCData = Colors; 
% changes some of the vertices of the contours to white and subcortical to
% black
Surfl.SurfData.FaceVertexCData(ind_lh,:) = ones(size(Surfl.SurfData.FaceVertexCData(ind_lh,:))); 
Surfl.SurfData.FaceVertexCData(ind_subctxl,:) = zeros(size(Surfl.SurfData.FaceVertexCData(ind_subctxl,:))); 

Surf = Compound_Surf([Surfr; Surfl]);

Surf_total{6,1} = Surfl; 
Surf_total{7,1} = Surfl; 
Surf_total{8,1} = Surf; 
Surf_total{9,1} = Surfr; 
Surf_total{10,1} = Surfr; 
[H] = Surface_Viewer_Dvlp(Surf_total,'colMap','winter');

% Surface_Viewer_Dvlp(Surf_total);
% saveFile = '/Users/jakub/Matlab/September_2018/Figuresv3/BrainPlots/'
% titleName = 'STD/SDprob ADults';
% [H] = test_Surface_Setter(Surf_total,TDdiffzscoreExtreme,SDdiffzscoreExtreme)
% suptitle(strcat('{\color{magenta} ',titleName,'}'))
%saveas(H,strcat(saveFile,'BrainPlots_TDSDprobAdults'),'jpg')
%saveas(H,strcat(saveFile,'BrainPlots_TDSDprobAdults'),'svg')
