%%%%%%%%%%%%%%%%%% SerialTrack-X-ray (ST-XR, in situ X-ray projection image particle tracking) %%%%%%%%%%%%%%%%%
%
% SerialTrack main file: run this file to execute the entire workflow
%
% ===================================================
% Dimension:            2D
% Particle rigidity:    hard 
% Tracking mode:        incremental
%
% ===================================================
% Author: Alex Landauer (NIST MML MMSD)
% Org. code from Jin Yang: https://github.com/FranckLab/SerialTrack
% Email: alex.landauer@gmail.com or landauer@nist.gov 
% Date: 2023.11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
close all; clear all; clc; clearvars -global
disp('************************************************');
disp('*** Welcome to SerialTrack Particle Tracking ***');
disp('************************************************');
addpath( './function/','./src/'); 

 
%% user defined parameters

%%%%% Problem dimension and units %%%%%
MPTPara.DIM = 2; % problem dimensionality, this is always "2" for ST-RX
MPTPara.xstep = 3.5665;  % unit: um/px, from the x-ray machine's parameters
MPTPara.tstep = 1;  % time step 

%%%%% Code mode %%%%%
MPTPara.mode = 'inc'; % 'inc': incremental mode, ST-XR is only set up for inc mode


%%%%% Particle rigidity %%%%%
MPTPara.parType = 'hard'; % {'hard': hard particle, ST-XR is only set up for inc mode

disp('************************************************');
disp(['Dimention: ',num2str(MPTPara.DIM)]);
disp(['Tracking mode: ',MPTPara.mode]);
disp(['Particle type: ',MPTPara.parType]);
disp('************************************************'); fprintf('\n');

%%%%% Image binary mask file %%%%%
MaskFileLoadingMode = 0; % {0: No mask file
                         %  1: Load only one mask file for all frames; 
                         %  2: Load one mask file for each single frame;
                         %  3: Load a MATLAB mat file for all frames;

if MaskFileLoadingMode == 3
    im_roi_mask_file_path = '.\im_roi.mat';  % TODO: Path of the mat file to be used as the mask file, if used
else
    im_roi_mask_file_path = '';  % If there is no mask mat file, leave it as empty;
end


%%%%% Particle detection and localization parameters %%%%%
%%%%% Bead parameters %%%%%
beadPara.thres = 0.7;           % Threshold for detecting particles
beadPara.circThresh = 0.2;      %circulatiry required to count as a particle
beadPara.smoothFac = 0.3;       %smoothing factor for active contours
beadPara.beadSize = 9;          % Estimated radius of a single particle [px]
beadPara.minSize = 2;           % Minimum area of a single particle [px^2]
beadPara.maxSize = 100;          % Maximum area of a single particle [px^2]
beadPara.winSize = [25, 25];      % Default [not used in 2D]
beadPara.dccd = [1,1];          % Default [not used in 2D]
beadPara.abc = [1,1];           % Default [not used in 2D]
beadPara.forloop = 1;           % Default [not used in 2D]
beadPara.randNoise = 1e-7;      % Default [not used in 2D]
beadPara.PSF = [];              % PSF function; Example: PSF = fspecial('disk', BeadPara.beadSize-1 ); % Disk blur
beadPara.distMissing = 250;       % Distance threshold to check whether particle has a match or not [px]
beadPara.color = 'black';       % Bead color: 'white' -or- 'black'


%%%%%%%%%%%%%%%%%%%%% SerialTrack particle tracking %%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Multiple particle tracking (MPT) Parameters %%%%%
MPTPara.f_o_s = 120;              % Size of search field: max(|u|,|v|) [px]
MPTPara.n_neighborsMax = 24;     % Max # of neighboring particles
MPTPara.n_neighborsMin = 1;      % Min # of neighboring particles
MPTPara.locSolver = 1;           % [not used] Local solver: 1-topology-based feature; 2-histogram-based feature first and then topology-based feature;
MPTPara.gbSolver = 3;            % Global step solver: 1-moving least square fitting; 2-global regularization; 3-ADMM iterations
MPTPara.outlrThres = 6;          % Threshold for removing outliers in TPT
MPTPara.maxIterNum = 20;         % Max ADMM iteration number
MPTPara.iterStopThres = 1e-4;    % ADMM iteration stopping threshold
MPTPara.strain_n_neighbors = 16;  % # of neighboring particles used in strain gauge
MPTPara.strain_f_o_s = 120;      % Size of virtual strain gauge [px]
MPTPara.usePrevResults = 0;      % Whether use previous results or not: 0-no; 1-yes;  


%%%% Postprocessing: Eulerian frame total strain %%%%%
MPTPara.edge_width = 45;         % ignores particle close to the edges [px]
MPTPara.grid_spacing = [3,3];    % grid spacing of the Eulerian mesh [px]
MPTPara.smoothness = 1e-2;       % Coefficient of regularization

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Execute SerialTrack particle tracking %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(MPTPara.mode,'inc')==1
    run_Serial_MPT_2D_hardpar_inc;
else
    disp('Only inc mode is available for SerailTrack-XR')
end
 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Postprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 run_post_process_eulerian_2d

