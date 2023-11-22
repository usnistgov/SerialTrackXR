% %%%%%%%%%%%%%%%%%% SerialTrack-XR (in situ X-ray projection image particle tracking) %%%%%%%%%%%%%%%%%
% Runscript for code "SerialTrack-XR", called by the "main" function
% ***********************************************
% Dimension:            2D
% Particle rigidity:    hard 
% Tracking mode:        incremental
% -----------------------------------------------
%
% -----------------------------------------------
% References
% [1] M Patel, SE Leggett, AK Landauer, IY Wong, C Franck. Rapid,
%     topology-based particle tracking for high-resolution measurements of
%     large complex 3D motion fields. Scientific Reports. 8:5581 (2018).
% [2] J Yang, L Hazlett, AK Landauer, C Franck. Augmented Lagrangian
%     Digital Volume Correlation (ALDVC). Experimental Mechanics (2020).
% [3] T Janke, R Schwarze, K Bauer. Part2Track: A MATLAB package for double
%     frame and time resolved Particle Tracking Velocimetry. 11, 100413, SoftwareX (2020).
% [4] J Heyman. TracTrac: a fast multi-object tracking algorithm for motion
%     estimation. Computers & Geosciences, vol 128, 11-18 (2019).
% [5] https://www.mathworks.com/matlabcentral/fileexchange/77347-gridded-interpolation-and-gradients-of-3d-scattered-data
% [6] https://www.mathworks.com/matlabcentral/fileexchange/61436-regularizend
% [7] Jin Yang, et al.. SerialTrack: ScalE and Rotation Invariant Augmented Lagrangian Particle Tracking. 
%     SoftwareX, Volume 19, 101204 (2022) https://www.sciencedirect.com/science/article/pii/S2352711022001224
% [8] Buyukozturk, S., Landauer, A., Summey, L. et al. High-Speed, 3D Volumetric Displacement and Strain Mapping 
%     in Soft Materials Using Light Field Microscopy. Exp Mech (2022). https://doi.org/10.1007/s11340-022-00885-z
% -----------------------------------------------
% Author: Alex Landauer (adapted from code by Jin Yang)
% Contact and support: alex.landauer@gmail.com or landauer@nist.gov
% Date: 2023.11
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
%%%%% Load 2D images %%%%%
[file_names,Img,MPTPara_temp] = funReadImage2; close all;
disp('%%%%%% Load images: Done! %%%%%%'); fprintf('\n');

%%%%% Update MPTPara %%%%%
MPTPara.gridxyROIRange = MPTPara_temp.gridxyROIRange;
MPTPara.LoadImgMethod = MPTPara_temp.LoadImgMethod;
MPTPara.ImgSize = MPTPara_temp.ImgSize;

%%%%% Load image mask file(s) %%%%%
try
    if MaskFileLoadingMode == 1
        %%%%% Load one image mask file for all frames %%%%%
        [file_mask_name,ImgMaskFiles] = funReadImageMask2; close all;
    elseif MaskFileLoadingMode == 2
        %%%%% Load multiple image mask files %%%%%
        [file_mask_name,ImgMaskFiles] = funReadImageMask2; close all;
    elseif MaskFileLoadingMode == 3
        try load(im_roi_mask_file_path); catch; end
        try MPTPara.ImgRefMask = im_roi'; % Load stored image roi if existed
        catch, MPTPara.ImgRefMask = ones(size(Img{1})); % Set up default image mask file
        end
    else
        MPTPara.ImgRefMask = ones(size(Img{1})); % Set up default image mask file
    end
    disp('%%%%%% Load image mask file: Done! %%%%%%'); fprintf('\n');
catch
    MPTPara.ImgRefMask = ones(size(Img{1})); % Set up default image mask file
end
  
  
   
%% ====== Detect and localize particles ======
%%%%% Particle detection parameters %%%%%
%%%%% Bead Parameters %%%%%
% beadPara.thres = 0.35 ;         % Threshold for detecting particles
% beadPara.beadSize = 0;          % Estimated radius of a single particle [px]
% beadPara.minSize = 2;           % Minimum radius of a single particle [px]
% beadPara.maxSize = 40;          % Maximum area of a single particle [px^2]
% beadPara.circThresh = 0.2;       %circularity thresh for defining aparticle
% beadPara.smoothFac = 0.3;       %smoothing factor for active contours
% beadPara.winSize = [5, 5];      % Default (not used)
% beadPara.dccd = [1,1];          % Default (not used)
% beadPara.abc = [1,1];           % Default (not used)
% beadPara.forloop = 1;           % Default (not used)
% beadPara.randNoise = 1e-7;      % Default (not used)
% beadPara.PSF = [];              % PSF function; Example: PSF = fspecial('disk', beadPara.beadSize-1 ); % Disk blur
% beadPara.distMissing = 2;       % Distance threshold to check whether particle has a match or not [px]
% beadPara.color = 'black';       % Foreground (particle) color: options, 'white' or 'black'


gridxy = MPTPara.gridxyROIRange;

fprintf('Brightfield correction?:  \n')
fprintf('     0: No;  \n')
fprintf('     1: Yes;  \n')
prompt = 'Input here: ';
brightfield_yn = input(prompt);

if brightfield_yn == 1
    disp('Select a brightfield correction image: ')
    [file,path] = uigetfile('*.tiff');
    if isequal(file,0)
       disp('User selected Cancel');
    else
       brightfield_correction_image_path = fullfile(path,file);
    end
    
    brightfield_image_ = imnlmfilt(double(imread(brightfield_correction_image_path)),...
        'ComparisonWindowSize',3,'SearchWindowSize',13);
    brightfield_image = brightfield_image_(gridxy.gridx(1):gridxy.gridx(2), gridxy.gridy(1):gridxy.gridy(2));
else
    brightfield_image = ones(size(gridxy.gridx(1):gridxy.gridx(2), gridxy.gridy(1):gridxy.gridy(2)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ImgSeqNum = 1; % First reference image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Apply the uploaded mask file %%%%%
try
    if (MaskFileLoadingMode==1) || (MaskFileLoadingMode==3)
        currImg = Img{ImgSeqNum}.*MPTPara.ImgRefMask;
    elseif MaskFileLoadingMode==2
        currImg = Img{ImgSeqNum}.*ImgMaskFiles{ImgSeqNum};
    else
        currImg = Img{ImgSeqNum};  
    end
catch
    currImg = Img{ImgSeqNum};  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

currImg = currImg(MPTPara.gridxyROIRange.gridx(1):MPTPara.gridxyROIRange.gridx(2), ...
                  MPTPara.gridxyROIRange.gridy(1):MPTPara.gridxyROIRange.gridy(2));

%%%%% If PSF is non-empty, SerialTrack performs deconvolution %%%%%
if ~isempty(beadPara.PSF)
    currImg = deconvlucy(currImg,beadPara.PSF,6);
    disp(['----- Deconvolution frame #',num2str(ImgSeqNum),' ------']);
end
%%%% figure, imshow(currImg);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Detect and localize particles %%%%%

currImg2_norm = currImg;

%pre-process to get initial threshold and size values
beadPara = funGetBeadPara(beadPara, MPTPara, currImg2_norm, brightfield_image);

%run the particle detection and localization
[x_, bw_img, beadPara] = funLocateParticlesAC_2d(currImg2_norm, beadPara, ImgSeqNum, brightfield_image);
x{1}{ImgSeqNum} = x_;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Store particle positions as "parCoordA" %%%%%
x{1}{ImgSeqNum} = x{1}{ImgSeqNum} + [MPTPara.gridxyROIRange.gridx(1)-1, MPTPara.gridxyROIRange.gridy(1)-1];
parCoordA = x{1}{ImgSeqNum};

%%%%% Remove bad parCoord outside the image area %%%%%
for tempi=1:2, parCoordA( parCoordA(:,tempi) > size(Img{ImgSeqNum},tempi), : ) = []; end
for tempi=1:2, parCoordA( parCoordA(:,tempi) < 1, : ) = []; end
 
%%%%% Plot %%%%%
figure, imshow(imread(file_names{1,1}))
hold on; plot( parCoordA(:,1), parCoordA(:,2), 'r.');
view(2); box on; axis equal; axis tight; set(gca,'fontsize',18);
title('Detected particles in ref image','fontweight','normal');
 
%%%%% Report detected beads # %%%%%
disp(['Detected particle # in ref image: ',num2str(size(parCoordA,1))]);
disp('%%%%%% Detect particles: Done! %%%%%%'); fprintf('\n');
 


%% %%%%% Initialization %%%%%
%%%%%  MPT Parameter %%%%%
%these are reasonable starting point "defaults":
% MPTPara.f_o_s = 30;              % Size of search field: max(|u|,|v|) [px]
% MPTPara.n_neighborsMax = 25;     % Max # of neighboring particles
% MPTPara.n_neighborsMin = 1;      % Min # of neighboring particles
% MPTPara.locSolver = 1;           % Local solver: 1-topology-based feature; 2-histogram-based feature first and then topology-based feature;
% MPTPara.gbSolver = 3;            % Global step solver: 1-moving least square fitting; 2-global regularization; 3-ADMM iterations
% MPTPara.smoothness = 1e-2;       % Coefficient of regularization
% MPTPara.outlrThres = 2;          % Threshold for removing outliers in TPT
% MPTPara.maxIterNum = 20;         % Max ADMM iteration number
% MPTPara.iterStopThres = 1e-2;    % ADMM iteration stopping threshold
% MPTPara.strain_n_neighbors = 20; % # of neighboring particles used in strain gauge
% MPTPara.strain_f_o_s = 60;       % Size of virtual strain gauge [px]
% MPTPara.usePrevResults = 0;      % Whether use previous results or not: 0-no; 1-yes;

%%%%%% To store results %%%%%
parCoord_prev = cell(length(Img)-1,1);      
parCoord_prev{1} = parCoordA;
track_A2B_prev = cell(length(Img)-1,1);     
track_B2A_prev = cell(length(Img)-1,1);
uv_B2A_prev = cell(length(Img)-1,1);

 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
track_ratio = nan*zeros(length(Img)-1,1);
DefType = 'exp'; defList = [2:1:length(Img)]';
fig_track_rat = figure; ax_TR = axes; hold on; plot(defList,track_ratio,'r^-.','linewidth',1);
title('Tracking ratio');
xlabel('Frame #'); ylabel('Incremental tracking ratio');
try axis([2,size(file_names,2),0,1]); catch, end
drawnow

resultDisp = cell(size(file_names,2)-1,1);
resultDefGrad = cell(size(file_names,2)-1,1);

plot_vectors = 0;
for ImgSeqNum = 2:length(Img)  % "ImgSeqNum" is the frame index
    
    disp(['====== Frame #',num2str(ImgSeqNum),' ======']);
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Apply the uploaded mask file %%%%%
    clear defImg;
    try 
        if (MaskFileLoadingMode==1) || (MaskFileLoadingMode==3)
            defImg = Img{ImgSeqNum}.*MPTPara.ImgRefMask;
        elseif MaskFileLoadingMode==2
            defImg = Img{1,ImgSeqNum}.*ImgMaskFiles{1,ImgSeqNum};
        else
            defImg = Img{ImgSeqNum};
        end
    catch
       defImg = Img{ImgSeqNum};
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%% SerialTrack particle tracking %%%%%
    [parCoordB_temp,uv_B2A_temp,resultDisp{ImgSeqNum-1},resultDefGrad{ImgSeqNum-1},track_A2B_temp,track_B2A_temp,matchRatio] = fun_SerialTrack_2D_HardPar( ...
        ImgSeqNum,defImg,brightfield_image,beadPara,MPTPara,parCoord_prev{ImgSeqNum-1},parCoord_prev(2:end),uv_B2A_prev,plot_vectors);
     
    %%%%% Store results %%%%%
    parCoord_prev{ImgSeqNum} = parCoordB_temp;
    uv_B2A_prev{ImgSeqNum-1} = uv_B2A_temp; % incremental displacement
    track_A2B_prev{ImgSeqNum-1} = track_A2B_temp;
    track_B2A_prev{ImgSeqNum-1} = track_B2A_temp;
     
    
    track_A2B = track_A2B_prev{ImgSeqNum-1}; 
    track_ratio(ImgSeqNum-1) = length(track_A2B(track_A2B>0))/length(track_A2B);
    
    plot(defList,track_ratio,'r^-.','linewidth',1);
    drawnow
    
end


%%%%% Save results %%%%%
disp('%%%%%% SerialTrack 2D hard particle tracking: Done! %%%%%%'); fprintf('\n');
results_file_name = [file_names{1,1}(1:end-4),'_tracking_results_STXR.mat'];
mkdir results
save(['./results/' results_file_name]);
 

 










