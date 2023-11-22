function beadPara = funGetBeadPara(beadPara,TPTpara,Img, brightfield_image)
% MATLAB script: funGetBeadPara.m
% ----------------------------------------------
%   This script is used to collect Bead Parameter settings from the user
%   Careful selection of parameters is important to correctly segment
%   particles from the background, by specifying the threshold and size
%   parameters
%
%   INPUT: BeadPara      Struct of bead parameters to use for segmentation
%          Img           An example input image (usually the first)
%          TPTpara       Struct of particle tracking parameters
%
%   OUTPUT: BeadPara    Updated bead parameters struct
%
% ----------------------------------------------
% Author: Alex Landauer
% Last time updated: 11/2023
% ==============================================

% Parameters

% first, get the threshold initial estimate
[BWseg,Img_] = funImBinarize(Img,beadPara,brightfield_image);

f1 = figure;
histogram(Img_(Img_ ~= 0))
drawnow
f2 = figure;
imshow(Img_,[]),colorbar,axis image
drawnow

beadPara.thres = input('Enter binarization threshold estimate: ');
try
    close(f1)
    close(f2)
catch
end
YN = 0;
while YN == 0
    f3 = figure;
    
    % Image thresholding
    [BWseg,Img_] = funImBinarize(Img,beadPara,brightfield_image);
    
    imshow(BWseg,[]),colorbar,axis image
    drawnow
    
    YN_ = input('Binarization okay? (Y/N) [N]: ','s');
    
    if strcmpi(YN_,'Y')
        YN = 1;
    else
        f1 = figure;
        histogram(Img_(Img_ ~= 0))
        drawnow
        f2 = figure;
        imshow(Img_,[]),colorbar,axis image
        drawnow
        beadPara.thres = input('Enter binarization threshold estimate from histogram: ');
    end
    
end
try
    close(f1)
catch
end
try
    close(f2)
catch
end
try
    close(f3)
catch
end


% second, using this threshold get the min and max sizes
% Parameters
% thres = beadPara.thres;    %Threshold value

% other params, set by user set later in this method
% minPixels = beadParameter.minSize;  %Minimum pixel count in blob for bead
% maxPixels = beadParameter.maxSize;  %Maximum pixel count in blob for bead

disp('%%%%%% Starting Binarization %%%%%%')

% Image thresholding
BWseg = funImBinarize(Img,beadPara,brightfield_image);

% Find bead blobs
CC = bwconncomp(BWseg);
numPixels = cellfun(@numel,CC.PixelIdxList);

%check the histogram of sizes and get user defined size limits
try
    nbins = sshist(numPixels);
catch
    nbins = 1;
end
nbins = max([15,nbins]);

h = figure;
histogram(numPixels,nbins)
minPixels = input('Enter min bead size: ');  %Minimum pixel count in blob for bead
maxPixels = input('Enter max bead size: ');  %Maximum pixel count in blob for bead
try
    close(h)
catch
end

beadPara.minSize = minPixels;  %Minimum pixel count in blob for bead
beadPara.maxSize = maxPixels;  %Maximum pixel count in blob for bead


