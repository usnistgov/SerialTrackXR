function [x, BWseg, beadPara] = funLocateParticlesAC_2d(Img, beadPara, img_num, brightfield_image)
% [x, BWseg, beadParameter] = funLocateParticlesAC_2d(I, beadParameter,img_num) locates particles in the image
%
% INPUTS
% -------------------------------------------------------------------------
%   I:                 Input volumetric image
%   beadParameter:     Parameters to detect particle position in images (see 'funSetUpBeadParams.m')
%   img_num:           Number of volume in the sequence
%   brightfield_image: background correction image
%
% OUTPUTS
% -------------------------------------------------------------------------
%   x:              pixel-level estimate of particle center in MxNxO format
%   BWseg:          Segmentation estimate of the image
%   beadParameter:  Parameters to detect particle position in images (see 'funSetUpBeadParams.m')
%
% ----------------------------------------------
% Author: Alex Landauer
% Last time updated: 11/2023
% ==============================================

% Parameters
thres = beadPara.thres;    %Threshold value [0,1]
minPixels = beadPara.minSize;    %Threshold value [0,1]
maxPixels = beadPara.maxSize;    %Threshold value [0,1]
smoothFac = beadPara.smoothFac; %active contour smoothing factor
circThresh = beadPara.circThresh; %circulatiry threshold to define a particle [0,1]

disp('%%%%%% Starting Binarization %%%%%%')

BWseg = funImBinarize(Img,beadPara, brightfield_image);

% Find bead blobs
CC = bwconncomp(BWseg);
numPixels = cellfun(@numel,CC.PixelIdxList);

beadBlob = numPixels>minPixels & numPixels<maxPixels;

% Find centroid of all the eligible blobs
% get region info
S_ = regionprops(CC, 'Centroid', 'Eccentricity', 'PixelList');

for ii = 1:length(S_)
    circularity(ii) = S_(ii).Eccentricity;
end

%get minimal region region props
s = regionprops(CC, 'Centroid');

%get centroid points of the eligible blobs
blobPts = round(double(struct2dataset(s)));
blobPts = blobPts(beadBlob'  & abs(circularity')>circThresh,:);
temp = blobPts;

% Convert to m,n,o coordinates
blobPts(:,1) = temp(:,2);
blobPts(:,2) = temp(:,1);
x = blobPts;

if size(blobPts,1) > 1
    disp(length(blobPts))
end

disp('%%%%%% Binarization complete! %%%%%%')

end

