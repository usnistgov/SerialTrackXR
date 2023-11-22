function [img_binary,Img_] = funImBinarize(Img,beadPara,brightfield_image)

%% preprocessing

[Img_] = imnlmfilt(Img,'ComparisonWindowSize',3,'SearchWindowSize',13);
% Img_(1,1) = (2^16)-1;
% Img_(1,2) = 1;
% Img_ = imflatfield(Img_,5);
Img_ = adapthisteq(rescale(Img_./brightfield_image),'ClipLimit',0.007);

% figure,imshow(Img_,[])

%% % Image thresholding
if strcmp(beadPara.color,'black')
    Img_ = (imcomplement(Img_));
    img_binary = imbinarize(Img_,beadPara.thres); %threshold-based to start
else
    Img_ = (Img_);
    img_binary = imbinarize(Img_,beadPara.thres); %threshold-based to start
end


% figure,imshow(img_binary,[])
% figure,imshow(Img_,[])

% % img_binary = edge(Img,'canny',[beadPara.thres/3,beadPara.thres]);
% %
% % % se_seg = strel('disk',1);
% % % img_binary = imdilate(img_binary, se_seg);
% %
% % img_binary = imfill(img_binary, 'holes');

% % % %dilate the result with a relatively large structuring element to
% % % %privide a good initial guess to active contours
% % % % se_bin = strel('disk',1);
% % % L = imdilate(img_binary, se_bin); %if using watershed, replace "vol_binary" with "L"
% % %
% % %
% % % %use active contour segmentation to refine binarization
% % % BWseg = activecontour(I,L,50,'Chan-Vese','SmoothFactor',smoothFac);
% % % % se_seg = strel('disk',1);
% % % % BWseg = imdilate(BWseg, se_seg); %dilate slightly to improve connectivity

% %distance transform the raw binarized image
% Di = bwdist(~img_binary);
% Di = -Di;
% % % Di(BW) = Inf;
% % % figure,imshow(Di,[])
% %now do a watershed to refine the binarized regions
% L = watershed(Di);
% L(~img_binary) = 0;
% L(L>1) = 1;

% % smoothFac = 0.01;
% % BWseg = activecontour(I,L,50,'Chan-Vese','ContractionBias',0.2,'SmoothFactor',smoothFac);
% % figure,imshow(BWseg,[])




