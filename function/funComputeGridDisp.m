function [gridPts,disp_grid,vecPts,u_vec] = funComputeGridDisp(x,u_tpt,track)
% This file is to run after TPT code to interpolate onto a grid.
% This will generate TPT post-processed displacement fields on a
% gridded single z-plane, all pushed back to the reference config grid
% ('engineering' displacements)
%
%--- INPUTS ---
%  x{time}{channel}: bead locations, we assume channel is {1}
%  u_tpt{time}{channel}: displacements
%  track{time}{channel}: tracking matches
%
%--- OUTPUTS ---
% gridPts : full set of grid points in matrix (ndgrid) form
% u_grid  : displacements at the gridPts points
% vecPts  : full set of points in col vec format (more concise)
% u_vec   : cell of vectors of displacements at the vecPts points
%
% Alex Landauer 2021-01-27

%% Parameters for regularization
winstepsize = [1,1]; % Mesh element size

%% ====== Generate uniform mesh ======
% Find domain borders
ymin = ceil(min(x{1}(:,1))); ymax = floor(max(x{1}(:,1)));
xmin = ceil(min(x{1}(:,2))); xmax = floor(max(x{1}(:,2)));

%set up grid
[Yq,Xq] = ndgrid(xmin:winstepsize(1):xmax,ymin:winstepsize(2):ymax);
gridPts{1} = Xq;
gridPts{2} = Yq;

brd_wd = 3; %pad size
[Yq_,Xq_] = ndgrid(xmin-brd_wd*winstepsize(1):winstepsize(1):xmax+brd_wd*winstepsize(1),...
    ymin-brd_wd*winstepsize(2):winstepsize(2):ymax+brd_wd*winstepsize(2));

% Combine coordinate points vectors
coordPts = [reshape(Xq_,prod(size(Xq_)),1),reshape(Yq_,prod(size(Yq_)),1)];

U = cell(length(u_tpt),1);
% ====== Interpolate nodal values ======
for t = 2:length(x)
    fprintf('Interpolating scattered points: %i of %i\n',t-1,length(u_tpt))
    u1temp = scatteredInterpolant(x{1}(track{t-1}>0,1),x{1}(track{t-1}>0,2),...
        u_tpt{t-1}(track{t-1}>0,1),'natural','none');
    v1temp = scatteredInterpolant(x{1}(track{t-1}>0,1),x{1}(track{t-1}>0,2),...
        u_tpt{t-1}(track{t-1}>0,2),'natural','none');
   
    %pad with zero and fill in any nans
    u_grid{t-1} = (padarray(u1temp(Xq,Yq),[brd_wd,brd_wd]));
    v_grid{t-1} = (padarray(v1temp(Xq,Yq),[brd_wd,brd_wd]));
    
    disp_grid{t-1}{1} = u_grid{t-1};
    disp_grid{t-1}{2} = v_grid{t-1};
    
    %Combine whole displacement vector
    U{t-1} = zeros(prod(size(u_grid{t-1}))*2,1);
    U{t-1}(1:2:end) = reshape(u_grid{t-1},prod(size(u_grid{t-1})),1);
    U{t-1}(2:2:end) = reshape(v_grid{t-1},prod(size(v_grid{t-1})),1);
    u_vec{t-1} = [reshape(u_grid{t-1},prod(size(u_grid{t-1})),1),...
        reshape(v_grid{t-1},prod(size(v_grid{t-1})),1)];
    
end

vecPts = coordPts;


% Plot solved denoised displacement field
% Plotdisp_show3(full(uhat),coordinatesFEM,elementsFEM);

