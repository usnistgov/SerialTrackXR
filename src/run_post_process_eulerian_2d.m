
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SerialTrack-XR post-processing
% ===================================================
% Author: Alex Landauer
% Org. code HR-VPTM: https://github.com/FranckLab/HR-VPTM
% Email: alex.landauer@gmail.com or landauer@nist.gov 
% Date: 2023.10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%To begin, load a results .mat file from "run_serial_MPT_2D_hardpar_inc" or
%have a run in the workspace from that function 

% set up vars
smoothness = MPTPara.smoothness;
grid_spacing = MPTPara.grid_spacing;

% Find domain borders

%automatically find data range
% MPTPara.xRange(1) = ceil(min(resultDisp{1}.parCoordB(:,1))); MPTPara.xRange(2) = floor(max(resultDisp{1}.parCoordB(:,1)));
% MPTPara.yRange(1) = ceil(min(resultDisp{1}.parCoordB(:,2))); MPTPara.yRange(2) = floor(max(resultDisp{1}.parCoordB(:,2)));

%here, we prefer to use the whole image
MPTPara.xRange(1) = 1;
MPTPara.xRange(2) = size(Img{1},1);
MPTPara.yRange(1) = 1;
MPTPara.yRange(2) = size(Img{1},2);


%define data range, pad edges to account for beads near edges moving outside
%the range of the particles in the first image
%But, for these experiments use the whole range - there are usually some
%edge effect near the top and bottom, but we get a more complete view of
%the specimen
xList = MPTPara.xstep.*(MPTPara.xRange(1):grid_spacing(1):MPTPara.xRange(2)); %MPTPara.xRange(1)-2*MPTPara.edge_width*MPTPara.xstep:grid_spacing(1):MPTPara.xRange(2)+2*MPTPara.edge_width*MPTPara.xstep;
yList = MPTPara.xstep.*(MPTPara.yRange(1):grid_spacing(2):MPTPara.yRange(2)); %MPTPara.yRange(1)-2*MPTPara.edge_width*MPTPara.xstep:grid_spacing(2):MPTPara.yRange(2)+2*MPTPara.edge_width*MPTPara.xstep;
[yGrid,xGrid] = meshgrid(yList,xList);


disp('%%%%% Interpolating tracking results on ref grid %%%%%'); fprintf('\n');
%interpolate scattered data onto eulerian grid
u_inc = cell(1,length(resultDisp));
H_inc = cell(1,length(resultDefGrad));
for ii = 1:length(resultDisp)
    
    if ~mod(ii,10)
        disp(ii)
    end
    
    %only keep beads that are within the valid range of coordinates
    coords_cur_ = MPTPara.xstep.*resultDisp{ii}.parCoordB;
    keep_coords = ~[coords_cur_(:,1) > MPTPara.xstep.*MPTPara.xRange(2)|...
                   coords_cur_(:,2) > MPTPara.xstep.*MPTPara.yRange(2)|...
                   coords_cur_(:,1) < MPTPara.xstep.*MPTPara.xRange(1)|...
                   coords_cur_(:,2) < MPTPara.xstep.*MPTPara.yRange(1)];
    coords_cur_disp = coords_cur_(keep_coords,:);
    
    disp_cur = MPTPara.xstep.*resultDisp{ii}.disp_A2B_parCoordB(keep_coords,:);
    
    %interpolate displacements onto the reference grid 
    [x_grid_ref,y_grid_ref,u_inc{ii}{1}] = ...
                           funScatter2Grid2D(coords_cur_disp(:,1),coords_cur_disp(:,2),disp_cur(:,1),grid_spacing,smoothness,xGrid,yGrid);
    [~,~,u_inc{ii}{2}] = funScatter2Grid2D(coords_cur_disp(:,1),coords_cur_disp(:,2),disp_cur(:,2),grid_spacing,smoothness,xGrid,yGrid);

    %repeat the process for the strain data
    coords_cur_ = MPTPara.xstep.*resultDefGrad{ii}.XY_refA;
    keep_coords = ~[coords_cur_(:,1) > MPTPara.xstep.*MPTPara.xRange(2)|...
                    coords_cur_(:,2) > MPTPara.xstep.*MPTPara.yRange(2)|...
                    coords_cur_(:,1) < MPTPara.xstep.*MPTPara.xRange(1)|...
                    coords_cur_(:,2) < MPTPara.xstep.*MPTPara.yRange(1)];
    coords_cur_strain = coords_cur_(keep_coords,:);
                                         
    %collect the deformation gradient tensor components for each particle
    F_cur_vec = resultDefGrad{ii}.F_A2B_refA;
    F_cur_ = zeros(length(coords_cur_),4);
    F_cur_(:,1) = F_cur_vec(1:4:end);
    F_cur_(:,2) = F_cur_vec(4:4:end);
    F_cur_(:,3) = F_cur_vec(2:4:end);
    F_cur_(:,4) = F_cur_vec(3:4:end);
    F_cur = F_cur_(keep_coords,:);
    
    %convert deformation gradient tensor components at each step to
    %imcremental deformations on the reference grid
    [~,~,H_inc{ii}{1,1}] = funScatter2Grid2D(coords_cur_strain(:,1),coords_cur_strain(:,2),F_cur(:,1),grid_spacing,smoothness,xGrid,yGrid);
    [~,~,H_inc{ii}{2,2}] = funScatter2Grid2D(coords_cur_strain(:,1),coords_cur_strain(:,2),F_cur(:,2),grid_spacing,smoothness,xGrid,yGrid);
    [~,~,H_inc{ii}{2,1}] = funScatter2Grid2D(coords_cur_strain(:,1),coords_cur_strain(:,2),F_cur(:,3),grid_spacing,smoothness,xGrid,yGrid);
    [~,~,H_inc{ii}{1,2}] = funScatter2Grid2D(coords_cur_strain(:,1),coords_cur_strain(:,2),F_cur(:,4),grid_spacing,smoothness,xGrid,yGrid);

    %get a mask for "nan" values in the field
    F_interp = scatteredInterpolant(coords_cur_strain(:,1),coords_cur_strain(:,2),F_cur(:,1),'nearest','none');
    F_for_mask = reshape(F_interp(xGrid(:),yGrid(:)), size(xGrid));
    nan_mask_F{ii}{1} = F_for_mask./F_for_mask;
    
end

%%
% compute the cumulative displacement field in the referential
% configuration
disp('%%%%% Cumulating displacements %%%%%'); fprintf('\n');
msh{1} = xGrid;
msh{2} = yGrid;
[u_total,nan_mask] = inc2cum(u_inc,grid_spacing(1),msh,'cubic');


%%
disp('%%%%% Cumulating deformation gradient %%%%%'); fprintf('\n');

% compute the cumulative deformation gradiant (F = e - I)
F_total = cell(length(H_inc)+1,1);
F_total{1} = cell(2);

%set upthe identity tensor
for j = 1:2
    for k = 1:2
        if k == j
            F_total{1}{j,k} = ones(size(H_inc{1}{j,k}));
        else
            F_total{1}{j,k} = zeros(size(H_inc{1}{j,k}));
        end
    end
end

for ii = 1:length(H_inc)
    if ~mod(ii,10)
        disp(ii)
    end
    for loc = 1:numel(H_inc{ii}{1,1})
        H_inc_pt(loc,1) = H_inc{ii}{1,1}(loc);
        H_inc_pt(loc,2) = H_inc{ii}{2,2}(loc);
        H_inc_pt(loc,3) = H_inc{ii}{2,1}(loc);
        H_inc_pt(loc,4) = H_inc{ii}{1,2}(loc);
        
        %Convert to stretch
        Fhat_inc(1,1,loc) = H_inc_pt(loc,1)+1;
        Fhat_inc(2,2,loc) = H_inc_pt(loc,2)+1;
        
        Fhat_inc(2,1,loc) = H_inc_pt(loc,3);
        Fhat_inc(1,2,loc) = H_inc_pt(loc,4);
        
        [i,j,k] = ind2sub(size(xGrid),loc);
        if ii == 1
            %F_total = Finc for first step
            for m = 1:2
                for n = 1:2
                    F_total{ii+1}{m,n}(i,j,k) = Fhat_inc(m,n,loc);
                end
            end
        else
            %current F
            for m = 1:2
                for n = 1:2
                    Fcur(m,n,loc) = F_total{ii}{m,n}(loc);
                end
            end
            
            %F_total = F0*F1*F2...Fn
            F_tot_ = Fcur(:,:,loc)*Fhat_inc(:,:,loc);
            
            for m = 1:2
                for n = 1:2
                    F_total{ii+1}{m,n}(i,j,k) = F_tot_(m,n);
                end
            end
        end
    end
    
end

%alternatively, compute the deforamtion gradient from the gridded
%displacement directly - usually this is more lossy though
% F_total_ = calculateFij(u_total,grid_spacing(1),[1,1,1],'optimal5');

%Calculate strains from deformation gradient
[E_total,e_total] = calculateEij_2d(F_total);

%% save the results
results_file_names_Eul = fullfile('results',['results_full_2D_EulTotal_',file_names{1}(1:end-4),'.mat']);
if ~exist('results','dir') 
   mkdir('results')
end

resultDispInc = resultDisp;
resultDefGradInc = resultDefGrad;
save(results_file_names_Eul);

disp('%%%%% Complete cumulated Eulerian results saved %%%%%'); fprintf('\n');

%% plotting options

%plot mean strain comps
bd_wd = 1;
clear mean_strain_* std_strain_* ste_strain_* mean_disp_* ste_disp_*

%compute mean and ste strain at each step
for ii = 2:length(E_total)
    N = length(resultDefGrad{ii-1}.XY_refA);
    mean_strain_11(ii-1,1) = mean(nan_mask_F{ii-1}{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd)...
        .*E_total{ii}{1,1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd),'all','omitnan');
    mean_strain_22(ii-1,1) = mean(nan_mask_F{ii-1}{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd)...
        .*E_total{ii}{2,2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd),'all','omitnan');
    mean_strain_12(ii-1,1) = mean(nan_mask_F{ii-1}{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd)...
        .*E_total{ii}{1,2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd),'all','omitnan');

    
    ste_strain_11(ii-1,1) = std(nan_mask_F{ii-1}{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd)...
        .*E_total{ii}{1,1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd),[],'all','omitnan')/sqrt(N);
    ste_strain_22(ii-1,1) = std(nan_mask_F{ii-1}{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd)...
        .*E_total{ii}{2,2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd),[],'all','omitnan')/sqrt(N);
    ste_strain_12(ii-1,1) = std(nan_mask_F{ii-1}{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd)...
        .*E_total{ii}{1,2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd),[],'all','omitnan')/sqrt(N);
end


% Compute mean and ste strain at each step
%   NOTE: this is meaningful for rigid body motion estimation, but not
%         particularly useful otherwise
for ii = 1:length(u_total)
    N = length(resultDefGrad{ii}.XY_refA);
    mean_disp_1(ii,1) = mean(nan_mask_F{ii}{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd)...
        .*u_total{ii}{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd),'all','omitnan');
    mean_disp_2(ii,1) = mean(nan_mask_F{ii}{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd)...
        .*u_total{ii}{2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd),'all','omitnan');
    
    ste_disp_1(ii,1) = std(nan_mask_F{ii}{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd)...
        .*u_total{ii}{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd),[],'all','omitnan')/sqrt(N);
    ste_disp_2(ii,1) = std(nan_mask_F{ii}{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd)...
        .*u_total{ii}{2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd),[],'all','omitnan')/sqrt(N);
    
end

% Now plot the mean strain vs step number with error bars from the standard
% error
step_num = 1:length(mean_strain_11);
figure

shadedErrorBar(step_num,100*mean_strain_11,100*ste_strain_11,'b-x',1)
hold on
shadedErrorBar(step_num,100*mean_strain_22,100*ste_strain_22,'g-*',1)
shadedErrorBar(step_num,100*mean_strain_12,100*ste_strain_12,'m-o',1)
xlabel('Step number')
ylabel('Lagrange strain [%]')
title('Mean strain comps; Strain noise floor shaded region; b=e11,r=e22,m=e12')

% Plot the mean displacement vs step number with error bars from the standard
% error
step_num =  1:length(mean_disp_1);

figure
shadedErrorBar(step_num,-mean_disp_1,ste_disp_1,'b-x',1)
hold on
shadedErrorBar(step_num,-mean_disp_2,ste_disp_2,'r-+',1)
xlabel('Nominal displacement [\mum]')
ylabel('Mean measured displacement [\mum]')
% axis([0.002*0,0.002*85,-.1,0.25])
title('Mean displacement w std err; b=u1,r=u2')


%% Plot contour maps for a specific frame
frame_num = 10; %step number
lvls = 20; % number of contour levels
bd_wd = 1; % cropping width on the edge

% plot the strain contours
figure,
subplot(1,3,1)
contourf(msh{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd),msh{2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd),...
    flipud(fliplr(nan_mask_F{frame_num}{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd)...
.*E_total{frame_num}{1,1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd))),lvls,'linestyle','none')
colorbar
axis image
colormap(cmrMap)
title('Lagrange strain, E11 [-]')

subplot(1,3,2)
contourf(msh{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd),msh{2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd),...
    flipud(fliplr(nan_mask_F{frame_num}{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd)...
.*E_total{frame_num}{2,2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd))),lvls,'linestyle','none')
colorbar
axis image
colormap(cmrMap)
title('Lagrange strain, E22 [-]')

subplot(1,3,3)
contourf(msh{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd),msh{2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd),...
    flipud(fliplr(nan_mask_F{frame_num}{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd)...
.*E_total{frame_num}{1,2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd))),lvls,'linestyle','none')
colorbar
axis image
colormap(cmrMap);
title('Lagrange strain, E12 [-]')

% plot the displacement contours

figure
subplot(1,2,1)
contourf(msh{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd),msh{2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd),...
    flipud(fliplr(nan_mask_F{frame_num}{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd)...
.*u_total{frame_num}{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd))),lvls,'linestyle','none')
colorbar
axis image
colormap(cmrMap)
title('Total displacement, u_1 [\mum]')

subplot(1,2,2)
contourf(msh{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd),msh{2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd),...
    flipud(fliplr(nan_mask_F{frame_num}{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd)...
.*u_total{frame_num}{2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd))),lvls,'linestyle','none')
colorbar
axis image
colormap(cmrMap)
title('Total displacement, u_2 [\mum]')


%% Plot contours with an overlay of the orginal specimen image

% set transparency of the overlay image
OrigDICImgTransparency = 0.5;

% add the orginal image to the frame
fig1 = figure; ax1=axes; 
h1 = imshow(rot90(Img{frame_num}),[],'InitialMagnification','fit');
axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2);  set(gca,'ydir','normal');
hold on; ax2=axes; 

% plot the contours
u_field_curr = nan_mask_F{ii}{1}.*u_total{frame_num}{2};
[h2, hCont] = contourf(msh{1}./MPTPara.xstep,msh{2}./MPTPara.xstep,fliplr(u_field_curr),lvls,'linestyle','none');%show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,strain_vonMises,'NoEdgeColor');
set(gca,'fontSize',18); view(2); box on; axis equal;  axis tight; colormap(cmrMap); clim auto 

% now adjust to make the overlay correct and lined up
linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'none'; ax1.FontName = 'Ariel';
ax1.XLabel.String = 'x_1 [px]';
ax1.YLabel.String = 'x_2 [px]';
ax1.XLabel.FontSize = 18;
ax1.YLabel.FontSize = 18;
ax1.Title.String = 'u_2 displacement, \mum';
ax1.Title.FontSize = 18;

%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(ax1.XTick*10)/10, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(ax1.YTick*10)/10, length(ax1.YTick) )' );
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'none';
  
drawnow;  % this is important, to ensure that FacePrims is ready
hFills = hCont.FacePrims;  % array of TriangleStrip objects
[hFills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
for idx = 1 : numel(hFills)
   hFills(idx).ColorData(4) = round(OrigDICImgTransparency*255);   % default=255
end


