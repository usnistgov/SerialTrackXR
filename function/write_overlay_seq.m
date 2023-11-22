%% Write out sequence of overlay images for a .gif of the deformation

contour_range = [-1,0.6];

for frame_num = 1:length(Img)

    fig1=figure; ax1=axes; 
    h1=imshow(rot90(imadjust(Img{frame_num})),[],'InitialMagnification','fit');
    axis on; axis equal; axis tight; box on; set(gca,'fontSize',12); view(2);  set(gca,'ydir','normal');
    hold on; ax2=axes; 
    
    E_field_curr = nan_mask_F{ii}{1}.*E_total{frame_num}{2,2};
    
    [h2, hCont] = contourf(msh{1}./MPTPara.xstep,msh{2}./MPTPara.xstep,fliplr(E_field_curr),lvls,'linestyle','none');%show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,strain_vonMises,'NoEdgeColor');
    set(gca,'fontSize',12); view(2); box on; axis equal;  axis tight; colormap(cmrMap); clim(contour_range); 
    
    linkaxes([ax1,ax2]);  % Link axes together
    ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
    colormap(ax1,'gray'); % Give each one its own colormap
    set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
    ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'none'; ax1.FontName = 'Ariel';
    ax1.XLabel.String = 'x_1 [px]';
    ax1.YLabel.String = 'x_2 [px]';
    ax1.XLabel.FontSize = 12;
    ax1.YLabel.FontSize = 12;
    ax1.Title.String = 'E22 strain';
    ax1.Title.FontSize = 12;
    
    %%%%% convert pixel unit to the physical world unit %%%%%
    xticklabels(ax1, num2cell(round(ax1.XTick*10)/10, length(ax1.XTick) )' );
    yticklabels(ax1, num2cell(round(ax1.YTick*10)/10, length(ax1.YTick) )' );
    cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'none';
      
    drawnow;  % this is important, to ensure that FacePrims is ready in the next line!
    hFills = hCont.FacePrims;  % array of TriangleStrip objects
    [hFills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
    for idx = 1 : numel(hFills)
       hFills(idx).ColorData(4) = round(OrigDICImgTransparency*255);   % default=255
    end

    saveas(fig1,['./results/imageseq2/',num2str(frame_num),'.tif'])
    close all

end
