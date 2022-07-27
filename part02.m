
%Part2; draw a line on the filament to generate a mask, this way every
%peak selected out side the mask will be ignored for analysis
% click along the filament and press Enter at the last point

imageData1b=imcrop(imageData1,rectfilament);
imageData2b=imcrop(imageData2,rectfilament);

    imcomp(imageData1b(:,:,1),imageData2b(:,:,1),'y','n');
    title("Draw a line on the filament to mask it")
    axis image;

current_fig_position =[100,100,1256,247.2];
set(gcf,'Position',current_fig_position)

lineSelect = impoly('Closed', false);

%ROI_area_trace_v4.p is required polygeom.p
%Change the mask size using ROI_length and ROI_width 

ROI_length = config.Mask_ROI_length; %def = 2;
ROI_width  = config.Mask_ROI_width;  %def=4;

Img_stack = cat(4, imageData1b, imageData2b);

ptraw = lineSelect.getPosition();
alongLineDist = [0; diag(squareform(pdist(ptraw)), 1)];
alongLineDist = cumsum(alongLineDist);
pt = interparc(round(max(alongLineDist)/ROI_length), ptraw(:,1), ptraw(:,2) , 'linear');

InterpPoints = reshape(pt(~isnan(pt)), [], 2);

points_here = InterpPoints;

xCh1=round(points_here(:,1));
yCh1=round(points_here(:,2));

set(gca, 'NextPlot', 'add');

if size(points_here, 1) > 1

    % Delete old line (assume red, dashed)
    delete(findobj('Parent', gca, 'color', 'r', 'LineStyle', ':'));
    ROI_line = plot(gca, pt(:,1), pt(:,2), ':', 'Color', 'r', 'LineWidth', 2);
   
    % Delete old start point (if it exists) and mark start point
    delete(findobj('Parent', gca, 'color', 'y', 'Marker', '*'));
    ROI_start = plot(gca, pt(1,1), pt(1,2), '*', 'Color', 'y', 'MarkerSize', 8);

    ROIstructHold = ROI_area_trace_v4(Img_stack(:,:,1,:), 1:2,...
        points_here, ROI_width, length(points_here), 1, 0);

    Segs_tot = ROIstructHold.Segs_tot;



    PixelIdxHold = cell(1, (size(Segs_tot, 1) - 1));
    for pol = 1:(size(Segs_tot, 1) - 1)

        % disp(pol)
        x_mask = [Segs_tot(pol, 4) Segs_tot(pol+1, 4) Segs_tot(pol+1, 6) ... 
            Segs_tot(pol, 6) Segs_tot(pol, 4)];
        y_mask = [Segs_tot(pol, 5) Segs_tot(pol+1, 5) Segs_tot(pol+1, 7) ...
            Segs_tot(pol, 7) Segs_tot(pol, 5)];

        maskHere = logical(poly2mask(x_mask, y_mask, size(Img_stack, 1), ...
            size(Img_stack, 2)));

        PixelIdxHold{pol} = find(maskHere);

    end

    FakeCC.Connectivity = 8;
    FakeCC.ImageSize = [size(Img_stack, 1), size(Img_stack, 2)];
    FakeCC.NumObjects = (size(Segs_tot, 1) - 1);
    FakeCC.PixelIdxList = PixelIdxHold;



    for chan = 1:size(Img_stack, 4);

        mean_mask = zeros((size(Segs_tot, 1) - 1), size(Img_stack, 3));
        std_mask = zeros((size(Segs_tot, 1) - 1), size(Img_stack, 3));

        for fk = 1:size(Img_stack, 3)

            regStruct = regionprops(FakeCC, Img_stack(:,:,fk,chan), 'MeanIntensity', 'PixelValues');

            for pol = 1:(size(Segs_tot, 1) - 1)

                mean_mask(pol, fk) = regStruct(pol).MeanIntensity;
                std_mask(pol, fk) = std(regStruct(pol).PixelValues);

            end


            ROI_struct(chan, fk).Mean_mask = mean_mask(:,fk);
            ROI_struct(chan, fk).std_mask = std_mask(:,fk);

            ROI_struct(chan, fk).x_center = ROIstructHold.x_center;
            ROI_struct(chan, fk).y_center = ROIstructHold.y_center;
            ROI_struct(chan, fk).center_dist = ROIstructHold.center_dist;
            ROI_struct(chan, fk).center_spacing = ROIstructHold.center_spacing;
            ROI_struct(chan, fk).Segs_tot = ROIstructHold.Segs_tot;
            ROI_struct(chan, fk).areas = ROIstructHold.areas;

        end

    end

    xx = [Segs_tot(:,4); flipud(Segs_tot(:,6)); Segs_tot(1,4)];
    yy = [Segs_tot(:,5); flipud(Segs_tot(:,7)); Segs_tot(1,5)];

    % Delete old line (assume red, solid)
    delete(findobj('Parent', gca, 'color', 'r', 'LineStyle', '-'));
    box_handles = plot(xx, yy, 'Color', 'r');

end

mask = poly2mask([Segs_tot(:,4); flipud(Segs_tot(:,6)); Segs_tot(1,4)], ...
    [Segs_tot(:,5); flipud(Segs_tot(:,7)); Segs_tot(1,5)], ...
    size(imageData1b, 1), size(imageData1b, 2));
mask2=mask;

pause(0.1)


imcomp(imageData1b,imageData2b,'y','n');
title('Select the rectangle around the background')
current_fig_position =[100,100,1256,247.2];
set(gcf,'Position',current_fig_position)
[Iback, rectnoise]=imcrop;

