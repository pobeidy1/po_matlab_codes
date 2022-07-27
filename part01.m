
imageData1 = imread(fullfile(pName, fName),1);
imageData1=double(imageData1);
imageData2 = imread(fullfile(pName, fName),2);
imageData2=double(imageData2);
imcomp(imageData1,imageData2,'y','n');

%current_fig_position = get(gcf,'Position') % get the position of the current figure
current_fig_position =[100,100,1256,247.2];
set(gcf,'Position',current_fig_position)
title('Select the rectangle around the filament and Enter to Crop','FontSize',14);
[Icrop, rectfilament]=imcrop;
close all force 

