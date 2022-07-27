# po_matlab_codes
 
%   D E S C R I P T I O N
%   ------------------------------------------------------------------------ 
% 
%     This code is prepared to calculate the number of Tpms on single 
%     actin filament in the following steps;
% 
%     a) Tropomyosin clusters were fitted with a 2D Gaussian function 
%     b) Signal is iteratively deconvolved from the overlapping clusters, 
%        until the background level is reached. 
%     c) The number of fluorophores in each cluster is calculated by
%        dividing the total integral intensity of the diffraction-limited
%        spot of a cluster by that of a single fluorophore. 
 

%   S E T T I N G;
%   ------------------------------------------------------------------------
%   Parameters in this section are defines as below and set manually to
%   desire value 
%   (a) peaks within a distance of 2σ pixels of each other, equivalent to 
%   the standard deviation, σ,  of the point-spread function (PSF), 
%   (b) small intensity peaks (with σ <1) remaining after this process 
%   were assigned to noise, as it is physically impossible to have PSF 
%   with σ <1 pixel, and excluded from further analysis (c) peaks with 
%   intensity value equivalent to a σ/2 or smaller and more than three 
%   pixels were excluded, (d) filaments with only one type of tropomyosin 
%   (usually small filament) or only one of each tropomyosins, were also 
%   excluded.

%   MATLAB V E R S I O N
%   ------------------------------------------------------------------------
%   This code was tested on Matlab 2022 version 

%   H O W  TO  T E S T 
%   ------------------------------------------------------------------------
%   Run the main file "main.m"
%   Select any file in "example_img "folder 
%   A new folder should be created where the outcome of the analysis is stored


%  C O N T A C T   I N F O
%   ------------------------------------------------------------------------
   Need more help? contact me at peyman.obeidy@gmail.com
