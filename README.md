**DESCRIPTION**

This code is designed to calculate the number of tropomyosin (Tpm) molecules present on a single actin filament in the following steps:

a) Tropomyosin clusters are fitted using a 2D Gaussian function.

b) Signal is deconvolved iteratively from the overlapping clusters until the background level is reached.

c) The number of fluorophores in each cluster is calculated by dividing the total integral intensity of the diffraction-limited spot of a cluster by that of a single fluorophore.

**SETTINGS**

The parameters in this section are defined and manually set to the desired value, as follows:

(a) Peaks within a distance of 2σ pixels of each other, equivalent to the standard deviation, σ, of the point-spread function (PSF).

(b) Small intensity peaks (with σ < 1) remaining after this process are assigned to noise, as it is physically impossible to have a PSF with σ < 1 pixel and excluded from further analysis.

(c) Peaks with intensity values equivalent to a σ/2 or smaller and more than three pixels are excluded.

(d) Filaments with only one type of tropomyosin (usually small filaments) or only one of each tropomyosin are also excluded.

**MATLAB VERSION**

This code was tested on Matlab 2022 version.

**HOW TO TEST**

To test this code:

a) Run the main file "main.m".

b) Select any file in the "example_img" folder.

c) A new folder will be created where the outcome of the analysis is stored.

**CONTACT INFO**

Need more help? Contact me at peyman.obeidy@gmail.com.



