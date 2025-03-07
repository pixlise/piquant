---
id: plot
title: Plot Files
---

 The plot file is a comma-separated-values file. The first line is the title of the plot. The second line contains the headers for each column. These headers will be used in the legend for the plot. The remaining lines have the energy (or channel number) in the first column and the energy-dependent spectrum information in subsequent columns. The number of columns is determined automatically and the spectrum data runs to the end of the file. This is an example of a plot file with most of the spectrum data removed. In the plot file, each component of the spectrum fit has the background added to it for better visual overlay with the measured spectrum.

PIQUANT X-ray Spectrum Energy (keV), meas, calc, bkg, sigma, residual, Fe_K, Ca_K, Na_K, Rh_L_inc
3.60761, 57576, 55720.1, 1292.67, 239.954, 1855.94, 1292.67, 55720.1, 1292.67, 1292.67
3.61764, 82441, 81193.8, 1300.26, 287.129, 1247.23, 1300.26, 81193.8, 1300.26, 1300.26
3.62767, 113485, 112983, 1307.53, 336.878, 502.055, 1307.53, 112983, 1307.53, 1307.53
3.6377, 150707, 149950, 1314.3, 388.213, 756.938, 1314.3, 149950, 1314.3, 1314.3, 1314.3, 1314.3
3.64773, 190301, 189685, 1320.82, 436.237, 616.016, 1320.82, 189685, 1320.82, 1320.82
3.65776, 229284, 228609, 1326.33, 478.838, 674.734, 1326.33, 228609, 1326.33, 1326.33
3.66779, 262717, 262442, 1331.33, 512.561, 274.719, 1331.33, 262442, 1331.33, 1331.33
3.67782, 285616, 286943, 1336.08, 534.432, -1327.44, 1336.08, 286943, 1336.08, 1336.08
3.68785, 298316, 298778, 1340.64, 546.185, -461.688, 1340.64, 298778, 1340.64, 1340.64
3.69788, 296816, 296263, 1345.11, 544.81, 552.875, 1345.11, 296263, 1345.11, 1345.11
3.70791, 280309, 279760, 1348.69, 529.444, 548.938, 1348.69, 279760, 1348.69, 1348.69 