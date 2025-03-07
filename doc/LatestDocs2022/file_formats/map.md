---
id: map
title: Map Files
---

 The map file is a comma-separated-values file with the location of each spectrum and the percents of each element in columns, with each spectrum on a separate line. The location information is precisely as read from the spectrum data file. The first line is a title from the spectrum file or the first spectrum file name. The second line is a list of the column headings including the element symbols for any element percent columns. The final columns may contain some diagnostic information, such as the goodness-of-fit parameter (reduced chi squared), to help distinguish any poorly-quantified spectra. This is an example of a simple map file produced by PIQUANT.

10-7-16-1 Fresh Surface X, Y, Z, Ca, Ti, Mn, Fe, Co, chisq
0.0000, 0.0000, -3.3500, 5.3299, 0.0225, 0.0401, 0.3760, 0.0376, 0.22
0.100000, 0.000000, -3.350000, 5.8467, 0.0116, 0.0397, 0.2188, 0.0105, 0.21
0.200000, 0.000000, -3.350000, 6.5944, 0.0006, 0.0484, 0.1590, 0.0024, 0.25
0.300000, 0.000000, -3.350000, 7.8825, 0.0003, 0.0665, 0.1776, 0.0010, 0.28
0.400000, 0.000000, -3.350000, 9.2958, 0.0000, 0.0850, 0.2114, 0.0158, 0.31
0.500000, 0.000000, -3.350000, 9.9761, 0.0027, 0.1363, 0.2460, 0.0157, 0.65
0.600000, 0.000000, -3.350000, 9.7767, 0.0064, 0.1133, 0.2620, 0.0029, 0.58
0.700000, 0.000000, -3.350000, 7.2178, 0.0100, 0.0690, 0.2177, 0.0011, 0.36
0.800000, 0.000000, -3.350000, 4.5039, 0.0268, 0.0352, 0.4258, 0.0000, 0.23
0.900000, 0.000000, -3.350000, 1.5232, 0.0164, 0.0145, 0.1541, 0.0000, 0.15
1.000000, 0.000000, -3.350000, 0.1047, 0.0030, 0.0000, 0.0074, 0.0000, 0.08 