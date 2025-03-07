---
id: calibration
title: Calibration Files
---

 ## Old Text Calibration File (.txt)

This file contained only the atomic numbers of elements for which ECFs were calculated from standards and the actual ECF values.

0 Lines that start with zero are skipped.
0 Number of element followed by atomic numbers
0 The immediately followed by a line of corresponding ECF values
0 (no comment allowed after the element list)
8 11 14 15 17 19 20 26 92 299
0.5686 0.5686 8.3526 2.2013 0.5686 0.2160 0.3353 0.4282 0.0000 0.0000
0 There are a few numbers at the end of the element and ECF lists.
0 They are ignored but must be present.

## New CSV Calibration File (.csv)

The new calibration file format is identical to the standards input file. Note the ECF and ECF error fields at the end of each element line. This information is added from the calibration results when the calibration file is written. The calibration file contains all of the information from the standards input file with the ECF fields added. The ECF values for each standard are included separately so that the ECFs used for a particular unknown can be adjusted as desired based on the complete information about each standard. Currently all standards are used and only the weights affect the calculation of ECFs for unknowns. The complete set of information is read in for possible future use in more sophisticated ECF calculations that may be individualized for each unknown. 