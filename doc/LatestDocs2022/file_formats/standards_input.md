---
id: standards_input
title: Standards Input Files
---

 ## Old Text Standards Input File (.txt)

This was an ad-hoc format that contained only the essential information for processing measured spectra from standard materials and producing element calibration factors. It contained just the spectrum file name and the percents of elements present in the standard. Each standard had a separate complete entry and the format was not very flexible. It is in the process of being replaced by a commas-separated-value format with much more information.

The file contained a single element fit list at the beginning. It was ignored but was required for proper reading of the file.

// Any line that starts with two slashes is a comment
// Calibration spectra and compositions for producing new calibration file

// Composition from "ReducedStandardsComposition.xls"

// The first non-comment line must be the element list for unknown analysis
// The number of elements followed by element symbols
// Blank or tab separated (not commas) and case sensitive
16 Na Mg Al Si P K Ca Ti V Cr Mn Fe Ni Cu Zn Sr

// Data from Standards_Emily_Oct2013
// With compositions from standardCompositions.txt of Dec. 2, 2013 3:20pm

BHVO-2G_100uA_28kV.mca
// The number of elements followed by element symbols and percents
// Blank or tab separated (not commas) and case sensitive
18 Na 1.6469 Mg 4.3599 Al 7.1449 Si 23.3250 P 0.2146 K 0.4317 Ca 8.1475 Ti 1.6362 V 0.0317 Cr 0.0280 Mn 0.1290 Fe 8.6029 Ni 0.0119 Cu 0.0127 Zn 0.0103 Sr 0.0389 Zr 0.0172 O 44.5883

## New CSV Standards Input File (.csv)

This format was developed to give much more control over individual standards, how they were fit, and how they are used to produce the element calibration factors (ECFs). It uses comma-separated values with few keywords and fixed fields for the elemental information. It is much more flexible and can be extended as necessary.

The new calibration file format is identical to this file. Note the ECF and ECF error fields at the end of each element line. If these fields are present in the standards input file they are ignored.

Comment, The information on each input line is given below
Comment, A more detailed description of the input values is given at the end of this file
Comment, This file is useable as-is
Comment, Each line consists of comma-separated values, either text or numbers
Comment, If a text value contains a comma, it must be enclosed in quotes "like this, for example"
Comment, Each line can start with either a keyword or an element symbol
Comment, Keywords are not case sensitive but element symbols ARE case sensitive

Comment, Each standard description starts with the Standard keyword, which resets all values
Comment, The standard keyword is optionally followed by one or more names for the standard
Standard,BHVO-2,"Basalt Hawaiian Volcanic Observatory"

Comment, The Comment keyword can appear anywhere and the rest of the line is ignored
Comment, except that comments are stored and written to the calibration file
Comment, Orig tot %: 100.2100

Comment, Each line that begins with an element symbol has fields in a fixed order
Comment, The fields are: Element symbol, Emission line [KLMN], Fit qualifier [XIFM], Type (see below), Amount, Uncertainty, Oxide ratio, and Weight
Comment, See the end of this file for a more complete description of these fields
Si , , , , 23.325%, 0.6a, 2, 1
Ti , , , , 1.6366%, 0.04a, 2, 1
Al , , , , 7.14475%, 0.2a, 1.5, 1
Fe , , , , 8.6029%, 0.2a, 1.5, 1
Mn , , , , 0.09990%, 0.004a, 1, 1
Mg , , , , 4.35988%, 0.12a, 1, 1
Ca , , , , 8.14739%, 0.2a, 1, 1
Na , , , , 1.6469%, 0.08a, 0.5, 1
K , , , , 0.43167%, 0.01a, 0.5, 1
P , , , , 0.11783%, 0.02a, 2.5, 1
V , , , , 317.0ppm, 11a, 0, 1
Cr , , , , 280.0ppm, 19a, 0, 1
Co , , , , 45.0ppm, 3a, 0, 1
Ni , , , , 119.0ppm, , 0, 1
Cu , , , , 127.0ppm, 7a, 0, 1
Zn , , , , 103.0ppm, 6a, 0, 1
Rb , , , , 9.8ppm, 1a, 0, 0
Sr , , , , 389.0ppm, 23a, 0, 1
Y , , , , 26.0ppm, 2a, 0, 1
Zr , , , , 172.0ppm, 11a, 0, 1

Comment, The Spectrum keyword gives the name of the measured spectrum file for this standard
Comment, This keyword triggers the fitting of the spectrum and writing the ECFs for this
Comment, standard to the calibration file.
SPECTRUM, USGS_BHVO2_He_28kV_20uA_1hr.mca

Comment, At this point, additional element lines can be entered, which will modify the
Comment, information already entered for this standard. Then another Spectrum line
Comment, can be entered to produce more calibration entries with the modified
Comment, information and a different measured spectrum.

Comment, Detailed description of the values on an element line (not case sensitive unless indicated)
Comment, Element symbol Case-sensitive, one or two letters, only the first one capitalized
Comment, Emission line [KLMN] Denotes which emission lines this entry applies to, allowed values are K, L, M, or N
Comment, (if no emission line is specified this applies to all lines for this element)
Comment, Fit qualifier [XIFM] Controls how the emission line will be treated in fitting the spectrum
Comment, Allowed values are X, I, F, and M
Comment, X means to exclude this emission line from the fitting
Comment, (it will not be calculated or used to fit the spectrum)
Comment, I means fit this line but ignore it, do not include it in the standard composition
Comment, F is reserved for future use
Comment, M is for matrix, it is treated the same as X
Comment, Type Spectrum component type, for future use in calibrating scatter peaks
Comment, Allowed values are inc (or Com) for incoherent or Compton scatter,
Comment, coh (or Ray) for coherent or Rayleigh scatter, and bkg for background
Comment, The default type is E for Element and is the only type used for now
Comment, Amount Amount of the element present in the standard composition, in percent
Comment, The amount can be followed by the percent sign (%)
Comment, The amount can also be entered as a fraction
Comment, in this case the number must be followed immediately by the letter f
Comment, The amount can also be entered as parts-per-million
Comment, in this case the number must be followed immediately by the letter p (can be ppm)
Comment, If amounts are entered for the same element on different lines, the last one read will be used.
Comment, Uncertainty Uncertainty in the amount of this element in the standard, default is relative percent
Comment, The uncertainty can be entered as an absolute value in the same units as the
Comment, amount but it must be followed immediately by the letter a (for absolute)
Comment, Oxide ratio Number of oxygen atoms attached to each atom of the element. Default is zero.
Comment, For example, for Al2O3 the oxide ratio would be 1.5
Comment, A value of -1 can be entered and will cause the default oxide ratio for this element to be used
Comment, Weight When computing the ECFs to be used to quantify an unknown, the ECF from this standard for this
Comment, element and emission line will be weighted by this factor. The default is unity.
Comment, Set this value to zero to include the element in the composition but not use the ECFs.
Comment, Zero is also useful for elements with small amounts whose peaks are too small for a good fit.
Comment, (The following two values are added when this information is written to the calibration file.
Comment, If they are included in the standards input file they will be ignored.)
Comment, ECF Element Calibration Factor
Comment, It is the fit coefficient that gives the ratio to actual peak size in the measured spectrum
Comment, relative to the calculated peak size using the fundamental parameters quantitative
Comment, calculation for this standard composition.
Comment, ECF sigma Uncertainty in the ECF as a relative percent based on the spectrum counting statistics
Comment, including any correlations in the fit.
Comment, Intensity Net peak intensity for this element emission line in the spectrum, from the fit.