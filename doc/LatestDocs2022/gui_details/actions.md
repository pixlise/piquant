---
id: actions
title: Actions
---

 The action buttons at the top of the interface allow the user to choose which step in the analysis process is to be performed. Each action corresponds to a radio button, where only one button can be selected at a time.

Each action button performs a particular step in the spectrum analysis process. In the following subsections details are given: what the step accomplishes, the inputs needed, and how to choose the inputs. The complete list of inputs is always displayed in the GUI but the ones that are not needed for the chosen action are grayed out and are not active. None are active when the GUI first appears and will remain that way until an action is chosen. Once an action is chosen, the active inputs must be given a value before the Go! button is pressed.

### Actions performed to obtain XRF spectra information using PIQUANT: 

![Actions hierarchy](../../../static/img/piquant/actions.png)

## Energy Calibration

This performs a calibration of the relation between channels in the pulse-height histogram and the energy of the incident X-rays. The algorithm takes as input the histogram and the elements corresponding to the two largest peaks in the histogram. The highest point in the histogram (the channel with the most counts) is found and taken as the first peak. The peak center is found via a 5-point average around the highest channel. An exclusion zone is then put into place around this peak from 15% above the peak channel to 10% below the peak channel. This is to prevent finding beta or other secondary peaks from the same element associated with this peak. The channel with the highest counts in the remaining part of the histogram is then found and treated similarly to find the center of the second peak. If either peak has fewer than 100 counts, an error message is generated and the calibration fails. The element fit control list is then inspected and the first two entries are used to associate the peaks with energy. The energy of the emission line from the elements in the list are used, and the K lines are used if no energy level is given. The order of the elements is not important, as the lower energy emission line will always be associated with the peak in the lower channel number. The energy of the starting channel and the energy per channel is computed and output to the screen.

For spectra with more than one detector, the energy calibration applies to the sum of the detectors. To obtain an energy calibration for a single detector use the –d option. Using –d,0 selects the first detector and the zero can be replaced by one for the second detector or any integer up to the number of detectors minus one. (See also Section 4.4.3.)

The energy calibration is also entered in the options text box. This will cause the new energy calibration to be used for any subsequent actions until a new energy calibration is run, the text is deleted from the options, or the GUI is restarted. This allows further processing of a histogram that has no energy calibration without having to enter the new energy calibration into the spectrum or configuration files. Note that the new energy calibration will override any energy calibration in the spectrum or configuration files. The new energy calibration remains in the command line options text box so that several spectra with the same energy calibration can be processed, or a calibration spectrum can be used to find the energy calibration and then other spectra taken at the same time can be processed using the new energy calibration.

For histograms with only one major peak, only one element can be entered in the element fit control list. In that case, the starting channel will be assumed to correspond to zero energy and only the energy per channel will be calculated.

### Inputs:

* Spectrum file
* Element fit control list (one or two elements), e.g., Ca Fe to designate the dominant calcium and iron Kα peaks found in BHVO2

### Outputs:

* Energy offset and slope to output text box and options text box

## Plot Spectrum

This action produces a plot of the spectrum read in from a spectrum file. The file format is detected automatically. The counts in each channel are plotted vs. energy or vs. channel number if no energy calibration is available. An energy calibration is sought in the options text box and the spectrum file in that order. (No configuration file is read in this case.)

This command can also be used to display plot files written earlier and any comma-separated-values file that has the title and headers on the first two lines in a form that matches plot files in addition to spectrum files.

### Inputs:

* Spectrum file or CSV plot file

### Outputs:	

* A plot file is written and displayed in the plot window

## Calculate Primary Spectrum

The fundamental parameters method is used to calculate the primary spectrum that is expected from the instrument configuration specified in the configuration file. The spectrum consists of the output of the X-ray tube (or other X-ray source) as modified by the filter material and optic (if any) in the primary beam. Because the spectrum consists of characteristic emission lines with very narrow energy width plus a continuum emission, the spectrum is calculated as if the primary beam were incident directly on the detector described in the configuration. (Do not try this in practice unless you have taken precautions to not overload and/or damage the detector, as the primary beam is typically too strong to enter the detector directly and will damage it without reducing its intensity somehow.) This approach yields a spectrum with the characteristic emission lines as peaks on the continuum background, which is easy to visualize in a plot. None of the materials in the emitted fluorescence nor the interactions in the target are included in the calculation, but the detector response is included (with escape peaks).

The resulting spectrum can thus be compared to a measurement of the primary spectrum if the beam is attenuated, or the solid angle reduced, sufficiently to get a measurement with an energy-dispersive detector. The calculation is multiplied by the live time in the configuration file, so set it to unity to get photons per second.

To obtain the actual spectrum that impinges on the target (rather than the spectrum that would be measured with a detector), set the detector active area to some very large number and its Be window thickness to zero. Note that the height of the characteristic emission lines will depend on the detector resolution and the continuum intensity will depend on the channel energy width, although the integral intensities for both will always be correct.

The plot file and its display in the plot window are the only outputs from this action, although the output text window has some useful ancillary information. Note also that an energy calibration must be available in the configuration file to perform this calculation, otherwise an error results and there is no plot output.

### Inputs:

* Configuration File

### Outputs:

* The calculated primary spectrum is written to the plot file and is displayed in the plot window

## Calculate Full Spectrum

The fundamental parameters method is used to calculate the spectrum expected from the instrument configuration specified in the configuration file and for the target composition of the first standard in the standards input file. The contributions from individual elements and energy levels are placed into different columns in the plot file, as are the background from scatter of the continuum and the characteristic peaks from the source via Compton and Rayleigh scatter. The full calculated spectrum obtained by summing these various contributions is also output.

The plot file and its display in the plot window are the only outputs from this action, although the output text window has some useful ancillary information including some information on PIXL L5 requirements for the X-ray Subsystem (XRS). Note that an energy calibration must be available to perform this calculation, otherwise an error results and there is no plot output.

### Inputs:	
* Configuration File
* Standards input file (reads only the first standard in the file)
### Outputs:	
* A plot file is written and displayed in the plot window

## Compare Measured to Calculated

This is identical to the full spectrum calculation except that a measured spectrum is read in and included in the plot. The energy calibration and number of channels are taken from the input spectrum.

### Inputs:

* Configuration File
* Standards input file (reads only the first standard in the file)
* Spectrum file

### Outputs:	

* A plot file is written and displayed in the plot window

## Fit One Standard with Plot

To improve quantitative elemental performance, PIQUANT uses standards with known composition to determine individual element calibration factors for each element. The calibration factors are found by fitting the standard using the same calculation and fitting routines as for unknown measurements but only retaining the fit coefficients and not modifying the composition of the standard. To achieve the best possible fit, there is considerable control over the fitting process.

When performing an actual calibration, many standards are processed at once via an automated sequence with the fit controls included in the standards input file. In this case the large number of plots are not produced. To make it easier to make the initial determination of the fit controls for a specific standard, this action will fit a single standard (the first one in the standards input file) and display a comprehensive plot of the fit quality. The display includes all of the contributions to the calculation as they were adjusted by the fit. The residual and the original measured spectrum are also displayed.

The text output has a list of all the components used in the fit along with their coefficients and an indication of the relative uncertainty of the coefficients assuming Poisson counting statistics for the measured spectrum. Note that this does not write a calibration file, it only produces the text and plot. Also, the measured spectrum file is read from the standards input file, not the spectrum file in the main window.

The element fit control list is used to determine which elements and which components are included in the fit to the measured spectrum. All elements in the standards composition will be included unless overridden by inputs in this list. This format and composition of this list are described in a separate section following (with the file selection boxes).

### Inputs:	

* Configuration file
* Standards input file (reads only the first standard in the file)
* Element fit control list (required)

### Outputs:	

* A plot file is written with the measured spectrum, fit, and all fit components
* The plot is displayed in the plot window

## Bulk Sum and Max Value

When a large number of short-dwell spectra are taken for a map, or several spectra are taken at different locations on a standard, it is often desirable to sum the spectra to produce a single spectrum with a longer dwell time. In such a sum, elements that appear in only one or a few individual spectra may get lost in the average. A spectrum that contains the maximum value for each channel across all the map spectra will highlight these rarely occurring elements.

The bulk sum and max value action takes a sequence of spectra as input and produces these two outputs – a sum across all the spectra and a spectrum with the maximum value for each channel. The spectra are aligned using the individual energy calibrations before summing. The sum is quantified after it is produced. The max value spectrum is only plotted and/or written to the plot file.

See the quantify action for more information on the processing of the sum, particularly the fit control list. See the map action for the sequence of spectrum names.

### Inputs:

* Configuration file
* Calibration file
* Spectrum file
* Element fit control list (required)

### Outputs:

* The output text box has entries with complete quantification and fit information.
* A plot file is written with the sum spectrum, maximum value spectrum, the fit, and all fit components from the quantification
* The plot is displayed in the plot window

## Calibrate

To produce the calibration file with the element calibration factors based on a wide set of standards, batch processing is used in the calibration action. The configuration file is used for the instrument description and applies to all the spectra processed together. The composition of each standard as well as the measured spectrum file are input via the standards input file.

The standards are processed as their information is read in from the standards input file. Any number of standards can be included and all will be processed once the Go! button is pressed. Any errors will prevent the writing of the calibration file, although processing of the standards will continue to locate any additional errors.

A calibration file will be written with the element calibration factors.

The element fit control list is used to determine which elements and which components are included in the fit to the measured spectrum. Any choices in this list will be applied to all spectra and will override any information in the standards input file. If nothing is specified in the standards input file then all elements in the given composition and all emission lines in the spectrum range will be included in the fit. This is usually not the best choice, so some entry in the element fit control list is required for this action. If the new standards input file format is used, any information in that file will be used but can be overridden using entries in the element fit control list.

### Inputs:

* Configuration file
* Standards input file
* Calibration file name (a new file will be created or an existing file overwritten)
* Element fit control list (optional)

### Outputs:

* A calibration file is written with the element calibration factors
* The output text box has an entry for each standard with complete calibration and fit information. (No plot or plot file is produced.)

## Quantify

This is used to determine the element abundances represented by the peaks in the measured spectrum of an unknown target. The instrument configuration used in the fundamental parameters calculation is input via the configuration file. The measured spectrum is input from the spectrum file name in the main window. The element calibration factors read in from the calibration file will be applied during the quantification.

Because the composition is unknown, the list of elements present in the target must be given in the element fit control list. If no energy level is chosen for the emission lines to be used to quantify the element, then the K lines will be used unless they are not in the spectrum range, in which case the L lines will be used. Likewise the M lines will be used if the K and L lines are not in the spectrum range, although this is rare. K, L, M or N lines can be chosen for quantification.

Element calibration factors are used to adjust the quantification results. Each factor is an average over the standards that contain the respective element. Any elements that are missing in all standards will use the average calibration factor calculated from the available elements. If the calibration file is the new format then the weights for each standard will be used in the average. If the calibration file is in the older text format, the average will be unweighted.

### Inputs:

* Configuration file
* Calibration file
* Spectrum file
* Element fit control list (required)

### Outputs:

* The output text box has entries with complete quantification and fit information.
* A plot file is written with the measured spectrum, fit, and all fit components
* The plot is displayed in the plot window

## Map

This action quantifies a sequence of spectra and produces a comma-separated-values file with the location of each spectrum (from the spectrum file) and the quantification results from each spectrum. The results for each spectrum appear on a single line in the output map file. This file can be used to create quantitative maps or line scans from the sequence of spectra. No plot is produced by this action, and it can take large amounts of time depending on how many spectra are in the sequence and how many elements are in the fit control list. (For 10,000 spectra it can take up to 10 hours; we are working to improve this.)

The sequence of spectra is specified by entering the file name of the first spectrum in the sequence. All spectra in numerical sequence and that can be found in the same folder or directory will be processed, starting with the given file name. Processing will stop as soon as a file in the sequence cannot be found. This can be modified using the –m option in the command line options to specify the number of files to process (e.g., –m,10 to process only the next 10 files; see Section 3.3 on command line options).

The spectrum files should have sequential numbers in the file names. Two formats for these sequence numbers are accepted. Once has the number at the end preceded by an underscore (e.g., spectrum_1.msa, spectrum_2.msa, …). The number must be immediately before the dot at the end of the file name. This is the name sequence produced by the LabView program that runs the breadboard and Stony Brook instruments. The second format is to include the number immediately after the letters Seq (case sensitive, e.g., abcSeq123spectrum.msa, abcSeq124spectrum.msa, …). The Seq and number digits can be anywhere in the file name but the characters Seq must not appear anywhere else and there must not be anything else between the Seq characters and the number. This is a historical format from the Mars Borehole instrument and is not presently used by PIXL. The sequencing code is modular and versatile, so other ways to indicate the sequence number in the file name can be implemented in future versions if necessary.

### Inputs:

* Configuration file
* Calibration file
* Spectrum file
* Element fit control list (required)

### Outputs:

* The output text box has entries with complete quantification and fit information for all the spectra
* A map file is written with the location and element percents for all spectra, one spectrum per line, and each element always in the same column
