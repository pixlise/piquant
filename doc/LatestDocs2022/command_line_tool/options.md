---
id: options
title: Command Line Options
---

 A few options can be added to the end of the argument list. They are denoted by a minus sign as the first character. Each option must not contain any blanks, although any number of options can be included and must be separated by one or more blanks (the system argument separator is relied upon for this task so its rules apply). No arguments will be inspected for options until all required arguments are found, so required file names that start with a minus sign are (in principle) possible, but not optional file name arguments. Numbers to be entered with the options must be separated from the option letter by a comma and subsequent numbers must be separated by commas with no blanks or tabs.

Options are under active development so check for what options are available on the current version. Options may be added without documentation in the future to aid in debugging, testing, code development, or other special purposes.

## List of Options

As of July 13th, 2021: 

* -e        Energy calibration
* -b        background parameters for overall or for low energy background
* -bh       background parameters for high energy background
* -bx       background parameters to control crossover from low to high energy background
* -T        Control for detector shelf adjustment factor and slope vs energy
* -d        Choose which of multiple detectors to include
* -m        Maximum number of spectrum files to read for map
* -q        Specify outputs to map file
* -f        Turn off adjustments to energy calibration in fits
* -g        Turn off adjustments to detector resolution in fits
* -c        Treat some elements as carbonates instead of oxides
* -t        Number of threads to create for map processing (for multiple CPUs)
* -s        Select standard from input file by number or name
* -w        Minimum weight in stds file for inclusion in evaluate output
* -u        Output evaluation file during Calibration or plot file during Evaluate
* -n        Normalization of element sum in percent (default is not mormalized)
* -Fe       Change iron default oxide ratio

## Energy Calibration Option, -e

An energy calibration option is available to override the calibration in the spectrum and/or configuration files, or to enter a calibration if one is not present in either of those files. It consists of a minus sign and the lower case letter e followed by two numbers separated by commas. The energy of the first channel and the energy increment per channel follow the letter e separated by commas. For example: –e,–5.25,10.3 (units are electron Volts and electron Volts per channel, and the energy per channel must be positive). At present there is only one energy calibration allowed, so it applies to all detectors if there is more than one in the spectrum file. In the future this may be handled differently as we gain more experience with multiple detectors.

## Background Control Options, -b, -bh, -bx

Options are included to override the defaults for background removal. Each consists of a minus sign and the lower case lettera b, bh, or bx.  -b controls the main background and will apply to the full spectrum unless the -bh or -bx options are used (see below).  For a SNIP background, the option is followed by up to 6 numbers separated by commas. The first number is the channel where the background algorithm will start processing, and the values in any channels below this will be ignored. The second number is the width in channels for the smoothing component of the background algorithm. The third number is the number of iterations to be performed by the background algorithm. These numbers can be used empirically by observing their effects on the net spectrum and the background curve displayed in the plots. To understand their action more fully see more on the SNIP background algorithm [Van Grieken and Markowicz, Handbook of X-Ray Spectrometry 2nd Ed., 2001]. For example, entering the option –b,25,12,24 would cause the background algorithm to start at channel 25, use a width of 12 channels, and perform 24 iterations. Note that all values are in channels, not energy, and must be integers. Any of the three values can be omitted and the default values will be used, but the commas must still be present if any values follow an omitted value.

If the start channel is zero, then the SNIP algorithm will start at the channel corresponding to the minimum energy in the configuration. If there is no energy calibration, or the minimum energy is zero, the SNIP algorithm will start at the lowest non-zero channel. If the width is zero it will be set to 12 and if the number of iterations is zero it will be set to 24.

There are now two extensions to the -b option. The first is as above but entering more than 5 parameters will cause a two-region SNIP. The 4th parameter is the start channel for the second region, the 5th parameter is the end channel for the second region, and the 6th parameter is the width for the second region. The number of iteration in the second region is the same as in the first region.

If the second parameter is negative, then a background is calculated using fundamental parameters and the sample composition. That background is used instead of the SNIP background. In this case, if the second parameter is zero (or negative) then the calculated background will be adjusted as part of the least squares fit. If the second parameter is positive, then the background will not be fit but it will be multiplied by the second parameter. This allows the calculated background to be adjusted manually for a better match.


In PIQUANT Version 3.1.2 and later, there are three background options:

* -b low-energy (or full-energy) background
* -bh high-energy background
* -bx low to high energy crossover (-bx,center,width in eV

-b or -bh parameters:
* -b,-1 use calculated background with scale factor = 1
* -b,-1,s scale calculated curve by scale factor s (see below)
* -b,0 use SNIP background with default parameters
* -b,0,s use SNIP with default parameters and scale factor
* -b,start,width,iterations,start2,end2,width2,scale
    * start = start channel
    * (start=zero => start 2 channels above the first non-zero channel)
    * width = smoothing width in channels
    * iterations = number of times to apply SNIP procedure
    * start2 = start of 2nd zone (with different width)
    * end2 = end of 2nd zone (return to original width)
    * width2 = width in 2nd zone
    * scale = scale factor as described below (=1 if omitted)

Scale factor s:

* s > 0 multiply background by a constant (default = 1)
* s = 0 find scale factor as part of full spectrum least-squares fit
* s < 0 use scale-under-peaks algorithm (described at right)

## Tail Options, -T

PIQUANT v.313 has the scale factor fixed for the calculated background and the -T option to control the extended detector shelf.

-T takes 1 to 3 parameters: the multiplicative factor, the slope, and the slope start. 

## Detector Selection Option, -d

For spectrum files that have more than one detector, the detector signals are summed prior to processing the spectrum. If only one detector should be processed, this option can select any of the detectors in the spectrum file and process only that signal. The option consists of a minus sign and the lower case letter d, followed by a comma and a single integer giving the number of the detector to be selected. The first detector is number zero. The numbers correspond to the columns in the spectrum file, starting with zero and going up to the number of columns minus one. An error results if the number is equal to or larger than the number of columns.

## Maximum Map or Bulk Sum Spectra, -m

For maps and bulk sum / max channel actions, the spectrum file gives a list of tiles (or the first file name in a sequence of spectra). Th -m option allows only a limited number of spectra in the sequence to be processed. Normally all the spectra that can be found with sequentially numbered file names will be processed. If this option is given processing stops after the indicated number of spectra are found. The format of the option is a minus sign and the lower case letter m, followed by a comma and a single integer giving the number of spectra to process. An example is –m,5 to process five spectra starting with the spectrum file given in the file selection box or command line argument. Four additional spectrum files will be processed (if found). If any spectrum file is not found processing stops with an error. See the section on Maps for more information about the spectrum file name sequencing.

## Map File Output Options, -q

For maps an output CSV file is generated with one line of data for each spectrum processed as part of the map. The information on these lines can be controlled via the –q option. The –q is followed by a comma and a set of letters strung together without any spaces or other separators (for example, -q,PIET to include percents, intensities, errors, and total counts). Each letter designates a particular type of information to include given in the list below. Some entries are one per element, and including that letter will include that entry for all of the elements in the element fit control list (in the order they appear in the element fit control list). Other entries are one per spectrum and will appear only once on each line. All of the information included in the map output file will be included on the header line (usually the second line in the file after the title line). Note that the letters are case sensitive. In general, upper case letters are used for values calculated as part of the quantitative analysis and lower case letters are used for auxiliary information from the spectrum file that is not used in PIQUANT but included in the map file for downstream processing. Output headers and values appear in the same order as their corresponding letters are entered in the –q option.

### List of Map Output Options

Options appear in map file in the order they appear in the -q input string. The default options are pPIETVXCRNFetsr.

* P        percents
* I        intensities
* E        fit relative errors
* L        fit coefficients
* K        Element Calibration Factors used for quantification
* G        Given percents for standards
* H        Element Calibration Factors used for quantification
* W        Matrix effect factor found during fundamental parameters calculation
* T        total counts
* X        reduced chi squared
* C        energy calibration
* R        detector resolution
* N        number of iterations
* F        file name
* S        Element sum %
* Q        sequence number in file name
* V        Live time
* M        Real time
* 7        Counts in 1 to 7.25 keV region
* x        X position     
* y        Y position   
* z        Z position  
* i        I position (in context image)   
* j        J position (in context image)
* s        SCLK       (spececraft clock time) 
* r        RTT        (Real Time Tracking token)
* d        DPC        (PIXL data product category number)
* p        PMC        PIXL motion counter
* e        Events        (DSPC total X-ray events in spectrum)
* t        Triggers      (DSPC total X-ray events detected including rejected pileup events)
* o        Overflows     (DSPC values above upper spectrum limit)
* u        Underflows    (DSPC values below lower spectrum limit)
* b        Baseline_samples   (number of DSPC baseline samples collected)
* a        Resets        (number of preamp resets detected)
* s        Saturates     (number of times DSPC ADC was saturated)
* l        Fast_livetime  (also called DSPC_Livetime, real time corrected for fast channel busy, not actual spectrum counting live time, does not account for rejected events)      
* n        Universal Sequence Number (USN) - location in PIXL non-volatile memory        
* U        Title

## Fit Control Options, -f and -g

These two options control the behavior for fitting the calculated spectrum components to the measured spectrum for unknowns and standards. Normally the energy calibration is adjusted to match the measured spectrum to the calculated energies so that the peaks align properly. The peak width (energy resolution) is also adjusted for the calculated spectrum to match the measured peak widths. Both of these adjustments assure that the intensities, taken from the calculated peaks, are as accurate as possible. Failing to match the energy calibration in particular will bias the calculation to smaller peak heights. However, some short integration spectra have enough noise that these adjustments cause unacceptable variations when fitting the spectra. For this reason the adjustments can be disabled using these options. –f disables the adjustment of the energy calibration. -g disables the adjustment of the detector resolution. All peaks are calculated using the resolution in the configuration file. Note that disabling the energy calibration adjustment may cause the resolution adjustment not to be done since the peak positions must match fairly well for the resolution adjustment to work properly. Likewise the peak energy shift can not be evaluated correctly if the energy mismatch is too great. The algorithm has built-in checks that prevent adjustment of either the energy or resolution if any of these checks fail.

* -f Disable adjustment of the energy calibration during spectrum fits
* -g Disable adjustment of the peak widths (detector resolution) during spectrum fits

## Carbonates Option, -c

If this option is included, the following elements will be treated as carbonates with the respective formulas: CaCO3, MgCO3, FeCO3, MnCO3, and SrCO3.

## Threads Option, -t

If this option is included, map processing only will be separated into the specified number of threads for processing on multiple CPUs.  This option has no effect for any other action.

## Normalization Option, -n

If this option is included, all quantifications will be forced to sum to the indicated percentage.  For normalization to 100%, use -n,100.

## Iron default oxide ratio Option, -Fe

This option is used to change the default oxide ration for iron (Fe) from its usual value of 1 (FeO).  For example, to treat all iron as Fe2O3 use -Fe,1.5 (oxide ratio is the number of exygen atoms per element atom).  Any oxide ratios entered in a standards input file take precedence over this value.

