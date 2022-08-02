// Copyright (c) 2018-2022 California Institute of Technology (“Caltech”) and
// University of Washington. U.S. Government sponsorship acknowledged.
// All rights reserved.
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of Caltech nor its operating division, the Jet Propulsion
//   Laboratory, nor the names of its contributors may be used to endorse or
//   promote products derived from this software without specific prior written
//   permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

//
//  PIQUANT_CommandLine.cpp
//  PIQUANT_CommandLine
//
//  Created by W. T. Elam on 1/14/2017.
//  Copyright (c) 2017 APL/UW. All rights reserved.

//  Version 2.alpha.79 released to Chris Heirwegh for testing March 22, 2017
//  Version 2.alpha.100 to Chris Heirwegh for testing May 1, 2017
//  Version 2.alpha.103 to Chris Heirwegh for testing May 15, 2017
//  Released as Version 2.00 on June 9, 2017
//  Released as Version 2.01 to Chris Heirwegh on June 27, 2017
//  Released as Version 2.02 to Chris Heirwegh on Sept. 30, 2017
//  Released first SEND_SDD_DATA histogram conversion tool on Nov. 29, 2017
//  Released as Version 2.10 on Dec. 15, 2017
//  Released as Version 2.2 on Jan. 8, 2018
//  Released as Version 2.21 on Jan. 31, 2018 to Chris and Joel
//  Released as Version 2.22 on Mar. 3, 2018 to Chris for testing -q option
//  Released as Version 2.23 on Mar. 7, 2018 to Chris for testing file list input
//  Released as Version 2.24 on Apr. 11, 2018 to Chris to fix nan resolution problem
//  Released as Version 2.25 on Apr. 18, 2018 to Matt to fix poor resolution fits to EM BHVO spectra
//  Released as Version 2.26 on Aug. 2, 2018 to Chris, Alan, and Les for testing
//  Released as Version 2.30 on Sep. 21, 2018 to PIXL team (for EM I&T)
//  Released as Version 2.41 on May 22, 2019 to Chris Heirwegh for use in PIXL FM Elemental Calibration
//  Released as Version 2.42 on June 9, 2019 to Chris Heirwegh for use in PIXL FM Elemental Calibration
//  Released as Version 2.43 on July 8, 2019 to Chris Heirwegh for use in PIXL FM Elemental Calibration
//  Released as Version 2.45 on Aug. 2, 2019 for use at Denver X-ray Conference QA II Workshop
//  Released as Version 2.46 on Nov. 4, 2019 to Chris Heirwegh for use in PIXL FM Elemental Calibration
//  Released as Version 2.47 on Nov. 7, 2019 to Chris Heirwegh for use in PIXL FM Elemental Calibration
//  Released as Version 2.5 on Jan. 10, 2020 to Chris Heirwegh and Joel Hurowitz
//          This will likely be the version that is transitioned to the JPL GitHub
//  Released as Version 3.0.1 on Jan. 6, 2021 to Chris Heirwegh for further testing on PIXL elemental calibration
//  Released as Version 3.0.2 on Jan. 8, 2021 to PIXLISE and PIXL GDS groups
//  Released as Version 3.0.4-V3_bug_fixes on Feb. 4, 2021 to Chris Heirwegh for further testing on PIXL elemental calibration with new bkg
//  Released as Version 3.1.0-V3_bug_fixes on Apr. 12, 2021 to Chris Heirwegh for final testing on PIXL FM elemental calibration data (UnityECFs)
//  Released as Version 3.1.2-V3_bkg_options on May 12, 2021 to Chris Heirwegh for final testing on PIXL FM elemental calibration data (UnityECFs)
//  Released as Version 3.1.3-V3_bkg_options on May 14, 2021 to Chris Heirwegh for final testing on PIXL FM elemental calibration data (UnityECFs)

#include <exception>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <time.h>
#include "parse_arguments.h"
#include "parse_element_list.h"
#include "read_spectrum_file.h"
#include "setupStandardsTXT.h"
#include "setupStandardsCSV.h"
#include "Element.h"
#include "borehole_read.h"
#include "read_EMSA_PIXL.h"
#include "map_spectrum_file_increment.h"
#include "energy_calibration.h"
#include "fpSetupConditions.h"
#include "quantBackground.h"
#include "quantCalculate.h"
#include "quantPrimarySpec.h"
#include "split_component.h"
#include "quantOpticResponse.h"
#include "quantStandard.h"
#include "quantUnknown.h"
#include "quantWriteResults.h"
#include "quantWriteCalibrationTXT.h"
#include "quantWriteCalibrationCSV.h"
#include "quantWriteMap.h"
#include "quantWritePlot.h"
#include "XraySpectrum.h"
#include "upper_trim.h"
#include "fpMain.h"
#include "XRFconditions.h"
#include "XRFconstants.h"
#include "XRFcontrols.h"
#include "XRFutilities.h"
#include "write_conditions.h"
#include "quantCombineSpectra.h"
#include "histogram_from_SDD_data.h"
#include "write_EDRhistogram_data.h"
#include "spectrumBulkSumMax.h"
#include "time_code.h"
#include "debug_stack.h"
#include "map_threading.h"
#include "setup_spectrum_parameters.h"
#include <mutex>
#include <thread>
#include <chrono>
#include "fpBeams.h"

#include "version.h"

#include "split_component.h"


//  Started development Jan. 14, 2017
//      This is a major re-write of the PIQUANT XRF analysis package
//      The primary reason is to allow separate fit coefficients for the K, L, M,and N
//          series of X-ray emission lines (rather than a single fit coefficient per element)
//      This will also change to a single subprocess with sub-commands.  This will
//          avoid duplication of code and allow a single GUI possible, as preferred by users.
//      The code will be reorganized to use the new XraySpectrum and XrayMaterial classes
//          to store the spectrum components and allow arbitrary composition of all
//          materials in the beam paths.  It will also facilitate future code to
//          handle multilayered targets.
//      These changes are also intended to make the package more flexible and easier to maintain.
//  Modified May 9, 2017 to append terminal output in case it is used as a log file for command line usage
//  Modified May 14, 2017
//      More digits in energy calibration slope for output text
//      Add new keyword to EMSA files for eMin (minimum energy for lines in spectrum)
//      Add an option to control where background removal starts, width, & # of iterations
//          (-b,start_ch,width_ch,iterations)  (include these in the Spectrum object as well as arguments struct)
//      Fix use of energy calibration in options for FIT_ONE_STD and CALCULATE sub-commands
//      Correct energy of ignored elements when energy calibration is changed by fit
//      Use an average of SNIP with and without LSQ (one is too large and one is too small, try average)
//  Modified May 26, 2017
//      Add new breadboard optic function (type 5)
//      Change linear least squares fit to prevent nan or inf returns
//      Make sure config energy calibration is used for compare
//      Change component description function to return string
//      Be sure one fit iteration is done after disabling negative coefficients
//  Modified June 7, 2017
//      Fix error in calculation of energy per channel and offset when tilt corrections were non-zero
//      Allow background to be a fitted component (change residual to meas vs calc, was net vs calc)
//      Fix convergence criteria (to be relative instead of absolute) - also increased max iterations back to 40
//  Modified June 9, 2017
//      Write date and time in ISO format to output header
//      Improve error processing and error messages for input files
//      Add background as a fitted component and fix exclusion of negative components in quantUnknown
//      Fix bug in writing components to plot file when no background present
//      Prevent crash in energy calibration when input spectrum has no channels
//  Modified June 26, 2017   V 2.01
//      Fix specifying and reading optic transmission file
//  Modified July 18, 2017
//      Consolidate code for selecting energy calibration and conditions to end of this file
//      Add map capability (read sequence of spectrum files and quantify each, writing line to map file)
//  Modified Sept. 13, 2017
//      Fix energy calibration for L lines
//      Update spectrum background when background components or coefficients are updated
//      Plot using only spectrum background, not background component (for consistency in plot)
//  Modified Sept. 23, 2017
//      Add bulk sum and maximum value sub-command for maps
//      Start working on detector selection and summing
//  Modified Sept. 29, 2017
//      Fix map spectrum file increment to handle "Seq" name format in addition to number at end after underscore
//      Write configuration info used in FP calculations (from conditionsVector with keywords and units)
//      Write warning message when calculated intensity is zero (would like to write which factor was zero, where to put it?)
//      Changed window type 3 from Brass to Carbon Fiber Reinforced Polymer (composition unknown, pure C for now)
//      Added option m for maximum number of map spectrum files to read (for testing and debugging, maybe later for multiprocessor scripts)
//      Check resolution change to be sure it's not too large before applying (also check for negative values)
//      Fix bug in new resolution calculation for one peak - forgot to take square root (this and above change in quantFitSpectrum)
//      Skip disabled components in loop to re-calculate ignore components  (both this and next change in quantUnknown & quantStandard)
//      Also disable ignore components with zero, negative, or nan calculated intensity
//      Check for very small components compared to largest and don't include in fit (to improve stability)
//      Two detectors (simple channel-by-channel sum in this version)
//      Return invalid arguments error messages in arguments.invalid_arguments instead of writing them to cout
//  Modified Oct. 27, 2017 (preparation for implementing CSV standards file)
//      Standards list struct moved to "XRFstandards.h" (include only one spectrum file name per entry)
//      Separate parser for individual element string with qualifiers
//      Add information for holding standards description to element list entries
//      Change setupStandardsTXT routine to use the element list entries from parse_element_list
//  Modified Nov. 1, 2017 (preparation for implementing CSV standards file)
//      Change all checks for file extension (.txt, .msa, etc.) to check only last few characters in name, and be case insensitive
//      Set up to read csv standards input files
//  Modified Nov. 8, 2017
//      Add ems sub-command to convert output of SEND_SDD_DATA command to EDR (csv) format
//  Modified Dec. 6, 2017
//      Finished testing ems sub-command with simulated iFSW data
//      More work on CSV standards input and processing plus CSV calibration files
//  Modified Dec. 8, 2017
//      Add SpectrumComponentType to element list and read from standards input CSV file
//          (so coherent and incoherent scatter could have calibration factors)
//      Separate input in standards CSV file into individual comma-separated fields
//  Modified Dec. 13, 02017
//      Add check for minimum energy to escape peaks in fpLineSpectrum (and all calls to it)
//  Modified Dec. 15, 2017
//      Add rebin of spectra before combining using individual energy calibrations
//      Don't plot other detectors if one detector is selected
//  Modified Dec. 23, 2017
//      Improve error messages in read_XIA_PIXL.cpp
//      Fix error combining single spectrum without energy calibration for plot only
//      Extend quantECFs to use standards list with weights if available
//      Error message and return if quant component has zero intensity, warning only for other components
//      Check for zero live time for any commends that require calculations
//  Modified Jan. 3, 2018
//      Move tab, single quote, double quote, blank, comma, and underscore definitions to XRFconstants.h
//      Ignore trailing spaces and tabs for element symbols
//      Fix energy start display in read_spectrum_file.cpp
//      Quantify bulk sum after summing
//  Modified Jan. 15, 2018
//      Fix headers in plot file for max value spectrum
//  Modified Jan. 17, 2018
//      Change XUNITS and YUNITS in MSA read to be case-insensitive
//      Process eV or keV XUNITS in MSA spectrum files
//  Modified Jan. 24, 2018
//      Add Si K alpha 3 and K alpha 4 satellite lines to escape peaks for Si detectors
//      Add total counts to information written about each spectrum read
//  Modified Jan. 26, 2018
//      Remove +- from uncertainty numbers in calibration file output
//      Add net peak intensity to calibration file output
//      Use average live time when combining multiple-detector spectra
//      Fixed total counts output in read_spectrum_file (was after new line)
//  Modified Jan. 31, 2018
//      Process new ##EVENTS and ##TRIGGERS keywords and make XIA livetime calculation
//      Add checks for less than or equal to zero for NPOINTS (0 OK for config files), XPERCHAN, LIVETIME, TRIGGERS, and EVENTS
//  Modified Mar. 2, 2018
//      Added -q option to specify outputs to map file (restructure quantWriteMap like quantWriteResults)
//      Use number of iterations stored in spectrum
//      Load quantification results into element list in quantWriteResults and use that for cal and map file outputs
//  Modified Mar. 7, 2018
//      Move oxideFormulaString and datetime to XRFutilities and add path/filename separation function
//      Recognize and process a list of spectrum files for map and bulk sum if spectrum file ends in ".txt"
//      Use utility to check file extensions
//      Add message about map file name and options for map command
//  Modified Apr. 11, 2018
//      Change ECF output to 4 decimal places (was zero decimal places)
//      Don't call quantWriteResults after quantUnknown if error
//      Modify detector resolution calculation to avoid sqrt of negative numbers
//  Modified Apr. 18, 2018
//      Increase acceptance for change to Fano factor to 40% to avoid fit problems in some noisy spectra  (quantFitSpectrum)
//  Modified July 18, 2018
//      Fix bug in XrayMaterial that failed to add oxygen when all oxide ratios were zero and one was changed
//      Also add flag to XrayMaterial constructor to accept element fractions as input, not oxide fractions (still handles oxide ratios)
//      Change quantUnknown to call fraction_input instead of fraction_oxide to adjust composition (necessary with above change)
//  Modified July 25, 2018
//      Write out some useful information if calculated intensity is zero or nan in several fp calculation functions
//      Check for zero ECFs from csv calibration file to prevent zero ECFs for any element
//      Relative intensities for K lines changed to use values from Scofield
//      Correct resolution of escape peaks to value at escape energy, not line energy
//  Modified July 27, 2018
//      Change -b command line arguments have any number of entries
//      Change background to 2-zone SNIP
//  Modified July 31, 2018
//      Fix quantWriteCalibrationTXT to use weight from standards input (if CSV standards file is used)
//          (Include same acceptance criteria based on min % and max rsd if old TXT standards file was used)
//      Use spectrum file name as name of standard since no names in TXT standards file
//      Implement linear energy calibration correction from Chris
//  Modified Aug. 1, 2018
//      Change end channel checks for background to avoid accepting small non-zero channel number when source kV is zero
//      Fix various minor problems with text outputs
//  Modified Sep. 19, 2018
//      Correct counts in range 1-7.25 keV (was total spectrum counts)
//      Fix bug in linear energy calibration correction (x-intercept was compared to offset, not energy_in)
//      Allow zero values for LIVETIME, REALTIME, TRIGGERS, and EVENTS (for testing MSA files without collecting a spectrum)
//      Make plot with all-zero measured spectrum work OK in quantWritePlot.cpp
//  Modified May 6, 2018
//      Fix write_EMSA_PIXL to write measured spectrum for bulk sum (not calculated spectrum, pass flag for bulk sum)
//      Add spectrum file name and PIQUANT version header to plot file, first line (configuration file for calculations)
//  Modified May 13, 2019
//      Crude correction to Rh L lines from X-ray tube to better fit Rh L scatter peaks in spectrum (PIXL flight tube)
//      Didn't match measured spectrum of SN44 flight tube, so removed it (must be something in optic)
//      In read_EMSA_PIXL, get rid of spectrum_hold struct, enter in XraySpectrum object when read in
//      In read_EMSA_PIXL, fix XIA live time calculation to include overflows and underflows
//      Get rid of fake RTT, USN, and PMC, and use real values from MSA file
//  Modified May 14, 2019
//      Radical change to quantWriteMap to enable writing all of the possible information, including aux info, using upper and lower case letters
//      Moved defaults for quant map outputs and background parameters to this file, just after allocation of arguments
//  Modified May 16, 2019
//      Add counts in 1-7.25 keV region to quant map outputs, move from quantWriteResults to XraySpectrum
//      Add capability to turn off adjustments to energy calibration and detector resolution in fits (-f and -g options)
//      Remove negative intensities from map and log output (leave coefficients as-is for now)
//      Add -s option to convolve Compton scatter components with a Gaussian (detector resolution)  brute force, very expensive in compute time
//  Modified May 21, 2019
//      Implement Evaluate action, changes to parse_element_list, parse_arguments, and quantWriteMap (plus this file)
//  Modified May 21, 2019
//      Include minor changes from Kyle Uckert to enable compilation on other machines (See e-mail today at 2:57pm)
//  Modified May 22, 2019
//      Fix last bugs in Evaluate (for now at least, output seems OK, use for FM calibration)
//  Modified May 21, 2019       (after Version 2.40)
//      Fix bug - no comma between energy start and energy per channel in map file
//  Modified May 25, 2019
//      Don't include ignore, exclude, or matrix elements in evaluate output list
//      Always put all elements in edge list for excitation energies, to avoid changes in calc int when element list changes (temporary fix)
//      Fix more things in secondary fluorescence, see comments in fpSecondary.cpp
//  Modified June 6, 2019
//      Handle EXCLUDE and MATRIX elements correctly (for standards and evaluate, and in unknowns but no way to enter it yet)
//      Implement negative oxide ratio entry means use default oxide ratio (from XrayMaterial.default_oxide_ratio)
//      Add output of oxide oxygen percent from oxide ratio of matrix elements (quantWriteResults)
//  Modified July 2, 2019
//      Allow for background to be adjusted manually instead of in least-squares fit, using -b option parameter
//  Modified July 3, 2019
//      Provide for use of FP calculated bkg instead of SNIP background (with or without fit)
//      Background is now always present as a component, but is excluded (and not used) if not fit
//      All of the background setup is taken care of in quantBackground (which returns 1 if FP calculated bkg is used)
//      Allow matrix element fraction or percent to be entered in element list (as C_M=23.7%)
//  Modified July 22, 2019
//      Change light element oxide ratios to zero (below Z=10)
//      Treat any elements that have no components in the spectrum as matrix elements, and print warning
//      Don't remove input percent from element list, let quantUnknown and quantWriteResults handle it (after no-line matrix elements added)
//      Change wording for oxygen from matrix element oxides
//  Modified Aug. 2, 2019
//      Prepared Version 2.45 for use at Denver X-ray Conference QA Workshop demo
//  Modified Nov. 4, 2019
//      Change map default options to those needed for generating element maps using beam location file (includes PMC in 1st column)
//      Sort master element list for Evaluate by atomic number
//      Calculate an uncertainty for each ECF using the standard deviation across the set of fit coefficients for the standards
//      Include ECF standard deviation in error output for quantification (err in map files, extra info in cout)
//      Add atomic number to calibration file as an extra column (to aid in sorting by element)
//      Calculate a weighted average certificate uncertainty associated with each ECF and include it in the ecf_error entry
//      Compute weighted average of ECF fit errors and use larger of that or ECF standard deviation as ECF error
//          See document "PIQUANT_Error_propagation_plan_11_04_2019.docx" from Chris Heirwegh Nov. 4, 2019 2:24 PM
//  Modified Nov. 6, 2019
//      Add some keyword revisions to read_XIA_PIXL because XIA changed file format
//  Modified Dec. 11, 2019
//      Make quantComponents and quantDefaults skip IGNORE elements (Prevents bug where IGNORE elements were reclassified as MATRIX)
//  Modified Dec. 16, 2019
//      Change map output options: seq# changed to Q, S added for element sum (to compare to 100%)
//      Add element qualifier O to force the element to be included in the evaluate list (with zeros if not in any standard in this run)
//  Modified Jan. 10, 2019  Check for energy <= zero in fpConvolve and skip that point to avoid nan result
//                          Changed to version 2.5.298 and built for release
//  Modified Jan. 22, 2019  Cosmetic changes to prepare to send to JPL GitHub
//  Modified Oct. 19, 2020  Include detector tails, Compton escape, and separate Rh L lines from X-ray tube
//  Modified Oct. 21, 2020  Add multiple background components (split for independent fitting) and include Compton escape as separate component
//  Modified Oct. 22, 2020  Add optic response sub-command and calculation framework (quantOpticResponse.cpp, based on quantStandard.cpp)
//  Modified Oct. 30, 2020  Changed XrayOptic to extrapolate below energy arrays to lowest transmission value (flat extrapolation)
//  Modified Nov. 2, 2020   Add PIXL FM optic as input from include file, with spline fit
//                          Select standard from standards file using -s option (by number or name)  [change -s for Compton convolve to -v]
//  Modified Nov. 5, 2020   Option -w to write evaluation file during calibrate (to plot ECFs, etc. during calibration development)
//  Modified Nov. 9, 2020   Change output to oxides instead of elements (debug problem with given % if entered as element)
//  Modified Nov. 30, 2020  Exclude from fit vector any components that have fit Boolean set to false
//                          Also add factor to compute coefficients of non-fit components from components used for quant (XRFcontrols.h)
//  Modified Dec. 4, 2020   Use SNIP background to fit anomalous background at low energies
//  Modified Dec. 7, 2020   Choose default primary spectrum number of channels to go past max energy for monochromatic source
//                          Added detector shelf calculation from photoelectron and Auger electron escape (active volume and front contact)
//  Modified Dec. 28, 2020  Fix bug where composition made big change and coefficients had big mis-match, update_calc gave bad values for detector shelf calcs
//                          Changed default bkg arguments for optic response from calc bkg to SNIP 2-zone
//                          Add optic response to plot output of primary spectrum calculation
//                          Improvements to optic response function, results into XrayOptic as optic # 7
//                          Use SNIP for low-energy background, disable shelf calculations
//  Modified Dec. 31, 2020  Rework standard selection for multi-standard commands, add plot file as command option for Evaluate
//  Modified Jan. 4, 2021   Include residual error from component into fit error for coefficient, re-work ECF uncertainty calculations
//                          Interpolate ECFs vs Z for missing ECFs
//  Modified Jan. 5, 2021   Disable detector tails (in XRFcontrols.h), minor changes to improve speed
//  Modified Jan. 7, 2021   Fix bug with formula_fraction in XrayMaterial
//                          Turn off sec fluor until quant accuracy gets better, for speed doing maps
//  Modified Jan. 13, 2021  Various bugs fixed, found with unit test data set on GitLab
//  Modified Feb. 2, 2021   Fix bug reading standards files with only commas on line, between standards
//                          Change bkg back to SNIP, with fit below ~6 keV and not fit above that (selectable via -a option)
//                          Avoid evaluating a standard with itself as a calibration standard
//  Modified Feb. 3, 2021   Lowered minimum weight to include in Evaluate calculations to 0.15 (in XRFcontrols.h)
//  Modified Feb. 24, 2021  Add "U" option to map out options to include spectrum aux info title (add standard names here in calibrate and evaluate)
//  Modified Feb. 26, 2021  Adjust the shape of the calculated background to generally match measured spectra using ramp multiplier
//                          Adjust the overall intensity to match each measured spectrum using an overall factor calculated via least squares with peaks blocked
//  Modified Mar. 8, 2021   Change tail Czero to be 0.75 below Si K edge and 0.92 above
//  Modified Mar. 9, 2021   Don't double-count any background components, exploring values for detector shelf and tail using pure compounds
//  Modified Mar. 12, 2021  Fix bug interpolating ECF when only one entry in interpolation list
//  Modified Mar. 15, 2021  Adjust coefficients for all components, not just non-fit (necessary for shelf calculation)
//                          Front contact shelf disabled to improve Na results (also use ECF of 1.3 for Na)
//  Modified Mar. 18, 2021  Compton escape disabled
//      Scholze & Procop shelf model, with adjustments vs photon energy and slope vs electron loss energy (to avoid clobbering Na peak)
//  Modified Mar. 18, 2021  Shelf overall factor of 30 and slope vs overall electron loss energy (see quantCalculate.cpp)
//  Modified Apr. 2, 2021   Major rearrangements in fpLineSpectrum to improve speed, move shelf calc from quantCalculate to fpLineSpectrum
//                          Remove bkg spline adjustment, tweak optic response to remove structure at 6.5 & 8.5 keV
//  Modified Apr. 6, 2021   Major revisions to components for background, now control bkg defaults in quantBackground
//  Modified Apr. 8, 2021   Ignore zero-energy components for done check (not an element component with a reasonable peak at a known energy)
//  Modified Apr. 11, 2021  Fix quantBckground to handle plot cmd properly (just use SNIP, don't insert components and update calc)
//  Modified Apr. 28, 2021  Add La and Lb1 component types to checkComponent (for extra Rh L lines, especially Lb1 fit separately)
//                          Modify adjustment factors for Rh L lines in XraySource (tubeLinesSewell) to better match BHVO spectrum from PIXL FM Elemental Calibration
//                          Remove option to have -b,0 perform SNIP bkg with default parameters (must now specify parameters to use SNIP)
//                          Add back crossover, but using SNIP at high energies and adj calc bkg at low wnergies (via option -a)
//                          Change default to this new crossover, since it works much better for trace elements (crossover = 7000 eV)
//                          Adjust SHELF_SLOPE to 1.45 (Na in SRM 694 is most sensitive to this)
//                          New error values after modifications made today, significant improvements
//  Modified May 1, 2021    Optic response tweaked by hand to get Sr thru Zr ECFs to unity and for better fit to Teflon bkg above 25 keV
//                          Tweaked 6 keV value to get Ca closer and better bkg under Cr      Errors updated
//  Modified May 10, 2021   Add -bh and -bx background options, eliminate -a option, implement more controls in quantBackground
//  Modified May 14, 2021   Fixed bug in scale_under_peaks with sigma_multiplier not being used
//                          Fixed logic error in quantCalculate where coefficient was reset to unity when it should be manual scale factor
//                          Add -T option to control detector shelf factor and slope (was fixed in fpLineSpectrum, moved to XrayDetector)
//  Modified May 25, 2021   Added symbol to grouped lines, for identification during debugging
//  Modified June 8, 2021   Change error reporting to absolute (was relative)  quantWriteResults.cpp
//  Modified June 9, 2021   Fix bug in energy per channel calculation that was disturbing convolution normalization  XraySpectrum.h
//                          Fix bug where absolute error was multiplied by 100 in map and eval files  quantWriteResults.h
//                          Add average geometry factor to bulk sum msa files
//                          Change Z range for mid-Z trace elements to eliminate Co (too much interference from Fe beta peak)
//                          Change interpolation range to better match actual error performance    quantWriteResults.cpp
//                          Change oxides and element reporting to match team wishes (e-mail from Joel 6/7/2021, 1:27 PM)
//                          Modify ECF interpolation vs Z so it doesn't extrapolate outside of available values (takes end value)
//  Modified June 10, 2021  Use single SNIP bkg for default in PLOT cmd without residual, etc.
//                          Include final decision on background defaults
//  Modified June 17, 2021  Peter Nemere found and fixed bug in detector shelf calculation - uninitialized variable
//  Modified June 19, 2021  Fixing above bug had repercussions for calibration, lots of flailing around to check shelf factor, slope, & tail parameters
//                          Final choice was shelf factor 1 & slope 0, tail unchanged, and front contact shelf enabled with thickness 150 um
//  Modified June 27, 2021  Add command line option to normalize element sum to 100% (or any value)
//  Modified July 9, 2021   Add command line option to change Fe oxide ratio (-Fe)
//  Modified July 10, 2021  Add simple pulse pileup calculation


//  Remaining FP anomalies as of June 2021
//      Rh L line intensities from X-ray tube (XraySource)   adjusted from database values
//          *** Be window thickness was incorrect when tube spectrum evaluated, so factors were adjusted a bit to better match BHVO ***
//      Optic response zero-energy value 2.3x first low energy value (to get ECFs to unity for Na thru Cl)
//      Calculated background from continuum scatter must be adjusted for each spectrum (overall factor, quantCalculate)
//      Detector tail Z-zero adjusted above and below Si edge (XrayDetector)
//      Rh L intensities in standard spectra vary a lot (diffraction?)
//      Compton escape disabled to improve Na results (quantCalculate)
//      Scholze & Procop shelf model, shape still not quite right, need more data on actual detector response
//          L line edge jumps not corrected for K edge jump when energy is above K edge (everything but K edges disabled for now)

using namespace std;


string getVersionString()
{
    ostringstream os;
    os << Piquant_VERSION_MAJOR << "." << Piquant_VERSION_MINOR << "." << Piquant_VERSION_PATCH << "-" << Piquant_VERSION_BRANCH;
    return os.str();
}

#define HEADER_STRING_1 "PIQUANT   Quantitative X-ray Fluorescence Analysis"
#define HEADER_STRING_2 "Written for PIXL, the Planetary Instrument for X-ray Lithochemistry"
#define HEADER_STRING_3 "   W. T. Elam   APL/UW"



int main (const int argc, const char * argv[])
{
try {
    installSegHandler();

	int result = 0;
	bool error = false;
	//  Get the start time of the process in clock ticks (execution time only)
	clock_t startTime, finishTime;
	startTime = clock();

//      Parse arguments, enum for sub-commands is defined in parse_arguments.h
    PIQUANT_SUBCOMMAND cmd;
    ARGUMENT_LIST arguments;

    //  Set default background arguments in quantbBckground

    result = parse_arguments( argc, argv, cmd, arguments );
    if( result < 0 && result >= -2020 ) return result;  //  catastrophic problem with argument list

    if( cmd == CALIBRATE || cmd == EVALUATE ) {
        //  Set default map output options for evaluate ( file name, given value, quant percent, quant error, and error relative to given)
        if( arguments.quant_map_outputs.length() == 0 ) arguments.quant_map_outputs = "GPEHKLF";
    }
    if( arguments.quant_map_outputs.length() == 0 ) {
        //  Change default map output options to match what is required for generating map files from beam location
        arguments.quant_map_outputs = "pPIETVXCRNFetsr";
        //  Set default quant map outputs (x, y, and z locations plus percents, original default outputs)
        //   OLD   arguments.quant_map_outputs = "xyzP";
    }
    if( arguments.iron_oxide_ratio >= 0 ) XrayMaterial::default_iron_oxide_ratio( arguments.iron_oxide_ratio );

    // If we're asked for the version, just print it and return
    if(cmd == PRINT_VERSION) {
        cout << getVersionString() << endl;
        return 0;
    }

    auto tm = time_code("PIQUANT");

    //  Redirect standard output to a file given in the arguments list
    bool cout_redirect = false;
    ofstream terminal_out_stream;
    streambuf *coutbuf = 0;    //  to save cout rdbuf to undo the redirect before exiting

    if( arguments.terminal_text_file.length() > 0 ) {
        terminal_out_stream.open( arguments.terminal_text_file.c_str(), ios::out | ios::app );
        if( ! terminal_out_stream ) {
            cout << "Can't open terminal output file, name is " << arguments.terminal_text_file << endl;
            error = true;
        } else {
            //  Redirect cout to terminal_out_stream
            coutbuf = std::cout.rdbuf(); //save old buf
            cout.rdbuf(terminal_out_stream.rdbuf()); //redirect std::cout
            cout_redirect = true;
        }
    }
    //  If problem is just invalid options, try to provide a message via the GUI
    if( !error && result < -2020 ) {
        if( arguments.invalid_arguments.length() > 0 ) {
            cout << "**** Invalid arguments: " << arguments.invalid_arguments << endl;
            cout << "No calculations performed." << endl;
        }
        return result;
    }
    if( error ) return -2020;
    ostream &termOutFile = cout;

    bool oxidesOutput = true;
//    bool oxidesOutput = false;

    //		write program header
	termOutFile << "-----------------------------------------------------------------" << endl;
	termOutFile << HEADER_STRING_1 << endl;
	termOutFile << HEADER_STRING_2 << endl;
	termOutFile << getVersionString() << HEADER_STRING_3 << endl;

	termOutFile << "Local time:  " << datetime() << endl;
	termOutFile << endl;

	termOutFile.setf( ios::fixed, ios::floatfield );
	termOutFile.precision(2);

    if( !arguments.fit_adjust_energy ) cout << "Adjustment of energy calibration during fits is disabled." << endl;
    if( !arguments.fit_adjust_width ) cout << "Adjustment of peak widths during fits is disabled." << endl;


    //**************************************************************************
    //      Read the configuration file
    //**************************************************************************

    //		storage to hold instrument parameters
    XRFconditionsInput condStruct_config;
    XRFconditions configConditions;    //  default instrument measurement conditions
    XraySpectrum configSpectrum;
    vector <float> bkg_split_energies;
    string configurationFileName( arguments.configuration_file );
    if( ( cmd == PRIMARY || cmd == CALCULATE || cmd == CALIBRATE || cmd == QUANTIFY || cmd == EVALUATE || cmd == MAP || cmd == COMPARE || cmd == FIT_ONE_STANDARD || cmd == OPTIC_RESPONSE  || cmd == BULK_SUM_MAX ) && ( ! error ) ) {
        if( configurationFileName.length() > 0 ) {
            if( check_file_extension( configurationFileName, "XSP" ) ) {
                // Open and read the borehole XRF file format (older version of EMSA format)
                float ev_start_cfg = 0;
                float ev_ch_cfg = 0;
                float live_time_cfg = 100;  //  For calculate, in case none in configuration file
                float x = 0, y = 0, z = 0;
                vector <float> spectrum_cfg;
                vector <string> spectrum_titles;
                //  Optic type is a number between 0 and 3
                result = borehole_read ( configurationFileName, condStruct_config.conditionsVector, spectrum_cfg, ev_start_cfg, ev_ch_cfg, live_time_cfg, spectrum_titles, x, y, z );
                if ( result != 0 ) {
                    termOutFile << "Can't read xsp configuration file, result = " << result << "  for file name " << configurationFileName << endl;
                    error = true;
                } else {
                    XraySpectrum temp_spec( spectrum_cfg, ev_start_cfg, ev_ch_cfg );
                    temp_spec.live_time( live_time_cfg );
                    configSpectrum = temp_spec;
                }
            } else if ( check_file_extension( configurationFileName, "MSA" ) ) {
                //      open and read the ISO 22029 2012 EMSA format file
                vector <XraySpectrum> spectrum_vec;
                string acq_date, acq_time, x_label, y_label, unit;
                //  Minimum energy not handled yet, assume 900 eV if this type of file is read
                result = read_EMSA_PIXL( configurationFileName, condStruct_config, spectrum_vec );
                if ( result != 0 ) {
                    termOutFile << "Can't read msa configuration file, result = " << result << "  for file name " << configurationFileName << endl;
                    if( result == -999999 ) {
                        termOutFile << "Invalid file format or missing required keyword." << endl;
                    } else {
                        termOutFile << "Error on line number = " << -result << "." << endl;
                    }
                    error = true;
                };
                if( spectrum_vec.size() > 0 ) configSpectrum = spectrum_vec[0];
            } else {
                termOutFile << "Can't read configuration file, unrecognized file type, for file name " << configurationFileName << endl;
                error = true;
            }
            if( ! error ) termOutFile << "Configuration read from file " << configurationFileName << endl;
        } else {
                termOutFile << "A configuration file is required for this sub-command." << endl;
                error = true;
        }
        //  Set up new instrument measurement conditions from configuration file
        if( ! error ) {
            result = fpSetupConditions ( condStruct_config, configConditions );
            if( result < 0 ) {
                termOutFile << "fpSetupConditions failed, result " << result << endl;
                termOutFile << "Error in parameter " << get_EMSA_keyword( -(result+100) ) << endl;
                error = true;
            }
        };
        //  Write the configuration information that will be used for calculations
//        result = write_conditions ( condStruct_config );
    }   //  Read the configuration file
    termOutFile << endl;

    //  Move detector shelf adjustment parameters from -T option (if any) to conditions vector
    if( arguments.detector_shelf_parameters.size() > 0 ) {
        condStruct_config.conditionsVector[DETECTOR_SHELF_FACTOR_INDEX] = arguments.detector_shelf_parameters[0];
        if( arguments.detector_shelf_parameters.size() > 1 )
            condStruct_config.conditionsVector[DETECTOR_SHELF_SLOPE_INDEX] = arguments.detector_shelf_parameters[1];
        if( arguments.detector_shelf_parameters.size() > 2 )
            condStruct_config.conditionsVector[DETECTOR_SHELF_SLOPE_START_INDEX] = arguments.detector_shelf_parameters[2];
     }

    //      Spectrum for use in all commands, put here to make available for plots
    XraySpectrum singleSpectrum = configSpectrum;
    //  Set the configuration conditions as the default
    XRFconditionsInput condStruct_spec;
    copy_conditions_struct( condStruct_config, condStruct_spec );

   //************************************************************************************
    //		Calculate the primary spectrum as if it were going straight into the detector
    //***********************************************************************************

    if( ( cmd == PRIMARY ) && ( ! error ) ) {
        // If an energy calibration was given in the argument list, it overrides
        if( arguments.eV_ch > 0 ) {
            singleSpectrum.calibration( arguments.eV_start, arguments.eV_ch );
            termOutFile << "Using energy calibration from option argument  ";
            termOutFile.precision(1);
            termOutFile << "  eV start = " << singleSpectrum.calibration().energyStart();
            termOutFile.precision(4);
            termOutFile << "  eV/ch = " << singleSpectrum.calibration().energyPerChannel();
            termOutFile << endl;
        }
        if( ! singleSpectrum.calibration().good() ) {
            termOutFile << "Bad energy calibration, can't calculate primary spectrum." << endl;
            error = true;
        } else {
            //  See if there is a spectrum in the configuration file, if so use it as the measured spectrum to compare to
            if( singleSpectrum.numberOfChannels() <= 0 ) {
                //  There is not a spectrum in the configuration file, so add one that is all zeros
                int nc = ( configConditions.source.kV() * 1000 - singleSpectrum.calibration().energyStart() )
                        / singleSpectrum.calibration().energyPerChannel();
                //  If the source is not an X-ray tube, go past source energy to be sure main peak is in plot
                if( ! configConditions.source.continuum() ) nc += nc/10;
                vector <float> temp_meas( nc, 0 );
                singleSpectrum.meas( temp_meas );
            }
            //  Put the configuration file name into the spectrum object, so it will go into the output file
            singleSpectrum.file_name( configurationFileName );
            //  Put the current date and time into the spectrum object
            singleSpectrum.aux_info_change().date = datetime().substr( 0, 11 );
            singleSpectrum.aux_info_change().time = datetime().substr( 12, 9 );
            termOutFile << "Calculating primary spectrum, live time " << singleSpectrum.live_time();
            termOutFile << ",   energy calibration (eV): ";
            termOutFile.precision(1);
            termOutFile << "  eV start = " << singleSpectrum.calibration().energyStart();
            termOutFile.precision(4);
            termOutFile << "  eV/ch = " << singleSpectrum.calibration().energyPerChannel();
            termOutFile << endl;
            result = quantPrimarySpec( configConditions, singleSpectrum );
            if( result < 0 ) {
                termOutFile << "quantPrimarySpec failed, result " << result << endl;
                if( result == -701 ) termOutFile << "Invalid number of channels: " << singleSpectrum.numberOfChannels() << endl;
                if( result == -705 ) termOutFile << "Invalid energy calibration." << endl;
                if( result == -706 ) termOutFile << "Invalid live time: " << singleSpectrum.live_time() << endl;
                error = true;
            };
        }
    }   //  if( ( cmd == PRIMARY ) && ( ! error ) ) {



    //***************************************************************************
    //		Setup list of standards with their XRF spectra and given compositions
    //***************************************************************************

    vector <StandardInformation> standards;
    int standard_index = -1;    //  Selection of standard from list in standards input file (from option -s)
    if( ( cmd == CALCULATE || cmd == CALIBRATE || cmd == EVALUATE || cmd == COMPARE || cmd == FIT_ONE_STANDARD || cmd == OPTIC_RESPONSE || cmd == EVALUATE ) && ( ! error ) ) {
        if( check_file_extension( arguments.standards_file, "TXT" ) ) {
            result = setupStandardsTXT( arguments.standards_file, termOutFile,
                    standards, MINIMUM_AMOUNT );
        } else if( check_file_extension( arguments.standards_file, "CSV" ) ) {
            result = setupStandardsCSV( arguments.standards_file, standards, MINIMUM_AMOUNT );
        } else {
            termOutFile << "Standards input files can only be .txt or .csv" << endl;
            result = -1;
        }
        if ( result != 0 ) {
            termOutFile << "Standards file read failed, result = " << result << endl;
            error = true;
        } else {
            termOutFile << "Standards file read OK, entries for " << standards.size() << " standards read in." << endl;
        };
        if( ( cmd == CALCULATE || cmd == COMPARE || cmd == FIT_ONE_STANDARD || cmd == OPTIC_RESPONSE || cmd == EVALUATE ) && ( ! error ) ) {
            //  Set up the selection from the list of standards
            if( standards.size() <= 0 ) {
                termOutFile << "No standards input, can't perform this action." << endl;
                error = true;
            } else if( arguments.standard_selected ) {
                standard_index = arguments.standard_selection;
                if( arguments.standard_name.length() > 0 ) {
                    unsigned int is;
                    for( is=0; is<standards.size(); is++ ) {
                        standard_index = -1;
                        unsigned int i_name;
                        for( i_name=0; i_name<standards[is].names.size(); i_name++ ) {
                            if( standards[is].names[i_name] == arguments.standard_name ) {
                                standard_index = is;
                                break;
                            }
                        }
                        if( standard_index >= 0 ) break;
                    }
                }
                if( standard_index < 0 || standard_index > standards.size()-1 ) {
                    termOutFile << "invalid standard selection: ";
                    if( arguments.standard_name.length() > 0 ) termOutFile << arguments.standard_name;
                    else termOutFile << arguments.standard_selection;
                    termOutFile << endl;
                    error = true;
                } else {
                    termOutFile << "Standard selected: ";
                    if( standards[standard_index].names.size() > 0 ) termOutFile << standards[standard_index].names[0];
                    termOutFile << "   (# " << standard_index << ")." << endl;
                }
            }
        }
        termOutFile << endl;
    }

    //  Default to first standard in standards list
    if( standard_index < 0 ) standard_index = 0;
    // Check that we have a valid standard_index
    if( (cmd == CALCULATE || cmd == COMPARE || cmd == OPTIC_RESPONSE) && (standard_index < 0 || standard_index >= standards.size()) ) {
        termOutFile << "No standard selected!" << endl;
        error = true;
    }

    //**************************************************************************
    //  Get element list from arguments
    //**************************************************************************

    vector <ElementListEntry> element_list;
    bool element_list_carbonates = arguments.carbonates;
    if( ( cmd == ENERGY_CAL || cmd == CALIBRATE || cmd == EVALUATE || cmd == QUANTIFY || cmd == MAP || cmd == FIT_ONE_STANDARD || cmd == OPTIC_RESPONSE || cmd == BULK_SUM_MAX ) && ( ! error ) ) {
        termOutFile << "Element list: " << arguments.element_list << endl;
        termOutFile << endl;
        error = parse_element_list( arguments.element_list, element_list, element_list_carbonates );
    }

    vector <XraySpectrum> spectrum_vec; //  Save this here so we can plot multiple detectors for PLOT sub-command

    //**************************************************************************
    //      Read a single spectrum and check it's energy calibration
    //**************************************************************************

    if( ( cmd == ENERGY_CAL || cmd == PLOT || cmd == QUANTIFY || cmd == COMPARE || cmd == OPTIC_RESPONSE ) && ( ! error ) ) {
        string spectrumPathName( arguments.spectrum_file );
        //      open and read the spectrum file
        spectrum_vec.clear();
        result = read_spectrum_file( termOutFile, spectrumPathName, spectrum_vec, condStruct_spec );
        if ( result != 0 ) {
            //  Error message already sent in read_spectrum_file
            error = true;
        };
        if( cmd != ENERGY_CAL ) {
            //  Set up energy calibration, background parameters, and measurement conditions
            setup_spectrum_parameters( arguments, configSpectrum.calibration(), spectrum_vec,
                    condStruct_config, condStruct_spec, termOutFile );
        }
        //  Combine the spectrum information from several detectors (or the selected detector) into the variable where they will be used
//        if( spectrum_vec.size() > 0 ) singleSpectrum = spectrum_vec[0];
        //      NB: quantCombineSpectra modifies the spectra in the input list to match them to a single energy axis
        //          for proper plotting
        result = quantCombineSpectra( spectrum_vec, singleSpectrum, arguments.detector_select );
        if( result < 0 && ! ( result = -3 && cmd == PLOT ) ) error = true;  //  Ignore can't combine error if only plotting
    }   //  if( cmd == ENERGY_CAL || cmd == PLOT || cmd == QUANTIFY )
    if( ( cmd == PLOT || cmd == QUANTIFY || cmd == COMPARE ) && ( ! error ) ) {
        //  Check for bad energy calibration and set to channels if just for plot
        if( ! singleSpectrum.calibration().good() ) {
            termOutFile << "*** Warning energy per channel is bad: " << singleSpectrum.calibration().energyPerChannel() << endl;
            if( cmd != PLOT ) error = true; //  Plot can be vs channels, all others are not possible without calibration
        }
        termOutFile << endl;
    }



    //*************************************************************************************
    //      Perform a calculation of the spectrum from the chosen standard
    //*************************************************************************************
    FPstorage fpStorageST; // If running single threaded, all the stuff called from here can just use this FPstorage

    if( ( cmd == CALCULATE || cmd == COMPARE ) && ( ! error ) ) {
        //  Put in the energy calibration from arguments for calculate, since we didn't read in a spectrum
        if( cmd == CALCULATE ) {
            if( arguments.eV_ch > 0 ) {
                singleSpectrum.calibration( arguments.eV_start, arguments.eV_ch );
                cout << "Using energy calibration from option argument  ";
                cout.precision(1);
                cout << "  eV start = " << singleSpectrum.calibration().energyStart();
                cout.precision(4);
                cout << "  eV/ch = " << singleSpectrum.calibration().energyPerChannel();
                cout << endl;
            }
            //  Put the configuration file name into the spectrum object, so it will go into the output file
            singleSpectrum.file_name( configurationFileName );
        }
        if( ! singleSpectrum.calibration().good() ) {
            termOutFile << "Bad energy calibration, can't calculate spectrum." << endl;
            error = true;
        } else {
            //  See if there is a spectrum in the configuration file, if so use it as the measured spectrum to compare to
            if( singleSpectrum.numberOfChannels() <= 0 ) {
                //  There is not a spectrum in the configuration file, so add one that is all zeros
                int nc = ( configConditions.source.kV() * 1000 - singleSpectrum.calibration().energyStart() )
                        / singleSpectrum.calibration().energyPerChannel();
                vector <float> temp_meas( nc, 0 );
                singleSpectrum.meas( temp_meas );
            }
//                standards[standard_index].mat.thickness( 1 * CM_MICRON );  //  **********   for calculation of metal films for cal target
//                standards[standard_index].mat.thickness( specimen_thickness * CM_MICRON );  //  **********   for calculation of cal target pucks vs thickness
                //  Set up components for the calculated spectrum
                vector <SpectrumComponent> components;
                vector <XrayLines> sourceLines;
                //  Load vector with emission lines from X-ray source
                configConditions.source.lines( sourceLines, configConditions.eMin );
                vector <XrayLines> pureLines;
                //  Load vector with pure element emission lines from specimen and set up FP calculations
                fpPrep(fpStorageST, standards[standard_index].mat, configConditions, pureLines );
                result = setupComponents( sourceLines, pureLines, components );
                if( result < 0 ) {
                    cout << "setupComponents failed, result is " << result << endl;
                    error = true;
                }
                //  Add components for detector tail, shelf, and Compton escape
                vector <XrayLines> dummy_lines;
                result = makeComponents( CONTINUUM, dummy_lines, components );
                if( Compton_escape_enable_flag ) result = makeComponents( DETECTOR_CE, dummy_lines, components, 0 );
                //  Add the components to the spectrum object
                int ic;
                for( ic=0; ic<components.size(); ic++ ) singleSpectrum.add_component( components[ic] );
                //  Setup convolution of Compton components (now brute force, so very expensive in compute time)
                singleSpectrum.convolve_Compton( arguments.convolve_Compton );
                // Calculate spectrum for this standard
                //standards[standard_index].mat.density( 2.2f );
                //standards[standard_index].mat.thickness( 0.2f );
                result = quantCalculate(fpStorageST, standards[standard_index].mat, configConditions, singleSpectrum );
                if ( result != 0 ) {
                    cout << "quantCalculate failed, result = " << result << endl;
                    if( result == -701 ) termOutFile << "Invalid number of channels: " << singleSpectrum.numberOfChannels() << endl;
                    if( result == -705 ) termOutFile << "Invalid energy calibration." << endl;
                    if( result == -706 ) termOutFile << "Invalid live time: " << singleSpectrum.live_time() << endl;
                    error = true;
                };
                //  Insert optic response (or primary spectrum) as a component to get it into the plot file (not added to calculation)
/*                SpectrumComponent put_in_optic;
                put_in_optic.type = OPTIC_TRANS;
                put_in_optic.enabled = false;   //  Keep out of calc
                put_in_optic.fit = false;
                put_in_optic.spectrum.resize( singleSpectrum.numberOfChannels() );
                int i;
                //for( i=0; i<put_in_optic.spectrum.size(); i++ ) put_in_optic.spectrum[i] = configConditions.optic.CheckTransmission( singleSpectrum.energy(i) );
                //  Change to include primary spectrum (with everything in beam path, but not detector effects) -- no emission lines
                put_in_optic.type = PRIMARY_CONTINUUM;
                for( i=0; i<put_in_optic.spectrum.size(); i++ ) {
                    vector <float> contEn(1);
                    vector <float> contInt(1);
                    contEn[0] = singleSpectrum.energy(i);
                    if( contEn[0] <= 0 ) contEn[0] = 1;
                    contInt[0] = configConditions.source.continuum( contEn[0] );
                    fpIncidentBeam ( configConditions, contEn, contInt );
                    put_in_optic.spectrum[i] = contInt[0];
                }
                singleSpectrum.add_component( put_in_optic );
            }   */
        float element_sum = 0;
        vector <ElementListEntry> element_list_std;
        if( ! error ) result = quantWriteResults( standards[standard_index].mat, configConditions.detector, element_list_std,
                                singleSpectrum, oxidesOutput, termOutFile, element_sum );
        }
        termOutFile << endl;
    }   //  if( ( cmd == CALCULATE ) && ( ! error ) )


    //*************************************************************************************
    //      Calculate or adjust the optic response using a blank standard and its spectrum
    //*************************************************************************************

    if( cmd == OPTIC_RESPONSE && ( ! error ) ) {
         if( singleSpectrum.live_time() <= 0 ) {
            termOutFile << "*** Error - live time is bad, can't use this standard for calibration. ***" << endl;
            error = true; //  Plot can be vs channels, all others are not possible without calibration
        }
        //  Set up new instrument measurement conditions
        XRFconditions calConditions;
        result = fpSetupConditions ( condStruct_spec, calConditions );
        if( result < 0 ) {
            termOutFile << "fpSetupConditions failed, result " << result;
            termOutFile << "   error in parameter with keyword " << get_EMSA_keyword( -(result+100) ) << endl;
            error = true;
            return -500 + result;
        };
        termOutFile << endl;
        result = quantOpticResponse( fpStorageST, standards[standard_index].mat, element_list, calConditions, singleSpectrum );
        if ( result < 0 ) {
            cout << "quantOpticResponse failed, result = " << result << endl;
            error = true;
        }
        vector <ElementListEntry> element_list_std;
        float element_sum = 0;
        if( ! error ) result = quantWriteResults( standards[standard_index].mat, calConditions.detector, element_list_std,
                                singleSpectrum, oxidesOutput, termOutFile, element_sum );
    }   //  if( cmd == OPTIC_RESPONSE && ( ! error ) )


    //*************************************************************************************
    //      Create a master element list for Calibrate and Evaluate commands
    //*************************************************************************************

    //  Generate a master list of elements to be written to the map file (with default entries)
    vector <ElementListEntry> element_list_eval_master;
    //  For debugging and evaluation of the calibration, write a synopsis file using the evaluate format
    //      Calibrate does not normally write this, see line 362 above, which sets empty map file name to a fixed value
    if( ( cmd == CALIBRATE || cmd == EVALUATE ) && ( ! error ) ) {
        //  First include the elements from the input list that have the "O" qualifier
        unsigned int i_list;
        for( i_list=0; i_list<element_list.size(); i_list++ ) {
            if( element_list[i_list].qualifier == OUTPUT ) {
                ElementListEntry temp = element_list[i_list];
                temp.percent = 0;   //  Use zero for percent of output-only elements instead of default -1
                add_element_list_entry( temp, element_list_eval_master, true );  //  Ignore qualifier for matches
            }
        }
        //  Now include all elements that appear in any of the standards
        unsigned int istd;
        for( istd=0; istd<standards.size(); istd++ ) {
            for( i_list=0; i_list<standards[istd].element_list.size(); i_list++ ) {
                if( standards[istd].element_list[i_list].qualifier == IGNORE
                    || standards[istd].element_list[i_list].qualifier == EXCLUDE
                    || standards[istd].element_list[i_list].qualifier == MATRIX ) continue;
                //  Only include elements that have weight greater then minimum (set in XRFcontrols.h)
                if( standards[istd].element_list[i_list].weight <= arguments.min_wgt_eval ) continue;
                //  Add only the element to master element list entry
                ElementListEntry element_list_temp;
                element_list_temp.element = standards[istd].element_list[i_list].element;
                //  add_element_list_entry will avoid duplicates
                add_element_list_entry( element_list_temp, element_list_eval_master, true );  //  Ignore qualifier for matches
            }
        }
        //  Sort the Evaluate master list into atomic number order (exchange sort)
        //      (only if there were no elements in the input list, if there were retain its order)
        unsigned int iel;
        if( element_list.size() <= 0 && element_list_eval_master.size() > 0 ) {
            for( iel=0; iel<element_list_eval_master.size()-1; iel++ ) {
                unsigned int iel_sort;
                for( iel_sort=iel+1; iel_sort<element_list_eval_master.size(); iel_sort++ ) {
                    if( element_list_eval_master[iel].element.Z() > element_list_eval_master[iel_sort].element.Z() ) {
                        ElementListEntry temp = element_list_eval_master[iel];
                        element_list_eval_master[iel] = element_list_eval_master[iel_sort];
                        element_list_eval_master[iel_sort] = temp;
                    }
                }
            }
        }
    }



    //*************************************************************************************
    //      Perform a calibration by processing the spectrum from each standard in the list
    //*************************************************************************************

    if( ( cmd == CALIBRATE || cmd == FIT_ONE_STANDARD ) && ( ! error ) ) {
        // The stream to the map file we're writing
        std::ofstream fout;
        bool wroteMapHeader = true;
        if( arguments.cal_eval_file.length() > 0 ) {
            fout.open(arguments.cal_eval_file);
            wroteMapHeader = false;
        }
        //  Loop over all of the standards in the list
        int istd;
        for( istd=0; istd<standards.size(); istd++ ) {
            if( ( cmd == FIT_ONE_STANDARD ) && istd != standard_index ) continue;    //  only fit one standard so we can plot it and check fit
            //  Changed to only one spectrum file for each entry in the list of standards
            spectrum_vec.clear();
            string cal_optic_file;
			XRFconditionsInput condStruct_cal;
			if( standards[istd].spectrumFileName.length() == 0 ) {
                termOutFile << "File name missing for standard number " << istd;
                if( standards[istd].names.size() > 0 ) termOutFile << "   name " << standards[istd].names[0];
                termOutFile << endl;
                error = true;
			} else {
                result = read_spectrum_file( termOutFile, standards[istd].spectrumFileName,
                        spectrum_vec, condStruct_cal );
                if ( result != 0 ) {
                    termOutFile << "read_spectrum_file failed, result = " << result << "   file " << standards[istd].spectrumFileName << endl;
                    error = true;
                    continue;
                };
            }
            //  Set up energy calibration, background parameters, and measurement conditions
            setup_spectrum_parameters( arguments, configSpectrum.calibration(), spectrum_vec,
                    condStruct_config, condStruct_cal, cout );
//                if( spectrum_vec.size() > 0 ) singleSpectrum = spectrum_vec[0];
            //  Combine the spectrum information from several detectors (or the selected detector) into the variable where they will be used
            //      NB: quantCombineSpectra modifies the spectra in the input list to match them to a single energy axis
            //          for proper plotting
            result = quantCombineSpectra( spectrum_vec, singleSpectrum, arguments.detector_select );
            if( result < 0 ) {
                error = true;
                break;
            };
            if( singleSpectrum.live_time() <= 0 ) {
                termOutFile << "*** Error - live time is bad, can't use this standard for calibration. ***" << endl;
                error = true; //  Plot can be vs channels, all others are not possible without calibration
                continue;
            }
            //  Add the standard name to the spectrum aux info titles (for output in map file if requested using "U")
            if( standards[istd].names.size() > 0 ) {
                singleSpectrum.aux_info_change().titles.clear();
                singleSpectrum.aux_info_change().titles.push_back( standards[istd].names[0] );
            }
            //  Set up new instrument measurement conditions
            XRFconditions calConditions;
            result = fpSetupConditions ( condStruct_cal, calConditions );
            if( result < 0 ) {
                termOutFile << "fpSetupConditions failed, result " << result;
                termOutFile << "   error in parameter with keyword " << get_EMSA_keyword( -(result+100) ) << endl;
                error = true;
                return -500 + result;
            };
            termOutFile << endl;
            //  Input element fit control list overrides the one from standards input
            vector <ElementListEntry> element_list_std;
            int i_list;
            //  Move element list from standards input file
            for( i_list=0; i_list<standards[istd].element_list.size(); i_list++ )
                add_element_list_entry( standards[istd].element_list[i_list], element_list_std );
            //  Append or replace entries using element list from arguments
            for( i_list=0; i_list<element_list.size(); i_list++ )
                add_element_list_entry( element_list[i_list], element_list_std );
            result = quantStandard(fpStorageST, standards[istd].mat, element_list_std, calConditions, singleSpectrum );
            if ( result < 0 ) {
                cout << "quantStandard failed, result = " << result << "   file " << standards[istd].spectrumFileName << endl;
                error = true;
            };
            //  Complete list of components in spectrum
/*            unsigned int ic;
            cout << "Complete list of components (" << singleSpectrum.numberOfComponents() << "):" << endl;
            for( ic=0; ic<singleSpectrum.numberOfComponents(); ic++ ) {
                cout << "        " << componentDescription( singleSpectrum.component(ic) );
                cout << "  enable " << singleSpectrum.component(ic).enabled;
                cout << "  fit " << singleSpectrum.component(ic).fit;
                cout << "  ignore " << singleSpectrum.component(ic).ignore;
                cout << "  incl " << singleSpectrum.component(ic).included;
                cout << "  bkg indx " << singleSpectrum.component(ic).bkg_index;
                cout << "  coeff " << singleSpectrum.component(ic).coefficient;
                cout << "  matrix " << singleSpectrum.component(ic).matrix;
                cout << "  non-fit fac " << singleSpectrum.component(ic).non_fit_factor;
                cout << "  quant " << singleSpectrum.component(ic).quant;
                cout << endl;
            }       */
            //  Write results to output file and put results in element list for map and calibration files
            //  This will fill the element list for the standards with the coefficients and other fit results
            if( error ) break;  //  Don't write results so error will be at end of putput
            float element_sum = 0;
            result = quantWriteResults( standards[istd].mat, calConditions.detector, element_list_std,
                                singleSpectrum, oxidesOutput, termOutFile, element_sum );
            //  Move expanded element list with results back into standards info list
            //      add_element_list_entry function will automatically select info from appropriate entry
            int iel;
            for( iel=0; iel<element_list_std.size(); iel++ ) {
                add_element_list_entry( element_list_std[iel], standards[istd].element_list );
            }
            //  Add spectrum (with coefficients but erased component spectra) to standards list
            if( cmd == CALIBRATE ) {
                singleSpectrum.clean(); //  No plot, don't need the component spectra, remove to save space
                standards[istd].spectrum = singleSpectrum;
            }
            //  For debugging and evaluation of the calibration, write a synopsis file using the evaluate format
            //      Calibrate does not normally write this, only then the -w option is invoked (which sets the file name)
            if( cmd == CALIBRATE && arguments.cal_eval_file.length() > 0 ) {
                //  Duplicate the master list (to preserve it with default values)
                vector <ElementListEntry> element_list_eval_copy;
                unsigned int iel;
                for( iel=0; iel<element_list_eval_master.size(); iel++ ) element_list_eval_copy.push_back( element_list_eval_master[iel] );
                //  Replace the element list entries with quantification results to the copy (for writing to the map file, preserving element order)
                for( iel=0; iel<element_list_std.size(); iel++ ) {
                    //  Add coefficient as ECF to this new element list (if weight above threshold)
                    element_list_std[iel].ecf = element_list_std[iel].coefficient;
                    if( element_list_std[iel].qualifier == IGNORE
                            || element_list_std[iel].qualifier == EXCLUDE
                            || element_list_std[iel].qualifier == MATRIX ) continue;
                    //  Only include elements that have weight greater then minimum (set in XRFcontrols.h)
                    if( element_list_std[iel].weight <= arguments.min_wgt_eval ) continue;
                    //  Replace the entry in the master list copy
                    int i_copy;
                    for( i_copy=0; i_copy<element_list_eval_copy.size(); i_copy++ ) {
                        if( ! ( element_list_eval_copy[i_copy].element == element_list_std[iel].element ) ) continue;
                        element_list_eval_copy[i_copy] = element_list_std[iel];
                        break;
                    }
                    //termOutFile << "Element list cal  std: " << standards[istd].names[0] << "   s " << XrayMaterial::formula_string( element_list_std[iel].element, element_list_std[iel].stoichiometry ) << "  ";
                    //termOutFile << "     p " << element_list_std[iel].percent;
                    //termOutFile << "     u " << element_list_std[iel].uncertainty;
                    //termOutFile << "     g " << element_list_std[iel].given;
                    //termOutFile << "     o " << element_list_std[iel].stoichiometry.formula_ratio;
                    //termOutFile << "     w " << element_list_std[iel].weight;
                    //termOutFile << "     k " << element_list_std[iel].ecf;
                    //termOutFile << endl;
                }

                // Write header if needed
                if(!wroteMapHeader)
                {
                    string eval_title = "Debug evaluation file for PIQUANT Calibrate sub-command";
                    result = quantWriteMapHeader(fout, eval_title, arguments.quant_map_outputs, element_list_eval_copy, oxidesOutput);
                    wroteMapHeader = true;

                    if ( result != 0 ) {
                        termOutFile << "quantWriteMapHeader failed, result = " << result << endl;
                        error = true;
                        break;
                    }
                }

                //  Write line to map file
                quantWriteMapRow(fout, arguments.quant_map_outputs, element_list_eval_copy,
                                calConditions.detector, singleSpectrum, element_sum);

            }   //  if( cmd == CALIBRATE && arguments.cal_eval_file.length() > 0 )
        }   //  end for( istd=0; istd<standards.size(); istd++ )
        //  Write the results to the calibration file
        //  Note that quantWriteResults must be called before this to populate the fit info for the standards
        termOutFile << endl;
        if( cmd == CALIBRATE && ! error ) {
            if( check_file_extension( arguments.calibration_file, "TXT" ) ) {
                result = quantWriteCalibrationTXT( standards, arguments.calibration_file, datetime() );
            } else if( check_file_extension( arguments.calibration_file, "CSV" ) ) {
                result = quantWriteCalibrationCSV( standards, arguments.calibration_file, datetime() );
            } else {
                termOutFile << "Calibration files can only be .txt or .csv" << endl;
                result = -1;
            }
            if( result < 0 ) {
                cout << "Write of calibration file failed, result = " << result << "   file " << arguments.calibration_file << endl;
                error = true;
            } else {
                termOutFile << "Calibration file written to " << arguments.calibration_file << endl;
            }
        } else if( cmd == CALIBRATE && error ) {
             cout << "Errors in processing standards, no calibration file written." << endl;
        }

    }   //  if( ( cmd == CALIBRATE || cmd == FIT_ONE_STANDARD ) && ( ! error ) )


    //**************************************************************************************
    //      Evaluate calibration by processing the spectrum from each standard as an unknown
    //**************************************************************************************

    if( ( cmd == EVALUATE ) && ( ! error ) ) {
        // The stream to the map file we're writing
        std::ofstream fout(arguments.map_file);
        bool wroteMapHeader = false;
        unsigned int istd;
        //  Loop over all of the standards in the list
        for( istd=0; istd<standards.size(); istd++ ) {
            if( arguments.standard_selected && standard_index != istd )continue;
            //  Read the spectrum file for this entry in the list of standards
            spectrum_vec.clear();
            string eval_optic_file;
			XRFconditionsInput condStruct_eval;
            result = read_spectrum_file( termOutFile, standards[istd].spectrumFileName,
                    spectrum_vec, condStruct_eval );
            if ( result != 0 ) {
                cout << "read_spectrum_file failed, result = " << result << "   file " << standards[istd].spectrumFileName << endl;
                error = true;
                continue;
            };
            //  Set up energy calibration, background parameters, and measurement conditions
            setup_spectrum_parameters( arguments, configSpectrum.calibration(), spectrum_vec,
                    condStruct_config, condStruct_eval, cout );
            //  Combine the spectrum information from several detectors (or the selected detector) into the variable where they will be used
            //      NB: quantCombineSpectra modifies the spectra in the input list to match them to a single energy axis
            //          for proper plotting
            result = quantCombineSpectra( spectrum_vec, singleSpectrum, arguments.detector_select );
            if( result < 0 ) {
                error = true;
                break;
            };
            //  Check live time
            if( singleSpectrum.live_time() <= 0 ) {
                termOutFile << "*** Error - live time is bad, can't use this standard for evaluation. ***" << endl;
                error = true; //  Plot can be vs channels, all others are not possible without calibration
                continue;
            }
            //Check energy calibration
            if( ! singleSpectrum.calibration().good() ) {
                termOutFile << "Bad energy calibration, can't quantify spectrum." << endl;
                error = true;
                break;
            }
            //  Set up new instrument measurement conditions
            XRFconditions evalConditions;
            result = fpSetupConditions ( condStruct_eval, evalConditions );
            if( result < 0 ) {
                termOutFile << "fpSetupConditions failed, result " << result;
                termOutFile << "   error in parameter with keyword " << get_EMSA_keyword( -(result+100) ) << endl;
                error = true;
                return -500 + result;
            };
            termOutFile << endl;
            //  Add the standard name to the spectrum aux info titles (for output in map file if requested using "U")
            if( standards[istd].names.size() > 0 ) {
                singleSpectrum.aux_info_change().titles.clear();
                singleSpectrum.aux_info_change().titles.push_back( standards[istd].names[0] );
            }
            //  Input element fit control list overrides the one from standards input
            vector <ElementListEntry> element_list_std_unk;
            unsigned int i_list;
            //  Move element list from standards input file
            for( i_list=0; i_list<standards[istd].element_list.size(); i_list++ ) {
                //  Only include elements that have positive weight (or are matrix elements)
                if( standards[istd].element_list[i_list].qualifier != MATRIX && standards[istd].element_list[i_list].weight <= arguments.min_wgt_eval ) continue;
                ElementListEntry element_list_std_unk_new = standards[istd].element_list[i_list];
                add_element_list_entry( element_list_std_unk_new, element_list_std_unk );
            }
            //  Append or replace entries using element list from arguments if FORCE is chosen
            for( i_list=0; i_list<element_list.size(); i_list++ ) {
                if( element_list[i_list].qualifier != FORCE ) continue;
                add_element_list_entry( element_list[i_list], element_list_std_unk );
            }
//            for( i_list=0; i_list<element_list_std_unk.size(); i_list++ )
//        termOutFile << "Evaluate  Std Unk Element List " << i_list << "  " << (element_list_std_unk[i_list].qualifier == MATRIX?"MATRIX":"Not Marix") << "  " << element_list_std_unk[i_list].element.symbol() << "   " << element_list_std_unk[i_list].percent << endl;
            //  Start with a fresh material and quantify as unknown using new element list from above
            XrayMaterial std_unknown;
            //  Put the names if the being evaluated into the spectrum object (so that it won't be used for making the ECFs)
            singleSpectrum.std_names( standards[istd].names );
            result = quantUnknown( std_unknown, element_list_std_unk, evalConditions, singleSpectrum, arguments.calibration_file, termOutFile );
            if ( result < 0 ) {
                termOutFile << "quantUnknown failed, result = " << result << "   file " << arguments.spectrum_file << endl;
                error = true;
            }
            //  Write results to output file and put results in element list for map file
            //  This will fill the element list for the standards with the coefficients and other fit results
            if( error ) break;  //  Don't write results so error will be at end of putput
            float element_sum = 0;
            result = quantWriteResults( std_unknown, evalConditions.detector, element_list_std_unk,
                                singleSpectrum, oxidesOutput, termOutFile, element_sum, true );
            //  Code to write plot file after evaluate when single standard chosen and -w option used to set file name
            if( arguments.standard_selected && arguments.cal_eval_file.length() > 0 ) {
                //  Write a CSV file of plotting information (if only one standard was included in the evaluation)
                //  **** Hard-coded file name, for debugging only ****
                termOutFile << endl;
                termOutFile << "Writing plot to file " << arguments.cal_eval_file;
                termOutFile << "      " << singleSpectrum.numberOfChannels() << " channels." << endl;
                result = quantWritePlot( singleSpectrum, arguments.cal_eval_file, cmd, arguments.detector_select, spectrum_vec, getVersionString() );
            }
            //  Duplicate the master list (to preserve it with default values)
            vector <ElementListEntry> element_list_eval_copy;
            unsigned int iel;
            for( iel=0; iel<element_list_eval_master.size(); iel++ ) element_list_eval_copy.push_back( element_list_eval_master[iel] );
            //  Replace the element list entries with quantification results to the copy (for writing to the map file, preserving element order)
            for( iel=0; iel<element_list_std_unk.size(); iel++ ) {
                if( element_list_std_unk[iel].given <= 0 ) continue;
                if( element_list_std_unk[iel].qualifier == IGNORE
                        || element_list_std_unk[iel].qualifier == EXCLUDE
                        || element_list_std_unk[iel].qualifier == MATRIX ) continue;
                //  Replace the entry in the master list copy
                int i_copy;
                for( i_copy=0; i_copy<element_list_eval_copy.size(); i_copy++ ) {
                    if( ! ( element_list_eval_copy[i_copy].element == element_list_std_unk[iel].element ) ) continue;
                    element_list_eval_copy[i_copy] = element_list_std_unk[iel];
                    break;
                }
                //termOutFile << "Element list std unk  s " << XrayMaterial::formula_string( element_list_std_unk[iel].element, element_list_std_unk[iel].stoichiometry ) << "  ";
                //termOutFile << "     lv " << element_list_std_unk[iel].quant_level;
                //termOutFile << "     q " << element_list_std_unk[iel].qualifier;
                //termOutFile << "     p " << element_list_std_unk[iel].percent;
                //termOutFile << "     u " << element_list_std_unk[iel].uncertainty;
                //termOutFile << "     o " << element_list_std_unk[iel].stoichiometry.formula_ratio;
                //termOutFile << "     w " << element_list_std_unk[iel].weight;
                //termOutFile << endl;
            }
/*            for( iel=0; iel<element_list_eval_copy.size(); iel++ ) {
                termOutFile << "Element list eval copy  s " << XrayMaterial::formula_string( element_list_eval_copy[iel].element, element_list_eval_copy[iel].stoichiometry ) << "  ";
                termOutFile << "     lv " << element_list_eval_copy[iel].quant_level;
                termOutFile << "     q " << element_list_eval_copy[iel].qualifier;
                termOutFile << "     p " << element_list_eval_copy[iel].percent;
                termOutFile << "     u " << element_list_eval_copy[iel].uncertainty;
                termOutFile << "     o " << element_list_eval_copy[iel].stoichiometry.formula_ratio;
                termOutFile << "     w " << element_list_eval_copy[iel].weight;
                termOutFile << endl;
            }
*/
            // Write header if needed
            if(!wroteMapHeader)
            {
                // TIMTIME: What title should go in here???
                result = quantWriteMapHeader(fout, "Insert Title Here", arguments.quant_map_outputs, element_list_eval_copy, oxidesOutput);
                wroteMapHeader = true;

                if ( result != 0 ) {
                    termOutFile << "quantWriteMapHeader failed, result = " << result << endl;
                    error = true;
                    break;
                }
            }

            //  Write line to map file
            quantWriteMapRow(fout, arguments.quant_map_outputs, element_list_eval_copy,
                            evalConditions.detector, singleSpectrum, element_sum);

            //  Add spectrum (with coefficients but erased component spectra) to standards list
            singleSpectrum.clean(); //  No plot, don't need the component spectra, remove to save space
            standards[istd].spectrum = singleSpectrum;
        }   //  Loop over all of the standards in the list

        if( ! error ) {
            cout << endl << endl;
            cout << "Map file with evaluate results written to " << arguments.map_file << endl;
            cout << "          quantitative output options " << arguments.quant_map_outputs << endl;
        }

    }   //  if( ( cmd == EVALUATE ) && ( ! error ) )


    //**************************************************************************
    //      Perform an energy calibration on the spectrum using the element list
    //**************************************************************************

    if( ( cmd == ENERGY_CAL ) && ( ! error ) ) {
        float ev_start, ev_ch;
        result = energy_calibrate( singleSpectrum.meas(), element_list, ev_start, ev_ch );
        if ( result < 0 ) {
            termOutFile << "Error finding energy calibration, result = " << result << endl;
            switch( result ) {
                case -1: termOutFile << "Not enough channels in spectrum." << endl; break;
                case -2: termOutFile << "Not enough counts in peak 1 (lower energy)." << endl; break;
                case -3: termOutFile << "Not enough counts in peak 2 (higher energy)." << endl; break;
                case -4: termOutFile << "First element has no emission lines in spectrum range." << endl; break;
                case -5: termOutFile << "Second element has no emission lines in spectrum range." << endl; break;
            }
            error = true;
        } else {
            singleSpectrum.calibration( ev_start, ev_ch );
        //  Write out new energy calibration
            termOutFile << endl;
            termOutFile << "Energy calibration ";
            if( result == 1 ) termOutFile << "(one peak) ";
            termOutFile.precision(1);
            termOutFile << "  eV start = " << ev_start;
            termOutFile.precision(4);
            termOutFile << "  eV/ch = " << ev_ch;
            //  Write energy calibration in an easily-searchable format
            //      so that GUI can find and save it for later subprocess calls
            termOutFile << "                                      ";
            termOutFile << "(-e," << ev_start <<"," << ev_ch << ")";
            termOutFile << endl;

        }
    }   //  if( cmd == ENERGY_CAL )


    //***************************************************************************
    //      Produce a quantitative element map
    //      Read a sequence of spectrum files and quantify each
    //      Write a line to the map file with element fractions for each spectrum
    //      Keep track of locations in the map space
    //***************************************************************************

    if( ( cmd == MAP || cmd == BULK_SUM_MAX ) && ( ! error ) ) {
        string map_spec_file = arguments.spectrum_file;
        //  Recognize and process a list of spectrum files for map and bulk sum if spectrum file ends in ".txt"
        enum FileListType
        {
            LST_INCREMENTING_FILES,
            LST_MSA_LIST,
            LST_PIXLISE_FILE
        };

        FileListType fileListType = LST_INCREMENTING_FILES;
        //  Save the path from the name of the file list and put it into the spectrum file names if they do not have a path included
        string map_spec_file_name;
        string map_spec_file_path;
        /*bool path_found =*/ extract_path( map_spec_file, map_spec_file_path, map_spec_file_name );
        if(check_file_extension( map_spec_file_name, "TXT")) {
            fileListType = LST_MSA_LIST;
        }
        else if(check_file_extension( map_spec_file_name, "PMCS")) {
            fileListType = LST_PIXLISE_FILE;
        }

        ifstream map_file_name_list;
        if( fileListType == LST_MSA_LIST ) {
            map_file_name_list.open( map_spec_file.c_str(), ios::in );
            if( !map_file_name_list ) {
                termOutFile << "Can't open list of spectrum file names from file " << map_spec_file_name << endl;
                error = true;
            }
        }

        ifstream map_file_pmc_list;
        string pixlise_file_name;
        string pmc_list_file_name = map_spec_file;
        if( fileListType == LST_PIXLISE_FILE ) {
            map_file_pmc_list.open( pmc_list_file_name.c_str(), ios::in );
            if( !map_file_pmc_list ) {
                termOutFile << "Can't open list of PMCs to process from file " << pmc_list_file_name << endl;
                error = true;
            }

            // Got the file open, get the first line, as this should be the name/path of the PIXLISE binary file
            getline( map_file_pmc_list, map_spec_file );

            if(map_spec_file.length() <= 4 || map_spec_file.substr(map_spec_file.length()-4) != ".bin") {
                termOutFile << "Did not find PIXLISE binary file name as first line of PMC list file " << pmc_list_file_name << ", read: " << map_spec_file << endl;
                error = true;
            } else {
                // Put the path in front
                map_spec_file = map_spec_file_path+map_spec_file;

                termOutFile << "Read PIXLISE dataset binary file name: " << map_spec_file << endl;
            }
        }

        int n_map_spectra = 0;
        //  We only need these for the BULK_SUM_MAX command, but they must be declared at this scope
        vector <float> bulk_sum;
        vector <float> max_value;
        float sum_live_time = 0;
        float sum_geometry = 0;
        int geometry_count = 0;
        //  Infinite loop prevention value
        int max_map_spectra = 1000000;
        //  For faster debugging, cut the map short
        if( arguments.max_map_arg > 0 ) max_map_spectra = arguments.max_map_arg;
        int sequence_number = -1;

        // If we're processing a map command, start the processing threads
        vector<std::thread> processThreads;
        if(cmd == MAP && !error)
        {
            setMapJobRunning(cmd == MAP);

            termOutFile << "Using " << arguments.map_threads << " threads to process maps." << endl;
            for(int c = 0; c < arguments.map_threads; c++)
            {
                processThreads.push_back(std::thread(processMapJob));
            }
        }

        string pmcLine;
        while( n_map_spectra < max_map_spectra && !error ) {    //  Just to prevent a really infinite loop
            //  Check for spectrum file names from list or from sequential names
            if( fileListType == LST_MSA_LIST ) {
                getline( map_file_name_list, map_spec_file );
                //  Check for end of file, blank line, or comment
                string line_check = upper_trim( map_spec_file );
                if( !map_file_name_list ) break;
                if( line_check.length() < 2 || line_check.substr( 0, 2 ) == COMMENT_STRING ) continue;
                //  Include the path from the name of the file list in the spectrum file name if it does not have a path included
                string maybe_path;
                if( ! extract_path( map_spec_file, maybe_path, map_spec_file_name ) ) {
                    maybe_path = map_spec_file_path + map_spec_file;
                    map_spec_file = maybe_path;
                }
            } else if(fileListType == LST_PIXLISE_FILE) {
                // Expecting this next line to contain a PMC to process, so read it as such
                getline(map_file_pmc_list, pmcLine);

                if(!map_file_pmc_list) {
                    break;
                }

                if(pmcLine.empty()) {
                    termOutFile << "Read invalid PMC from \"" << map_spec_file << "\": \"" << pmcLine << "\"" << endl;
                    error = true;
                    break;
                }
            } else {
                if( n_map_spectra > 0 ) {
                    result = map_spectrum_file_increment( map_spec_file, sequence_number );
                    if( result < 0 ) break;
                    //  See if there is actually a file with this name (if not we are done)
                    ifstream map_spec_stream( map_spec_file.c_str() );
                    if( ! map_spec_stream ) break;
                    map_spec_stream.close();
                } else {
                    //  Call increment function to check and get sequence number for first spectrum
                    string dummy = map_spec_file;
                    result = map_spectrum_file_increment( dummy, sequence_number );
                    sequence_number--;
                }
                if ( result != 0 ) {
                    termOutFile << "No sequence number found in first spectrum file name." << endl;
                    error = true;
                    break;
                }
            }

            if(cmd == BULK_SUM_MAX) {
                result = spectrumBulkSumMax(map_spec_file,

                    condStruct_config,

                    arguments,
                    oxidesOutput,

                    configSpectrum,
                    n_map_spectra,

                    sequence_number,
                    bulk_sum,
                    max_value,
                    sum_live_time,
                    sum_geometry,
                    geometry_count,

                    singleSpectrum,
                    error);
            }
            else // cmd == MAP
            {
                // Let it run, results end up in output queue
                queueMapSpectrum(map_spec_file,

                    condStruct_config,

                    arguments,
                    oxidesOutput,

                    configSpectrum,
                    n_map_spectra,

                    element_list,

                    sequence_number,
                    pmcLine);
            }

            if(result < -1)
            {
                return result;
            }

            if(result == -1)
            {
                break;
            }

            if(result == 0 && error == false)
            {
                n_map_spectra++;

                if( n_map_spectra > max_map_spectra ) {
                    termOutFile << "Maximum number of map spectra exceeded: " << max_map_spectra << endl;
                    break;
                }
            }
        }   //  while( n_map_spectra < max_map_spectra )
        if(cmd == BULK_SUM_MAX ) {
            //  Finish this command by putting results in singleSpectrum for plotting
            singleSpectrum.meas( bulk_sum );
            singleSpectrum.max_value( max_value );
            singleSpectrum.live_time( sum_live_time );
            if( geometry_count > 0 ) {
                float avg_geometry = sum_geometry/geometry_count;
                termOutFile << "Bulk Sum  geometry factor " << avg_geometry << "  (from " << sum_geometry << " spectrum files)." << endl;
                singleSpectrum.geometry( avg_geometry );
                condStruct_spec.conditionsVector[GEOMETRY_INDEX] = avg_geometry;
            }
        }

	    if(cmd == MAP)
        {
            // Wait for all threads to finish
            termOutFile << "Waiting for process threads to finish..." << endl;
            setMapJobRunning(false);
            for(int c = 0; c < processThreads.size(); c++)
            {
                processThreads[c].join();
            }

            termOutFile << "Threads finished, outputting..." << endl;

            termOutFile << endl << endl;

            outputMapFile(termOutFile, arguments, element_list, oxidesOutput);
        }

    }   //  if( ( cmd == MAP || cmd == BULK_SUM_MAX ) )

    //*************************************************************************************
    //      Quantify the composition of an unknown by analyzing the measured XRF spectrum
    //*************************************************************************************

    //  This must come after BULK_SUM_MAX processing to quantify the resultant sum
    if( ( cmd == QUANTIFY || cmd == BULK_SUM_MAX ) && ( ! error ) ) {
        termOutFile << endl;
        if( ! singleSpectrum.calibration().good() ) {
            termOutFile << "Bad energy calibration, can't quantify spectrum." << endl;
            error = true;
        } else if( singleSpectrum.live_time() <= 0 ) {
            termOutFile << "*** Error - live time is bad, can't quantify this spectrum. ***" << endl;
            error = true;
        } else {
            //  Set up new instrument measurement conditions
            XRFconditions unkConditions;
            result = fpSetupConditions ( condStruct_spec, unkConditions );
            if( result < 0 ) {
                termOutFile << "fpSetupConditions failed, result " << result << endl;
                error = true;
                return -500 + result;
            };
            XrayMaterial unknown;
            result = quantUnknown( unknown, element_list, unkConditions, singleSpectrum, arguments.calibration_file, cout );
            if ( result < 0 ) {
                cout << "quantUnknown failed, result = " << result << "   file " << arguments.spectrum_file << endl;
                error = true;
            };
            //  Write results to output file and put results in element list for map and calibration files
            float element_sum = 0;
            //  Normalize result if argument is not zero
            if( arguments.normalization > 0 ) unknown.normalize( arguments.normalization / 100 );
            if( ! error ) result = quantWriteResults( unknown, unkConditions.detector, element_list,
                                singleSpectrum, oxidesOutput, termOutFile, element_sum );
        }
    }   //  if( ( cmd == QUANTIFY || cmd == BULK_SUM_MAX ) && ( ! error ) )



    //*******************************************************************************************
    //      Write a CSV file with plotting information or an MSA file with spectrum information
    //*******************************************************************************************

    if( ( cmd == CALCULATE || cmd == PLOT || cmd == PRIMARY || cmd == COMPARE || cmd == FIT_ONE_STANDARD || cmd == OPTIC_RESPONSE  || cmd == QUANTIFY|| cmd == BULK_SUM_MAX  )
                && ( ! error ) && ( arguments.plot_file.length() > 0 ) ) {
        //  For plot only, show the background (but only if there is an energy calibration to set peak width or it is set by an option)
        if( cmd == PLOT ) {
            if( singleSpectrum.calibration().good() ) {
                //  Use the energy calibration to set peak width
                if( cmd == PLOT ) quantBackground( configConditions, singleSpectrum, true );
            } else {
                //  See if there was a background width set as an option in the command line arguments
                vector <float> bkg_params;
                singleSpectrum.get_bkg_parameters( bkg_params );
                int width_chan = 0;
                if( bkg_params.size() > 1 ) width_chan = bkg_params[1];
                if( width_chan > 0 ) quantBackground( configConditions, singleSpectrum, true );
            }
		}
		//  Check the plot file name for an MSA extension
        if ( check_file_extension( arguments.plot_file, "MSA" ) ) {
            //  Write an MSA file with the spectrum information
            termOutFile << endl;
            termOutFile << "Writing calculated or measured spectrum to file " << arguments.plot_file;
            termOutFile << "      " << singleSpectrum.numberOfChannels() << " channels." << endl;
            //  For bulk sum, write measured spectrum (sum), not calculated spectrum from quantification
            bool bulk_sum = ( cmd == BULK_SUM_MAX );
            result = write_EMSA_PIXL( singleSpectrum, arguments.plot_file, bulk_sum );
        } else {
            //  Write a CSV file of plotting information
            termOutFile << endl;
            termOutFile << "Writing plot to file " << arguments.plot_file;
            termOutFile << "      " << singleSpectrum.numberOfChannels() << " channels." << endl;
            result = quantWritePlot( singleSpectrum, arguments.plot_file, cmd, arguments.detector_select, spectrum_vec, getVersionString() );
        }
     }


    //**************************************************************************
    //      Convert output of SEND_SDD_DATA command to EDR (csv) format
    //**************************************************************************

    if( ( cmd == EM_SDD_DATA ) && ( ! error ) ) {
        spectrum_vec.clear();
        result = histogram_from_SDD_data( arguments.spectrum_file,
                arguments.map_file, spectrum_vec );
        if ( result < 0 ) {
            termOutFile << "Reading input SDD data file (SDF contents) failed, result = " << result << endl;
            error = true;
        } else {
            termOutFile << "Writing " << spectrum_vec.size() << " histograms to CSV file in EDR format, two to a line, file name " << arguments.map_file << endl;
            int i_spec;
            //  Write histograms to EDR format CSV file, two to a line
            for( i_spec=0; i_spec<spectrum_vec.size()/2; i_spec+=2 ) {
                //  Write line to map file
                result = write_EDRhistogram_data( i_spec, spectrum_vec[i_spec],
                                    spectrum_vec[i_spec+1], arguments.map_file );
                if ( result != 0 ) {
                    termOutFile << "write_EDRhistogram_data failed, result = " << result << endl;
                    switch( result ) {
                        case -1: termOutFile << "Output file name has zero length." << endl; break;
                        case -2: termOutFile << "Could not open output file." << endl; break;
                    }
                    error = true;
                    break;
                };
            }
        }
    }   //  if( cmd == EM_SDD_DATA )



	finishTime = clock();
	double duration = (double) ( finishTime - startTime ) / CLOCKS_PER_SEC;
	termOutFile << endl << "Execution finished";
	termOutFile << ", CPU time " << duration << " secs." << endl;
    termOutFile << endl;

    if( cout_redirect ) {
        cout.rdbuf(coutbuf); //reset to standard output again
        terminal_out_stream.close();
    }

    if( error ) return -1;
    return 0;

    } catch( string s ) {
        cerr << "Bad error: " << s << endl;
        return -1;
    }

};
