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

#include "quantUnknown.h"
#include "quantWriteCalibrationTXT.h"
#include "quantBackground.h"
#include "quantComponents.h"
#include "quantIgnore.h"
#include "quantCalculate.h"
#include "quantFitSpectrum.h"
#include "setupStandardsCSV.h"
#include "quantECFs.h"
#include "fpMain.h"
#include "fpLineSpectrum.h"
#include "Lfit.h"
#include "Fit.h"
#include "upper_trim.h"
#include "differentiate.h"
#include "XRFconstants.h"
#include "XRFcontrols.h"
#include "XRFutilities.h"


using namespace std;

//      Quantify the composition of an unknown material
//          by analyzing the measured XRF spectrum
//      Fit the spectrum to calculated components
//      Return calculated spectrum components and best-fit coefficients

//  Written Mar. 29, 2017
//      Adapted from quantStandard.cpp
//  Modified June 9, 2017
//      Add background as a fitted component
//      Fix exclusion of negative components
//  Modified June 27, 2017
//      Write error message when element list for unknown is empty
//  Modified Sept. 29, 2017
//      Skip disabled components in loop to re-calculate ignore components
//      Also disable ignore components with zero, negative, or nan calculated intensity
//  Modified Dec. 13, 02017
//      Add check for minimum energy to escape peaks (passed to fpLineSpectrum)
//  Modified Mar. 2, 2018
//      Store iterations in spectrum
//      Store ECFs used for quantification in element list for output to map file
//      Use utility to check file extensions
//  Modified Apr. 12, 2018
//      Remove fitElements vector and use unkElements vector already defined referencing same target
//          (This should have no effect on anything)
//      Disable components if coefficient is negative and number of iterations is >= minimum
//          (was min - 1, now matches condition for setting fraction to zero)
//  Modified July 18, 2018 to call fraction_input instead of fraction_oxide to adjust composition
//  Modified June 6, 2019
//      Handle MATRIX elements correctly (for unknowns and for evaluate)
//      Remove new_element_list and new_fraction_list (not used anywhere)
//  Modified July 2, 2019
//      Allow for background to be adjusted manually instead of in least-squares fit, using -b option parameter
//  Modified July 3, 2019
//      Provide for use of calculated instead of SNIP background (with or without fit)
//      Background is now always present as a component, but is excluded (and not used) if not fit
//      All of the background setup is taken care of in quantBackground (which returns 1 if calculated bkg is used)
//  Modified July 22, 2019
//      Change quantDefaults to treat any elements that have no components in the spectrum as matrix elements, and print warning
//          (Changed quantDefaults argument list from elements to element list; moved loading matrix elements into XrayMaterial to after component setup)
//  Modified Nov. 4, 2019
//      Include standard deviation for ECFs in element list (from quantECFs)
//  Modified Dec. 16, 2019
//      Add element qualifier OUTPUT ("O") to force the element to be included in the evaluate list (with zeros if not in any standard in this run)
//      Elements in input list must be ignored here if they have that qualifier
//      Re-arrange how elements are added to only add if qualifier matches (NO_QUALIFIER or FORCE to add quantified element, MATRIX to add matrix element)
//  Modified Nov. 4, 2020
//      Allow for multiple background regions for fit, change how bkg components added (in quantBackground)
//  Modified Nov. 30, 2020
//      Exclude from fit vector any components that have fit Boolean set to false
//      Also add factor to compute coefficients of non-fit components from components used for quant
//  Modified Dec. 15, 2020
//      Modified for elements to be included as carbonates instead of oxides
//  Modified Dec. 29, 2020  Use SNIP for low-energy background, disable shelf calculations
//  Modified Feb. 2, 2021   Avoid evaluating a standard with itself as a calibration standard
//  Modified Feb. 26, 2021  Add DETECTOR_SHELF and DETECTOR_CE spectrum component types
//  Modified Apr. 5, 2021   Remove DETECTOR_SHELF spectrum component, shelf now in fpLineSpectrum
//                          Disable calibration file (all ECFs are unity, don't open or read cal file)
//  Modified Apr. 6, 2021   Add CONTINUUM component type and sort out how to handle background and Compton escape (remove Det shelf component)
//  Modified July 10, 2021  Add simple pulse pileup calculation - change return for fpLineSpectrum (note ignore peaks not included in pileup)


int quantUnknown( XrayMaterial &unknown, vector <ElementListEntry> &element_list,
        XRFconditions &conditions, XraySpectrum &unkSpectrum, const std::string &calFileName, std::ostream &logger ) {
//		check input parameters
	if( ! unkSpectrum.calibration().good() ) return -520;
	if( unkSpectrum.live_time() <= 0 ) return -521;
	int nChan = unkSpectrum.numberOfChannels();

	int result;

	//  Old-style txt calibration file with only elements and ECFs
    vector <Element> cal_element_list;
    vector <float> cal_factor_list;
    //  Standards list with ECFs and weights from a csv calibration file
    vector <StandardInformation> cal_standards;
    //  Read in ECF list or standards information from calibration file
    if( check_file_extension( calFileName, "TXT" ) ) {
        int ne_in = quantReadCalibrationTXT( calFileName, cal_element_list, cal_factor_list, logger );
        if( ne_in > 0 ) {
            logger << "Calibration file read OK, " << cal_element_list.size() << " element calibration factors." << endl;
        } else {
            logger << "No element calibration factors read in from file." << endl;
        }
    } else if( check_file_extension( calFileName, "CSV" ) ) {
        result = setupStandardsCSV( calFileName, cal_standards, MINIMUM_AMOUNT );
        if ( result != 0 ) {
            logger << "Calibration file read failed, result = " << result << endl;
        } else {
            logger << "Calibration file read OK, entries for " << cal_standards.size() << " standards read in." << endl;
            logger << endl;
            //  Disable any standards in the list that match the one being evaluated
            const vector <string> &eval_names = unkSpectrum.std_names();
            if( eval_names.size() > 0 ) {
                unsigned int is;
                int stds_count = 0;
                for( is=0; is<cal_standards.size(); is++ ) {
                    bool name_match_found = false;
                    unsigned int in;
                    for( in=0; in<cal_standards[is].names.size(); in++ ) {
                        unsigned int in_eval;
                        for( in_eval=0; in_eval<eval_names.size(); in_eval++ ) {
                            if( eval_names[in_eval] == cal_standards[is].names[in] ) name_match_found = true;
                            //cout << "Std disable " << eval_names[in_eval] << "     " << cal_standards[is].names[in] << "   " << name_match_found << endl;
                        }
                    }
                    if( name_match_found ) {
                        cal_standards[is].disable = true;
                        logger << "Standard    " << ( cal_standards[is].names.size()>0?cal_standards[is].names[0]:"");
                        logger << " (# " << is << ") is disabled for this evaluation." << endl;
                    } else {
                        stds_count++;
                        cal_standards[is].disable = false;
                    }
                }
                if( stds_count == 0 ) {
                    logger << "Error - no calibration standards for " << eval_names[0] << " during Evaluate." << endl;
                }
            }
        };
    } else {
        logger << "Calibration files can only be .txt or .csv" << endl;
        result = -1;
    }


    //  Set up the list of elements in the unknown with initial guess at fractions
    int ie;
    float trial_fraction = 1.0f / element_list.size();
    for( ie=0; ie<element_list.size(); ie++ ) {
        if( element_list[ie].qualifier == NO_QUALIFIER || element_list[ie].qualifier == FORCE ) {
            unknown.add_element( element_list[ie].element, trial_fraction, element_list[ie].stoichiometry );
        }
    }
    const vector <Element> &unk_elements = unknown.original_element_list();
    if( unk_elements.size() <= 0 ) {
        logger << "No elements specified for unknown quantification." << endl;
        return -580;
    }

    //  Find the calibration factors for the elements in the specimen
    vector <float> unk_factors_list;
    //  Empty vector to use default ECF calculation (not composition-specific, maybe in the future)
    vector <float> unk_fractions_dummy;
    //  Added return of ECF standard deviations            Nov. 4, 2019
    vector <float> unk_ECF_SDs;
    //  Standards list will be used if populated, otherwise element & ECF lists
    result = quantECFs( cal_standards, cal_element_list, cal_factor_list, unk_elements, unk_fractions_dummy, unk_factors_list, unk_ECF_SDs, logger );
    //  Put ECFs into element list for output to results and map file
    for( ie=0; ie<element_list.size(); ie++ ) {
        int iu;
        for( iu=0; iu<unk_elements.size(); iu++ ) {
            if( ! ( element_list[ie].element == unk_elements[iu] ) ) continue;
            element_list[ie].ecf = unk_factors_list[iu];
            element_list[ie].ecf_sigma = unk_ECF_SDs[iu];
            break;
        }
    }

    //  Set up components for the calculated spectrum
    vector <SpectrumComponent> components;
    vector <XrayLines> pureLines;
    if( Compton_escape_enable_flag ) result = makeComponents( DETECTOR_CE, pureLines, components, 0 );

    // Include components for any elements to be included in fit but ignored in composition
    vector <XrayLines> ignoreLines;
    result = quantIgnore( element_list, conditions, unkSpectrum, ignoreLines );
    if( result < 0 ) {
        logger << "quantIgnore failed to set up components for ignored elements, result is " << result << endl;
        return -540 + result;
    }

    vector <XrayLines> sourceLines;
    //  Load vector with emission lines from X-ray source
    conditions.source.lines( sourceLines, conditions.eMin );
    //  Load vector with pure element emission lines from specimen and set up FP calculations
    FPstorage fpStorage;
    fpPrep(fpStorage, unknown, conditions, pureLines );
    int i;
    for( i=0; i<pureLines.size(); i++ ) pureLines[i].commonFactor( unkSpectrum.live_time() );
    //  Copy list of pure element lines and remove any matrix elements before setting up spectrum components
    vector <XrayLines> pureLines_nonMatrix;
    for( i=0; i<pureLines.size(); i++ ) {
        bool is_matrix = false;
        for( ie=0; ie<element_list.size(); ie++ ) {
            if( ! ( pureLines[i].edge().element() == element_list[ie].element ) ) continue;
            if( element_list[ie].qualifier == MATRIX ) is_matrix = true;
            break;
        }
        if( ! is_matrix ) pureLines_nonMatrix.push_back( pureLines[i] );
    }
    //  Set up components for everything except background (bkg will be handled in quantBackground.cpp)
    result = setupComponents( sourceLines, pureLines_nonMatrix, components );
    if( result < 0 ) {
        logger << "setupComponents failed, result is " << result << endl;
        return -540 + result;
    }
    //  Use the element list to choose components to quantify
    result = quantComponents( element_list, components );
    if( result < 0 ) {
        logger << "quantComponents failed, result is " << result << endl;
        return -550 + result;
    }
    //  See if there is a component to quantify each element and pick a default if not
    result = quantDefaults( element_list, components );
    if( result < 0 ) {
        logger << "quantDefaults failed, result is " << result << endl;
        return -560 + result;
    }
    //  Put in extra components for debugging extra intensity in tube scatter peaks from L line
    for( i=0; i<sourceLines.size(); i++ ) {
        if( sourceLines[i].edge().index() == L3 ) {
            vector <XrayLines> temp_lines;
            temp_lines.push_back( sourceLines[i] );
            result = makeComponents( La, temp_lines, components );
            if( result < 0 ) {
                cout << "makeComponents failed for extra La line, result is " << result << endl;
                return -760 + result;
            }
        } else if( sourceLines[i].edge().index() == L2 ) {
            vector <XrayLines> temp_lines;
            temp_lines.push_back( sourceLines[i] );
            result = makeComponents( Lb1, temp_lines, components );
            if( result < 0 ) {
                cout << "makeComponents failed for extra Lb1 line, result is " << result << endl;
                return -770 + result;
            }
        }
    }
    //  Add the components to the spectrum object
    int ic;
    for( ic=0; ic<components.size(); ic++ ) {
        //  Leave out Compton lines from tube L edges (fit with extra La and Lb1 lines above)
        if( components[ic].type == COMPTON && components[ic].level == L ) continue;
//        if( components[ic].type == RAYLEIGH && components[ic].level == L ) continue;
        if( components[ic].type == La ) continue;
//        if( components[ic].type == Lb1 ) continue;
        //  Check for non-fit components and set factor for coefficients to track quant components
        if( components[ic].type == ELEMENT ) {
            unsigned int ic_q;
            bool different_quant_found = false;
            for( ic_q=0; ic_q<components.size(); ic_q++ ) {
                //  Skip this component and any components with a different element
                if( ic_q == ic ) continue;
                if( ! ( components[ic_q].element == components[ic].element ) ) continue;
                //  Skip any components with the same element and level (must be a different level for quantification)
                if( components[ic_q].type == ELEMENT && components[ic_q].level == components[ic].level ) continue;
                if( ! components[ic_q].quant ) continue;
                different_quant_found = true;
            }
            if( different_quant_found ) {
                components[ic].fit = false;   //  Don't fit this component
                //  Relate its coefficient to the element's quant component
                if( components[ic].level == L ) components[ic].non_fit_factor = COEFF_RATIO_L_K;    //  defined in XRFcontrols.h
                if( components[ic].level == M ) components[ic].non_fit_factor = COEFF_RATIO_M_L;    //  defined in XRFcontrols.h
            }
        }
        unkSpectrum.add_component( components[ic] );
    }

    //  Load the matrix elements into the XrayMaterial object
    for( ie=0; ie<element_list.size(); ie++ ) {
        if( element_list[ie].qualifier == MATRIX ) {  //
            unknown.add_element( element_list[ie].element, element_list[ie].percent/100, element_list[ie].stoichiometry );
            unknown.uncertainty( element_list[ie].element, element_list[ie].uncertainty / 100 );
        }
    }
    //  Set up FP calculations with final element list (including matrix)
    vector <XrayLines> pureLines_matrix;    //  This gets ignored since there are no matrix elements with useful emission lines
    fpPrep (fpStorage, unknown, conditions, pureLines_matrix );

    //  Check the bkg parameters and calculate the background using the SNIP digital filter (if selected)
    //  Checks for calculated background and for bkg fit (and adds components if fit)
    result = quantBackground( conditions, unkSpectrum );
    if ( result < 0 ) {
        logger << "quantBackground failed, result = " << result << endl;
        return -530 + result;
    };

    //  Fit the components to the measured spectrum
    int iterations = 0;
    bool done = false;
    while( iterations < MAX_ITERATIONS && ( ! done )) { //  MAX_ITERATIONS Defined in XRFcontrols.h
        iterations++;
        // Calculate spectrum for this unknown, updating component spectra
        result = quantCalculate(fpStorage, unknown, conditions, unkSpectrum );
        if ( result != 0 ) {
            logger << "quantCalculate failed, result = " << result << endl;
            return -560 + result;
        };
        //  Also re-calculate the ignored elements since the energy calibration may have been adjusted
        for( ic=0; ic<unkSpectrum.numberOfComponents(); ic++ ) {
            SpectrumComponent updated_component = unkSpectrum.component( ic );
            if( updated_component.type != ELEMENT ) continue;
            if( ! updated_component.ignore ) continue;
            if( ! updated_component.enabled ) continue;
            int i;
            for( i=0; i<updated_component.spectrum.size(); i++ ) updated_component.spectrum[i] = 0;
            int il;
            for ( il=0; il<ignoreLines.size(); il++ ) {
                if ( ignoreLines[il].numberOfLines() <= 0 ) continue;
        //			get approximate energy for detector resolution and bkg noise threshold
                float en = ignoreLines[il].energy(0);
        //			get background noise for threshold
                int k = unkSpectrum.channel( en );
                float threshold = 1;
                if ( k >= 0 && k < nChan && unkSpectrum.bkg()[k] > 0 ) threshold = 0.1f * sqrt( unkSpectrum.bkg()[k] );
                vector <LineGroup> dummy;
                fpLineSpectrum( ignoreLines[il], conditions.detector, threshold, unkSpectrum.calibration(), conditions.eMin, dummy, updated_component );
            };
            //  Check for zero (or nan) and disable (also write message)
            //  Only if not already disabled to avoid many messages (check is at top of loop)
            float sum = 0;
            for( i=0; i<updated_component.spectrum.size(); i++ ) sum += updated_component.spectrum[i];
            if( sum <= 0 || isnan( sum ) ) {
                logger << "*** Warning - calculated intensity is zero (or negative or nan) for";
                logger << " ignored component " << componentDescription( updated_component );
                logger << " (it is being disabled).   " << sum << endl;
                unkSpectrum.disable( ic );
            }
            //  Put the new calculation into the XraySpectrum object
            unkSpectrum.update_component( updated_component );
        }
        result = quantFitSpectrum( conditions, unkSpectrum, logger );
        if ( result < 0 ) {
            logger << "quantFitSpectrum failed, result = " << result << endl;
            return -570 + result;
        } else if( result == 0 ) {
            done = true;
        };
        if( iterations < MINIMUM_ITERATIONS ) done = false;   //  Defined in XRFcontrols.h
        //  Improve composition using fit
        bool all_zero = true;
        for( ie=0; ie<unk_elements.size(); ie++ ) {
            float fraction = unknown.fraction_input( unk_elements[ie] );
            float fraction_save = fraction;
            float coeff = unkSpectrum.coefficient( unk_elements[ie] );
            if( coeff > 0 ) {
                fraction *= coeff / unk_factors_list[ie];
                //  Calculate an adjusted coefficient to better match now fraction (for use in shelf calcs)
                float adj_coeff = unk_factors_list[ie];
                if( fraction > 0 ) adj_coeff = coeff * fraction_save / fraction;
                unkSpectrum.adjusted_coefficient( unk_elements[ie], adj_coeff );
            } else if( coeff == COEFFICIENT_NO_COMPONENT ) {
                continue;   //  No spectrum component for this element, don't adjust it
            } else {
                //  Need to do fit at least once after disabling negative coeff
                if( iterations <= MINIMUM_ITERATIONS - 1 ) {
                    fraction = NEGLIGIBLE_FRACTION;   //  Defined in XRFcontrols.h
                    unkSpectrum.adjusted_coefficient( unk_elements[ie], MINIMUM );
                } else {
                    fraction = 0;
                }
            }
            if( isinf(fraction) ) logger << "Frac nan " << unk_elements[ie].symbol() << "  " << fraction_save << "  " << fraction << "  " << coeff << "  " << unk_factors_list[ie] << endl;
            unknown.fraction( unk_elements[ie], fraction );
            if( fraction > 0 ) all_zero = false;
        }
        //  The following lines make the fit a non-negative least squares algorithm following Lawson and Hanson (1974)
        for( ic=0; ic<unkSpectrum.numberOfComponents(); ic++ ) {
            //  Need to do fit at least once after disabling negative coeff
            if( unkSpectrum.component( ic ).coefficient < 0 && iterations >= MINIMUM_ITERATIONS ) {
                unkSpectrum.disable( ic );
            }
        }
        if( all_zero ) unknown.fraction( unk_elements[0], MINIMUM );
//        if( iterations >= 1 ) done = true;
    }

    unkSpectrum.iterations( iterations );
	return iterations;

};
