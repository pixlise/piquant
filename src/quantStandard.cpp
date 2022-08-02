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

#include "quantStandard.h"
#include "quantBackground.h"
#include "quantComponents.h"
#include "quantIgnore.h"
#include "quantCalculate.h"
#include "quantFitSpectrum.h"
#include "fpMain.h"
#include "fpLineSpectrum.h"
#include "fpConvolve.h"
#include "Lfit.h"
#include "Fit.h"
#include "differentiate.h"
#include "XRFconstants.h"
#include "XRFcontrols.h"


using namespace std;

//		Process a standard material of given composition
//      Fit the spectrum to calculated components
//      Return calculated spectrum components and best-fit coefficients

//  Written Mar. 16, 2017
//      Based on stdQuant.cpp of Aug. 10, 2016
//      And specElementFit.cpp of June 27, 2016
//  Modified June 7, 2017
//      Add background as a fitted component
//  Modified Sept. 29, 2017
//      Skip disabled components in loop to re-calculate ignore components
//      Also disable ignore components with zero, negative, or nan calculated intensity
//  Modified Dec. 13, 02017
//      Add check for minimum energy to escape peaks (passed to fpLineSpectrum)
//  Modified Mar. 2, 2018
//      Store iterations in spectrum
//  Modified June 6, 2019
//      Handle MATRIX elements correctly (remove from pure lines before setting up spectrum components)
//  Modified July 2, 2019
//      Allow for background to be adjusted manually instead of in least-squares fit, using -b option parameter
//  Modified July 3, 2019
//      Provide for use of calculated instead of SNIP background (with or without fit)
//      Background is now always present as a component, but is excluded (and not used) if not fit
//      All of the background setup is taken care of in quantBackground (which returns 1 if calculated bkg is used)
//  Modified July 22, 2019
//      Change quantDefaults to treat any elements that have no components in the spectrum as matrix elements, and print warning
//          (Changed quantDefaults argument list from elements to element list)
//  Modified Nov. 4, 2020
//      Allow for multiple background regions for fit, change how bkg components added (in quantBackground)
//  Modified Nov. 30, 2020
//      Exclude from fit vector any components that have fit Boolean set to false
//      Also add factor to compute coefficients of non-fit components from components used for quant
//  Modified Dec. 4, 2020   Use SNIP background to fit anomalous background at low energies
//  Modified Dec. 29, 2020  Use SNIP for low-energy background, disable shelf calculations
//  Modified Apr. 6, 2021   Add CONTINUUM component type and sort out how to handle background and Compton escape (remove Det shelf component)
//  Modified July 10, 2021  Add simple pulse pileup calculation - change return for fpLineSpectrum (note ignore peaks not included in pileup)



int quantStandard(FPstorage &fpStorage, const XrayMaterial &standard, vector <ElementListEntry> element_list,
        XRFconditions &conditions, XraySpectrum &stdSpectrum ) {
//		check input parameters
	if( ! stdSpectrum.calibration().good() ) return -520;
	if( stdSpectrum.live_time() <= 0 ) return -521;
	int nChan = stdSpectrum.numberOfChannels();

	int result;

    //  Set up components for the calculated spectrum
    vector <SpectrumComponent> components;
    vector <XrayLines> pureLines;
    if( Compton_escape_enable_flag ) result = makeComponents( DETECTOR_CE, pureLines, components, 0 );

    // Include components for any elements to be included in fit but ignored in composition
    vector <XrayLines> ignoreLines;
    result = quantIgnore( element_list, conditions, stdSpectrum, ignoreLines );
    if( result < 0 ) {
        cout << "quantIgnore failed to set up components for ignored elements, result is " << result << endl;
        return -540 + result;
    }

    vector <XrayLines> sourceLines;
    //  Load vector with emission lines from X-ray source
    conditions.source.lines( sourceLines, conditions.eMin );
    //  Load vector with pure element emission lines from specimen and set up FP calculations
    fpPrep(fpStorage, standard, conditions, pureLines );
    //  Copy list of pure element lines and remove any matrix elements before setting up spectrum components
    vector <XrayLines> pureLines_nonMatrix;
    unsigned int i;
    for( i=0; i<pureLines.size(); i++ ) {
        bool is_matrix = false;
        unsigned int ie;
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
        cout << "setupComponents failed, result is " << result << endl;
        return -540 + result;
    }

    //  Use the element list to choose components to quantify or remove
    result = quantComponents( element_list, components );
    if( result < 0 ) {
        cout << "quantComponents failed, result is " << result << endl;
        return -550 + result;
    }
    //  See if there is a component to quantify each element and pick a default if not
    result = quantDefaults( element_list, components );
    if( result < 0 ) {
        cout << "quantDefaults failed, result is " << result << endl;
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
    unsigned int ic;
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
        stdSpectrum.add_component( components[ic] );
    }

    //  Check the bkg parameters and calculate the background using the SNIP digital filter (if selected)
    //  Checks for calculated background and for bkg fit (and adds components if fit)
    result = quantBackground( conditions, stdSpectrum );
    if ( result < 0 ) {
        cout << "quantBackground failed, result = " << result << endl;
        return -530 + result;
    };

    //  Fit the components to the measured spectrum (without changing the composition)
    int iterations = 0;
    bool done = false;
    while( iterations < MAX_ITERATIONS && ( ! done )) {
        iterations++;
        // Calculate spectrum for this standard, updating component spectra
        result = quantCalculate(fpStorage, standard, conditions, stdSpectrum );
        if ( result != 0 ) {
            cout << "quantCalculate failed, result = " << result << endl;
            return -570 + result;
        };
        //  Also re-calculate the ignored elements since the energy calibration may have been adjusted
        for( ic=0; ic<stdSpectrum.numberOfComponents(); ic++ ) {
            SpectrumComponent updated_component = stdSpectrum.component( ic );
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
                int k = stdSpectrum.channel( en );
                float threshold = 1;
                if ( k >= 0 && k < nChan && stdSpectrum.bkg()[k] > 0 ) threshold = 0.1f * sqrt( stdSpectrum.bkg()[k] );
                vector <LineGroup> dummy;
                fpLineSpectrum( ignoreLines[il], conditions.detector, threshold, stdSpectrum.calibration(), conditions.eMin, dummy, updated_component );
            };
            //  Check for zero (or nan) and disable (also write message)
            //  Only if not already disabled to avoid many messages (check is at top of loop)
            float sum = 0;
            for( i=0; i<updated_component.spectrum.size(); i++ ) sum += updated_component.spectrum[i];
            if( sum <= 0 || isnan( sum ) ) {
                cout << "*** Warning - calculated intensity is zero (or negative or nan) for";
                cout << " ignored component " << componentDescription( updated_component );
                cout << " (it is being disabled).   " << sum << endl;
                stdSpectrum.disable( ic );
            }
            //  Put the new calculation into the XraySpectrum object
            stdSpectrum.update_component( updated_component );
        }
        result = quantFitSpectrum( conditions, stdSpectrum, cout );
        if ( result < 0 ) {
            cout << "quantFitSpectrum failed, result = " << result << endl;
            return -580 + result;
        } else if( result == 0 ) {
            done = true;
        };
//        cout << "New fit     iter " << iterations << "    chi sq " << stdSpectrum.chisq() << endl;
//        cout << endl;
        if( iterations < MINIMUM_ITERATIONS ) done = false;   //  Defined in XRFcontrols.h
        //  The following lines make the fit a non-negative least squares algorithm following Lawson and Hanson (1974)
        //  Check for negative coefficients and disable if more than a few iterations
        //      (Don't disable immediately because first fits may be very far off)
        for( ic=0; ic<stdSpectrum.numberOfComponents(); ic++ ) {
            //  Need to do fit at least once after disabling negative coeff
            if( stdSpectrum.component( ic ).coefficient < 0 && iterations >= MINIMUM_ITERATIONS - 1 ) {
                stdSpectrum.disable( ic );
            }
        }
    }

    stdSpectrum.iterations( iterations );
	return iterations;

};
