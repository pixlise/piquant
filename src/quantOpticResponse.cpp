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

#include "quantOpticResponse.h"
#include "quantBackground.h"
#include "quantComponents.h"
#include "quantIgnore.h"
#include "quantCalculate.h"
#include "quantFitSpectrum.h"
#include "fpMain.h"
#include "fpLineSpectrum.h"
#include "Lfit.h"
#include "Fit.h"
#include "differentiate.h"
#include "XRFconstants.h"
#include "XRFcontrols.h"
#include "spline.h"
#include "fpBeams.h"


using namespace std;

//		Process a standard material of given composition
//      and a measured spectrum to obtain or adjust
//      the response curve of the optic in the primary beam

//  Written Oct. 22, 2020
//      Based on quantStandard.cpp of July 22, 2019
//  Modified Apr. 28, 2021  Lots of changes in other parts of program, update this to work with them
//                          Make this more compatible with unity ECFs approach
//                          Change zero-energy val;ue multiplier with respect to first low energy value
//  Modified May 10, 2021   Put optic energies in spectrum object, not pass to quantCalculate (needed for new bkg options)
//  Modified July 10, 2021  Add simple pulse pileup calculation - change return for fpLineSpectrum (note these peaks not included in pileup)

#define ZERO_EN_OPTIC_MULTIPLIER 2.3f

int quantOpticResponse(FPstorage &fpStorage, const XrayMaterial &standard, vector <ElementListEntry> element_list,
        XRFconditions &conditions, XraySpectrum &stdSpectrum ) {
//		check input parameters
	if( ! stdSpectrum.calibration().good() ) return -520;
	if( stdSpectrum.live_time() <= 0 ) return -521;
	int nChan = stdSpectrum.numberOfChannels();

	int result;
    unsigned int i;

    //  Check the bkg parameters and calculate the background using the SNIP digital filter (no_calc = true)
    result = quantBackground( conditions, stdSpectrum );
    if ( result < 0 ) {
        cout << "quantBackground failed, result = " << result << endl;
        return -530 + result;
    };

    //  Set up the new optic response function

    //  Energies for optic response (for least-squares fit in separate regions)const int n_regions = 11;   //      (1.84 is Si K edge, discontinuity possible)
//    float split_en_def = { 0, 1.84, 1.84, 4, 6, 8, 10, 12, 15, 20, 30 };
//    const int n_regions = split_en_def.size();   //      (1.84 is Si K edge, discontinuity possible)
//    float split_en_def[n_regions] = { 0, 2, 4, 6, 8, 10, 12, 15, 20, 25, 30 };
//    const int n_regions = 9;   //      (1.84 is Si K edge, discontinuity possible)
 //   float split_en_def[n_regions] = { 2, 4, 6, 8, 10, 12, 15, 20, 30 };   //  Best so far, Nov. 1, 2020 1:43pm
    vector <float> optic_energies = { 0, 4, 6, 8, 10, 12, 14, 16, 18, 20, 25, 30 };
    const int n_regions = optic_energies.size();
    int jj;
    for( jj=0; jj<n_regions; jj++ ) optic_energies[jj] = optic_energies[jj] * 1000;
    //  Set up spline approximation to optic response
    vector <float> optic_values( n_regions, 1 );
    vector <float> optic_derivatives( n_regions, 0 );
    cout.precision(4);

    //  Calculate the X-ray tube output and primary beam intensity at the optic energies
    //      This is for algorithm deployment, not needed in production
    vector <float> tube_output(n_regions);
    for( jj=0; jj<n_regions; jj++ ) tube_output[jj] = conditions.source.continuum( optic_energies[jj] );
    vector <float> primary_beam(n_regions);
    for( jj=0; jj<n_regions; jj++ ) primary_beam[jj] = tube_output[jj]; //  continuum only
    fpIncidentBeam ( conditions, optic_energies, primary_beam );
    //	result is per keV, so multiply by channel width in keV to get counts in each channel
	for( jj=0; jj<n_regions; jj++ ) primary_beam[jj] *= stdSpectrum.calibration().energyPerChannel() / 1000;
    cout.precision(0);
    cout << endl;
    cout << "   Optic energies   ";
    for( i=0; i<optic_energies.size(); i++ ) cout << ",  " << optic_energies[i];
    cout << endl;
    cout << "   Tube output   ";
    for( i=0; i<tube_output.size(); i++ ) cout << ",  " << tube_output[i];
    cout << endl;
    cout << "   Primary beam   ";
    for( i=0; i<primary_beam.size(); i++ ) cout << ",  " << primary_beam[i];
    cout << endl;
    cout.precision(4);


    //  Fit to the SNIP background to get an initial adjustment to the response
    //      before including any peaks in the fit

    //  Make a fake spectrum with only background
    XraySpectrum bkg_spec( stdSpectrum.bkg(), stdSpectrum.calibration().energyStart(), stdSpectrum.calibration().energyPerChannel() );
    bkg_spec.live_time( stdSpectrum.live_time() );
    bkg_spec.real_time( stdSpectrum.real_time() );
    //  Set up components for the calculated spectrum
    vector <SpectrumComponent> bkg_components;
    vector <XrayLines> pureLines;
    //  Set up components for background only in several split components (empty lines vector)
    result = makeComponents( CONTINUUM, pureLines, bkg_components, n_regions );
    if( result < 0 ) {
        cout << "makeComponents failed, result is " << result << endl;
        return -540 + result;
    }
    //  Add the components to the spectrum object
    int ic;
    for( ic=0; ic<bkg_components.size(); ic++ ) {
        bkg_components[ic].plot = true;
        bkg_spec.add_component( bkg_components[ic] );
    }
    bkg_spec.put_bkg_split( optic_energies );
    //  Fit the components to the SNIP background iteratively
    const int MAX_ITERATIONS_SNIP = 1;
    int iterations_SNIP = 0;
    bool done_SNIP = false;
    while( iterations_SNIP < MAX_ITERATIONS_SNIP && ( ! done_SNIP )) {
        iterations_SNIP++;
        //  Load vector with pure element emission lines from specimen and set up FP calculations
        fpPrep(fpStorage, standard, conditions, pureLines );
        // Calculate spectrum for this standard, updating component spectra
        result = quantCalculate(fpStorage, standard, conditions, bkg_spec );
        if ( result != 0 ) {
            cout << "quantCalculate failed, result = " << result << endl;
            return -570 + result;
        };
        result = quantFitSpectrum( conditions, bkg_spec, cout );
        if ( result < 0 ) {
            cout << "quantFitSpectrum failed, result = " << result << endl;
            return -580 + result;
        } else if( result == 0 ) {
            done_SNIP = true;
        };
        cout.precision(4);
        //  Keep track of largest coefficient to adjust low energy end (which is not fit well under Rh L peaks)
        float max_value = 0;
        //  Use coefficients from fit to modify response
        cout << "Initial fit  to SNIP background    iter " << iterations_SNIP << "  chi sq " << bkg_spec.chisq() << endl;
        cout << "   Fit coefficients   ";
        for( ic=0; ic<bkg_spec.numberOfComponents(); ic++ ) {
            SpectrumComponent updated_component = bkg_spec.component( ic );
            if( updated_component.type != CONTINUUM ) continue;
            if( updated_component.ignore ) continue;
            if( ! updated_component.enabled ) continue;
            float coeff = updated_component.coefficient;
            cout << ",  " << coeff;
            if( updated_component.bkg_index >= 0 && updated_component.bkg_index < n_regions ) {
                if( coeff > 0 ) {
                    optic_values[updated_component.bkg_index] *= coeff;
                    if( optic_values[updated_component.bkg_index] > max_value ) max_value = optic_values[updated_component.bkg_index];
                } else {
                    //  For negative fit coefficient, reduce optic response in this region but don't make it negative
                    optic_values[updated_component.bkg_index] *= 0.3f;
                }
            }
        }
        //  Modify response at low energy to avoid underestimation due to Rl L peaks
//        optic_values[0] = max_value / 2;
//        optic_values[0] = optic_values[1] / 2;
        //  Modify response at low energy to fix calculated intensity for Na thru Cl
        optic_values[0] = optic_values[1] * ZERO_EN_OPTIC_MULTIPLIER;
        cout << endl;
        cout << "   Optic values   ";
        for( i=0; i<optic_values.size(); i++ ) cout << ",  " << optic_values[i];
        cout << endl;
        //  Perform a spline fit to the calculated response
        float initial_slope = ( optic_values[1] - optic_values[0] ) / ( optic_energies[1] - optic_energies[0] );
        spline( optic_energies, optic_values, initial_slope, 0, optic_derivatives);
        cout << "   Optic deriv       ";
        cout.setf( ios::scientific, ios::floatfield );
        for( i=0; i<optic_derivatives.size(); i++ ) cout << ",  " << optic_derivatives[i];
        cout.setf( ios::fixed, ios::floatfield );
        cout << endl;
        //  Make a new optic with the calculated response
        XrayOptic new_optic( optic_energies, optic_values, optic_derivatives );
        //  Put the new optic into the measurement conditions
        conditions.optic = new_optic;
    }
    //  Redo the calculation with the latest optic response function but without a fit
    fpPrep(fpStorage, standard, conditions, pureLines );
    // Calculate spectrum for this standard, updating component spectra
    result = quantCalculate(fpStorage, standard, conditions, bkg_spec );
    if ( result != 0 ) {
        cout << "quantCalculate failed, result = " << result << endl;
        return -570 + result;
    };
    //  Put the optic response function into the calculation slot to plot it
    vector <float> optic_resp_calc( bkg_spec.numberOfChannels() );
    for( i=0; i<optic_resp_calc.size(); i++ ) optic_resp_calc[i] = conditions.optic.CheckTransmission( bkg_spec.energy(i) );
    bkg_spec.calc( optic_resp_calc );
    bkg_spec.iterations( iterations_SNIP );
    //  Replace the standard spectrum with the SNIP spectrum (for checking this intermediate result)
//    cout << "*** Stopping after SNIP fit to optic response. ***" << endl;
    stdSpectrum = bkg_spec;
    return iterations_SNIP;

    //  Now fit the calculated spectrum to the full measured spectrum
    //      and use the background fit coefficients to adjust the optic response

   //  Set up components for the calculated spectrum
    vector <SpectrumComponent> components;
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
    //  Set up components for everything except background
    result = setupComponents( sourceLines, pureLines_nonMatrix, components );
    if( result < 0 ) {
        cout << "setupComponents failed, result is " << result << endl;
        return -540 + result;
    }
    //  Set up components for background fit (with multiple regions for background)
    result = makeComponents( CONTINUUM, pureLines_nonMatrix, components, n_regions );
    if( result < 0 ) {
        cout << "makeComponents failed, result is " << result << endl;
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
    //  Add the components to the spectrum object - update the background since it never changes
//    int ic;
    for( ic=0; ic<components.size(); ic++ ) {
        //  Leave out Compton lines from tube L edges (fit with extra La and Lb1 lines above)
        if( components[ic].type == COMPTON && components[ic].level == L ) continue;
//        if( components[ic].type == RAYLEIGH && components[ic].level == L ) continue;
        if( components[ic].type == La ) continue;
//        if( components[ic].type == Lb1 ) continue;
        components[ic].plot = true;
        stdSpectrum.add_component( components[ic] );
    }
    stdSpectrum.put_bkg_split( optic_energies );
    //  Fit the components to the measured spectrum (without changing the composition)
    int iterations = 0;
    bool done = false;
    const int MAX_ITERATIONS_MEAS = 1;
    while( iterations < MAX_ITERATIONS_MEAS && ( ! done )) {
        iterations++;
        // Calculate spectrum for this standard, updating component spectra
        result = quantCalculate(fpStorage, standard, conditions, stdSpectrum );
        if ( result != 0 ) {
            cout << "quantCalculate failed, result = " << result << endl;
            return -570 + result;
        };
        //  Also re-calculate the ignored elements
        result = quantIgnore( element_list, conditions, stdSpectrum, ignoreLines );
        if( result < 0 ) {
            cout << "quantIgnore failed to set up components for ignored elements, result is " << result << endl;
            return -540 + result;
        }
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
        //  Keep track of largest coefficient to adjust low energy end (which is not fit well under Rh L peaks)
        float max_value = 0;
        //  Use coefficients from fit to modify response, check for convergence
        cout << "New fit     iter " << iterations << "    chi sq " << stdSpectrum.chisq() << endl;
        cout << "   Fit coefficients   ";
        for( ic=0; ic<stdSpectrum.numberOfComponents(); ic++ ) {
            SpectrumComponent updated_component = stdSpectrum.component( ic );
            if( updated_component.type != CONTINUUM ) continue;
            if( updated_component.ignore ) continue;
            if( ! updated_component.enabled ) continue;
            float coeff = updated_component.coefficient;
            cout << ",  " << coeff;
            if( updated_component.bkg_index >= 0 && updated_component.bkg_index < n_regions ) {
                if( coeff > 0 ) {
                    optic_values[updated_component.bkg_index] *= coeff;
                    if( optic_values[updated_component.bkg_index] > max_value ) max_value = optic_values[updated_component.bkg_index];
                } else {
                    //  For negative fit coefficient, reduce optic response in this region but don't make it negative
                    optic_values[updated_component.bkg_index] *= 0.3f;
                }
           }
       }
        //  Modify response at low energy to avoid underestimation due to Rl L peaks
//        optic_values[0] = max_value / 2;
//        optic_values[0] = optic_values[1] / 2;
        //  Modify response at low energy to fix calculated intensity for Na thru Cl
        optic_values[0] = optic_values[1] * ZERO_EN_OPTIC_MULTIPLIER;
        cout << endl;
        cout << "   Optic values       ";
        for( i=0; i<optic_values.size(); i++ ) cout << ",  " << optic_values[i];
        cout << endl;
        //  Perform a spline fit to the calculated response
        float initial_slope = ( optic_values[1] - optic_values[0] ) / ( optic_energies[1] - optic_energies[0] );
        spline( optic_energies, optic_values, initial_slope, 0, optic_derivatives);
        cout << "   Optic deriv       ";
        cout.setf( ios::scientific, ios::floatfield );
        for( i=0; i<optic_derivatives.size(); i++ ) cout << ",  " << optic_derivatives[i];
        cout.setf( ios::fixed, ios::floatfield );
        cout << endl;
        //  Make a new optic with the calculated response
        XrayOptic new_optic( optic_energies, optic_values, optic_derivatives );
        //  Put the new optic into the measurement conditions
        conditions.optic = new_optic;
        //  Re-initialize FP calculations
        fpPrep(fpStorage, standard, conditions, pureLines );
        if( iterations < MINIMUM_ITERATIONS ) done = false;   //  Defined in XRFcontrols.h
        //  Don't disable negative components since they are taken care of in optic value adjustments
    }

    //  Add the optic response (without updating calculation)
    SpectrumComponent optic_component;
    optic_component.type = OPTIC_TRANS;
    optic_component.fit = false;
    optic_component.enabled = false;    //  Make sure it does not get included in the calculation, plot only
    optic_component.plot = true;
    optic_component.bkg = true; //  This makes it plot without adding the background
    optic_component.spectrum.resize( stdSpectrum.numberOfChannels(), 0 );
    unsigned int is;
    for( is=0; is<stdSpectrum.numberOfChannels(); is++ ) {
        float en = stdSpectrum.energy( is );
        optic_component.spectrum[is] = conditions.optic.CheckTransmission( en );
    }
    stdSpectrum.add_component( optic_component );
    stdSpectrum.iterations( iterations );
	return iterations;

};
