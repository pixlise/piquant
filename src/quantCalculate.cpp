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

#include "quantCalculate.h"
#include "fpLineSpectrum.h"
#include "fpMain.h"
#include "fpConvolve.h"
#include "Element.h"
#include "XrayEdge.h"
#include "XrayLines.h"
#include "XRFconstants.h"
#include "XRFcontrols.h"
#include "quantComponents.h"
#include "split_component.h"
#include "scale_under_peaks.h"


using namespace std;

//		Perform fundamental parameters calculation for a standard sample
//			and calculate theoretical XRF spectrum

//  Modified June 1, 2014
//      Add Si and Fe calculated intensity printout in stdCalcSpec.cpp
//  Re-written Feb. 2, 2017 from stdCalcSpec.cpp
//      Use XrayMaterial class for specimen composition and thickness
//      Use new conditions structure and setup for fp calculations
//      Use XraySpectrum class for calculated spectrum (energy calibration + fit and bkg)
//  Modified Feb. 10, 2017
//      Use SpectrumComponent to calculate parts of spectrum (defined in fpComponents)
//  Modified Sept. 29, 2017
//`     Warn if calculated intensity is zero, negative, or nan and disable (only if not already disabled)
//  Modified Dec. 23, 2017
//      Error message and return if quant component has zero intensity, warning only for other components
//  Modified Sep. 22, 2020
//      Add simple Compton escape calculation to test algorithm
//  Modified Oct. 21, 2020
//      Put Compton escape calculation in separate component if that type component is available
//      Split up background into multiple components (to be fit separately)
//  Modified Dec. 4, 2020   Use SNIP background to fit anomalous background at low energies
//  Modified Dec. 9, 2020   Add detector electron loss shelf calculation
//  Modified Dec. 28, 2020  Reset coefficients to unity before shelf calculations
//                          Use SNIP for low-energy background, disable shelf calculations
//  Modified Feb. 26, 2021  Adjust the shape of the calculated background to generally match measured spectra using ramp multiplier
//                          Adjust the overall intensity to match this measured spectrum using an overall factor calculated via least squares with peaks blocked from sums
//  Modified Feb. 26, 2021  Adjust the shape of the calculated background using spline fit to Teflon scatter (with new unity ECF optic)
//  Modified Mar. 17, 2021  Try truncated shelf from PIXE empirical detector models
//  Modified Mar. 18, 2021  Abandon truncated shelf, back to Scholze & Procop shelf model, add adjustments vs photon energy and slope vs electron loss energy (non-flat shelf)
//  Modified Mar. 18, 2021  Shelf overall factor of 30 and overall slope vs electron loss energy
//  Modified Apr. 2, 2021   Remove shelf calculations (moved to fpLineSpectrum)
//  Modified Apr. 6, 2021   Include calculated continuum bkg as a separate component
//  Modified Apr.11, 2021   Use coefficient to adjust calc bkg, to see how far off cal is from measurement
//  Modified Apr.28, 2021   UFix bug that changed calc bkg coeff when not adjusted (replacing multiplier in bkg parameters)
//  Modified May 10, 2021   Use scale_under_peaks function when scale factor is negative in arguments
//  Modified May 14, 2021   Fixed logic error where coefficient was reset to unity when it should be manual scale factor
//  Modified July 10, 2021  Add simple pulse pileup calculation


const vector<float> X_BkgAdj;
const vector<float> Y_BkgAdj;
const vector<float> D_BkgAdj;
//const vector<float> X_BkgAdj = {        0,  4000,      6000,    8000,   10000,   12000,   14000,   16000,   18000,   20000,   25000,   30000 };
//const vector<float> Y_BkgAdj = {        0,  0.9934,  0.99,  0.9929,  0.99,  0.9909,  1.0121,  1.0320,  1.2285,  1.3594,  0.8755,  0.0880 };
//const vector<float> D_BkgAdj = {         6.4504e-08,  -1.2901e-07,  1.2809e-07,  -9.5202e-08,  4.2017e-08,  1.2952e-08,  -3.5783e-08,  8.3143e-08,  -3.2014e-08,  -5.3319e-08,  -3.2588e-08,  1.1080e-07 };
//const vector<float> Y_BkgAdj = {        1,     1,       1,       1,       1,     1,       1,       1,       1,       1,     1,   1 };
//const vector<float> D_BkgAdj( X_BkgAdj.size(), 0 );
//const vector<float> X_BkgAdj = {        0,  2000,   4000,    6000,    8000,   10000,   12000,   14000,   16000,   18000,   20000,   25000,   30000 };
//const vector<float> Y_BkgAdj = {        10,     10,      1,    0.80,  0.9929,  0.9483,  0.9609,  1.0121,  1.0320,  1.2285,  1.3594,  0.8755,  0.8755 };
//const vector<float> D_BkgAdj( X_BkgAdj.size(), 0 );


int quantCalculate(const FPstorage &fpStorage, const XrayMaterial &specimen, const XRFconditions &conditions_in,
            XraySpectrum &spectrum ) {
//		check input parameters
	if( spectrum.numberOfChannels() <= 0 ) return -701;
	if( ! spectrum.calibration().good() ) return -705;
	if( spectrum.live_time() <= 0 ) return -706;
	int nChan = spectrum.numberOfChannels();
    int i;

//			generate calculated emission line intensities for all elements using this sample composition
//cout << "Starting fp calc." << endl;
	vector <XrayLines> sampleLines(0);
	fpCalc (fpStorage, specimen, conditions_in, sampleLines );
//			correct for spectrum live time
    for( i=0; i<sampleLines.size(); i++ ) sampleLines[i].commonFactor( spectrum.live_time() );

    //  See if the background should be calculated by looking for a continuum component (also find Compton escape if any)
    unsigned int ic;
    int i_bkg_component = -1;
    int index_ce = -1;
    int index_pileup = -1;
    float sigma_mult = 0;
    for( ic=0; ic<spectrum.numberOfComponents(); ic++ ) {
        SpectrumComponent updated_component = spectrum.component( ic );
        if( index_ce < 0 && updated_component.type == DETECTOR_CE ) index_ce = ic;
        if( index_pileup < 0 && updated_component.type == PILEUP ) index_pileup = ic;
        if( updated_component.type != CONTINUUM ) continue;
        if( i_bkg_component < 0 ) {
            i_bkg_component = ic; //  Save for use below if only one component
            sigma_mult = updated_component.scale_under;
        }
    }
   //  Get any background crossover parameters stored with the spectrum and see if this is a defauly bkg calculation
    vector <float> bkg_split_energies;
    spectrum.get_bkg_split( bkg_split_energies );

    //  List of grouped lines for pulse pileup calculation
    vector <LineGroup> simple_pileup_list;

    //		calculate continuum background (if desired)
    if( i_bkg_component >= 0 ) {
//cout << "Starting background calculation." << endl;
        vector <float> temp_bkg( spectrum.numberOfChannels() );
        fpContScat(fpStorage, spectrum.calibration(), specimen, conditions_in, temp_bkg );
        //			correct for spectrum live time
        for( i=0; i<temp_bkg.size(); i++ ) temp_bkg[i] *= spectrum.live_time();
        if( spectrum.convolve_Compton() ) fpConvolve( conditions_in.detector, spectrum.calibration(), temp_bkg );
        //  Adjust the shape of the calculated background using spline fit to Teflon scatter (with new unity ECF optic)
        if( X_BkgAdj.size() > 0 ) for( i=0; i<temp_bkg.size(); i++ ) temp_bkg[i] *= splint( X_BkgAdj, Y_BkgAdj, D_BkgAdj, spectrum.energy(i) );
        float bkg_factor = 1;
        //  Adjust the overall intensity to match measured spectrum if desired (returns unity if measured spectrum is zero size)
        if( sigma_mult > 0 ) bkg_factor = scale_under_peaks( temp_bkg, spectrum.meas(), spectrum.sigma(), sigma_mult );
        if( bkg_split_energies.size() > 0 ) {
            //  Split up calculated background if desired (used for optic response and crossover background)
            for( ic=0; ic<spectrum.numberOfComponents(); ic++ ) {
                SpectrumComponent updated_component = spectrum.component( ic );
                if( updated_component.type != CONTINUUM ) continue;
                updated_component.spectrum.resize( temp_bkg.size(), 0 );
                for( i=0; i<temp_bkg.size(); i++ ) {
                    float e = spectrum.energy( i );
                    float split = split_weight( e, bkg_split_energies, updated_component.bkg_index );
                    updated_component.spectrum[i] = temp_bkg[i] * split;
                }
                //  Put the new calculation into the XraySpectrum object
                spectrum.update_component( updated_component );
                //  Put factor from adjustment above into coefficient
                if( sigma_mult > 0 ) spectrum.update_coefficient( ic, bkg_factor );
            }
        } else {
            //  Put the single-component calculated background into the XraySpectrum object
            SpectrumComponent updated_component = spectrum.component( i_bkg_component );
            updated_component.spectrum.resize( temp_bkg.size(), 0 );
            for( i=0; i<temp_bkg.size(); i++ ) updated_component.spectrum[i] = temp_bkg[i];
            spectrum.update_component( updated_component );
            //  Put factor from adjustment above into coefficient
            if( sigma_mult > 0 ) spectrum.update_coefficient( i_bkg_component, bkg_factor );
        }
    }
    //  Get the background into the spectrum as a sum of the bkg components
    spectrum.update_calc();


//cout << "Starting spectrum peak calculation." << endl;
//		calculate the peaks from characteristic emission lines
    //  Calculate the contribution to the spectrum from each component in the XraySpectrum object
    //  Only process components corresponding to elements in the specimen (not ignore or other extra components)
    vector <Element> specimen_elements = specimen.element_list();
    for( ic=0; ic<spectrum.numberOfComponents(); ic++ ) {
        SpectrumComponent updated_component = spectrum.component( ic );
        if( updated_component.type != ELEMENT ) continue;
        if( ! updated_component.enabled ) continue;
        bool found = false;
        int ie;
        for( ie=0; ie<specimen_elements.size(); ie++ ) {
            if( specimen_elements[ie] == updated_component.element ) {
                found = true;
                break;
            }
        }
        if( ! found ) continue;
        updated_component.spectrum.resize( spectrum.numberOfChannels() ,0 );
        //  resize does not necessarily set all values to zero (if vector is already the right size)
        int i;
        for( i=0; i<updated_component.spectrum.size(); i++ ) updated_component.spectrum[i] = 0;
        int il;
        for ( il=0; il<sampleLines.size(); il++ ) {
            if ( sampleLines[il].numberOfLines() <= 0 ) continue;
    //			get approximate energy for detector resolution and bkg noise threshold
            float en = sampleLines[il].energy(0);
    //			get background noise for threshold
            int k = spectrum.channel( en );
            float threshold = 1;
            if ( k >= 0 && k < nChan && spectrum.bkg()[k] > 0 ) threshold = 0.1f * sqrt( spectrum.bkg()[k] );
            fpLineSpectrum( sampleLines[il], conditions_in.detector, threshold, spectrum.calibration(), conditions_in.eMin, simple_pileup_list, updated_component );
        };
        //  Check for zero (or nan) and disable (also write message)
        //  Only if not already disabled to avoid many messages (check is at top of loop)
        float sum = 0;
        for( i=0; i<updated_component.spectrum.size(); i++ ) sum += updated_component.spectrum[i];
//        cout << "fpLineSpectrum " << componentDescription( updated_component ) << "  " << sum << endl;
        if( updated_component.quant && ( sum <= 0 || isnan( sum ) ) ) {
            cout << "*** Error - calculated intensity is zero (or negative or nan) for";
            cout << " component " << componentDescription( updated_component ) << "  " << sum << endl;
            return -710;
        } else if( sum <= 0 || isnan( sum ) ) {
            cout << "*** Warning - calculated intensity is zero (or negative or nan) for";
            cout << " component " << componentDescription( updated_component );
            cout << " (it is being disabled).   " << sum << endl;
            spectrum.disable( ic );
        }
        //  Put the new calculation into the XraySpectrum object
        spectrum.update_component( updated_component );
    }

//cout << "Starting scatter peak calculation." << endl;


//		calculate the peaks from Rayleigh scatter of the source characteristic lines
	vector <XrayLines> scatterLines(0);
    fpRayleigh(fpStorage, specimen, conditions_in, scatterLines );
//			correct for spectrum live time
    for( i=0; i<scatterLines.size(); i++ ) scatterLines[i].commonFactor( spectrum.live_time() );

//		rearrange the source lines (to debug extra La or Lb1 intensity)
    bool enable_La = false;
    bool enable_Lb1 = false;
    //  Check to see if components for these extra lines were included
    for( ic=0; ic<spectrum.numberOfComponents(); ic++ ) {
        if( spectrum.component( ic ).type == La ) enable_La = true;
        if( spectrum.component( ic ).type == Lb1 ) enable_Lb1 = true;
    }
    XrayLines scatterLines_La;
    XrayLines scatterLines_Lb1;
    for( i=0; i<scatterLines.size(); i++ ) {
        if( enable_La && scatterLines[i].edge().index() == L3 ) {
            scatterLines_La = scatterLines[i];
            int li;
            for( li=0; li<scatterLines_La.numberOfLines(); li++ ) {
                //  Separate L alpha lines (L3-M4,5)
                if( scatterLines[i].edgeSource(li).index() == M4 || scatterLines[i].edgeSource(li).index() == M5 ) {
                    scatterLines[i].factor( li, 0 );
                } else {
                    scatterLines_La.factor( li, 0 );
                }
            }
        } else if( enable_Lb1 && scatterLines[i].edge().index() == L2 ) {
            scatterLines_Lb1 = scatterLines[i];
            int li;
            for( li=0; li<scatterLines_Lb1.numberOfLines(); li++ ) {
                //  Separate L beta 1 line (L2-M4)
                if( scatterLines[i].edgeSource(li).index() == M4 ) {
                    scatterLines[i].factor( li, 0 );
                } else {
                    scatterLines_Lb1.factor( li, 0 );
                }
            }
        }
    }

//cout << "Starting Rayleigh peak calculation." << endl;
    //  Calculate the contribution to the spectrum from each Rayleigh scatter component in the XraySpectrum object
    //int ic;
    for( ic=0; ic<spectrum.numberOfComponents(); ic++ ) {
        SpectrumComponent updated_component = spectrum.component( ic );
        if( updated_component.type != RAYLEIGH ) continue;
        updated_component.spectrum.resize( spectrum.numberOfChannels() ,0 );
        //  resize does not necessarily set all values to zero (if vector is already the right size)
        int i;
        for( i=0; i<updated_component.spectrum.size(); i++ ) updated_component.spectrum[i] = 0;
        int il;
        for ( il=0; il<scatterLines.size(); il++ ) {
            if ( scatterLines[il].numberOfLines() <= 0 ) continue;
    //			get approximate energy for detector resolution and bkg noise threshold
            float en = scatterLines[il].energy(0);
    //			get background noise for threshold
            int k = spectrum.channel( en );
            float threshold = 1;
            if ( k >= 0 && k < nChan && spectrum.bkg()[k] > 0 ) threshold = 0.1f * sqrt( spectrum.bkg()[k] );
            fpLineSpectrum( scatterLines[il], conditions_in.detector, threshold, spectrum.calibration(), conditions_in.eMin, simple_pileup_list, updated_component );
        }
        //  Put the new calculation into the XraySpectrum object
        spectrum.update_component( updated_component );
    }

//cout << "Starting Compton peak calculation." << endl;
	//		calculate the contribution from Compton scatter of the source characteristic lines
    //int ic;
    for( ic=0; ic<spectrum.numberOfComponents(); ic++ ) {
        SpectrumComponent updated_component = spectrum.component( ic );
        if( updated_component.type != COMPTON ) continue;
        updated_component.spectrum.resize( spectrum.numberOfChannels() ,0 );
        //  resize does not necessarily set all values to zero (if vector is already the right size)
        int i;
        for( i=0; i<updated_component.spectrum.size(); i++ ) updated_component.spectrum[i] = 0;
        fpCompton(fpStorage, spectrum.calibration(), specimen,
					conditions_in, updated_component );
        //	correct for spectrum live time
        for( i=0; i<nChan; i++ ) updated_component.spectrum[i] *= spectrum.live_time();
        if( spectrum.convolve_Compton() ) fpConvolve( conditions_in.detector, spectrum.calibration(), updated_component.spectrum );
        spectrum.update_component( updated_component );
    }

//cout << "Starting extra peak calculation." << endl;
	//		calculate the contribution from the extra individual source lines separated above (to debug extra La or Lb1 intensity)
    //int ic;
    for( ic=0; ic<spectrum.numberOfComponents(); ic++ ) {
        SpectrumComponent updated_component = spectrum.component( ic );
        if( updated_component.type == La ) {
            updated_component.spectrum.resize( spectrum.numberOfChannels() ,0 );
            //  resize does not necessarily set all values to zero (if vector is already the right size)
            for( i=0; i<updated_component.spectrum.size(); i++ ) updated_component.spectrum[i] = 0;
            if ( scatterLines_La.numberOfLines() <= 0 ) continue;
    //			get approximate energy for detector resolution and bkg noise threshold
            float en = scatterLines_La.energy(0);
    //			get background noise for threshold
            int k = spectrum.channel( en );
            float threshold = 1;
            if ( k >= 0 && k < nChan && spectrum.bkg()[k] > 0 ) threshold = 0.1f * sqrt( spectrum.bkg()[k] );
            fpLineSpectrum( scatterLines_La, conditions_in.detector, threshold, spectrum.calibration(), conditions_in.eMin, simple_pileup_list, updated_component );
            spectrum.update_component( updated_component );
        } else if( updated_component.type == Lb1 ) {
            updated_component.spectrum.resize( spectrum.numberOfChannels() ,0 );
            //  resize does not necessarily set all values to zero (if vector is already the right size)
            for( i=0; i<updated_component.spectrum.size(); i++ ) updated_component.spectrum[i] = 0;

            if ( scatterLines_Lb1.numberOfLines() <= 0 ) continue;
    //			get approximate energy for detector resolution and bkg noise threshold
            float en = scatterLines_Lb1.energy(0);
    //			get background noise for threshold
            int k = spectrum.channel( en );
            float threshold = 1;
            if ( k >= 0 && k < nChan && spectrum.bkg()[k] > 0 ) threshold = 0.1f * sqrt( spectrum.bkg()[k] );
            fpLineSpectrum( scatterLines_Lb1, conditions_in.detector, threshold, spectrum.calibration(), conditions_in.eMin, simple_pileup_list, updated_component );
            spectrum.update_component( updated_component );
        }
    }

    //  Detector shelf calculations from Compton escape (electron loss shelf calculated in fpLineSpectrum

//cout << "Starting Compton escape calculation." << endl;
    spectrum.adjust_coefficients(); //  Adjust coefficients to better match new composition for shelf calculation
    spectrum.update_calc(); //  Get all photons incident on detector into shelf calculation, from the new calculations above
    vector <float> ce_calc( spectrum.numberOfChannels(), 0 );
    if( index_ce >= 0 ) {
        //  Calculate Compton escape shelf at low energies
        //  Loop over channels in the full calculation and add Compton escape contribution
        unsigned int i_ce;
        for( i_ce=0; i_ce<ce_calc.size(); i_ce++ ) {
            float spec_energy = spectrum.energy( i_ce );
            if( spec_energy < conditions_in.eMin ) continue;
            //  Check if Compton escape is possible for this channel (or any higher channels)
            float min_ce_energy = conditions_in.detector.ce_minimum( spec_energy );
            if( min_ce_energy > conditions_in.source.kV() * 1000 ) break;
            unsigned int min_ce_channel = spectrum.channel( min_ce_energy );
            if( min_ce_channel >= spectrum.calc().size() - 1 ) break;
            unsigned int is;
            for( is=min_ce_channel; is<spectrum.calc().size(); is++ ) {
                float meas_intensity = spectrum.calc()[is];
                if( meas_intensity <= 0 ) continue;
                float inc_energy = spectrum.energy( is );
                //  Find original intensity incident on detector by dividing by response at this energy
                float det_resp = conditions_in.detector.response( inc_energy );
                if( det_resp <= 0 ) continue;
                float incoming_int =  meas_intensity / det_resp;
                //  Compton escape for this spectrum channel from the incident energy
                float ce_intensity = incoming_int * conditions_in.detector.ce_fraction( inc_energy, spec_energy );
                //  Add the Compton escape intensity to the background channel
                ce_calc[i_ce] += ce_intensity;
                //  Include the same intensity in the channel that was the source of the Compton escape, for debugging
                //ce_calc[is] += ce_intensity;
            }
        }
        fpConvolve( conditions_in.detector, spectrum.calibration(), ce_calc );
        //  Add Compton escape to its component
        SpectrumComponent bkg_component = spectrum.component( index_ce );
        bkg_component.spectrum.resize( spectrum.numberOfChannels() );
        for( i_ce=0; i_ce<ce_calc.size(); i_ce++ ) bkg_component.spectrum[i_ce] = ce_calc[i_ce];
        spectrum.update_component( bkg_component );
    }
    spectrum.update_calc();

    //  Pulse pileup calculation using line group list from fpLineSpectrum

    vector <float> pileup_calc( spectrum.numberOfChannels(), 0 );
    if( index_pileup >= 0 ) {
//cout << "Starting pulse pileup calculation." << endl;
        //  Simple pileup calculation is just product of intensities times pulse resolving time divided by live time
        float resolving_time = conditions_in.detector.pileup_time();
        float pileup_factor = resolving_time / spectrum.live_time();
        //  Loop over line group list in nested loops to get all combinations
        unsigned int ig;
        for( ig=0; ig<simple_pileup_list.size(); ig++ ) {
            float line1_energy = simple_pileup_list[ig].energy;
            float line1_intensity = simple_pileup_list[ig].intensity;
            if( line1_energy < conditions_in.eMin ) continue;
            if( line1_intensity <= 0 ) continue;
            unsigned int ig2;
            for( ig2=0; ig2<simple_pileup_list.size(); ig2++ ) {
                float line2_energy = simple_pileup_list[ig2].energy;
                float line2_intensity = simple_pileup_list[ig2].intensity;
                if( line2_energy < conditions_in.eMin ) continue;
                if( line2_intensity <= 0 ) continue;
                float pileup_energy = line1_energy + line2_energy;
                float pileup_intensity = line1_intensity * line2_intensity * pileup_factor;
                unsigned int pileup_ch1 = spectrum.channel( pileup_energy );
                if( pileup_ch1 + 1 >= spectrum.numberOfChannels() ) continue;
                unsigned int pileup_ch2 = pileup_ch1 + 1;
                //  Make sure the pileup energy is between the two channels
                if( spectrum.energy( pileup_ch1 ) > pileup_energy ) {
                    pileup_ch2 = pileup_ch1;
                    pileup_ch1 = pileup_ch2 - 1;
                }
                //  Place the pileup intensity in two channels proportionally to get peak in right place
                float ch1_energy = spectrum.energy( pileup_ch1 );
                float ch2_energy = spectrum.energy( pileup_ch2 );
                float delta_energy = ch2_energy - ch1_energy;
                if( delta_energy == 0 ) continue;
                float split_intensity = pileup_intensity / delta_energy;
                pileup_calc[pileup_ch1] += split_intensity * ( pileup_energy - ch1_energy );
                pileup_calc[pileup_ch2] += split_intensity * ( ch2_energy - pileup_energy );
            }
        }
        //  Broaden the peaks by the appropriate Gaussian
        fpConvolve( conditions_in.detector, spectrum.calibration(), pileup_calc );
        //  Add the result to the pileup component
        SpectrumComponent pileup_component = spectrum.component( index_pileup );
        pileup_component.spectrum.resize( spectrum.numberOfChannels() );
        for( ig=0; ig<pileup_calc.size(); ig++ ) pileup_component.spectrum[ig] = pileup_calc[ig];
        spectrum.update_component( pileup_component );
    }
    spectrum.update_calc();

	return 0;
}
