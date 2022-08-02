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

#include <iostream>
#include "Element.h"
#include "XrayLines.h"
#include "XRFconstants.h"
#include "XRFcontrols.h"
#include "fpBeams.h"
#include "fpLineSpectrum.h"
#include "fpConvolve.h"
#include "quantPrimarySpec.h"

using namespace std;

//		Perform fundamental parameters calculation of the
//		predicted measured primary spectrum from an
//		X-ray source plus anything in the primary beam

// Modified from stdCalcSpec.cpp of June 1, 2015
//	Plus code from fpExcitation.cpp, fpContScat.cpp, and fpLineScat.cpp
//  Re-written Feb. 15, 2017 from primaryCalcSpec.cpp and based on quantCalculate.cpp
//      For PIQUANT Version 2, using new conditions and XraySpectrum class
//  Modified April 12, 2017 to fix bugs in incident beam corrections for lines
//  Modified Dec. 13, 02017
//      Add check for minimum energy to escape peaks (passed to fpLineSpectrum)
//  Modified Dec. 7, 2020
//      Added detector shelf calculation from photoelectron and Auger electron escape (active volume and front contact)
//  Modified Dec. 29, 2020
//      Add optic response to plot output of primary spectrum calculation
//  Modified Apr. 6, 2021   Add CONTINUUM component type and sort out how to handle background and Compton escape (remove Det shelf component)
//                          Major rearrangements in fpLineSpectrum to improve speed, move shelf calc from quantCalculate to fpLineSpectrum
//  Modified July 10, 2021  Add simple pulse pileup calculation - change return for fpLineSpectrum (note these peaks not included in pileup)


int quantPrimarySpec( const XRFconditions &conditions_in, XraySpectrum &primary_spectrum ) {
//		check input parameters
	if( primary_spectrum.numberOfChannels() <= 0 ) return -701;
	if( ! primary_spectrum.calibration().good() ) return -705;
	if( primary_spectrum.live_time() <= 0 ) return -706;
	int nChan = primary_spectrum.numberOfChannels();

    //**************************************************************************
    //      calculate contribution to spectrum from characteristic lines
    //**************************************************************************

//		get list with intensities of tube characteristic lines
	vector <XrayLines> sourceLines;
	conditions_in.source.lines( sourceLines, conditions_in.eMin );
//    apply incident beam corrections and detector response
	int edgeIndex;
	for ( edgeIndex=0; edgeIndex<sourceLines.size(); edgeIndex++ ) {
		int lineIndex;
		for ( lineIndex=0; lineIndex<sourceLines[edgeIndex].numberOfLines(); lineIndex++ ) {
//				calculate intensity for each line
			float lineEn = sourceLines[edgeIndex].energy( lineIndex );
			float lineInt = sourceLines[edgeIndex].factor( lineIndex );
//			    apply incident beam corrections
            lineInt *= fpIncidentBeam ( lineEn, conditions_in );
//		        apply detector response correction
            lineInt *= conditions_in.detector.response( lineEn );
//              save corrected intensity factor in line object
            sourceLines[edgeIndex].factor( lineIndex, lineInt );
        };
	};
//			correct for spectrum live time
    int i;
    for( i=0; i<sourceLines.size(); i++ ) sourceLines[i].commonFactor( primary_spectrum.live_time() );
/*    //  Save strongest line intensity for detector shelf calculation
    float line_en;
    float line_int = 0;
    for( i=0; i<sourceLines.size(); i++ ) {
        unsigned int il;
        for( il=0; il<sourceLines[i].numberOfLines(); il++ ) {
            if( sourceLines[i].intensity(il) > line_int ) {
                line_en = sourceLines[i].energy(il);
                line_int = sourceLines[i].intensity(il);
            }
        }
    }
*/    //   Put all of the necessary components into the spectrum object
    vector <SpectrumComponent> components;
    //  Put in components from list of source emission lines
    int result = makeComponents( PRIMARY_LINES, sourceLines, components );
    if( result < 0 ) {
        cout << "Primary lines makeComponents failed, result is " << result << endl;
        return result;
    }
	//cout << "Starting calculation for characteristic lines" << endl;
    //  Calculate the contribution to the spectrum from each component in the XraySpectrum object
    int ic;
    for( ic=0; ic<components.size(); ic++ ) {
        if( components[ic].type != PRIMARY_LINES ) continue;
        components[ic].spectrum.resize( nChan ,0 );
        int il;
        for ( il=0; il<sourceLines.size(); il++ ) {
            if ( sourceLines[il].numberOfLines() <= 0 ) continue;
    //			get approximate energy for detector resolution and bkg noise threshold
            float threshold = 1;
            vector <LineGroup> dummy;
            fpLineSpectrum( sourceLines[il], conditions_in.detector, threshold, primary_spectrum.calibration(), conditions_in.eMin, dummy, components[ic] );
        };
        //  Put the new calculation into the XraySpectrum object
        primary_spectrum.add_component( components[ic] );
    }

    //**************************************************************************
    //      calculate contribution to spectrum from continuum (if any)
    //**************************************************************************

	//cout << "Starting continuum calculation." << endl;
	if( conditions_in.source.continuum() ) {
        //  Put in component for source continuum
        vector <SpectrumComponent> continuum_component;
        vector <XrayLines> dummy_lines;
        result = makeComponents( PRIMARY_CONTINUUM, dummy_lines, continuum_component );
        if( result < 0 || continuum_component.size() < 1 ) {
            cout << "Primary continuum makeComponents failed, result is " << result << endl;
            return result;
        }
        //  Put in a loop in case continuum is broken up into more that one component in the future
        for( ic=0; ic<continuum_component.size(); ic++ ) {
            if( continuum_component[ic].type != PRIMARY_CONTINUUM ) continue;
            continuum_component[ic].spectrum.resize( nChan );
            int iChan;
            for( iChan=0; iChan<nChan; iChan++ ) {
                float contEn = primary_spectrum.energy( iChan );
                //			find continuum intensity at desired energy
                float contInt = conditions_in.source.continuum( contEn );
    //			    apply incident beam corrections
                contInt *= fpIncidentBeam ( contEn, conditions_in );
    //		        apply detector response correction
                float detResp = conditions_in.detector.response( contEn );
                contInt *= detResp;
                //			result is per keV, so multiply by channel width in keV to get counts in each channel
                contInt *= primary_spectrum.calibration().energyPerChannel( iChan ) / 1000;
                contInt *= primary_spectrum.live_time();
                continuum_component[ic].spectrum[iChan] = contInt;
            };
            //  Convolve the continuum with the detector broadening
            //  This was moved to after the Compton escape calculation
//            fpConvolve( conditions_in.detector, primary_spectrum.calibration(), continuum_component[ic].spectrum );
            //  Put the new calculation into the XraySpectrum object
            primary_spectrum.add_component( continuum_component[ic] );
        };
    };  //  if( conditions_in.source.continuum() )

    primary_spectrum.update_calc();

    //  Find the primary_spectrum component for the continuum
    for( ic=0; ic<primary_spectrum.numberOfComponents(); ic++ ) {
        SpectrumComponent ce_component = primary_spectrum.component( ic );
        if( ce_component.type != DETECTOR_CE ) continue;
        //  Calculate Compton escape shelf at low energies
        //  Loop over channels in the continuum and add all Compton escape contributions there
        unsigned int i_ce;
        for( i_ce=0; i_ce<ce_component.spectrum.size(); i_ce++ ) {
            float spec_energy = primary_spectrum.energy( i_ce );
            if( spec_energy < conditions_in.eMin ) continue;
            //  Check if Compton escape is possible for this channel (or any higher channels)
            float min_ce_energy = conditions_in.detector.ce_minimum( spec_energy );
            if( min_ce_energy > conditions_in.source.kV() * 1000 ) break;
            unsigned int min_ce_channel = primary_spectrum.channel( min_ce_energy );
            if( min_ce_channel >= primary_spectrum.calc().size() - 1 ) break;
            unsigned int is;
            for( is=min_ce_channel; is<primary_spectrum.calc().size(); is++ ) {
                float meas_intensity = primary_spectrum.calc()[is];
                if( meas_intensity <= 0 ) continue;
                float inc_energy = primary_spectrum.energy( is );
                //  Find original intensity incident on detector by dividing by response at this energy
                float det_resp = conditions_in.detector.response( inc_energy );
                if( det_resp <= 0 ) continue;
                float eV_ch = primary_spectrum.calibration().energyPerChannel( is );
                //  Get original intensity incident on detector
                float incoming_int =  meas_intensity / det_resp / eV_ch;
                //  Compton escape for this spectrum channel from the incident energy
                float ce_intensity = incoming_int * conditions_in.detector.ce_fraction( inc_energy, spec_energy );
                //  Add the Compton escape intensity to the background channel
                ce_component.spectrum[i_ce] += ce_intensity;
            }
        }
        //  Convolve the continuum with the detector broadening
        fpConvolve( conditions_in.detector, primary_spectrum.calibration(), ce_component.spectrum );
        primary_spectrum.update_component( ce_component );
        break;
    }
    //  Put the new background into the calculation
    primary_spectrum.update_calc();
    //  Add the optic response (disable so it won't be included in calculation)
    SpectrumComponent optic_component;
    optic_component.type = OPTIC_TRANS;
    optic_component.fit = false;
    optic_component.enabled = false;
    optic_component.spectrum.resize( primary_spectrum.numberOfChannels(), 0 );
    unsigned int is;
    for( is=0; is<primary_spectrum.numberOfChannels(); is++ ) {
        float en = primary_spectrum.energy( is );
        optic_component.spectrum[is] = conditions_in.optic.CheckTransmission( en );
    }
    primary_spectrum.add_component( optic_component );


	return 0;

};

