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

#include "quantIgnore.h"
#include "quantComponents.h"
#include "fpMain.h"
#include "fpLineSpectrum.h"
#include "XRFconstants.h"
#include "XRFcontrols.h"


using namespace std;

//  Process any element to be included in the fit but ignored in the composition
//  Added May 14, 2017
//  Modified Dec. 13, 02017
//      Add check for minimum energy to escape peaks (passed to fpLineSpectrum)
//  Modified July 10, 2021  Add simple pulse pileup calculation - change return for fpLineSpectrum (note ignore peaks not included in pileup)

int quantIgnore( const std::vector <ElementListEntry> element_list, XRFconditions &conditions,
                        XraySpectrum &spectrum, std::vector <XrayLines> &ignoreLines ) {
//		check input parameters
	if( ! spectrum.calibration().good() ) return -520;
	if( spectrum.live_time() <= 0 ) return -521;

    //  Set up components for the spectrum
    vector <SpectrumComponent> components;
    vector <XrayLines> pureLines;

    FPstorage fpStorage;

    //  Include components for any elements to be ignored
    int ie;
    for( ie=0; ie<element_list.size(); ie++ ) {
        if( element_list[ie].qualifier != IGNORE ) continue;
        //  Single-element material
        XrayMaterial temp_mat( element_list[ie].element );
        //  Use intensity of lines from pure element
        fpPrep(fpStorage, temp_mat, conditions, pureLines );
        int i;
        for( i=0; i<pureLines.size(); i++ ) pureLines[i].commonFactor( spectrum.live_time() );
        int result = makeComponents( ELEMENT, pureLines, components );
        if( result < 0 ) {
            cout << "makeComponents failed, result is " << result;
            cout << "  for ignored element " << element_list[ie].element.symbol() << endl;
            return -540 + result;
        }
        //  Make the initial calculation of the component spectrum (must be re-calculated if energy calibration changes during fit)
        int ic;
        for( ic=0; ic<components.size(); ic++ ) {
            components[ic].quant = false;
            components[ic].ignore = true;
            components[ic].spectrum.resize( spectrum.numberOfChannels() );
            int il = 0;
            for( il=0; il<pureLines.size(); il++ ) {
                vector <LineGroup> dummy;
                fpLineSpectrum( pureLines[il], conditions.detector, 1,
                   spectrum.calibration(), conditions.eMin, dummy, components[ic] );
            }
        }
        //  Return the XrayLines objects for use in re-calculating the ignore components
        int il;
        for( il=0; il<pureLines.size(); il++ ) ignoreLines.push_back( pureLines[il] );
    }

    int ic;
    for( ic=0; ic<components.size(); ic++ ) spectrum.add_component( components[ic] );

	return 0;

};
