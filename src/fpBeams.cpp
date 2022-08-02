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
//  fpBeams.cpp
//  PIQUANT
//
//  Created by W. T. Elam on 1/28/2017.
//  Copyright (c) 2017 APL/UW. All rights reserved.
//
#include "Element.h"
#include "XRFconstants.h"
#include "fpBeams.h"


//  Calculates transmission for things in the incident and emergent X-ray beams
//	Written Jan. 28, 2017 from code in fpMain.cpp, mostly fpContScat function
//      Extensive re-write to use new XRFconditions and XrayMaterial class
//      also add extra information from ISO standard version of EMSA format


using namespace std;

float fpIncidentBeam( const float energy, const XRFconditions &conditions ) {
    float factor = 1;
    //			apply incident beam filtration
    if( conditions.filter.thickness() > 0 )
        factor *= conditions.filter.transmission( energy );
    // Apply Optic transmission function
    if ( !conditions.optic.DefaultCheck() )
        factor *= conditions.optic.CheckTransmission( energy );
    //  Apply effect of dust on optic
    if( conditions.dust_on_optic.thickness() > 0 )
        factor *= conditions.dust_on_optic.transmission( energy );
    //  Apply effect of atmosphere path
    if( conditions.incidentPath.thickness() > 0 )
        factor *= conditions.incidentPath.transmission( energy );
    //  Account for source solid angle
    if( conditions.solidAngleSource > 0 )
        factor *= conditions.solidAngleSource;  //  Note that this is not steradians, converted in fpSetupConditions.cpp
    //  Apply effect of dust on specimen (to incoming beam)
    if( conditions.dust_on_specimen.thickness() > 0 )
        factor *= conditions.dust_on_specimen.transmission( energy, conditions.excitCosecant );
    //  Apply effect of window in front of specimen (to incoming beam)
    if( conditions.window.thickness() > 0 )
        factor *= conditions.window.transmission( energy, conditions.excitCosecant );
    return factor;
};

float fpEmergentBeam( const float energy, const XRFconditions &conditions ) {
    float factor = 1;
    //  Account for geometric factor
    if( conditions.geometryFactor > 0 )
        factor *= conditions.geometryFactor;
    //  Apply effect of dust on specimen (to outgoing beam)
    if( conditions.dust_on_specimen.thickness() > 0 )
        factor *= conditions.dust_on_specimen.transmission( energy, conditions.emergCosecant );
    //  Apply effect of window in front of specimen (to outgoing beam)
    if( conditions.window.thickness() > 0 )
        factor *= conditions.window.transmission( energy, conditions.emergCosecant );
    //  Apply effect of atmosphere path
    if( conditions.emergentPath.thickness() > 0 )
        factor *= conditions.emergentPath.transmission( energy );
    //  Apply effect of dust on detector
    if( conditions.dust_on_detector.thickness() > 0 )
        factor *= conditions.dust_on_detector.transmission( energy );
    //  Account for detector solid angle
    if( conditions.solidAngleDetector > 0 )
        factor *= conditions.solidAngleDetector;  //  Note that this is not steradians, converted in fpSetupConditions.cpp
    return factor;
};

void fpIncidentBeam( const XRFconditions &conditions, const std::vector <float> &energies,
                    std::vector <float> &intensities ) {
    int i;
    for( i=0; i<energies.size(); i++ ) intensities[i] *= fpIncidentBeam( energies[i], conditions );
    return;
};

void fpEmergentBeam( const XRFconditions &conditions, const std::vector <float> &energies,
                    std::vector <float> &intensities ) {
    int i;
    for( i=0; i<energies.size(); i++ ) intensities[i] *= fpEmergentBeam( energies[i], conditions );
    return;
};

