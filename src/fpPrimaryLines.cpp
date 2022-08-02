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
#include "fpBeams.h"
#include "fpPrimaryLines.h"

using namespace std;

//		Perform fundamental parameters calculation of the
//		predicted measured primary spectrum from an
//		X-ray source plus anything in the primary beam

// Modified from stdCalcSpec.cpp of June 1, 2015
//	Plus code from fpExcitation.cpp, fpContScat.cpp, and fpLineScat.cpp
//  Re-written Feb. 15, 2017  from parts of primaryCalcSpec.cpp
//      But mostly using new fpRayleigh in fpMain.cpp
//      For PIQUANT Version 2, using new conditions and fpBeams functions
//      This just sets the factors in the XrayLines objects for the calculated intensity values

void fpPrimaryLines( const XRFconditions &conditions_in, std::vector <XrayLines> &sourceLines ) {

    // Calculate contribution to spectrum from characteristic source emission lines

	//cout << "Starting calculation for characteristic lines" << endl;
//		get list with intensities of tube characteristic lines
	conditions_in.source.lines( sourceLines, conditions_in.eMin );
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
			float detResp = conditions_in.detector.response( lineEn );
			lineInt *= detResp;
//              put the calculated Rayleigh intensity into the factor for this emission line
			sourceLines[edgeIndex].factor( lineIndex, lineInt );
        };
	};

	return;

};
