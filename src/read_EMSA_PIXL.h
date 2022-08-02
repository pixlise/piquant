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

#ifndef read_EMSA_PIXL_h
#define read_EMSA_PIXL_h

#include <vector>
#include <string>
#include "XRFconditions.h"
#include "XraySpectrum.h"


//	reads X-ray Fluorescence spectrum files in EMSA/MAS format (ISO 22029 2012)
//	This version is explicitly for the Planetary Instrument for X-ray Lithochemistry (PIXL)
//	It has many user-defined keywords (following the standard) and is set up for two parallel detectors


int read_EMSA_PIXL ( const std::string spectrumFileName, 		XRFconditionsInput &conditionsStructEMSA,
		std::vector <XraySpectrum> &spectrum_vector);

//  Only writes spectrum information, not configuration information
int write_EMSA_PIXL ( const XraySpectrum &spectrum,
                const std::string spectrumFileName,
                const bool meas_only = false );

//  Utility function to allow keywords to be included with error messages
const string get_EMSA_keyword( const int index );

//		Utility function to return units for each conditions parameter read using an EMSA keyword (or user keyword)
//      If the condition is an enumerated choice, providing the value will return more useful text
const string get_EMSA_units( const int index, const int value = -999 );

//  integer returns
//		0	file read OK
//		-999999 	can't open file or file read error
//		any other value is negative of line number where error was found

#endif
