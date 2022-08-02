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

#ifndef borehole_read_h
#define borehole_read_h

#include <vector>
#include <string>

//	reads asc files from Borehole XRF software    Aug. 2007   W. T. Elam  APL/UW
//  Modified Jan. 18, 2017 to use vector for conditionsArray and remove using in header
//  Modified Sept. 14, 2017 to read location information and title strings


int borehole_read ( const std::string spectrumFileName, std::vector <float> &conditionsArray,
				   std::vector <float> &spectrum, float &ev_start, float &ev_ch, float &live_time,
				   std::vector <std::string> &title_strings, float &x, float &y, float &z );

//  integer returns
//		0	file read OK
//		-1	can't open file
//		-2	invalid format
//		-3	invalid version
//		-4	XUNITS not eV
//		-5	YUNITS not COUNTS
//		-6	not an XRF file
//		-7	NPOINTS keyword not found in file
//		-8	End of File encountered while reading spectrum data

#endif
