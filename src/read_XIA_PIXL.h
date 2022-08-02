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

#ifndef read_XIA_PIXL_h
#define read_XIA_PIXL_h

#include <vector>
#include <string>
#include "XraySpectrum.h"


//	reads X-ray Fluorescence spectrum files written by ProSpect for Ketek DPP and XIA MicroDXP


int read_XIA_PIXL ( const std::string spectrumFileName, XraySpectrum &spectrum_data_out,
				std::string &spec_acq_date, std::string &spec_title,
				std::string &spec_sample, std::string &spec_unitID );

//  integer returns
//		0	file read OK
//		-1	can't open file
//		-2	invalid version (not File Version = 2 on first line of file)
//		-3	not Ketek or XIA format (missing value for MCA data keyword)
//		-4	not Ketek or XIA format (MCA data does not include "ProSpect")
//		-5	XUNITS not eV
//		-6	YUNITS not COUNTS
//		-7	not an XRF file
//		-8	NPOINTS keyword not found in file
//		-9	Not Y column format (XY not handled yet)
//		-10	Too many columns (maximum 2 for Y format) or zero columns
//      -11 Not enough values with #XPERCHAN keyword (must be one per column)
//      -12 Not enough values with #OFFSET keyword (must be one per column)
//      -13 Not enough values with #LIVETIME keyword (must be one per column)
//		-20	End of file encountered while reading spectrum data
//		-21	Bad data value encountered while reading spectrum data
//		-22	Not enough data entries on the line while reading spectrum data
//		-23	Number of channels not correct on line before spectrum starts
#endif
