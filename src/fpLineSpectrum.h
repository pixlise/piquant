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

#ifndef fpLineSpectrum_h
#define fpLineSpectrum_h

#include <vector>
#include "XrayLines.h"
#include "XraySpectrum.h"
#include "XrayDetector.h"

//  Info for consolidated emission lines, used for tails, shelf, and pileup effects
struct LineGroup {
    float energy = 0;
    float intensity = 0;
    int number = 0;
    float tail_sum = 0;
    string symbol;
};


//	Generates calculated spectrum from Xray Lines object that has been loaded with intensity factor
//		Calculation is Gaussian with fwhm as input, integral matches line intensity
//		Calculated spectrum is counts in each channel
//	Added check for zero or negative energy at low channels     Dec. 12, 2011

void fpLineSpectrum( const XrayLines &lines_in, const XrayDetector detector, const float threshold_in,
				   const XrayEnergyCal &cal_in, const float eMin, std::vector <LineGroup> &pileup_list,
				   SpectrumComponent &component_out );

#endif
