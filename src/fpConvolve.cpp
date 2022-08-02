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
#include <vector>
#include <math.h>
#include "fpConvolve.h"
#include "XrayDetector.h"
#include "XRFcontrols.h"

using namespace std;

//	Convolves calculated spectrum of any component with a Gaussian
//      Width from detector resolution at each channel energy     May 16, 2019
//  Brute force convolution, expensive but accurate, Gaussian info from fpLineSpectrum.cpp
//  Modified    Jan. 10, 2019   Check for energy <= zero and skip that point to avoid nan result
//  Modified Apr. 2, 2021       Added check for zero spectrum value in loop, to skip as many calculations as possible


void fpConvolve( const XrayDetector detector, const XrayEnergyCal cal_in, vector <float> &spectrum_out ) {
	int ns = spectrum_out.size();
	if ( ns <= 0 ) return;
	vector <float> convolve_result( ns );
    unsigned int j;
//		calculated spectrum
	for ( j=0; j<ns; j++ ) {
        convolve_result[j] = 0;
		float el = cal_in.energy( float(j) );
		if( el <= 0 ) continue;
        float fwhm_in = detector.resolution( el );
        float alpha = FWHM_SIGMA * FWHM_SIGMA / (fwhm_in*fwhm_in);	//	to get Gaussian exponent for correct fwhm 4*ln(2)
		float thresh = 0;
//		thresh = threshold_in / intensity;
//		if( threshold_in > 0.1f * max ) thresh = 0.1f;
        thresh = 1e-7f;
		if ( thresh > 1 ) continue;	//	line is below noise level
		int kMin = 0;
		int kMax = ns;
		if ( thresh > 0 ) {
			thresh = sqrt( - log ( thresh ) / alpha );
			kMin = cal_in.channel( el - thresh ) - 2 - j;
			kMax = cal_in.channel( el + thresh ) + 2 - j;
		};
		float norm = 1 / ( fwhm_in * GAUSSIAN_INTEGRAL ) //   Gaussian integral is sqrt(PI/4ln2)*fwhm
			* cal_in.energyPerChannel( j );	//	to get counts per channel
        int k;
		for ( k=kMin; k<kMax; k++ ) {
            int j_shift = j + k;
            if( j_shift < 0 || j_shift >= ns ) continue;
            if( spectrum_out[j_shift] == 0 ) continue;
			float en = cal_in.energy( j_shift );
			if( en <= 0 ) continue;
			float diff = en - el;
			convolve_result[j] += norm * spectrum_out[j_shift] * exp ( -alpha * diff*diff );
		};

	};
	for ( j=0; j<ns; j++ ) spectrum_out[j] = convolve_result[j];
    return;
};
