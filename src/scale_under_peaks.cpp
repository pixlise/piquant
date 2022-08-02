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

#include <vector>
#include <math.h>
#include "scale_under_peaks.h"

using namespace std;

//  Modified May 14, 2021   Fixed bug with sigma_multiplier not being used

float scale_under_peaks( const std::vector <float> &curve,  const std::vector <float> &meas, const std::vector <float> &sigma,
            const float sigma_multiplier ) {
    //  Calculates a scale factor via least-squares match of the input curve to the measured spectrum.
    //      Then all channels more than sigma_multiplier times the statistical error (sigma) are eliminated
    //      from the fit and the process repeated.  This continues until no more channels are eliminated
    //      (usually only a few iterations, 10 maximum; S=-5 works well).
    float scale_factor = 1;
    unsigned int ns = curve.size();
    if( ns < meas.size() ) ns = meas.size();
    if( ns < sigma.size() ) ns = sigma.size();
    if( ns == 0 ) return scale_factor;
    vector <bool> included( ns, true );
    //  Use a for loop to keep compute time short and avoid any infinite loops
    int loop_count;
    for( loop_count=0; loop_count<10; loop_count++ ) {
        //  Use overall factor calculated via least squares to adjust background
        float sum_numerator = 0;
        float sum_denominator = 0;
        int count_included = 0;
        unsigned int i;
        for( i=0; i<ns; i++ ) {
            if( !included[i] ) continue;
            sum_numerator += meas[i] * scale_factor * curve[i];
            sum_denominator += scale_factor * curve[i] * scale_factor * curve[i];
            count_included++;
        }
        if( sum_denominator == 0 ) break;
        float bkg_adjustment = sum_numerator / sum_denominator;
        if( fabs( 1 - bkg_adjustment ) < 0.001f ) break;
        scale_factor *= bkg_adjustment;
        //  Eliminate everything more than 5-sigma above background so it is always under any peaks
        for( i=0; i<ns; i++ )
            if( meas[i] > scale_factor * curve[i] + sigma_multiplier * sigma[i] ) included[i] = false;
        //cout << "bkg adj loop " << loop_count << "   sum m " << sum_numerator << "   sum c " << sum_denominator << "   f " << scale_factor << "  # incl " << count_included << endl;
    }
	return scale_factor;
}

