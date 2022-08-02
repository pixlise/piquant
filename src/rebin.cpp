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

#include "rebin.h"
#include <math.h>
#include <cstdlib>

//  Modified May 6, 2019
//      Rebin bug from Lauren (e-mail Nov. 13, 2018 1;35PM)
//      There was a bug that only shows up when the spectrum does not actually need to be re-binned.
//      I think line 22 should be changed to:
//      	float hi_old = x_old[n_old-1] + ( x_old[n_old-1] - x_old[n_old-2] ) / 2;
//      This sets hi_old to the higher bound of the highest bin rather than the lower bound of the highest bin, otherwise it falls apart at line 59.  This seems to solve it.


using namespace std;

int rebin( const vector <float> &x_old, const vector <float> &y_old,
		const vector <float> &x_new, vector <float> &y_new ) {
	int n_old = x_old.size();
	if ( n_old < 2 ) return -1;
	if ( y_old.size() < n_old ) return -1;
	int n_new = x_new.size();
	if ( n_new < 2 ) return -2;
	if ( x_old[1] <= x_old[0] ) return -3;
	if ( x_new[1] <= x_new[0] ) return -4;
	y_new.resize( n_new, 0 );

//		bin ends are half way between successive x values
	int i;

//		find lo and hi limits of old data (ends of first and last bins)
	float lo_old = x_old[0] - ( x_old[1] - x_old[0] ) / 2;
//	float hi_old = x_old[n_old-1] + ( x_old[n_old-2] - x_old[n_old-1] ) / 2;
	float hi_old = x_old[n_old-1] + ( x_old[n_old-1] - x_old[n_old-2] ) / 2;

//		loop over all new bins
	int k;
	for ( k=0; k<n_new; k++ ) {

//			find lo and hi limits of new bin
		float lo_k;
		if ( k > 0 ) {
			lo_k = ( x_new[k-1] + x_new[k] ) / 2;
		} else {
			lo_k = x_new[k] - ( x_new[k+1] - x_new[k] ) / 2;
		};
		float hi_k;
		if ( k < n_new-1 ) {
			hi_k = ( x_new[k] + x_new[k+1] ) / 2;
		} else {
			hi_k = x_new[k] + ( x_new[k] - x_new[k-1] ) / 2;
		};

//			find index in x_old for nearest point to lo_k
		int i_lo = -1;	//	before start of old data
		if ( lo_k > hi_old ) i_lo = n_old;	//	after end of old data
		if ( lo_k >= lo_old && lo_k <= hi_old ) {
//				lo_k is within old data bounds, find it with binary search
			int lo = 0;
			int hi = n_old - 1;
			int i_search = 0;
			while ( abs(hi-lo) > 1 ) {
				i_search = ( lo + hi ) / 2;
				if ( x_old[i_search] < lo_k ) lo = i_search;
				if ( x_old[i_search] >= lo_k ) hi = i_search;
			};
			if ( fabs( x_old[lo] - lo_k ) < fabs( x_old[hi] - lo_k ) ) i_lo = lo; else i_lo = hi;
//				find limits of bin in x_old containing lo_k
			float lo_i_lo = lo_old;
			if ( i_lo > 0 ) lo_i_lo = ( x_old[i_lo-1] + x_old[i_lo] ) / 2;
			float hi_i_lo = hi_old;
			if ( i_lo < n_new-1 ) hi_i_lo = ( x_old[i_lo] + x_old[i_lo+1] ) / 2;
//				calculate overlap of new bin with old bin
			float overlap_lo = ( hi_i_lo - lo_k ) / ( hi_i_lo - lo_i_lo );
//				add counts from this old bin to new bin
			y_new[k] += y_old[i_lo] * overlap_lo;
		};


//			find index in x_old for nearest point to hi_k
		int i_hi = -1;	//	before start of old data
		if ( hi_k > hi_old ) i_hi = n_old;	//	after end of old data
		if ( hi_k > lo_old && hi_k < hi_old ) {
//				hi_k is with old data bounds, find it with binary search
			int lo = 0;
			int hi = n_old - 1;
			int i_search = 0;
			while ( abs(hi-lo) > 1 ) {
				i_search = ( lo + hi ) / 2;
				if ( x_old[i_search] < hi_k ) lo = i_search;
				if ( x_old[i_search] >= hi_k ) hi = i_search;
			};
			if ( fabs( x_old[lo] - hi_k ) < fabs( x_old[hi] - hi_k ) ) i_hi = lo; else i_hi = hi;
//				find limits of bin in x_old containing hi_k
			float lo_i_hi = lo_old;
			if ( i_hi > 0 ) lo_i_hi = ( x_old[i_hi-1] + x_old[i_hi] ) / 2;
			float hi_i_hi = hi_old;
			if ( i_hi < n_new-1 ) hi_i_hi = ( x_old[i_hi] + x_old[i_hi+1] ) / 2;
//				calculate overlap of new bin with old bin
			float overlap_hi = ( hi_k - lo_i_hi ) / ( hi_i_hi - lo_i_hi );
//				add counts from this old bin to new bin
			y_new[k] += y_old[i_hi] * overlap_hi;
//				if new bin is entirely within old bin, calculate overlap diferently
			if ( i_hi == i_lo ) y_new[k] = y_old[i_hi] * ( hi_k - lo_k ) / ( hi_i_hi - lo_i_hi ) ;
		};

//			if there are any bins between i_lo and i_hi, include them
		if ( i_hi - i_lo > 1 ) for(i=i_lo+1; i<=i_hi-1; i++ ) y_new[k] += y_old[i];

	};	//		end loop over all new bins

	return 0;

};

