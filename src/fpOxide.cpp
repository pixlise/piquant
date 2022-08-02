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

#include "fpOxide.h"

using namespace std;

bool fpOxide  ( vector <Element> &elements, vector <float> &x, const vector <float> &oxideRatios ) {
// modifies element list and fractions to account for oxygen in oxide components of sample
//	oxide ratios are atomic ratios of oxygen for each element
//     Copyright 2001  W. T. Elam

//		see if there are any oxide components
	bool oxide = false;
	int i;
	for ( i=0; i<elements.size(); i++ ) {
		if ( oxideRatios[i] > 0.0 ) oxide = true;
	};
	if ( ! oxide ) return oxide;

//		add oxygen to list of elements 
	Element oxygen ( 8 );
	elements.push_back ( oxygen );
	x.push_back ( 0.0 );
	int oIndex = x.size()-1;
	for ( i=0; i<oIndex; i++ ) {
		if ( oxideRatios[i] > 0.0 ) {
			float oa = oxideRatios[i] * oxygen.atomicWeight() / elements[i].atomicWeight();
			x[oIndex] += x[i] * oa / ( 1.0 + oa );
			x[i] /= ( 1.0 + oa );
		};
	};
	return oxide;
};
