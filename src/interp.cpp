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

#include "interp.h"
#include <string>
#include <cmath>

//  Adapted from "Numerical Recipes in C"
//  Modified to check for n <= 0     April 17, 2009   WTE

using namespace std;

float interp(const float x, const float xa[], const float ya[], const int n)
//	linear interpolation of arrays (not vectors)
{
	int klo,khi,k;
	float h,b,a;

	if( n <= 0 ) return 0;
	klo=0;
	khi=n-1;
//		handle reverse order
	if ( xa[klo] > xa[khi] ) {
		klo=n-1;
		khi=0;
	};
	while (abs(khi-klo) > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) { throw string("interp: Bad xa input"); };
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	 return (a*ya[klo]+b*ya[khi]);
};


float interp(const float x, const vector <float> xa, const vector <float> ya)
//	linear interpolation of vectors
{
	int n = xa.size();
	if( n <= 0 ) return 0;
	if ( n != ya.size() ) { throw string("interp: xa and ya size mismatch"); };
	int klo,khi,k;
	float h,b,a;

	klo=0;
	khi=n-1;
//		handle reverse order
	if ( xa[klo] > xa[khi] ) {
		klo=n-1;
		khi=0;
	};
	while (abs(khi-klo) > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	};
	h=xa[khi]-xa[klo];
	if (h == 0.0) { throw string("interp: Bad xa input"); };
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	 return (a*ya[klo]+b*ya[khi]);
};
