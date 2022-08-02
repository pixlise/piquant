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

/*
 *  upper_trim.cpp
 *  APS1IDanalysis
 *
 *  Created by W. T. Elam on 4/8/11.
 *
 */

//  Modified May 1, 2017 to ignore trailing tabs
//  Modified Jan. 3, 2017
//      Move tab, single quote, double quote, blank, comma, and underscore definitions to XRFconstants.h

#include "upper_trim.h"
#include "XRFconstants.h"

using namespace std;

string upper_trim( const string inStr ) {
	string outStr;
	string lower;
	lower = "abcdefghijklmnopqrstuvwxyz";
	string upper;
	upper = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	string tab( TAB_CHARACTER );
	int i;
//		convert to upper case and keep track of trailing blanks
	int trail = -1;
	for( i=0; i<inStr.length(); i++ ) {
		if( inStr.substr( i, 1 ) == BLANK_CHARACTER || inStr.substr( i, 1 ) == tab ) {
			if( trail < 0 ) trail = i;
		} else {
			trail = -1;
		};
		int j = lower.find( inStr.substr(i,1) );
		if( j < 0 || j >= lower.length() ) {
			outStr += inStr.substr(i,1);
		} else {
			outStr += upper.substr( j, 1 );
		};
	};
//		trim trailing blanks
	if( trail > 0 ) outStr.erase( trail, outStr.length() - trail );
	return outStr;
};
