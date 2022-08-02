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
#include <ios>
#include <iosfwd>
#include "XrayXsectTable.h"
//#include "xrayxsct_data.h"
#include "xrayxsct_data_cmh.h"
#include "nofx.h"

using namespace std;

//	Elam Ravel Sieber database class

//	Modified Oct. 1, 2015 to get data from compiled array, not binary file
//  Modified Nov. 9, 2016 to use attenuation cross-sections as modified by Chris Heirwegh
//      (These replace the XCOM cross-sections with the Chantler values for selected elements)
//      (only change here is name of include file with database data)
//      Na 0.1 keV to Na K binding edge
//      Mg 0.1 keV to Mg K binding edge
//      Al 0.1 keV to Al K binding edge
//      Si 0.1 keV to Si K binding edge
//      O at E points 1.0, 1.0, 1.5 and 2.0.
//	Chantler, C.T., Olsen, K., Dragoset, R.A., Chang, J., Kishore, A.R., Kotochigova, S.A., and Zucker, D.S. (2005), X-Ray Form Factor, Attenuation and Scattering Tables (version 2.1). [Online] Available: http://physics.nist.gov/ffast [year, month day]. National Institute of Standards and Technology, Gaithersburg, MD.
//	Originally published as Chantler, C.T., J. Phys. Chem. Ref. Data 29(4), 597-1048 (2000); and Chantler, C.T., J. Phys. Chem. Ref. Data 24, 71-643 (1995).
//  Modified July 25, 2019 to prevent na returns when energy is too low for tables

XrayXsectTable::XrayXsectTable(const Element& el) {
	int i;

//		check to see if atomic number is within range
	if ( el.Z()<1 || el.Z()>maxZ ) {
		throw string("XrayXsectTable: element not in database" );
	} else {
		thisElement = el;
		int thisZ = el.Z();
		int db_pointer = db_index[thisZ];
//			get number of entries for coherent and incoherent tables and resize vectors
		numberCoherent = db_numberCoherent[thisZ];
		energiesCoh.resize(numberCoherent);
		sigmaCoh.resize(numberCoherent);
		splineCoh.resize(numberCoherent);
		sigmaIncoh.resize(numberCoherent);
		splineIncoh.resize(numberCoherent);
//			read data into vectors
		float f;
		for (i=0; i<energiesCoh.size(); i++ ) {
			int thisPointer = db_pointer + i * 5;
			f = db_data[thisPointer+0];
			energiesCoh[i] = f;
			f = db_data[thisPointer+1];
			sigmaCoh[i] = f;
			f = db_data[thisPointer+2];
			splineCoh[i] = f;
			f = db_data[thisPointer+3];
			sigmaIncoh[i] = f;
			f = db_data[thisPointer+4];
			splineIncoh[i] = f;
		};
//			get number of entries for photoabsorption table and resize vectors
		numberPhoto = db_numberPhoto[thisZ];
		energiesPhoto.resize(numberPhoto);
		sigmaPhoto.resize(numberPhoto);
		splinePhoto.resize(numberPhoto);
		db_pointer = db_index[thisZ] + 5 * numberCoherent;
//			read data into vectors
		for (i=0; i<energiesPhoto.size(); i++ ) {
			int thisPointer = db_pointer + i * 3;
			f = db_data[thisPointer+0];
			energiesPhoto[i] = f;
			f = db_data[thisPointer+1];
			sigmaPhoto[i] = f;
			f = db_data[thisPointer+2];
			splinePhoto[i] = f;
		};

	};
};

//	default constructor for initializing vectors
XrayXsectTable::XrayXsectTable() {
	numberCoherent = 0;
	numberPhoto = 0;
	energiesCoh.resize(0);
	sigmaCoh.resize(0);
	splineCoh.resize(0);
	sigmaIncoh.resize(0);
	splineIncoh.resize(0);
	energiesPhoto.resize(0);
	sigmaPhoto.resize(0);
	splinePhoto.resize(0);
};


float XrayXsectTable::photoCalc(const float energy) const {
	if ( energy <= 0.0 ) return 0.0;
	int i = nofx ( energiesPhoto, log(energy) );
//		if the requested energy exactly matches an edge discontintuity pair,
//			 return the cross section value below the edge
	if ( log(energy) != energiesPhoto[i] ) {
		return exp ( splint(energiesPhoto, sigmaPhoto, splinePhoto, log(energy) ) ) ;
	} else {
		bool pair = false;
		if ( i>0 && energiesPhoto[i] == energiesPhoto[i-1] ) pair = true;
		if ( i<energiesPhoto.size()-1 && energiesPhoto[i] == energiesPhoto[i+1]) pair = true;
		float value = 0;
		if ( pair ) {
			value = exp ( splint(energiesPhoto, sigmaPhoto, splinePhoto, log(energy)-1.0e-6 ) ) ;
		} else {
			value = exp ( splint(energiesPhoto, sigmaPhoto, splinePhoto, log(energy) ) ) ;
		};
		if( isnan(value) ) return 0;
		return value;
	};
};
