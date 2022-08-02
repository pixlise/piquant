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

#include "fpPrimary.h"

//  Modified July 25, 2018
//      Write out some useful information if calculated intensity is zero or nan


using namespace std;

float fpPrimary  ( XrayLines &line, const vector <float> &elementAbs, const float ci,
				const vector <float> &excitEnergies, const vector <float> &excitIntensities,
			 	const float muSi, const vector <float> &sampleIncAbs, const float sinPsi1,
			 	const float sinPsi2, const float q, const float massThickness ) {
//		calculates primary fluorescence of an x-ray emission line using the
//			fundamental parameters equation
//     Copyright 2001  W. T. Elam
	float a = sinPsi1 / sinPsi2;
	float rk = line.edge().jump();
	float esubi = line.edge().yield() * ( rk - 1 ) / rk;
	float amu = a * muSi;
	float ee = line.edge().energy();
//		integrate over incident intensity
//			this assumes that incident intensities have already been multiplied by
//			the appropriate energy intervals and any integration coefficients
//			and that the energies are ordered from largest to smallest
	int i;
	float integral = 0.0;
	float cnt = 0;
	for ( i=0; i<excitEnergies.size(); i++ ) {
//			stop if incident energy is below absorption edge energy
		if ( excitEnergies[i] < ee ) break;
//			calculate absorption for this element at incident energy
		float temp = ( elementAbs[i] * excitIntensities[i] ) / ( sampleIncAbs[i] + amu );
		if( massThickness > 0 ) {
			float expArg = ( sampleIncAbs[i] + amu ) * massThickness / sinPsi1;
			if( expArg < EXP_FLOAT_TEST ) temp *= ( 1 - exp( - expArg ) );
		};
		integral += temp;
		cnt += excitIntensities[i];
	};
//		line relative intensity will be taken care of by XrayLines intensity member function
    cout.precision(4);
//    if( line.edge().element().Z() == 26 ) cout << "fpPri " << line.edge().element().symbol() << "_" << line.edge().symbol() << "  " << ee << "  " << a << "  " << esubi << "  " << amu << "  " << q << "  " << ci << "  " << integral << "  " << cnt << endl;
    if( isnan( integral ) ) cout << "fpPri error " << ee << "  " << a << "  " << esubi << "  " << a << "  " << muSi << "  " << q << "  " << ci << "  " << integral << endl;
	return q * esubi * ci * integral;
};
