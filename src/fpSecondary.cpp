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

#include "fpSecondary.h"

//  Modified May 25, 2019
//      Fix some things in this calculation, to match equations (originally implemented with mistakes)
//      Original lines are commented out and marked with today's data

using namespace std;

float fpSecondary  ( XrayLines &line, const float eiAbs, const float ci,
			 	const XrayLines &exLine, const int exLineIndex, const vector <float> &ejAbs,
			 	const float cj, const vector <float> &excitEnergies, const vector <float> &excitIntensities,
			 	const float muSi, const float muSj, const vector <float> &sampleIncAbs, const float sinPsi1,
			 	const float sinPsi2, const float q, const float massThickness ) {
//		calculates secondary fluorescence of an x-ray emission line excited by
//			an intermediate line using the fundamental parameters equation
//     Copyright 2001  W. T. Elam
	float a = sinPsi1 / sinPsi2;
	float rk = line.edge().jump();
	float esubi = line.edge().yield() * ( rk - 1 ) / rk;
	float rkj = exLine.edge().jump();
	float esubj = exLine.edge().yield() * ( rkj - 1 ) / rkj;
	esubj *= exLine.relative(exLineIndex);
	float amu = a * muSi;
	float ee = exLine.edge().energy();
//		beta term is independent of incident energy
	float beta = muSi / sinPsi2 / muSj;
//	float betaTerm = log ( 1.0 + beta ) / beta;                 //  Modified May 25, 2019
	float betaTerm = log ( 1.0 + beta ) / ( muSi / sinPsi2 );
//		integrate over incident intensity
//			this assumes that incident intensities have already been multiplied by
//			the appropriate energy intervals and any integration coefficients
//			and that the energies are ordered from largest to smallest
	int i;
	float integral = 0.0;
	for ( i=0; i<excitEnergies.size(); i++ ) {
//			stop if incident energy is below absorption edge energy
		if ( excitEnergies[i] < ee ) break;
//			calculate alpha term
		float alpha = sampleIncAbs[i] / sinPsi1 / muSj;
		float lzero = ( log ( 1.0 + alpha ) / ( sampleIncAbs[i] / sinPsi1 ) ) + betaTerm;
//		float lzero = ( log ( 1.0 + alpha ) / alpha ) + betaTerm;       //  Modified May 25, 2019
		float temp = lzero * ( ejAbs[i] * excitIntensities[i] ) / ( sampleIncAbs[i] + amu );
//			Rough approximation for very thin films
//			This is not correct, it should be calculated via the Mantler equations
		if( massThickness > 0 ) {
			float expArg = ( sampleIncAbs[i] + amu ) * massThickness / sinPsi1;
			if( expArg < THIN_SEC_FLUOR_TEST ) temp = 0;
		};
		integral += temp;
	};
//		line relative intensity will be taken care of by XrayLines intensity member function
//	cout << q << "  " << esubi << "  " << ci << "  " << esubj << "  " << cj << "  " << muij << "  " << integral << "  " << muSj << endl;
//	return 0.5 * q * esubi * ci * esubj * cj * eiAbs * integral / muSj; //  Modified May 25, 2019
	return 0.5 * q * esubi * ci * esubj * cj * eiAbs * integral;
};
