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

#include <math.h>
#include <iostream>
#include "Sewell_tube_calc.h"
#include "XRFconstants.h"

//  X-ray tube spectrum calculation based on electron probe model of Sewell
//      Reference D. A. Sewell, G. Love, and V. D. Scott, J. Phys. D: Appl. Phys. 20(1987) 1567-1573.
//  Written June 23, 2020

using namespace std;

float Sewell_j( float z ) {
	return 13.5f * z;
};

float Sewell_eta( float z, float energy ) {
//		backscattered electron fraction
//		Ref Myklebust and Newbury, Electron Probe Quantitation,
//			Heinrich and Newbury, eds. (Plenum, New York) 1991, p 177.
	float z2 = z * z;
	float z3 = z2 * z;
	float hz = ( -1112.8f + 30.29f * z - 0.15498f * z2 ) * 1.0e-4f;
	float eta = -52.3791f;
	eta += 150.48371f * z;
	eta += -1.67373f * z2;
	eta += 0.00716f * z3;
	eta *= 1.0e-4f;
	return eta * ( 1 + hz * log ( energy / 20000 ) );
};

float Sewell_r( float u0, float eta, float tilt ) {
//		backscatter correction factor
	float one_u0 = 1 / u0;
	float eta_tilted = eta;
	if ( tilt != 0 ) eta_tilted = 0.891f * pow ( eta / 0.891f, cos ( tilt ) );
	float exponent = 1.5f - 1.5f * pow ( u0, -0.25f );
	float i = 0.3f * ( -one_u0 + exp ( exponent ) );
	exponent = 1 - 2.3f * pow ( u0, -4 );
	float g = ( 0.368f - 0.075f * log ( u0 ) ) * exp ( exponent );
	float rhs = i + eta_tilted * g;
	return 1 - eta_tilted * pow ( rhs, 1.667f );
};

float Sewell_S_lines( float u0, float j, float ee, float z_a ) {
//		stopping power used for tube lines
	if ( u0 <= 1 || ee <= 0 ) return 0;
	float s = 1;
	float uRatio = ( sqrt ( u0 ) - 1 ) / ( u0 - 1 );
	s += 16.05f * sqrt ( j / ee ) * pow ( uRatio, 1.07f );
	s /= z_a;
	if ( s <= 0 ) s = MINIMUM;
//	s *= 1.47e-6f  * ( u0 * log ( u0 ) + ( 1 - u0 ) ) / ( ee * ee * 1e-6f );
	s *= 1.47e-6f  * ( u0 * log ( u0 ) + ( 1 - u0 ) );
	s *= 1.85e17f;
	return s;
};

float Sewell_S_continuum( float u0, float j, float energy, float z_a ) {
//		stopping power used for continuum
	if ( u0 < 1 || energy <= 0 ) return 0;
	float s = 1;
	float uRatio = ( sqrt ( u0 ) - 1 ) / ( u0 - 1 );
	s += 16.05f * sqrt ( j / energy ) * pow ( uRatio, 1.07 );
	s /= z_a;
	if ( s <= 0 ) s = MINIMUM;
//	s *= ( u0 * log ( u0 ) + ( 1 - u0 ) ) / ( energy * energy );
//	return 6.37e17 * s;
	s *= ( u0 * log ( u0 ) + ( 1 - u0 ) );
	return 1.96e10 * s;
};

float Sewell_stopping( float energy, float j, float z_a ) {
	float v = energy / j;
	float j1000 = j / 1000;
	float den = 1.18e-5f * sqrt ( v ) + 1.47e-6f * v;
	float r = - ( z_a / j1000 ) / den;
	return r;
};

float Sewell_Q( float energy, float ec ) {
//		electron ionization (Bethe form)
	float u = energy / ec;
	float logU = log ( u );
	float ec_kV = ec / 1000;
	float q = logU / ( u * ec_kV * ec_kV );
	return q;
};

float Sewell_pz( float j, float energy, float eta, float u0,
			   float z_a, float tilt ) {
	float num1 = 1.1e-5f * pow ( j, 1.1f ) * pow ( energy, 1.2f );
	float num2 = 3e-6f * pow ( j, 0.13f ) * pow ( energy, 1.75f )
		* pow ( energy, -0.0008f * eta * u0 );
	float den1 = 1.1f + 6.5f * eta + 3.5f * j - 3 * pow ( eta, 1.5f );
	float logu0 = log ( u0 );;
	float num = ( num1 + num2 ) * logu0;
	float den = ( den1 * logu0 + 1 + 0.08f / eta ) * z_a;
	float pz = num / den;
	if ( tilt != 0 ) {
		float alpha = 0.708f;
		float fac = 1 - alpha + alpha * cos ( tilt );
		pz *= fac;
	};
	return pz;
};

float Sewell_h( float u0, float z, float eta, float tilt ) {
	float x = 1.225f - 1.25f * eta;
	float a1 = 2.2f + 1.88e-3f * z;
	float a3 = 0.01f + 7.19e-3f * z;
	float a2 = ( a1 - 1 ) * exp ( a3 );
	float h90 = a1 - a2 * exp ( -a3 * pow ( u0, x ) );
	float h_tilt = h90;
	if ( tilt != 0 ) {
		float u023 = pow ( u0, -2.0f/3.0f );
		float c_tilt = cos ( tilt );
		float fac = 0.44f + 0.56f * u023 + 0.56f * ( 1 - u023 ) * c_tilt * c_tilt;
		h_tilt = h90 * fac;
	};
	return h_tilt;
};

float Sewell_pz_m( float pz, float u0, float z, float tilt ) {
	float term2 = ( 0.662f + 0.443f * pow ( u0, 0.2f ) ) / sqrt ( z ) ;
	float f =  0.29 + term2;
//		tilt correction
	if ( tilt != 0 ) {
		float c_tilt = cos ( tilt );
		float f_tilt = 1.18f * c_tilt * c_tilt - 0.18f;
		if ( f_tilt < 0.01f ) f_tilt = 0.01f;
		f *= f_tilt;
	};
	if ( f <= 0 ) f = MINIMUM;
	return pz * f;
};

float Sewell_pz_r( float pz, float pz_m, float h ) {
//		X-ray range in Scott / Love quadrilateral Phi_rho_Z absorption correction
//		solve quadratic equation pz = ( pz_m**2 + h*pz_r**2 + h*pz_m*pz_r ) / 3*( pz_m + h*pz_r )
	float m_3b = pz_m - 3 * pz;
	float b = h * m_3b;
	float t_2 = b*b - 4 * h * ( pz_m * m_3b );
	if ( t_2 < 0 ) return pz_m;
	t_2 = sqrt ( t_2 );
	float r = ( - b - t_2 ) / ( 2 * h );
	if ( r <= pz_m ) r = ( - b + t_2 ) / ( 2 * h );
	return r;
};

float Sewell_phi_pz( float pz, float pz_m, float pz_r, float h ) {
//	quadrilateral model
	float phiZero = 1;
	if ( pz <= 0 ) return phiZero;
	if ( pz > pz_r ) return 0;
	float phi;
	if ( pz < pz_m ) {
		phi = phiZero * ( 1 + ( h - 1 ) * pz / pz_m );
	} else {
		phi = phiZero * h * ( pz_r - pz ) / ( pz_r - pz_m );
	};
	return phi;
};

float Sewell_f( float chi, float pz_m, float pz_r, float pz, float h ) {
	float n11 = - exp ( - chi * pz_m );
	float n21 = - n11;
	float f = n21;
	float n12 = h * exp ( - chi * pz_r );
	float n13 = chi * ( pz_r - pz_m );
	float n22 = pz_r - h * pz_r;
	float n23 = h * pz_r - pz_r;
	float t_1 = n11 + n12 + n13 - h + 1;
	float d1 = pz_r - pz_m;
	float d2 = ( pz_m + h * pz_r ) * chi*chi;
	f = 2 * ( t_1 + ( n21 * n22 + n23 ) / pz_m ) / ( d1 * d2 );
	return f;
};

float Reed_f( float cj, float tau_ij, float sigma_j, float ri, float omega_j, float ai,
				float aj, float ui, float uj, float chi_i, float sigmaLenard ) {
//		characteristic fluorescence correction
//			Ref S. J. B. Reed, "Electron Microprobe Analysis", 2nd Ed.
//			(Cambridge Univ. Press, Cambridge) 1997. ISBN 0-521-59944-X
	float jump_i = ( ri - 1 ) / ri;
	float bigUfactor_j = uj * log ( uj ) - uj + 1;
	float bigUfactor_i = ui * log ( ui ) - ui + 1;
	float f = cj / 2;
	f *= ( tau_ij / sigma_j ) * jump_i * omega_j * ( ai / aj );
	f *= bigUfactor_j / bigUfactor_i;
	float u = chi_i / sigma_j;
	float corr = log ( 1 + u ) / u;
	float v = sigmaLenard / sigma_j;
	corr += log ( 1 + v ) / v;
	f *= corr;
	return f;
};
