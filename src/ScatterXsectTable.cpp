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
#include "interp.h"
#include "XrayEdge.h"
#include "XRFconstants.h"
#include "ScatterXsectTable.h"
#include "scatxsct_data.h"

using namespace std;

//	Elam Ravel Sieber database class

//  Modified Oct. 1, 2015 to use data from include file (compiled in)
//      instead of binary file read at execution

ScatterXsectTable::ScatterXsectTable(const Element& el) {
	int i;

//		check to see if atomic number is within range
	if ( el.Z()<1 || el.Z()>maxZ ) {
		throw string("ScatterXsectTable: element not in database" );
	} else {
		thisElement = el;
        int thisZ = thisElement.Z();
		int db_pointer = db_index[thisZ];
//			get number of entries for coherent table and resize vectors
		numberIncoherent = db_numberIncoherent[thisZ];
		xIncoh.resize(numberIncoherent);
		sofx.resize(numberIncoherent);
//			move data into vectors
		float f;
		for (i=0; i<xIncoh.size(); i++ ) {
			int thisPointer = db_pointer + i * 2;
			f = db_data[thisPointer+0];
			xIncoh[i] = f;
			f = db_data[thisPointer+1];
			sofx[i] = f;
		};
//			get number of entries for incoherent table and resize vectors
		numberCoherent = db_numberCoherent[thisZ];
		xCoh.resize(numberCoherent);
		fofx.resize(numberCoherent);
		db_pointer = db_index[thisZ] + 2 * numberIncoherent;
//			move data into vectors
		for (i=0; i<xCoh.size(); i++ ) {
			int thisPointer = db_pointer + i * 2;
			f = db_data[thisPointer+0];
			xCoh[i] = f;
			f = db_data[thisPointer+1];
			fofx[i] = f;
		};
	};
};

//	default constructor for initializing vectors
ScatterXsectTable::ScatterXsectTable() {
	numberCoherent = 0;
	numberIncoherent = 0;
};

float ScatterXsectTable::eC ( const float energy_in, const float angle_in ) {
//	calculates Compton energy
	float cosTheta = cos ( angle_in );
	float oneMinusCosTheta = 1.0f - cosTheta;
	float alpha = energy_in / ME;
	return energy_in / ( 1 + alpha * oneMinusCosTheta );
};

float ScatterXsectTable::eComptonUp ( const float energy_in, const float angle_in ){
//	calculates incidnet energy for given Compton-scattered energy
	float cosTheta = cos ( angle_in );
	float oneMinusCosTheta = 1.0f - cosTheta;
	float alpha = energy_in / ME;
	return energy_in / ( 1 - alpha * oneMinusCosTheta );
};


float ScatterXsectTable::cohCalc(const float energy_in, const float angle_in) const {
	if ( energy_in <= 0.0f ) return 0.0;
	float cosTheta = cos ( angle_in );
	float cos2Theta = cosTheta * cosTheta;
	float x = sin ( angle_in / 2.0 ) * energy_in / (HC*1000.0f);		//	sin(theta/2) / lambda, momentum transfer variable
	float onePlusCos2Theta = 1.0 + cos2Theta;
	float f = interp ( x, xCoh, fofx );
	float sigmaRayleigh = 0.5 * RE2 * f * f * onePlusCos2Theta;
	sigmaRayleigh *= AVOGADRO / thisElement.atomicWeight();				//	differential cross section at theta
	return sigmaRayleigh;
};

//			single-diferential Compton cross sections (vs incident energy and angle)
//			calculated using the incoherent scatter function
float ScatterXsectTable::incohCalc(const float energy_in, const float angle_in) const {
	if ( energy_in <= 0.0f ) return 0.0;
	float cosTheta = cos ( angle_in );
	float cos2Theta = cosTheta * cosTheta;
	float oneMinusCosTheta = 1.0f - cosTheta;
	float alpha = energy_in / ME;
	float den = 1.0f + alpha * oneMinusCosTheta;
	float h = ( 1.0f + cos2Theta + alpha * alpha * oneMinusCosTheta * oneMinusCosTheta
		/ ( 1.0f + alpha * oneMinusCosTheta )     )     / ( den * den );
//			sin(theta/2) / lambda, momentum transfer variable
	float x = sin ( angle_in / 2.0 ) * energy_in / (HC*1000.0f);
//		RE2 is classical electron radius (squared) in cm2 (defined in file FPconstants.h)
	float sigmaCompton = 0.5f * RE2 * h * interp ( x, xIncoh, sofx );
//		convert from cm2 per atom to cm2/gm (AVOGADRO defined in file FPconstants.h)
	sigmaCompton *= AVOGADRO / thisElement.atomicWeight();
	return sigmaCompton;
};

//	doubly-differential Compton cross sections (vs incident energy, angle, and scattered energy)
//	calculated using the relativistic impulse approximation and analytic Compton profiles
//		(See references in header.)
float ScatterXsectTable::incohCalc(const float energy_in, const float angle_in, const float ePrime_in) const {

	float cosTheta = cos ( angle_in );
	float ec = eC ( energy_in, angle_in );
//		get all of the occupied electron orbitals (identified by corresponding x-ray absorption edges)
	vector <EdgeIndex> edgeIndices;
	vector <XrayEdge> edges;
	int ns = XrayEdge::numberOccupied( edgeIndices, thisElement );
	int jLoop;
	for ( jLoop=0; jLoop<ns; jLoop++ ) {
		XrayEdge edge ( thisElement, edgeIndices[jLoop] );
		edges.push_back ( edge );
	};
//cout << "edges " << thisElement.symbol() << " " << edges.size() << endl;
	float q = sqrt ( energy_in*energy_in + ePrime_in*ePrime_in - 2.0 * energy_in * ePrime_in * cosTheta );
	float pz = energy_in * ( ePrime_in - ec ) / ( ec * q );
//		ME is electron rest energy in eV (defined in file FPconstants.h)
	float bigR = (energy_in/ME) * (  sqrt( 1.0 + pz*pz ) + ( energy_in - ePrime_in * cosTheta ) * pz / q );
	float bigRinv = 1.0/bigR;
	float bigRprime = bigR - (ePrime_in/ME) * ( energy_in/ec - 1 );
	float bigRprimeInv = 1.0/bigRprime;
	float invFac = bigRinv - bigRprimeInv;
	float bigX = bigR*bigRprimeInv + bigRprime*bigRinv + 2 * invFac + invFac*invFac;
	float jAll = 0.0;
	for ( jLoop=0; jLoop<edges.size(); jLoop++ ) {
//			account for kinematic limit
		if ( edges[jLoop].energy() > energy_in-ePrime_in ) continue;
		float jZero = edges[jLoop].jzero();
//			convert from hydrogenic units (m e**2 / hBar, as in Biggs et. al.) to units of mc (as in Brusa et. al.)
//				using the fine structure constant, alpha
		jZero *= ALPHA_INV;
		float jZeroArg = 1.0 + 2.0 * jZero * fabs(pz);
		jAll += edges[jLoop].occupancy() * jZero * jZeroArg * exp ( 0.5 - 0.5 * jZeroArg * jZeroArg );
	};
//		RE2 is classical electron radius (squared) in cm2 (defined in file FPconstants.h)
	float ddcs = 0.5 * RE2 * (ePrime_in/energy_in)*(1.0/q) * sqrt( 1.0 + pz*pz ) * bigX * jAll;
//		convert from cm2 per atom to cm2/gm (AVOGADRO defined in file FPconstants.h)
	return ddcs * ( AVOGADRO / thisElement.atomicWeight() );

};
