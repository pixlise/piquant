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
#include <vector>
#include <algorithm>
#include "Element.h"
#include "XrayEdge.h"
#include "XrayLines.h"
#include "XraySource.h"
#include "nofx.h"
#include "fpExcitation.h"

//  Modified, Sep. 30, 2013    W. T. Elam and Nick Yang
//  For optic, include center energy +/- bandwidth in list of energies
//      so that narrow band optics will not be missed

using namespace std;

void fpExcitation  ( const vector <XrayEdge> &sampleEdges,  const XraySource &source, const float opticCenter[],
		const float minEnergy_in, vector <float> &excitEnergies, vector <float> &excitIntensities ) {
// sets up vectors of excitation energies and intensities
//		for fundamental parameters x-ray fluorescence calculation
//	output is in photons per second per steradian per milliAmp
//		(continuum intensities are multiplied by energy intervals for easy integration)
// excitEnergies and excitIntensities are loaded and returned
//     Copyright 2001  W. T. Elam
	int i;	//	loop counter
//		temporary vectors to hold information
//			(output vectors will be sized and loaded at end)
	vector <float> energies;
	vector <float> intensities;
//		peak inensity to use in eliminating low intensities to save calculations
	float peakIntensity = 0.0;
//		minimum energy from conditions or source, whichever is larger, and maximum energy from source
	float minEnergy = minEnergy_in;
	if( source.minEnergy() > minEnergy ) minEnergy = source.minEnergy();
// if source is an x-ray tube, start with continuum
	if ( source.continuum() ) {
//			build list of anode x-ray edges using XraySource member function
		vector <XrayEdge> anodeEdges;
		source.edges ( anodeEdges );
//			include each edge as a pair of energies, plus min and max energy
//			add energies in any order, sort later
//			anode edge pairs
		for ( i=0; i<anodeEdges.size(); i++ ) {
			if( anodeEdges[i].energy() < minEnergy || anodeEdges[i].energy() > source.voltage() ) continue;
			energies.push_back( anodeEdges[i].energy() - 0.5 );
			energies.push_back( anodeEdges[i].energy() + 0.5 );
		};
//			list of sample edges contains only those excited by this source
		for ( i=0; i<sampleEdges.size(); i++ ) {
			if( sampleEdges[i].energy() < minEnergy || sampleEdges[i].energy() > source.voltage() ) continue;
			energies.push_back(  sampleEdges[i].energy() - 0.5 );
			energies.push_back( sampleEdges[i].energy() + 0.5 );
		};
		if( opticCenter[0] > 0 ) {
			energies.push_back( opticCenter[0] );
			energies.push_back( opticCenter[1] );
		}
		energies.push_back( minEnergy );
		energies.push_back( source.voltage() );
//			sort list of energies from largest to smallest
		sort( energies.begin(), energies.end() );
		reverse( energies.begin(), energies.end() );
//			add list of corresponding intensities
		intensities.resize ( energies.size() );

//			calculate source continuum intensities for each energy
		for ( i = 0; i < energies.size(); i++ ) {
//				source intensity for continuum is photons per sec. per steradian per milliAmp per keV
			intensities[i] = source.continuum ( energies[i] );
		};
//				(note that we can't use a for loop since vector size and indexing are changing)
//				(leave the highest energy in place)
		i = 1;
		float peakIntensity = 0.0;
//			keep track of the energy just below the last energy with some intensity
		float lowestEnergy = energies[energies.size()-1];
		bool saveLowest = false;
		while ( i < energies.size() ) {
			if ( intensities[i] <= 1.0e-6*peakIntensity )	{
				if ( saveLowest) {
					saveLowest = false;
					lowestEnergy = energies[i];
				};
				energies.erase ( energies.begin() + i );
				intensities.erase ( intensities.begin() + i );
			} else {
				peakIntensity = ( intensities[i]>peakIntensity ? intensities[i] : peakIntensity );
				saveLowest = true;
				i++;
			};
		};
//			terminate the list with the energy just below the last energy with some intensity
		energies.push_back ( lowestEnergy );
		intensities.push_back ( 0.0 );

//			prepare list for integration, using trapeziodal rule

//			check each pair of points to see if the integral is significantly improved
//				by adding another energy point between the existing points
		bool listCheck = true;
		while ( listCheck ) {
//				set flag to say everything is OK,
//					change if any points are added so check will be made again
			listCheck = false;
			for ( i = 0; i < energies.size()-1; i++ ) {
//					treat energies between edge energies as separate regions, since things
//						change discontinuously at edges (either sample or anode)
				if ( fabs(energies[i] - energies[i+1]) < 1.1 ) continue;
//					calculate integral over region without any intermediate points
//						convert eV to keV
				float integral = ( 0.001 * fabs ( energies[i+1] - energies[i] ) )
					* 0.5 * ( intensities[i+1] + intensities[i] );
//					 compare this to the integral with one additional point in the middle
				float middleEnergy  = 0.5 * ( energies[i+1] + energies[i] );
				float middleIntensity = source.continuum ( middleEnergy );
				float testIntegral = ( 0.001 * fabs ( middleEnergy - energies[i] ) )
					* 0.5 * ( middleIntensity + intensities[i] )
					+ ( 0.001 * fabs ( energies[i+1] - middleEnergy ) )
					* 0.5 * ( intensities[i+1] + middleIntensity );
//						add the middle point to the list if difference in integrals is > 1%
				float norm = ( testIntegral == 0.0 ? integral : testIntegral );
				if ( norm == 0.0 ) continue;
				if ( fabs( (testIntegral-integral)/norm ) > 0.005 ) {
//						insert works at location before interator position
					energies.insert ( energies.begin()+i+1, middleEnergy);
					intensities.insert ( intensities.begin()+i+1, middleIntensity);
//						set flag to check everything again
					listCheck = true;
				};
			};
		};

//			now pre-multiply list by energy intervals and integration coefficients
//				so integral can be performed just by multiplying by list intensity values
		for ( i = 0; i < energies.size(); i++ ) {
			if ( i==0 ) {
//					start of list
				intensities[i] *= 0.5 * 0.001 * fabs ( energies[i+1] - energies[i] );
			} else {
				if ( i==energies.size()-1 ) {
//						end of list
					intensities[i] *= 0.001 * fabs ( energies[i] - energies[i-1]) * 0.5;
				} else {
//						intermediate point in region, integral has terms from both sides
//							(convert eV to keV)
					intensities[i] *= 0.001 * fabs ( energies[i+1] - energies[i-1] ) * 0.5;
				};
			};
		};

//		end of continuum
	};

//		add characteristic line intensities (no integration necessary for these,
//			since energy width is very small)
// 		get list of XraySource emission lines using XraySource member function
	vector <XrayLines> sourceLines;
	source.lines ( sourceLines, minEnergy_in );
	for ( i=0; i<sourceLines.size(); i++ ) {
		int lineIndex;
		for ( lineIndex=0; lineIndex<sourceLines[i].numberOfLines(); lineIndex++ ) {
//			insert each line in order by energy, so that excitation conditions
//				for each sample edge will be calculated correctly
			float lineEnergy = sourceLines[i].energy(lineIndex);
//				find nearest energy in list
			int location = nofx ( energies, lineEnergy );
//				be sure location points to energy just past line energy
			if ( energies.size() > 0 && energies[location] > lineEnergy ) location += 1;
//				sourceLines already has intensity information
			float inten = sourceLines[i].intensity(lineIndex);
			if ( inten <= 1.0e-6*peakIntensity ) continue;
			peakIntensity = ( inten>peakIntensity ? inten : peakIntensity );
			energies.insert ( energies.begin() + location, lineEnergy );
			intensities.insert ( intensities.begin() + location, inten );
		};
	};
//		end of characteristic lines

//		load output vectors and return
	excitEnergies.resize ( energies.size() );
	for ( i=0; i<energies.size(); i++ ) excitEnergies[i] = energies[i];
	excitIntensities.resize ( intensities.size() );
	for ( i=0; i<intensities.size(); i++ ) excitIntensities[i] = intensities[i];
//	for ( i=0; i<intensities.size(); i++ ) {
//        cout << "fpExcit  " << excitIntensities[i] << "  " << intensities[i] << endl;
//    }
	return ;
};
