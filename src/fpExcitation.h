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

#ifndef fpExcitation_h
#define fpExcitation_h

#include <vector>
#include "Element.h"
#include "XrayEdge.h"
#include "XraySource.h"

// added tube current  Nov. 30, 2011
//  Modified, Sep. 30, 2013    W. T. Elam and Nick Yang
//  For optic, include center energy +/- bandwidth in list of energies
//      so that narrow band optics will not be missed
//  opticCenter is two-element array, float input, ignored if opticCenter[0] <= 0

using namespace std;

void fpExcitation  ( const vector <XrayEdge> &sampleEdges,  const XraySource &source, const float opticCenter[], 
		const float minEnergy_in, vector <float> &excitEnergies, vector <float> &excitIntensities ) ;
// sets up vectors of fpExcitation energies and intensities 
//		for fundamental parameters x-ray fluorescence calculation
// excitEnergies and excitIntensities are loaded and returned
//     Copyright 2001  W. T. Elam

#endif
