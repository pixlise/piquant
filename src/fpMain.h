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

#ifndef fpMain_h
#define fpMain_h

#include <vector>
#include "XRFconditions.h"
#include "Element.h"
#include "XrayEdge.h"
#include "XraySpectrum.h"
#include "XrayXsectTable.h"
#include "ScatterXsectTable.h"

//  Re-written Feb. 2, 2017
//      Use XrayMaterial class for specimen composition, thickness, and X-ray parameters
//      Use new conditions structure and setup for fp calculations

struct FPstorage {
	std::vector <Element> sampleElements;
	std::vector <XrayEdge> sampleEdges;
	std::vector <int> elementIndices;
	std::vector <float> excitEnergies;
	std::vector <float> excitIntensities;
	float sinExcit;
	float sinEmerg;
	float geometry;
	std::vector <XrayLines> pureLines;
};

string FPstorage_toString(const FPstorage &storage);


//		constructs a list of line energies with corresponding Elements for peak ID
//		(lines within detector FWHM are combined as weighted average)
void fpIDlist ( const float detRes, const float eMin, const float eMax,
			   std::vector <float> &energies, std::vector <Element> &elements );

//		prepare info for FP calculations of an element list and return pure element intensities
void fpPrep(FPstorage &storage, const XrayMaterial sample, const XRFconditions &conditions_in,
            std::vector <XrayLines> &pureLines );

//		perform FP calculations for a specific sample composition
void fpCalc(const FPstorage &storage, const XrayMaterial sample, const XRFconditions &conditions_in,
            std::vector <XrayLines> &sampleLines );

void fpContScat(const FPstorage &storage, const XrayEnergyCal &cal_in, const XrayMaterial &sample,
				const XRFconditions &conditions_in, std::vector <float> &continuumSpec );

void fpRayleigh(const FPstorage &storage, const XrayMaterial &sample, const XRFconditions &conditions_in,
            std::vector <XrayLines> &scatterLines );

void fpCompton(const FPstorage &storage, const XrayEnergyCal &cal_in, const XrayMaterial &sample,
				const XRFconditions &conditions_in, SpectrumComponent &component_out );

#endif
