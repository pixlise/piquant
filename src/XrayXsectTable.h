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

#ifndef XrayXsectTable_h
#define XrayXsectTable_h

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <math.h>
#include "Element.h"
#include "spline.h"

using namespace std;

//	Elam Ravel Sieber database class

class XrayXsectTable {
//	contains and evaluates a table of x-ray scattering cross sections 
//		for a particular element
//		member functions return cross section values at a single energy
//		class contains entire table for all energies
//		data read from database file (currently hard-wired in implementation file)

public:
	XrayXsectTable(const Element& el);
	XrayXsectTable();
//		equality operators
	bool operator == ( const XrayXsectTable& test ) const {
		return ( thisElement == test.element() ); };
	float coherent(const float energy) const { return cohCalc(energy); };
	float incoherent(const float energy) const { return incohCalc(energy); };
	float photo(const float energy) const { return photoCalc(energy); };
	float total(const float energy) const {
		return cohCalc(energy) + incohCalc(energy) + photoCalc(energy) ;
	};
	float inelastic(const float energy) const {
		return incohCalc(energy) + photoCalc(energy) ;
	};
	const Element& element() const { return thisElement; };
	float Rayleigh(const float energy) const { return cohCalc(energy); };
	float Compton(const float energy) const { return incohCalc(energy); };
//		access functions for private data and table values 
	static const float maxEnergy() { return 1000000.0; };
	static const float minEnergy() { return 100.0; };
	int cohCount() const { return numberCoherent; };
	int photoCount() const { return numberPhoto; };
	float coherentEnergy(const int index) const { return energiesCoh.at(index); };
	float coherentValue(const int index) const { return sigmaCoh.at(index); };
	float coherentSpline(const int index) const { return splineCoh.at(index); };
	float incoherentEnergy(const int index) const { return energiesCoh.at(index); };
	float incoherentValue(const int index) const { return sigmaIncoh.at(index); };
	float incoherentSpline(const int index) const { return splineIncoh.at(index); };
	float photoEnergy(const int index) const { return energiesPhoto.at(index); };
	float photoValue(const int index) const { return sigmaPhoto.at(index); };
	float photoSpline(const int index) const { return splinePhoto.at(index); };

private:
	float cohCalc(const float energy) const {
		return energy > 0.0 ? exp ( splint(energiesCoh, sigmaCoh, splineCoh, log(energy) ) ) :0.0;
	};
	float incohCalc(const float energy) const { 
		return energy > 0.0 ? exp ( splint(energiesCoh, sigmaIncoh, splineIncoh, log(energy) ) ) :0.0;
	};	
//			photoCalc must check for exact match to edge energy and always return xsect below edge
	float photoCalc(const float energy) const ;

//individual data
	Element thisElement;
	int numberCoherent;
	vector<float> energiesCoh;
	vector<float> sigmaCoh;
	vector<float> splineCoh;
	vector<float> sigmaIncoh;
	vector<float> splineIncoh;
	int numberPhoto;
	vector<float> energiesPhoto;
	vector<float> sigmaPhoto;
	vector<float> splinePhoto; 
};
#endif
