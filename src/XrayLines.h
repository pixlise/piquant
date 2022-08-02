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

#ifndef XrayLines_h
#define XrayLines_h

#include <vector>
#include <string>
#include "XrayEdge.h"

struct XrayLinesInfo {
    XrayEdge edgeOccupied;
    std::string lineIUPAC;
    float lineEnergy;
    float lineIntensity;
    float lineFactor;
    float matrixFactor = 1;
};

std::string XrayLinesInfo_toString(const XrayLinesInfo &);


//	Elam Ravel Sieber database class

//	Modified Sept. 20, 2010 to accomodate synchrotron sources in XRFanalysis DLL   WTE
//	Modified Aug. 20, 2012 to accomodate special fits of a subset of lines         WTE
//	Modified Feb. 20, 2017 to include common factor for all lines (for live time)
//	Modified Nov. 24, 2020 to include matrix effect factor

class XrayLines {
//		creates a list of x-ray emission lines from a vacancy
//			in the energy level associated with an absorption edge
//		keeps energies and relative intensities in memory to save search and
//			interpolation time during repetitive calculations
public:
//		default constructor just allocates variables, for vector allocation
	XrayLines ( );
//		constructor for single-entry line
	XrayLines ( const float energy );
//      modified constructor to produce subset of available lines for special fits   Aug 20, 2012
	XrayLines ( const XrayEdge& newEdge, const float separation = 0,
               const float energy_low_limit = 0, const float energy_high_limit = 1000000 );
	XrayLines ( const XrayLines& newLines, const float separation = 0 );
	int operator == ( const XrayLines& testLines ) const {
		return ( edgeVacant == testLines.edge() ); };
//			don't check index to save time (may have to add if too many confusing errors)
	const float energy ( const int index ) const { return lineList[index].lineEnergy; };
	const float energy ( const std::string& symbol ) const;
	std::string symbolIUPAC ( const int index ) const { return lineList.at(index).lineIUPAC; };
	const std::string& symbolSiegbahn ( const int index ) const;
	const int index ( const std::string& symbol ) const;
//			intensity gives relative intensity modified by user factor (emitted intensities, for example)
	const float intensity ( const int index ) const { return lineList[index].lineIntensity*lineList[index].lineFactor*commonFactor_value; };
	const float relative ( const int index ) const { return lineList[index].lineIntensity; };
//          Line width from width of energy levels
	const float width ( const int index ) const;
//			factor allows setting and retrieving user factors
	float factor ( const int index ) const { return lineList[index].lineFactor; };
	float factor ( const int index, const float newFactor ) { lineList[index].lineFactor = newFactor; return newFactor; };
	float matrix ( const int index ) const { return lineList[index].matrixFactor; };
	float matrix ( const int index, const float newFactor ) { lineList[index].matrixFactor = newFactor; return newFactor; };
//			common factor allows setting and retrieving user factor for all lines (put in for live time)
	const float commonFactor() const { return commonFactor_value; };
	void commonFactor( const float newFactor ) { commonFactor_value = newFactor; return; };
//		private data access functions
	const int numberOfLines () const { return lineCount; };
	const XrayEdge& edge() const { return edgeVacant; };
	const XrayEdge& edgeSource ( const int index ) const { return lineList.at(index).edgeOccupied; };

    std::string toString() const;

private:
	XrayEdge edgeVacant;
	int lineCount = 0;
	float commonFactor_value = 1;
    std::vector <XrayLinesInfo> lineList;
//		the following constants are used only to translate IUPAC to Siegbahn notation
	static const std::string LINE_NAMES_Siegbahn[];
	static const std::string LINE_NAMES_IUPAC[];
	float intSymbol ( std::string& symbol, const XrayEdge& upper, const XrayEdge& lower );
	float linePolyCalc ( const float z, const float array[] ) ;
	static const std::string blank;
    void merge_peaks( const float separation, std::vector <XrayLinesInfo> &lineList_in );
};

#endif
