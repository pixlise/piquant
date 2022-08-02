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

#ifndef XraySource_h
#define XraySource_h
#include "XRFconstants.h"
#include "Element.h"
#include "XrayEdge.h"
#include "XrayLines.h"
#include "XrayXsectTable.h"
#include "XrayMaterial.h"

using namespace std;

enum XrayIsotope {Cd109=1};

class XraySource {
//     Copyright 2001  W. T. Elam
//		Modified March 27, 2007 for end-window X-ray tubes
//		Modified June 18, 2009 to add Cd-109 radioisotope source
//		Modified Sept. 20, 2010 for synchrotron sources (requires modified XrayLines class)
//		Modified Dec. 1, 2011 tube current added (milliAmps)
//      Modified May 10, 2015 to return tube current
//      Modified May 21, 2020 to add constructor from XrayLines objects and continuum arrays
//                                  so tube spectrum can be used from a different calculation
//                            move Be window absorption out of header file for Ebel continuum

public:
//		must have default constructor to declare arrays
//		Default constructor
	XraySource ( );
//		X-ray tube
	XraySource( Element anodeIn, float kVIn, float incAngleIn,
			float takeoffAngleIn, float mmBeIn, float anodeThicknessIn = 0, float tubeCurrent_in = 1 );
	XraySource( XrayMaterial anodeIn, float kVIn, float incAngleIn,
			float takeoffAngleIn, float mmBeIn, float anodeThicknessIn = 0, float tubeCurrent_in = 1 );
//		X-ray source from some other calculation
//          assumes line intensites are ph/sec/sr/mA and continuum intensities are ph/sec/keV/sr/mA
	XraySource( std::vector <XrayLines> lines_in, std::vector <float> continuum_energies_in,
			std::vector <float> continuum_intensities_in, float mmBeIn, float kVIn = 0, float tubeCurrent_in = 1 );
//		radioisotope source with activity in Becquerels (Bq, or disentegrations per second, 3.7e10 Bq = 1 Curie)
	XraySource( XrayIsotope isotopeIn, float activityIn, float mmBeIn );
//		synchrotron source: single energy with intensity input in photons per second
	XraySource( const float energyIn, const float activityIn, float mmBeIn );
//		Member functions
	float continuum ( const float energy ) const {
        if( energy <= 0 ) return 0;
        else if( tube ) return continuum_Ebel( energy );
//        else if( tube ) return continuum_Sewell( energy );
        else if( external ) return continuum_ext( energy );
        else return 0;	};
    int lines ( vector <XrayLines> &linesOut, const float eMin ) const {
        if( tube && endWindow ) return tubeLinesEbel ( linesOut, eMin );
        else if( tube && !endWindow ) return tubeLinesSewell ( linesOut, eMin );
        else if( sr ) return srLines( linesOut, eMin );
        else if( external ) return extLines( linesOut, eMin );
        else if( activity > 0 ) return isotopeLines ( linesOut, eMin );
        else { linesOut.clear(); return 0; }  };
	int edges ( vector <XrayEdge> &edgesOut ) const { return tubeEdges ( edgesOut ); };
//		set and/or change parameters
	void setVoltage ( const float voltageIn ) { tubeVoltage = voltageIn; };
	void setkV ( const float kVIn ) { tubeVoltage = kVIn * 1000.0; };  //  Convert eV to keV
//			convert degrees to radians and store sine of electron incidence angle
	void setIncAngle ( const float incAngleIn ) { incSin = sin(incAngleIn*RADDEG); inc_angle = incAngleIn*RADDEG; };
//			convert degrees to radians and store sine of takeoff angle
	void setTakeoffAngle ( const float takeoffAngleIn ) { takeoffSin = sin(takeoffAngleIn*RADDEG); takeoff_angle = takeoffAngleIn*RADDEG; };
//			store Be window thickness in cm2/gm = density [gm/cm3] * thickness [mm] / 10.0 [cm/mm]
	void setmmBe ( const float mmBeIn ) { window.thickness( mmBeIn/10.0 );  return; };  //  convert mm to cm
//		private data access functions
	const bool continuum ( ) const { return ( tube  || external ); };
	const XrayMaterial &anode ( ) const { return target; };
	const float minEnergy ( ) const { return XrayXsectTable::minEnergy(); };
//		this is maximum energy available from sources other than x-ray tubes
	const float voltage ( ) const { return tubeVoltage; };
	const float kV ( ) const { return tubeVoltage / 1000.0; };  //  Convert eV to keV
	const float current ( ) const { return tubeCurrent; };
	const float incAngle ( ) const { return DEGRAD*asin(incSin); };
	const float takeoffAngle ( ) const { return DEGRAD*asin(takeoffSin); };
	const float mmBe ( ) const { return window.thickness() * 10.0; };   //  Convert cm to mm

	string toString() const;

private:
	bool tube;
	bool sr;
	bool external;
	XrayMaterial target;
	float tubeVoltage;
	float incSin;
	float takeoffSin;
	float inc_angle;
	float takeoff_angle;
	float tubeCurrent;		//	added Dec. 1, 2011 (milliAmps)
	float activity;	//	added June 18, 2009 for radioisotope source
	bool endWindow;
	XrayMaterial window;
	//  These 3 added May 21, 2020, for use with constructor that supplies these from an external source
	vector <XrayLines> tube_lines_ext;   //  line intensites are ph/sec/sr/mA
	vector <float> continuum_energies;
	vector <float> continuum_intensities;   //  continuum intensities are ph/sec/keV/sr/mA


//		functions for radioisotope source calculations
	int isotopeLines ( vector <XrayLines> &lines, const float eMin ) const ;
//		functions for synchrotron source calculations
	int srLines( vector <XrayLines> &lines, const float eMin ) const ;
//		functions used to calculate x-ray tube output
	float continuum_Ebel ( const float energy ) const;
	float absCorr_Ebel ( const float energy ) const ;
	float rhozbar_Ebel ( const float energy ) const ;
//		functions to calculate characteristic lines and intensities
	int tubeLinesEbel ( vector <XrayLines> &lines, const float eMin ) const;
	int tubeEdges ( vector <XrayEdge> &edges ) const;
	int extLines( vector <XrayLines> &lines, const float eMin ) const;
	float continuum_ext ( const float energy ) const;
	int tubeLinesSewell ( vector <XrayLines> &lines, const float eMin ) const;
//	float continuum_Sewell ( const float energy ) const;

};

#endif
