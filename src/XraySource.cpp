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
#include <math.h>
#include "XRFconstants.h"
#include "XraySource.h"
#include "Element.h"
#include "XrayEdge.h"
#include "XrayLines.h"
#include "XrayXsectTable.h"
#include "interp.h"
#include "Sewell_tube_calc.h"
#include <sstream>

#define Rh_adjust_min_Z 43
#define Rh_adjust_max_Z 46
#define Rh_Lb6_factor 25.0f          //  2.922 keV
#define Rh_Lb2_15_factor 0.7f       //  3.002 keV
#define Rh_Ll_factor 1.0f           //  2.375 keV
#define Rh_Ln_factor 1.0f           //  2.517 keV
#define Rh_Lg1_factor 1.0f          //  3.144 keV
#define Rh_Lb34_factor 1.0f         //  2.900 keV
#define Rh_Lg23_factor 1.0f         //  3.362 keV
//  Factors prior to April 28 adjustment to match BHVO spectrum from Elemental Calibration
//#define Rh_Lb6_factor 50.0f          //  2.922 keV
//#define Rh_Lb2_15_factor 0.7f       //  3.002 keV
//#define Rh_Ll_factor 0.0f           //  2.375 keV
//#define Rh_Ln_factor 1.0f           //  2.517 keV
//#define Rh_Lg1_factor 0.25f          //  3.144 keV
//#define Rh_Lb34_factor 0.0f         //  2.900 keV
//#define Rh_Lg23_factor 0.3f         //  3.362 keV


using namespace std;

//	this class currently represents an x-ray tube, but may include radioactive isotopes in the future
//	calculations of the x-ray tube intensity are from:
//			Horst Ebel, "X-ray Tube Spectra", X-RAY SPECTROMETRY 28, 255-266 (1999).
//		Modified March 27, 2007 for end-window X-ray tubes    WTE
//		Modified June 18, 2009 to add Cd-109 radioisotope source   WTE

//			initialize static constants in private data
//	Element XraySource::window("Be");
//	XrayXsectTable XraySource::beAbs(window);
//     Copyright 2001  W. T. Elam

//		tube current added Dec. 1, 2011
//  Modified Feb. 11, 2017 to exclude emission lines less than the minimum energy for spectrum analysis (eMin)
//  Modified May 13, 2019
//      Crude correction to Rh L lines from X-ray tube to better fit Rh L scatter peaks in spectrum (PIXL flight tube)
//      Didn't match measured spectrum of SN44 flight tube, so removed it (must be something in optic)
//  Modified May 21, 2020 to add constructor from XrayLines objects and continuum arrays
//                                  so tube spectrum can be used from a different calculation
//  Modified June 22, 2020  tube spectrum calculation using prescription of Sewell [Love-Scott with 3-point phi(rho*z) and tilt]
//           June 23, 2020  Best results is using Ebel 1999 for continuum and Sewell for lines with Coster_Kronig (reduction only)
//                          Constants:  Ebel continuum 1.35e9    Sewell S 1.85e17   match Rh K, L, and Continuum from PIXL FM tube SN044 28kV
//  Modified Aug. 14, 2020  Major revisions to use XrayMaterial class for multi-element anode and window
//  Modified Apr. 28, 2021  Modify adjustment factors for Rh L lines to better match BHVO spectrum from PIXL FM Elemental Calibration


//		default constructor
XraySource::XraySource ( ) {
	tube = false;
	tubeVoltage = 0.0;
	setIncAngle ( 90.0 );
	setTakeoffAngle ( 90.0 );
	setmmBe ( 0.0 );
	tubeCurrent = 0;
	endWindow = false;
	activity = 0;
	sr = false;
	external = false;
};

//		X-ray tube source

XraySource::XraySource( Element anodeIn, float kVIn, float incAngleIn,
			float takeoffAngleIn, float mmBeIn, float anodeThicknessIn, float tubeCurrent_in ) {
//		standard x-ray tube constructor
	activity = 0;
	sr = false;
	external = false;
	tube = true;
	XrayMaterial temp_anode( anodeIn );
	target = temp_anode;
	Element Be( 4 );
	XrayMaterial temp_window( Be );
	window = temp_window;
	tubeVoltage = kVIn * 1000.0;
	setIncAngle ( incAngleIn );
	if( takeoffAngleIn < 0 ) {
		setTakeoffAngle ( -takeoffAngleIn );
		endWindow = true;
	} else {
		setTakeoffAngle ( takeoffAngleIn );
		endWindow = false;
	};
	setmmBe ( mmBeIn );
	tubeCurrent = tubeCurrent_in;
//		convert microns to cm
	target.thickness( CM_MICRON * anodeThicknessIn );
};

XraySource::XraySource( XrayMaterial anodeIn, float kVIn, float incAngleIn,
			float takeoffAngleIn, float mmBeIn, float anodeThicknessIn, float tubeCurrent_in ) {
//		standard x-ray tube constructor
	activity = 0;
	sr = false;
	external = false;
	tube = true;
	target = anodeIn;
	Element Be( 4 );
	XrayMaterial temp_window( Be );
	window = temp_window;
	tubeVoltage = kVIn * 1000.0;
	setIncAngle ( incAngleIn );
	if( takeoffAngleIn < 0 ) {
		setTakeoffAngle ( -takeoffAngleIn );
		endWindow = true;
	} else {
		setTakeoffAngle ( takeoffAngleIn );
		endWindow = false;
	};
	setmmBe ( mmBeIn );
	tubeCurrent = tubeCurrent_in;
//		convert microns to cm
	target.thickness( CM_MICRON * anodeThicknessIn );
};

XraySource::XraySource( std::vector <XrayLines> lines_in, std::vector <float> continuum_energies_in,
			std::vector <float> continuum_intensities_in, float mmBeIn, float kVIn, float tubeCurrent_in ) {
//		use X-ray line intensities and continuum intensities from an external source
	activity = 0;
	sr = false;
	tube = false;
	external = true;
	Element Be( 4 );
	XrayMaterial temp_window( Be );
	window = temp_window;
	//  Copy X-ray emission lines into object storage
	tube_lines_ext.resize( lines_in.size() );
	//  Also find maximum energy and anode element
	float max_energy = 0;
	unsigned int k;
	for( k=0; k<lines_in.size(); k++ ) {
        if( max_energy == 0 ) {
            XrayMaterial temp_anode( lines_in[k].edge().element() );
            target = temp_anode;
        }
        int l;
        for( l=0; l<lines_in[k].numberOfLines(); l++ ) if( lines_in[k].energy(l) > max_energy ) max_energy = lines_in[k].energy(l);
        tube_lines_ext[k] = lines_in[k];
    }
	//  Copy X-ray emission lines and continuum arrays into object storage
	continuum_energies.resize( continuum_energies_in.size() );
	continuum_intensities.resize( continuum_energies_in.size() );
	for( k=0; k<continuum_energies_in.size(); k++ ) {
        continuum_energies[k] = continuum_energies_in[k];
        continuum_intensities[k] = continuum_intensities_in[k];
        if( continuum_energies_in[k] > max_energy ) max_energy = continuum_energies_in[k];
	}
    //  Set tube voltage to input value or maximum energy found in input vectors
	if( kVIn > 0 ) tubeVoltage = kVIn * 1000.0;
	else tubeVoltage = max_energy;
	setmmBe ( mmBeIn );
	tubeCurrent = tubeCurrent_in;
};

//		Radioisotope source

XraySource::XraySource( XrayIsotope isotopeIn, float activityIn, float mmBeIn ) {
//		default setttings for variables asociated with X-ray tube
	tube = false;
	sr = false;
	external = false;
	Element Be( 4 );
	XrayMaterial temp_window( Be );
	window = temp_window;
	setIncAngle ( 90.0 );
	setTakeoffAngle ( 90.0 );
	tubeCurrent = 1;
	endWindow = false;
//		allow for the possibility of a Be window on the isotope source
	setmmBe ( mmBeIn );
	if( activityIn > 0 ) activity = activityIn; else activity = 0;
//		set the target variable to the isotope element
	switch( isotopeIn ) {
		case Cd109:	{
			Element isotopeElement( 47 );	//	Ag <= Cd after electron capture
            XrayMaterial temp_anode( isotopeElement );
            target = temp_anode;
			tubeVoltage = 25515;	//	Ag K edge + 1 eV
			break;
		};
		default: {
			activity = 0;
			tubeVoltage = 0;
		};
	};
};

//		synchrotron source: single energy with intensity input in photons per second
XraySource::XraySource( const float energyIn, const float activityIn, float mmBeIn ) {
	sr = true;
	tube = false;
	external = false;
	setIncAngle ( 90.0 );
	setTakeoffAngle ( 90.0 );
	tubeCurrent = 1;
	endWindow = false;
//		allow for the possibility of a Be window on the source
	Element Be( 4 );
	XrayMaterial temp_window( Be );
	window = temp_window;
	setmmBe ( mmBeIn );
	if( activityIn > 0 ) activity = activityIn; else activity = 0;
	if( energyIn > 0 ) tubeVoltage = energyIn;
};


//		functions for radioisotope source calculations

int XraySource::isotopeLines ( vector <XrayLines> &lines, const float eMin ) const {
//		calculate the intensity of emission lines in photons per second per steradian
//			for a radioisotope source of the specified activity in Bq
//		resizes and loads vector of XrayLines and sets intensity factors in XrayLines class

//		generate vector of XrayEdge objects for target element
//		include only levels (edges) which can be excited (use utility function)
//		loads vector of XrayEdge objects which can be excited by decay
	lines.clear();
	if( activity <= 0 ) return 0;
	vector <XrayEdge> edges;
//		now calculate the intensity of each line, using radioisotope activity
//		put the calculated intensity into the factor included in the class data for that purpose
	unsigned int edgeIndex;
	for ( edgeIndex=0; edgeIndex<edges.size(); edgeIndex++ ) {
//			generate Xray Lines object which gives lines emitted by vacancy at this edge
		XrayLines thisLine ( edges[edgeIndex] );
		int nlEdge = thisLine.numberOfLines();
//			skip it if no lines with any intensity
		if ( nlEdge <= 0 ) continue;
        float max_line_energy = 0;
		int lineIndex;
		for ( lineIndex=0; lineIndex<nlEdge; lineIndex++ ) {
			float lineEnergy = thisLine.energy(lineIndex);
			if( lineEnergy > max_line_energy ) max_line_energy = lineEnergy;
//				ignore self-absorption for now
			if( edges[edgeIndex].index() == K1 ) {
//					K edge - activity directly excites these lines
				float intensityFactor = activity
					* window.transmission(lineEnergy);
				thisLine.factor(lineIndex, intensityFactor );
			} else {
//					excitation by cascade from K emissions
			};
		};
        //  ignore if less than minimum energy for spectrum analysis
        if( max_line_energy < eMin ) continue;
//			insert into output vector
		lines.insert ( lines.end(), thisLine );
	};
	return lines.size();
};


//		functions for synchrotron source calculations

int XraySource::srLines( vector <XrayLines> &lines, const float eMin ) const {
//		instantiate modified XrayLines object using given energy
	XrayLines srLine( tubeVoltage );
//		put intensity into factors and keep track of max energy
    float max_line_energy = 0;
	int i;
	for( i=0; i<srLine.numberOfLines(); i++ ) {
        srLine.factor( i, activity );
        if( srLine.energy(i) > max_line_energy ) max_line_energy = srLine.energy(i);
    }
    //  ignore if less than minimum energy for spectrum analysis
    if( max_line_energy < eMin ) return lines.size();
//		insert into output vector
	lines.insert ( lines.end(), srLine );
	return lines.size();
};




//	functions used to calculate continuum intensity at a given energy

float XraySource::continuum_Ebel ( const float energy ) const {
	if ( ! tube ) return 0.0;
	if ( energy >= tubeVoltage ) return 0.0;
	float CONST = 1.35e9;	// photons per sec. per steradian per milliAmp per keV
	float u0 = tubeVoltage / energy;
	float x = 1.109 - 0.00435 * target.avgZ() + 0.00175 * (tubeVoltage/1000.0);
	float sigma;
	sigma = target.avgZ() * pow ( u0-1, x );
	float corr;
	corr = absCorr_Ebel ( energy );
	return CONST * 1.0f * sigma * corr * window.transmission(energy) * tubeCurrent;
};



float XraySource::absCorr_Ebel ( const float energy ) const {
//		calculate target absorption correction
//			This version based on work of Love and Scott
//			equidistribution versus depth (quadrangle function vs depth)
	float tau = target.photo(energy);
	float rhoZbar = rhozbar_Ebel ( energy );
	float tau_term = 2.0 * tau * rhoZbar * incSin / takeoffSin;
	float fp;
	if( ! endWindow ) {
//			side window X-ray tube, use Ebel's expression
		fp = ( 1.0 - exp ( -tau_term ) ) / tau_term;
	} else {
//			end window X-ray tube, use modified expression
		float rhoZanode_term = target.mass_thickness() * tau / takeoffSin;
		if( rhoZanode_term > tau_term ) {
			fp = exp( -rhoZanode_term ) * ( exp ( tau_term ) - 1 ) / tau_term;
		} else {
			fp = ( 1.0 - exp ( -rhoZanode_term ) ) / tau_term;
		};
	};
	return fp;
};


float XraySource::rhozbar_Ebel ( const float energy ) const {
// average depth of generation of x-rays per Love and Scott distribution
	float v = tubeVoltage/1000.0;
	float u0 = tubeVoltage / energy;
	float m;
	m = 0.1382 - ( 0.9211 / sqrtf ( float(target.avgZ()) ) );
	float logZ = logf ( float(target.avgZ()) );
	float eta ;
	eta = pow ( v, m ) * ( 0.1904 - 0.2236 * logZ + 0.1292 * logZ * logZ
		- 0.0149 * logZ * logZ * logZ );
	float rhozm;
	rhozm = target.avgAoverZ();
	float j = 0.00135 * target.avgZ();
	rhozm *= ( 0.787e-5 * sqrt ( float(j) ) * pow ( v, 1.5f) + 0.735e-6f * v * v );
	float rhoZbar_ratio;
	rhoZbar_ratio = ( 0.49269 - 1.0987 * eta + 0.78557 * eta * eta ) * log ( u0 );
	rhoZbar_ratio /= ( 0.70256 - 1.09865 * eta + 1.0046 * eta * eta + log ( u0) );
	return rhoZbar_ratio * rhozm;
};


//	functions used to calculate characteristic line intensities

int XraySource::tubeEdges ( vector <XrayEdge> &edges ) const {
//		loads vector of XrayEdge objects which are excited by electrons of
//			specified voltage hitting target element
//		this is also used for radioisotopes by making an appropriate voltage setting
    edges.clear();
    //  Get list of elements in the target material
    const vector <Element> &elements = target.element_list();
    unsigned int ie;
    for( ie=0; ie<elements.size(); ie++ ) {
        int ne;
        vector <EdgeIndex> edgeList;
        XrayEdge::numberOfEdges ( edgeList, elements[ie], tubeVoltage );
        ne = edgeList.size();
        int edgeIndex;
        for ( edgeIndex=0; edgeIndex<ne; edgeIndex++ ) {
            XrayEdge thisEdge ( elements[ie], edgeList[edgeIndex] );
            edges.insert ( edges.end(), thisEdge );
        };
    };
	return edges.size();
};

int XraySource::tubeLinesEbel ( vector <XrayLines> &lines, const float eMin ) const {
//		calculate the intensity of emission lines in photons per second per steradian per milliAmp
//			for electron beam striking single-element target, as in X-ray tube
//		resizes and loads vector of XrayLines and sets intensity factors in XrayLines class

	if ( ! tube ) {
		lines.clear();
		return 0;
	};
//		generate vector of XrayEdge objects for target element
//		include only levels (edges) which can be excited (use utility function)
	int ne;
	vector <XrayEdge> lineEdges(0);
	ne = tubeEdges( lineEdges );
	int edgeIndex;
//		calculate some things which depend only on target and geometry
//	float CONSTKL = 6e13f;		// photons per second per steradian per milliAmp
	float CONSTKL = 8.11e13f;		// photons per second per steradian per milliAmp
	float z = target.avgZ();
	float j = 13.5 * z;
//		for each edge, calculate the number of vacancies produced and
//			multiply by other parameters which don't depend on which line
	vector <float> vacancies(ne);
	for ( edgeIndex=0; edgeIndex<ne; edgeIndex++ ) {
		float ec = lineEdges[edgeIndex].energy();
		float zk = lineEdges[edgeIndex].degeneracy();
		float u0 = tubeVoltage / ec;
//			B sub k depends on whether it is a K or L edge
		float bk;
		float corr_fac = 1;
		switch ( lineEdges[edgeIndex].level() ) {
			case K:	bk = 0.35f;	corr_fac = 1.0f;    break;
			case L:	bk = 0.25f;	corr_fac = 1.0f;    break;
			default:	bk = 0.2f;	break;			//   wild guess
		};
		float log_u0 = log ( u0 );
		float oneOverS = zk * bk * ( u0 * log_u0 + 1 - u0 ) / z;
		float bigFraction = sqrt(u0) * log_u0 + 2.0 * ( 1.0 - sqrt(u0) );
		bigFraction /= u0 * log_u0 + 1.0 - u0;
		oneOverS *= 1 + 16.05  * sqrt ( j / ec ) * bigFraction;
		float r = 1 - 0.008151 * z + 3.613e-5 * z*z;
		r += 0.009583 * z * exp ( -u0 ) + 0.001141 * ec / 1000.0;
		vacancies[edgeIndex] = CONSTKL * corr_fac * oneOverS * r;
	};
//	Coster_Kronig transitions from primary to secondary edges
/*	for ( edgeIndex=0; edgeIndex<ne; edgeIndex++ ) {
		int secEdgeIndex;
		for ( secEdgeIndex=0; secEdgeIndex<ne; secEdgeIndex++ ) {
            if( secEdgeIndex == edgeIndex ) continue;
			float cktemp = lineEdges[edgeIndex].cktotal( lineEdges[secEdgeIndex] );
			if ( cktemp <= 0.0 ) continue;
//			vacancies[secEdgeIndex] += cktemp * vacancies[edgeIndex];
			vacancies[edgeIndex] -= cktemp * vacancies[edgeIndex];
        }
	}
*/
//		now calculate the intensity of each line, using above info
//			equation from Ebel does not consider Coster-Kronig transition probabilities
//		put the calculated intensity into the factor included in the class data for that purpose
	for ( edgeIndex=0; edgeIndex<ne; edgeIndex++ ) {
//			generate Xray Lines object which gives lines emitted by vacancy at this edge
		XrayLines thisLine ( lineEdges[edgeIndex] );
		int nlEdge = thisLine.numberOfLines();
//			skip it if no lines with any intensity
		if ( nlEdge <= 0 ) continue;
        float max_line_energy = 0;
		int lineIndex;
		for ( lineIndex=0; lineIndex<nlEdge; lineIndex++ ) {
			float lineEnergy = thisLine.energy(lineIndex);
			if( lineEnergy > max_line_energy ) max_line_energy = lineEnergy;
            float intensityFactor = 0;
            //  Ignore if less than minimum energy for spectrum analysis
            if( lineEnergy > eMin ) {
                intensityFactor = vacancies[edgeIndex]  * thisLine.edge().yield()
                    * absCorr_Ebel ( lineEnergy )
                    * window.transmission(lineEnergy)
                    * tubeCurrent;
            }
            //  Fix database problem for 4d transition metals (especially Rh)
            if( Rh_adjust_min_Z <= thisLine.edge().element().Z() && thisLine.edge().element().Z() <= Rh_adjust_max_Z ) {
                if( thisLine.symbolIUPAC( lineIndex ) == "L3-N4,5" ) intensityFactor *= Rh_Lb2_15_factor;
            }
			thisLine.factor(lineIndex, intensityFactor );
//            cout << "tube lines Ebel " << thisLine.symbolSiegbahn(lineIndex) << "   " << thisLine.symbolIUPAC(lineIndex) << "   " << thisLine.energy(lineIndex)/1000 << "   " << thisLine.intensity(lineIndex)/1e8 << endl;
		};
        //  ignore if less than minimum energy for spectrum analysis
        if( max_line_energy < eMin ) continue;
//			insert into output vector
		lines.insert ( lines.end(), thisLine );
	};
	return lines.size();
};


//      X-ray tube with calculated spectrum read in from a file (lines and continuum)

int XraySource::extLines( vector <XrayLines> &lines, const float eMin ) const {
    lines.clear();
	if ( ! external ) return 0;
    //  Copy X-ray emission lines from object storage to argument
    float max_line_energy = 0;
    unsigned int k;
    for( k=0; k<tube_lines_ext.size(); k++ ) {
        XrayLines thisLine = tube_lines_ext[k];
        //  Skip it if no lines with any intensity
		if ( thisLine.numberOfLines() <= 0 ) continue;
        int lineIndex;
        for ( lineIndex=0; lineIndex<thisLine.numberOfLines(); lineIndex++ ) {
            float lineEnergy = thisLine.energy(lineIndex);
            if( lineEnergy > max_line_energy ) max_line_energy = lineEnergy;
            //  Include beryllium window absorption and emission current (assumes intensities are in ph/sec/sr/mA)
            float intensityFactor = thisLine.factor( lineIndex );
                intensityFactor *= window.transmission(lineEnergy) * tubeCurrent;
            thisLine.factor(lineIndex, intensityFactor );
        }
        //  ignore if less than minimum energy for spectrum analysis
        if( max_line_energy < eMin ) continue;
        lines.push_back( thisLine );
    }
	return lines.size();
};

float XraySource::continuum_ext ( const float energy ) const {
	if ( ! external ) return 0.0;
	if ( energy >= tubeVoltage ) return 0.0;
    //  Interpolate continuum intensity at desired energy
	float inten =  interp( energy, continuum_energies, continuum_intensities );
	//  Be window and emission current (assumes continuum intensities are ph/sec/keV/sr/mA
	return inten * window.transmission(energy) * tubeCurrent;
};


//      X-ray tube with calculated spectrum read in from a file (lines and continuum)

int XraySource::tubeLinesSewell ( vector <XrayLines> &lines, const float eMin ) const {
//		Calculate the intensity of emission lines in photons per second per steradian per milliAmp
//			for electron beam striking single-element target, as in X-ray tube
//		Resizes and loads vector of XrayLines and sets intensity factors in XrayLines class

	if ( ! tube ) {
		lines.clear();
		return 0;
	};
	if( endWindow ) return -1;  //  Don't have correct absorption correction yet
//		Generate vector of XrayEdge objects for target element
//		Include only levels (edges) which can be excited (use utility function)
	int ne;
	vector <XrayEdge> lineEdges(0);
	ne = tubeEdges( lineEdges );
	vector <float> vacancies(ne);
	int edgeIndex;
//		Calculate some things which depend only on target and geometry
	float z = target.avgZ();
    float z_a = target.avgZoverA();
    float j = Sewell_j( z );   //  Note j must be averaged as log(j) for multi-element targets
    float eta = Sewell_eta( z, tubeVoltage );
    float tilt = PI/2 - inc_angle;
    if( tilt < 0 ) return -2;
//		For each edge, calculate the number of vacancies produced and
//			multiply by other parameters which don't depend on which line
	for ( edgeIndex=0; edgeIndex<ne; edgeIndex++ ) {
		float ec = lineEdges[edgeIndex].energy();
		float u0 = tubeVoltage / ec;
		if( u0 <= 1 ) continue;
        float R = Sewell_r( u0, eta, tilt );
        float S = Sewell_S_lines( u0, j/1000, ec/1000, z_a );
        float occ = lineEdges[edgeIndex].occupancy();
        vacancies[edgeIndex] = occ * R * S;
    }
//	Coster_Kronig transitions from primary to secondary edges
	for ( edgeIndex=0; edgeIndex<ne; edgeIndex++ ) {
		int secEdgeIndex;
		for ( secEdgeIndex=0; secEdgeIndex<ne; secEdgeIndex++ ) {
            if( secEdgeIndex == edgeIndex ) continue;
			float cktemp = lineEdges[edgeIndex].cktotal( lineEdges[secEdgeIndex] );
			if ( cktemp <= 0.0 ) continue;
//			vacancies[secEdgeIndex] += cktemp * vacancies[edgeIndex];
			vacancies[edgeIndex] -= cktemp * vacancies[edgeIndex];
        }
	}
        //		Now calculate the intensity of each line
        //		Put the calculated intensity into the factor included in the class data for that purpose
//			Generate Xray Lines object which gives lines emitted by vacancy at this edge
	for ( edgeIndex=0; edgeIndex<ne; edgeIndex++ ) {
		XrayLines thisLine ( lineEdges[edgeIndex] );
		int nlEdge = thisLine.numberOfLines();
//			skip it if no lines with any intensity
		if ( nlEdge <= 0 ) continue;
		float ec = lineEdges[edgeIndex].energy();
		float u0 = tubeVoltage / ec;
        float yield = lineEdges[edgeIndex].yield();
//			Calculate model parameters for depth distribution (phi-rho-z function)
        float pz = Sewell_pz( j/1000, tubeVoltage/1000, eta, u0, z_a, tilt );
        float h = Sewell_h( u0, z, eta, tilt );
        float pz_m = Sewell_pz_m(pz, u0, z, tilt );
        float pz_r = Sewell_pz_r( pz, pz_m, h );
        float max_line_energy = 0;
		int lineIndex;
		for ( lineIndex=0; lineIndex<nlEdge; lineIndex++ ) {
			float lineEnergy = thisLine.energy(lineIndex);
			if( lineEnergy > max_line_energy ) max_line_energy = lineEnergy;
            float intensityFactor = 0;
            //  Ignore if less than minimum energy for spectrum analysis
            if( lineEnergy > eMin ) {
                float sigma = target.photo(lineEnergy);
                float chi = sigma / takeoffSin;
                float f = Sewell_f( chi, pz_m, pz_r, pz, h );
                intensityFactor = vacancies[edgeIndex] * yield * f * window.transmission(lineEnergy) * tubeCurrent;
            }
            //  Fix database problem for 4d transition metals (especially Rh)
            if( Rh_adjust_min_Z <= thisLine.edge().element().Z() && thisLine.edge().element().Z() <= Rh_adjust_max_Z ) {
                if( thisLine.symbolIUPAC( lineIndex ) == "L3-N1" ) intensityFactor *= Rh_Lb6_factor;   //  Lb6
                if( thisLine.symbolIUPAC( lineIndex ) == "L3-N4,5" ) intensityFactor *= Rh_Lb2_15_factor; //  Lb2,15
                if( thisLine.symbolIUPAC( lineIndex ) == "L3-M1" ) intensityFactor *= Rh_Ll_factor; //  Ll
                if( thisLine.symbolIUPAC( lineIndex ) == "L2-M1" ) intensityFactor *= Rh_Ln_factor; //  Ln
                if( thisLine.symbolIUPAC( lineIndex ) == "L2-N4" ) intensityFactor *= Rh_Lg1_factor; //  Lg1
                if( thisLine.symbolIUPAC( lineIndex ) == "L1-M2" || thisLine.symbolIUPAC( lineIndex ) == "L1-M3" ) intensityFactor *= Rh_Lb34_factor; //  Lb3,4
                if( thisLine.symbolIUPAC( lineIndex ) == "L1-N2" || thisLine.symbolIUPAC( lineIndex ) == "L1-N3" ) intensityFactor *= Rh_Lg23_factor; //  Lg2,3
            }
            //  Adjustment to get calculated XRF intensities below 3 keV to match measurements of pure compounds
            //if( thisLine.edge().element().Z() == 45 && thisLine.edge().level() == L ) intensityFactor *= 1.85f;
			thisLine.factor(lineIndex, intensityFactor );
		};
        //  Ignore entire set if less than minimum energy for spectrum analysis
        if( max_line_energy < eMin ) continue;
//			Insert into output vector
		lines.insert ( lines.end(), thisLine );
	};
	return lines.size();
};


/*
float XraySource::continuum_Sewell( const float energy ) const {
	if ( ! tube ) return 0.0;
	if ( energy >= tubeVoltage ) return 0.0;
	if( endWindow ) return 0.0;  //  Don't have correct absorption correction yet
	float CONST = 1.35e9;	// photons per sec. per steradian per milliAmp per keV
	float u0 = tubeVoltage / energy;
	float z = target.avgZ();
    float z_a = target.avgZoverA();
    float j = Sewell_j( z );   //  Note j must be averaged as log(j) for multi-element targets
    float eta = Sewell_eta( z, tubeVoltage );
    float tilt = PI/2 - inc_angle;
    if( tilt < 0 ) return -2;
    float R = Sewell_r( u0, eta, tilt );
    float S = Sewell_S_continuum( u0, j, energy, z_a );
    float pz = Sewell_pz( j/1000, tubeVoltage/1000, eta, u0, z_a, tilt );
    float h = Sewell_h( u0, z, eta, tilt );
    float pz_m = Sewell_pz_m(pz, u0, z, tilt );
    float pz_r = Sewell_pz_r( pz, pz_m, h );
	float sigma = target.photo(energy);
    float chi = sigma / takeoffSin;
//    float f = Sewell_f( chi, pz_m, pz_r, pz, h );
//    float xEbel = 1;
	float xEbel = 1.109 - 0.00435 * target.avgZ() + 0.00175 * (tubeVoltage/1000.0);
    float crossover = 1.0f * j;
    // photons per sec. per steradian per milliAmp per keV
    float const_continuum = 1.35e9; //  Ebel 1999
//    float const_continuum = 2.06e9; //  sigmaEbel * R * f
//    float const_continuum = 1.83e9;
    float sigmaJ = const_continuum * z * pow ( tubeVoltage / crossover - 1, xEbel );
    float sigmaEbel = const_continuum * z * pow ( u0-1, xEbel );
//    sigmaEbel = ( 1 - exp ( - energy / crossover ) ) * sigmaEbel + exp ( - energy / crossover ) * sigmaJ;
//	float continuum = S * R * f;   //  almost OK at peak, too low at high energy, much worse with / ( energy * energy )
//	float continuum = S * f;   //    way off, too big at peak and too low at high energy
//	float continuum = sigmaEbel * R * f; //  shifts peak ~1keV too high
//	float continuum = sigmaEbel * f;  //  too small at peak, good at high energy, maybe a little better with xEbel = 1 but barely
    float   f = absCorr_Ebel ( energy );
	float continuum = sigmaEbel * f;  //  Ebel 1999 works very best for continuum
//    if( fabs( energy - 10000 ) < 100 ) cout << "Sewell cont " << energy << "  " << continuum << "  " << f << "  " << R << "  " << S << endl;
	return continuum * window.transmission(energy) * tubeCurrent;
};
*/

string XraySource::toString() const
{
    ostringstream os;
    os << "XraySource: " << endl;
    os << "  tube=" << tube << endl;
    os << "  sr=" << sr << endl;
    os << "  target=" << endl << target.toString() << endl;
    os << "  tubeVoltage=" << tubeVoltage << endl;
    os << "  incSin=" << incSin << endl;
    os << "  takeoffSin=" << takeoffSin << endl;
    os << "  tubeCurrent=" << tubeCurrent << endl;
    os << "  activity=" << activity << endl;
    os << "  targetAbs=" << "TODO" << endl;
    os << "  endWindow=" << endWindow << endl;
    os << "  beAbs=" << "TODO" << endl;
    return os.str();
}

