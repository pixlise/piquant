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

//
//  fpSetupConditions.cpp
//  PIQUANT
//
//  Created by W. T. Elam on 1/27/2017.
//  Copyright (c) 2017 APL/UW. All rights reserved.
//
#include "Element.h"
#include "XRFconstants.h"
#include "XrayDetector.h"
#include "read_tube_spectrum.h"
#include "parse_element_list.h"
#include "fpSetupConditions.h"

//  Sets up XRF measurement conditions from input array of float values
//  Default values for solid angles and geometry are set here (to unity)
//	Written Jan. 27, 2017 from XRFanalysisDLL.cpp, condQuant function
//      Extensive re-write to use XrayMaterials class for most elements in beam
//      also add extra information from ISO standard version of EMSA format
//  Modified May 26, 2017 To add optic type 5 for new breadboard
//  Modified June 26, 2017 to trap error for bad optic file and return error code
//  Modified Sept. 29, 2017
//      Changed window type 3 from Brass to Carbon Fiber Reinforced Polymer (composition unknown, pure C for now)
//  Modified May 22, 2020 to put conditions vector and optic file name in struct
//                          add file name for X-ray tube spectrum input from external calculation
//  Modified Nov. 2, 2020   Add PIXL FM optic type, number 7
//  Modified Apr. 28, 2020  Add PIXL FM optic type, number 8 (calculated with correct Be window for X-ray tube)
//  Modified May 14, 2021   Move shelf factor and slope to XrayDetector and control via -T option
//  Modified May 14, 2021   Move shelf factor and slope to XrayDetector and control via -T option (via conditions struct)


using namespace std;

int fpSetupConditions ( const XRFconditionsInput &cond_in, XRFconditions &conditions_out ) {

//		X-ray source
	int anodeZ = cond_in.conditionsVector[ANODE_Z_INDEX];
	if( cond_in.tube_file_name.length() > 0 ) {
        //  X-ray tube source with intensities read in from file (from external calculation)
        //  Read file of intensities and return error if file can't be read
        vector <XrayLines> tube_lines_in;
        vector <float> continuum_en_in;
        vector <float> continuum_int_in;
        string title;
        float kV = 0;
        int result = read_tube_spectrum( cond_in.tube_file_name, tube_lines_in, kV,
                            continuum_en_in, continuum_int_in, title );
        if( result < 0 ) return -100 - (XRF_PARAMETER_TUBE_FILE);
        if( title.length() > 0 ) conditions_out.tube_file_title = title;
		if( cond_in.conditionsVector[KV_INDEX] < 0 ) return -100 - (KV_INDEX);
		if( cond_in.conditionsVector[KV_INDEX] > 0 ) kV = cond_in.conditionsVector[KV_INDEX];
		if( cond_in.conditionsVector[TUBE_BE_WINDOW_INDEX] < 0 ) return -100 - (TUBE_BE_WINDOW_INDEX);
		float mmBe = cond_in.conditionsVector[TUBE_BE_WINDOW_INDEX];
		float tubeCurrent = 0;
		if( cond_in.conditionsVector[TUBE_CURRENT_INDEX] <= 0 ) return -100 - (TUBE_CURRENT_INDEX);
        tubeCurrent = cond_in.conditionsVector[TUBE_CURRENT_INDEX];    //  milliAmps
		XraySource tube_ext_source( tube_lines_in, continuum_en_in, continuum_int_in, mmBe, kV, tubeCurrent );
		conditions_out.source = tube_ext_source;
	} else if( Element::check_Z( anodeZ ) ) {
		Element anode(anodeZ);
//			x-ray tube kilovolts, takeoff angle, and incidence angle
		if( cond_in.conditionsVector[KV_INDEX] <= 0 ) return -100 - (KV_INDEX);
		float kV = cond_in.conditionsVector[KV_INDEX];
		if( cond_in.conditionsVector[TUBE_INC_ANGLE_INDEX] <= 0 || cond_in.conditionsVector[TUBE_INC_ANGLE_INDEX] > 90 ) return -100 - (TUBE_INC_ANGLE_INDEX);
		float incAngle = cond_in.conditionsVector[TUBE_INC_ANGLE_INDEX];
		if( cond_in.conditionsVector[TUBE_TAKEOFF_ANGLE_INDEX] < -90 || cond_in.conditionsVector[TUBE_TAKEOFF_ANGLE_INDEX] > 90 ) return -100 - (TUBE_TAKEOFF_ANGLE_INDEX);
		if( cond_in.conditionsVector[TUBE_TAKEOFF_ANGLE_INDEX] == 0 ) return -100 - (TUBE_TAKEOFF_ANGLE_INDEX);
		float takeoffAngle = cond_in.conditionsVector[TUBE_TAKEOFF_ANGLE_INDEX];
//			Be window thickness in millimeters
		if( cond_in.conditionsVector[TUBE_BE_WINDOW_INDEX] < 0 ) return -100 - (TUBE_BE_WINDOW_INDEX);
		float mmBe = cond_in.conditionsVector[TUBE_BE_WINDOW_INDEX];
		float tubeCurrent = 0;
		if( cond_in.conditionsVector[TUBE_CURRENT_INDEX] < 0 ) return -100 - (TUBE_CURRENT_INDEX);
		if( cond_in.conditionsVector[TUBE_CURRENT_INDEX] > 0 ) {
			tubeCurrent = cond_in.conditionsVector[TUBE_CURRENT_INDEX];    //  milliAmps
		} else {
//			put in default value of 20 microAmps for compatibility with previous versions
			tubeCurrent = 0.020f;
		};
//			Instantiate x-ray source object (this one is an x-ray tube)
//			Negative takeoff angle to indicate end window X-ray tube with 1.2 micron thick anode
//				(anode thickness determined from scatter calculations for NRIXS experiments, changed May 28, 2009   W. T. Elam)
		XraySource source(anode, kV, incAngle, takeoffAngle, mmBe, 1.2f, tubeCurrent );
		conditions_out.source = source;
	} else if( anodeZ == -1 ) {
//			Monochromatic synchrotron source - set up as a fake characteristic line with one energy
		float energy = cond_in.conditionsVector[KV_INDEX] * 1000;
		float intensity = cond_in.conditionsVector[TUBE_INC_ANGLE_INDEX];
		float mmBe = cond_in.conditionsVector[TUBE_BE_WINDOW_INDEX];
		XraySource source( energy, intensity, mmBe );
		conditions_out.source = source;
	} else if( anodeZ == 0 && cond_in.anode_element_list.length() > 0 ) {
        //  Anode is an element list, parse it and make an XrayMaterial for the anode
        vector <ElementListEntry> list_of_anode_elements;
        bool dummy; //  Required for carbonates in unknowns, not used here
        parse_element_list( cond_in.anode_element_list, list_of_anode_elements, dummy );
        vector <Element> anode_elements;
        vector <float> anode_fractions;
        int i;
        for( i=0; i<list_of_anode_elements.size(); i++ ) {
            anode_elements.push_back( list_of_anode_elements[i].element );
            anode_fractions.push_back(  list_of_anode_elements[i].percent / 100 );
        }
        XrayMaterial anode_mat( anode_elements, anode_fractions );
//			x-ray tube kilovolts, takeoff angle, and incidence angle
		if( cond_in.conditionsVector[KV_INDEX] <= 0 ) return -100 - (KV_INDEX);
		float kV = cond_in.conditionsVector[KV_INDEX];
		if( cond_in.conditionsVector[TUBE_INC_ANGLE_INDEX] <= 0 || cond_in.conditionsVector[TUBE_INC_ANGLE_INDEX] > 90 ) return -100 - (TUBE_INC_ANGLE_INDEX);
		float incAngle = cond_in.conditionsVector[TUBE_INC_ANGLE_INDEX];
		if( cond_in.conditionsVector[TUBE_TAKEOFF_ANGLE_INDEX] < -90 || cond_in.conditionsVector[TUBE_TAKEOFF_ANGLE_INDEX] > 90 ) return -100 - (TUBE_TAKEOFF_ANGLE_INDEX);
		if( cond_in.conditionsVector[TUBE_TAKEOFF_ANGLE_INDEX] == 0 ) return -100 - (TUBE_TAKEOFF_ANGLE_INDEX);
		float takeoffAngle = cond_in.conditionsVector[TUBE_TAKEOFF_ANGLE_INDEX];
//			Be window thickness in millimeters
		if( cond_in.conditionsVector[TUBE_BE_WINDOW_INDEX] < 0 ) return -100 - (TUBE_BE_WINDOW_INDEX);
		float mmBe = cond_in.conditionsVector[TUBE_BE_WINDOW_INDEX];
		float tubeCurrent = 0;
		if( cond_in.conditionsVector[TUBE_CURRENT_INDEX] < 0 ) return -100 - (TUBE_CURRENT_INDEX);
		if( cond_in.conditionsVector[TUBE_CURRENT_INDEX] > 0 ) {
			tubeCurrent = cond_in.conditionsVector[TUBE_CURRENT_INDEX];    //  milliAmps
		} else {
//			put in default value of 20 microAmps for compatibility with previous versions
			tubeCurrent = 0.020f;
		};
//			Instantiate x-ray source object (this one is an x-ray tube)
//			Negative takeoff angle to indicate end window X-ray tube with 1.2 micron thick anode
//				(anode thickness determined from scatter calculations for NRIXS experiments, changed May 28, 2009   W. T. Elam)
		XraySource source(anode_mat, kV, incAngle, takeoffAngle, mmBe, 1.2f, tubeCurrent );
		conditions_out.source = source;
	} else {
        return -100 - (ANODE_Z_INDEX);
    }
//		Instantiate x-ray source object => Cd-109 radioisotope
//			Activity 10 microCurie = 3.7e6 Bq
//			Leave Be window same as Ag anode X-ray tube      June 25, 2009   W. T. Elam)
//	XraySource source( Cd109, 3.7e6, mmBe );
//	conditions_out.source = source;

//		incident beam filter (if any)
	int filterZ = cond_in.conditionsVector[FILTER_Z_INDEX];
	if( ! Element::check_Z( filterZ ) ) return -100 - (FILTER_Z_INDEX);
	if( cond_in.conditionsVector[FILTER_THICK_INDEX] < 0 ) return -100 - (FILTER_THICK_INDEX);
	float filterThickness = cond_in.conditionsVector[FILTER_THICK_INDEX] * CM_MICRON;   //  Filter thickness is in microns
	Element filter_element( filterZ );
	XrayMaterial filter( filter_element );
	filter.thickness( filterThickness );
	conditions_out.filter = filter;

//		Xray Optic
	// NY Edit, Added an opt
	int opticType = cond_in.conditionsVector[TEST_OPTIC_TYPE_INDEX];
	//cout << "Optical Type is (DLL):\t" << opticType << endl;
	switch( opticType ) {
		case 0:
		case 1:
		{
			//  No optic, use default constructor
			break;
		}
		case 2:
		{
			XrayOptic tempOptic(24000, 4000, 0.001, BOXCAR);
			conditions_out.optic = tempOptic;
			break;
		}
		case 3:
		{
			XrayOptic tempOptic(0, 0, 1, PIXL);    //  changed to zeros for center and bandwidth WTE Sep. 30, 2013
			conditions_out.optic = tempOptic;
			break;
        }
		case 4:
		{
            //  Use try and catch in case file can't be opened
            try {
                XrayOptic tempOptic( cond_in.optic_file_name );    //  changed to zeros for center and bandwidth WTE Sep. 30, 2013
                conditions_out.optic = tempOptic;
            } catch (exception& e) {
                return -100 - (XRF_PARAMETER_OPTIC_FILE);
            }
			break;
        }
		case 5:
		{
			XrayOptic tempOptic(0, 0, 1, NEW_BB);    //  added May 26, 2017 with efficiency curve for new breadboard from Chris Heirwegh
			conditions_out.optic = tempOptic;
			break;
        }
		case 7:
		{
			XrayOptic tempOptic(0, 0, 1, PIXL_FM_OPTIC_OLD);    //  added Nov. 2, 2020 with efficiency curve for PIXL FM (in work)
			conditions_out.optic = tempOptic;
			break;
        }
		case 8:
		{
			XrayOptic tempOptic(0, 0, 1, PIXL_FM_OPTIC);    //  added Apr. 28, 2021 with efficiency curve for PIXL FM (calculated with correct Be window for X-ray tube)
			conditions_out.optic = tempOptic;
			break;
        }
		default: 	return -100 - (TEST_OPTIC_TYPE_INDEX);	//  Integer input not found, error return

	};
	// NY Edit end

//      dust on optic

    //**************************************************************************
    //  Dust composition and thickness from J.L. Campbell et al.,
    //  Nuclear Instruments and Methods in Physics Research B 323 (2014) 49–58.
    //  Page 57, Section 12. Conclusions  (see text, note sulfur is not in Table 3)
    //  The film analyzed on Sol 34 (thickness ~100 nm) can be described by a mixture of
    //  approximately 10 wt.% MgO, 62.5 wt.% Fe2O3 , 3.9 wt.% Na2O, 3.6 wt.% Cl
    //  and 20 wt.% SO3
    //**************************************************************************
/*
    const int n_dust = 5;
    //                              MgO Fe2O3 NaO  Cl   SO3
    const int z_dust[n_dust] =    { 12, 26,   11,  17,  16 };
    const float pct_dust[n_dust] = { 10, 62.5, 3.9, 3.6, 20 };
    const float dust_thickness = 100;   //  nanometers
    float frac_dust[n_dust];
    int id;
    for( id=0; id<n_dust; id++ ) frac_dust[id] = pct_dust[id] / 100;
    XrayMaterial dust( n_dust, z_dust, frac_dust, true );   //  oxides = true
    dust.thickness( dust_thickness * NM_CM );  //  convert nanometers to centimeters
    conditions_out.dust_on_optic = dust;
    conditions_out.dust_on_specimen = dust;
    conditions_out.dust_on_detector = dust;
*/

//		incident atmosphere path
//		atmosphere descriptions - atomic number, weight fractions, and atoms in gas molecule
//				Earth and Mars atmosphere compositions from CRC Handbook of Chemistry and Physics,
//				77th Ed., David R. Lide, Ed., CRC Press (Boca Raton, 1996),
//				ISBN 0-8493-0477-6, p 14-3.     Note that these are VOLUME FRACTIONS
//				Other data from http://nssdc.gsfc.nasa.gov/planetary/factsheet/
//			See spreadsheet EarthMarsAtmosphere.xls for conversions
	int pathType = cond_in.conditionsVector[PATH_TYPE_INDEX];
	switch( pathType ) {
		case VACUUM:	break;  //  Use default constructor
		case HELIUM:    {
                Element he( "He" );
                XrayMaterial path( he );
                path.density( he.atomicWeight() / GAS_MOLE_VOLUME );	//	gm/cm3 at STP
                conditions_out.incidentPath = path;
                break;
            }
		case MARS:	    {
//		            Mars 95.5% CO2, 2.7% N2, 0.2% O2, 1.6% Ar by volume, pressure about 7 milliBars
                const int Mars_N =4;
                const int Mars_Z[Mars_N] = { 6, 7, 8, 18 };	//	C, N, O, Ar
                const float Mars_F[Mars_N] = { 0.265f, 0.017f, 0.707f, 0.011f };
                const float Mars_density = 0.00002f;	//	gm/cm3  corrected, 0.02 kg/m3  Dec. 2, 2013
                XrayMaterial path( Mars_N, Mars_Z, Mars_F );
                path.density( Mars_density );	//	gm/cm3 at STP
                conditions_out.incidentPath = path;
                break;
            }
		case HE_MARS:	return -100 - (PATH_TYPE_INDEX);    break;
		case AIR:	//  These two are the same
		case EARTH:		    {
//		            Earth 78.1% N2, 20.9% O2, 0.9% Ar by volume, pressure one standard atmosphere
                const int Earth_N = 3;
                const int Earth_Z[Earth_N] = { 7, 8, 18 };	//	N2, O2, Ar
                const float Earth_F[Earth_N] = { 0.758f, 0.232f, 0.01f };
                const float Earth_density = 0.00122f;	//	gm/cm3
                XrayMaterial path( Earth_N, Earth_Z, Earth_F );
                path.density( Earth_density );	//	gm/cm3 at STP
                conditions_out.incidentPath = path;
                break;
            }
		default: break;  //  Use default constructor (vacuum)
	};
	if( cond_in.conditionsVector[INC_PATH_LENGTH_INDEX] < 0 ) return -100 - (INC_PATH_LENGTH_INDEX);
	float incPathLength = cond_in.conditionsVector[INC_PATH_LENGTH_INDEX];
	conditions_out.incidentPath.thickness( incPathLength );

//		source geometry
	float solidAngleSource = cond_in.conditionsVector[SOURCE_SOLID_ANGLE_INDEX];
	if ( solidAngleSource < 0 ) return -100 - (SOURCE_SOLID_ANGLE_INDEX);
	if( solidAngleSource == 0 ) solidAngleSource= FOUR_PI;
	conditions_out.solidAngleSource = solidAngleSource / FOUR_PI;   //  Note that this is just a factor in conditions, not steradians
	float excitAngle = cond_in.conditionsVector[EXCIT_ANGLE_INDEX];
	if ( excitAngle < 0 || excitAngle > 90.0 ) return -100 - (EXCIT_ANGLE_INDEX);
	conditions_out.excitAngle = excitAngle;
	float sinExcit = sin ( conditions_out.excitAngle * RADDEG );
	if ( sinExcit < 1.0e-6f ) sinExcit = 1.0e-6f;
	conditions_out.excitCosecant = 1 / sinExcit;

//      specimen window and dust
//      see above for dust info
	int winType = cond_in.conditionsVector[WINDOW_TYPE_INDEX];
	switch( winType ) {
		case NO_WINDOW:	break;      //  Use default constructor
		case B4C:	{
                const int B4C_N =2;
                const int B4C_Z[B4C_N] = { 5, 6 };
                const float B4C_F[B4C_N] = { 0.7826f, 0.2174f };
            //	const float B4C_den = 1.2f;
            //	Fe added to B4C window based on measured spectra April 26, 2006
            //	const int B4C_N =3;
            //	const int B4C_Z[B4C_N] = { 5, 6, 26 };
            //	const float B4C_F[B4C_N] = { 0.7826f, 0.2174f, .025f };
                const float B4C_den = 3.0f;	//	modified based on measured spectra   April 26, 2006
                XrayMaterial window( B4C_N, B4C_Z, B4C_F );
                window.density( B4C_den );	//	gm/cm3 at STP
                conditions_out.window = window;
                break;
            }
		case PLASTIC:	{
                const int Plas_N = 3;
                const int Plas_Z[Plas_N] = { 1, 6, 8 };
                const float Plas_F[Plas_N] = { 0.1f, 0.7f, 0.2f };
                const float Plas_den = 1.2f;
                XrayMaterial window( Plas_N, Plas_Z, Plas_F );
                window.density( Plas_den );	//	gm/cm3 at STP
                conditions_out.window = window;
                break;
            }
		case CFRP:	{
                //	CFRP (Carbon fiber reinforced polymer)
                //	No composition available yet, so use pure carbon
                const int CFRP_N = 1;
                const int CFRP_Z[CFRP_N] = { 6 };
                const float CFRP_F[CFRP_N] = { 1 };
                const float CFRP_den = 2.3f;
                XrayMaterial window( CFRP_N, CFRP_Z, CFRP_F );
                window.density( CFRP_den );	//	gm/cm3 at STP
//                const int Brass_N = 2;
//                const int Brass_Z[Brass_N] = { 29, 30 };
//                const float Brass_F[Brass_N] = { 0.70f, 0.30f };
//                const float Brass_den = 8.0f;
//                XrayMaterial window( Brass_N, Brass_Z, Brass_F );
//                window.density( Brass_den );	//	gm/cm3 at STP
                conditions_out.window = window;
                break;
            }
		case ZR:	{	    //	zirconium
                Element zr( "Zr" );
                XrayMaterial window( zr );
                conditions_out.window = window;
                break;
            }
		case AL:	{	    //	Aluminum
                Element al( "Al" );
                XrayMaterial window( al );
                conditions_out.window = window;
                break;
            }
		case NYLON:	{
                //		http://physics.nist.gov/cgi-bin/Star/compos.pl?matno=210
                const int Nylon_N = 4;
                const int Nylon_Z[Nylon_N] = { 1, 6, 7, 8 };
                const float Nylon_F[Nylon_N] = { 0.107062f, 0.680449f, 0.099189f, 0.113300f };
                const float Nylon_den = 1.14f;
                XrayMaterial window( Nylon_N, Nylon_Z, Nylon_F );
                window.density( Nylon_den );	//	gm/cm3 at STP
                conditions_out.window = window;
                break;
            }
		case NYLON_ZR:	{
                //		for MTXRF resin from Oli via Grundl, ZrO2 loaded version
                const int NylonZr_N = 5;
                const int NylonZr_Z[NylonZr_N] = { 1, 6, 7, 8, 40 };
                const float NylonZr_F[NylonZr_N] = { 0.107062f, 0.675449f, 0.099189f, 0.113300f, 0.005f };
                const float NylonZr_den = 1.14f;
                XrayMaterial window( NylonZr_N, NylonZr_Z, NylonZr_F );
                window.density( NylonZr_den );	//	gm/cm3 at STP
                conditions_out.window = window;
                break;
            }
		case AL2O3:	{	    //	Aluminum oxide
                //      Alumina, added for PIXL X-ray tube ceramic body   Al2O3
                Element al( "Al" );
                XrayMaterial window( al, true );    //  oxides = true
                const float Al2O3_den = 3.965f; //  CRC Handbook 51st Ed. 1970  p B-64.
                window.density( Al2O3_den );
                conditions_out.window = window;
                break;
            }
		default: return -100 - (WINDOW_TYPE_INDEX);	break;
	};
	if( cond_in.conditionsVector[WINDOW_THICK_INDEX] < 0 ) return -100 - (WINDOW_THICK_INDEX);
	float winThick = cond_in.conditionsVector[WINDOW_THICK_INDEX] * CM_MICRON;  //  Window thickness is microns
	conditions_out.window.thickness( winThick );

//      geometry factor
	float geometryFactor = cond_in.conditionsVector[GEOMETRY_INDEX];
	if ( geometryFactor < 0 ) return -100 - (GEOMETRY_INDEX);
	if( geometryFactor == 0 ) geometryFactor= 1;
	conditions_out.geometryFactor = geometryFactor;

//      emergent beam geometry
	float emergAngle = cond_in.conditionsVector[EMERG_ANGLE_INDEX];
	if ( emergAngle < 0 || emergAngle > 90.0 ) return -100 - (EMERG_ANGLE_INDEX);
	conditions_out.emergAngle = emergAngle;
	float sinEmerg = sin ( conditions_out.emergAngle * RADDEG );
	if ( sinEmerg < 1.0e-6f ) sinEmerg = 1.0e-6f;
	conditions_out.emergCosecant = 1 / sinEmerg;

//      emergent beam atmosphere path
    conditions_out.emergentPath = conditions_out.incidentPath;
	if( cond_in.conditionsVector[EMERG_PATH_LENGTH_INDEX] < 0 ) return -100 - (EMERG_PATH_LENGTH_INDEX);
	float emerPathLength = cond_in.conditionsVector[EMERG_PATH_LENGTH_INDEX];
	conditions_out.emergentPath.thickness( emerPathLength );

//		set up detector using type from conditions array
	float solidAngleDetector = cond_in.conditionsVector[DET_SOLID_ANGLE_INDEX];
	if ( solidAngleDetector < 0 ) return -100 - (DET_SOLID_ANGLE_INDEX);
	if( solidAngleDetector == 0 ) solidAngleDetector= FOUR_PI;
	conditions_out.solidAngleDetector = solidAngleDetector / FOUR_PI;   //  Note that this is just a factor in conditions, not steradians
	int detType = cond_in.conditionsVector[DETECTOR_TYPE_INDEX];
	if( detType <= NO_DETECTOR || detType >= BAD_DETECTOR ) return -100 - (DETECTOR_TYPE_INDEX);
	float detRes = cond_in.conditionsVector[DET_RESOLUTION_INDEX];
	if( detRes < 0  ) return -100 - (DET_RESOLUTION_INDEX);
	float detWinThick = cond_in.conditionsVector[DET_BE_WINDOW_INDEX];
	if( detWinThick < 0 ) return -100 - (DET_BE_WINDOW_INDEX);
	float detActiveThick = cond_in.conditionsVector[DET_ACTIVE_THICK_INDEX];
	if( detActiveThick < 0 ) return -100 - (DET_ACTIVE_THICK_INDEX);
	DetectorType det_type_enum = NO_DETECTOR;
	switch( detType ) { //  Necessary to convert from int to enum type for constructor call
        case SI_PIN: det_type_enum = SI_PIN;    break;
        case SI_SDD: det_type_enum = SI_SDD;    break;
        case CD_TE: det_type_enum = CD_TE;    break;
        case HP_GE: det_type_enum = HP_GE;    break;
	};
	//  Detector window thickness is in microns for constructor, but in cm in configuration file
	//  Detector active layer thickness is in mm for constructor, but in cm in configuration file
    XrayDetector det_temp( detRes, detWinThick / CM_MICRON, 0, detActiveThick / CM_MM, det_type_enum );
    //  Move shelf adjustment factor and slope into new detector if non-zero
    if( cond_in.conditionsVector[DETECTOR_SHELF_FACTOR_INDEX] > 0 )
        det_temp.set_shelf_factor( cond_in.conditionsVector[DETECTOR_SHELF_FACTOR_INDEX] );
    if( cond_in.conditionsVector[DETECTOR_SHELF_SLOPE_INDEX] > 0 )
        det_temp.set_shelf_slope( cond_in.conditionsVector[DETECTOR_SHELF_SLOPE_INDEX] );
    if( cond_in.conditionsVector[DETECTOR_SHELF_SLOPE_START_INDEX] > 0 )
        det_temp.set_shelf_slope_start( cond_in.conditionsVector[DETECTOR_SHELF_SLOPE_START_INDEX] );
     conditions_out.detector = det_temp;

	float eMin = cond_in.conditionsVector[MINIMUM_ENERGY_INDEX];
	if( eMin < 0 || eMin > cond_in.conditionsVector[KV_INDEX] * 1000 ) return -100 - (MINIMUM_ENERGY_INDEX);
	if( eMin > 0 ) conditions_out.eMin = eMin;

	return 0;
};


