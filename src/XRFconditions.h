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

#ifndef XRFconditions_h
#define XRFconditions_h

#include <vector>
#include <string>
#include "Element.h"
#include "XraySource.h"
#include "XrayOptic.h"
#include "XrayMaterial.h"
#include "XrayDetector.h"


// structure to hold XRF instrument parameters   Jan. 22, 2006 W. T. Elam  APL/UW
//	added X-ray tube current     Nov. 30, 2011
//  move tube current to XraySource class    May 10, 2015
//  change AmpTekDet to XrayDetector class   May 11, 2016
//  extensive re-write Jan. 27, 2017
//      to use XrayMaterial class for most elements in beam
//      and to include information to help decode conditions vector
//          (from XRFparameterIndexEMSA.h)
//      also add extra information from ISO standard version of EMSA format,
//          like azimuth angle and separate solid angles for source and detector
//  Modified May 14, 2017 to add default value for eMin
//  Modified July 31, 2018 to add linear energy calibration correction
//  Modified May 22, 2020 to put conditions vector and optic file name in struct
//                          add file name for X-ray tube spectrum input from external calculation
//  Modified May 14, 2021   Move shelf factor and slope to XrayDetector and control via -T option

enum XrayAtmosphere { VACUUM = 1, HELIUM, MARS, HE_MARS, AIR, EARTH };

enum XrayWindowMaterials { NO_WINDOW = 0, B4C, PLASTIC, CFRP, ZR, AL, NYLON, NYLON_ZR, AL2O3 };

struct XRFconditions {
	XraySource source;
	XrayMaterial filter;
	XrayOptic optic;			// NY edit 9-13-2013
	XrayMaterial dust_on_optic; //  Added Dec. 22, 2016 for error budget
	XrayMaterial incidentPath;
	float solidAngleSource;
	float excitAngle;
	float excitCosecant;
	float geometryFactor;
	XrayMaterial dust_on_specimen;  //  treat like window
	XrayMaterial window;
	float emergAngle;
	float emergCosecant;
	XrayMaterial emergentPath;
	XrayMaterial dust_on_detector;  //  treat just before detector response
	float solidAngleDetector;
	XrayDetector detector;
	float eMin = 900;   //  Default value for compatibility with previous configurations
    std::string tube_file_title;
};

std::string toString(const XRFconditions &cond);


//		Indices for parameters in conditions array
// Added X-ray tube current     Nov. 30, 2011 (also added parentheses)
// Tabbed 'break;' over to align in a column					 NY 9-9-2013
// Added test_optic_param to allow specification of optic types, NY 9-10-2013

#define XRF_PARAMETER_FIRST			-1
#define ANODE_Z_INDEX				(XRF_PARAMETER_FIRST+1)
#define KV_INDEX					(XRF_PARAMETER_FIRST+2)
#define TUBE_INC_ANGLE_INDEX		(XRF_PARAMETER_FIRST+3)
#define TUBE_TAKEOFF_ANGLE_INDEX 	(XRF_PARAMETER_FIRST+4)
#define TUBE_BE_WINDOW_INDEX		(XRF_PARAMETER_FIRST+5)
#define TUBE_CURRENT_INDEX			(XRF_PARAMETER_FIRST+6)
#define FILTER_Z_INDEX				(XRF_PARAMETER_FIRST+7)
#define FILTER_THICK_INDEX			(XRF_PARAMETER_FIRST+8)
#define EXCIT_ANGLE_INDEX			(XRF_PARAMETER_FIRST+9)
#define EMERG_ANGLE_INDEX			(XRF_PARAMETER_FIRST+10)
#define AZIMUTH_ANGLE_INDEX			(XRF_PARAMETER_FIRST+11)
#define XTILT_ANGLE_INDEX			(XRF_PARAMETER_FIRST+12)
#define YTILT_ANGLE_INDEX			(XRF_PARAMETER_FIRST+13)
#define X_POSITION_INDEX			(XRF_PARAMETER_FIRST+14)
#define Y_POSITION_INDEX			(XRF_PARAMETER_FIRST+15)
#define Z_POSITION_INDEX			(XRF_PARAMETER_FIRST+16)
#define SOURCE_SOLID_ANGLE_INDEX	(XRF_PARAMETER_FIRST+17)
#define DET_SOLID_ANGLE_INDEX		(XRF_PARAMETER_FIRST+18)
#define GEOMETRY_INDEX			    (XRF_PARAMETER_FIRST+19)
#define PATH_TYPE_INDEX				(XRF_PARAMETER_FIRST+20)
#define INC_PATH_LENGTH_INDEX		(XRF_PARAMETER_FIRST+21)
#define EMERG_PATH_LENGTH_INDEX		(XRF_PARAMETER_FIRST+22)
#define WINDOW_TYPE_INDEX			(XRF_PARAMETER_FIRST+23)
#define WINDOW_THICK_INDEX			(XRF_PARAMETER_FIRST+24)
#define DETECTOR_TYPE_INDEX			(XRF_PARAMETER_FIRST+25)
#define DET_RESOLUTION_INDEX		(XRF_PARAMETER_FIRST+26)
#define DET_BE_WINDOW_INDEX			(XRF_PARAMETER_FIRST+27)
#define DET_ACTIVE_THICK_INDEX		(XRF_PARAMETER_FIRST+28)
#define TEST_OPTIC_TYPE_INDEX		(XRF_PARAMETER_FIRST+29)		// NY Edit 9-10-2013  Kept for compatibility, zero => no optic
#define MINIMUM_ENERGY_INDEX		(XRF_PARAMETER_FIRST+30)
#define ENERGY_CORRECTION_SLOPE_INDEX		(XRF_PARAMETER_FIRST+31)        //  July 31, 2018
#define ENERGY_CORRECTION_OFFSET_INDEX		(XRF_PARAMETER_FIRST+32)
#define DETECTOR_SHELF_FACTOR_INDEX		(XRF_PARAMETER_FIRST+33)    //  Added May 14, 2021
#define DETECTOR_SHELF_SLOPE_INDEX		(XRF_PARAMETER_FIRST+34)
#define DETECTOR_SHELF_SLOPE_START_INDEX		(XRF_PARAMETER_FIRST+35)
#define XRF_PARAMETER_LAST			(XRF_PARAMETER_FIRST+36)

//  These are used in fpSetupConditions and read_EMSA_PIXL to find keywords for error messages when optic file or tube file are not readaable
#define XRF_PARAMETER_OPTIC_FILE			(XRF_PARAMETER_LAST+10)
#define XRF_PARAMETER_TUBE_FILE			(XRF_PARAMETER_LAST+11)


struct XRFconditionsInput {
    XRFconditionsInput() {
        conditionsVector.resize( XRF_PARAMETER_LAST, 0 );
    }

    std::vector <float> conditionsVector;
    std::string optic_file_name;
    std::string tube_file_name;
    std::string anode_element_list;
};

#endif
