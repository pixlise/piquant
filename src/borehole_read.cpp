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

#include <exception>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdlib.h>
#include "borehole_read.h"
#include "XRFconditions.h"
#include "upper_trim.h"

//  Added parameter input for optic type   W. T. Elam APL/UW Sep. 30, 2013
//  (Already included in XRFparameterIndex.h and XRFconditions.h by Nick Yang Sep. 2013)
// Modified Dec. 8, 2016
//  Remove training CR if found (file with Windows line endings read on Linux or Mac)
// Modified Jan. 18, 2017
//  Use vector for conditionsArray and remove "using namespace std" from header
// Modified Feb. 1, 2017
//  Change SOLID_ANGLE_INDEX to GEOMETRY_INDEX for better compatibility with new EMSA format
//  Put combined solid angle from old format files into Geometry factor, since source and detector solid angles are separate in new format
// Modified Sept. 14, 2017
//  Read and return location information and title strings

	using namespace std;

void parseEMASkeyword( const string strIn, const string delim, string &sKeyword, string &sValue );

int borehole_read ( const string spectrumFileName, vector <float> &conditionsArray,
				   vector <float> &spectrum, float &ev_start, float &ev_ch, float &live_time,
				   vector <string> &title_strings, float &x, float &y, float &z ) {

    // Check that conditions array is large enough
    if(conditionsArray.size() < XRF_PARAMETER_LAST) {
        cout << "borehole_read: conditionsArray size " << conditionsArray.size() << " less than expected: " << XRF_PARAMETER_LAST << endl;
        return -1;
    }

//		set up keywords for conditions
vector <string> paramName( XRF_PARAMETER_LAST + 1 );
paramName[ANODE_Z_INDEX] = "anode_z";
paramName[KV_INDEX] = "beamkv";
paramName[TUBE_INC_ANGLE_INDEX] = "tube_inc_angle";
paramName[TUBE_TAKEOFF_ANGLE_INDEX] = "tube_takeoff_angle";
paramName[TUBE_BE_WINDOW_INDEX] = "tube_be_window";
paramName[TUBE_CURRENT_INDEX] = "tube_current";
paramName[FILTER_Z_INDEX] = "filter_z";
paramName[FILTER_THICK_INDEX] = "filter_thick";
paramName[EXCIT_ANGLE_INDEX] = "excit_angle";
paramName[EMERG_ANGLE_INDEX] = "emerg_angle";
paramName[SOURCE_SOLID_ANGLE_INDEX] = "solid_angle";
paramName[PATH_TYPE_INDEX] = "path_type";
paramName[INC_PATH_LENGTH_INDEX] = "inc_path_length";
paramName[EMERG_PATH_LENGTH_INDEX] = "emerg_path_length";
paramName[WINDOW_TYPE_INDEX] = "window_type";
paramName[WINDOW_THICK_INDEX] = "window_thick";
paramName[DETECTOR_TYPE_INDEX] = "detector_type";
paramName[TEST_OPTIC_TYPE_INDEX] = "optic_type";
paramName[MINIMUM_ENERGY_INDEX] = "minimum_energy";


//		open input file
	ifstream inputFile(spectrumFileName.c_str(), ios::in);
	if ( !inputFile ) {
		return -1;
	};

//		read file, interpret data, and place in structure

    string strRead, sKeyword, sValue;
	int numChannels;
    bool bNumChanFound;
    bNumChanFound = false;
	while (sKeyword != "#SPECTRUM") {

        //Process EMSA keywords
		getline( inputFile, strRead );
		if ( !inputFile ) {
			return -2;
		};
		//  Get rid of trailing CR if file is Windows line endings on Linux or Mac
		if( (int) (strRead.data()[strRead.length()-1]) == 13 ) strRead.erase( strRead.length()-1,1);
		parseEMASkeyword(strRead, ":", sKeyword, sValue);
        sKeyword = upper_trim(sKeyword);
//			work-around for early typographical error in spectrum save
		if( sKeyword == "SOLID_ANLGE" ) sKeyword = "SOLID_ANGLE";
		if (sKeyword == "#FORMAT") {
            if (sValue != "EMSA/MAS Spectral Data File") {
                return -2;
            };
        } else if (sKeyword == "#VERSION") {
            if (sValue != "1.0") {
                return -3;
            };
        } else if (sKeyword == "#TITLE") {
            title_strings.push_back( sValue );
        } else if (sKeyword == "#DATE") {
 //           txtAcqDate.TEXT = sValue
        } else if (sKeyword == "#TIME") {
 //           txtAcqTime.TEXT = sValue
        } else if (sKeyword == "#OWNER") {
 //           txtOperator.TEXT = sValue
        } else if (sKeyword == "#NPOINTS") {
			istringstream val( sValue );
            val >> numChannels;
            bNumChanFound = true;
        } else if (sKeyword == "#XUNITS") {
            if (sValue != "eV") {
                return -4;
            };
        } else if (sKeyword == "#YUNITS") {
            if (sValue != "COUNTS") {
                return -5;
            };
        } else if (sKeyword == "#XPERCHAN") {
			istringstream val( sValue );
            val >> ev_ch;
        } else if (sKeyword == "#OFFSET") {
 			istringstream val( sValue );
			val >> ev_start;
        } else if (sKeyword == "#LIVETIME") {
 			istringstream val( sValue );
			val >> live_time;
        } else if (sKeyword == "#XPOSITION") {
 			istringstream val( sValue );
			val >> x;
        } else if (sKeyword == "#YPOSITION") {
 			istringstream val( sValue );
			val >> y;
        } else if (sKeyword == "#ZPOSITION") {
 			istringstream val( sValue );
			val >> z;
       } else if (sKeyword == "#SIGNALTYPE") {
            if (sValue != "XRF") {
                return -6;
            };
		};

//			Process Conditions keywords (not part of EMAS format)
		int i;
		for( i=0; i<XRF_PARAMETER_LAST + 1; i++ ) {
			string upper_test = upper_trim( paramName[i] );
            if (sKeyword == upper_test) {
 				istringstream val( sValue );
 				//  Conversions necessary for backward compatibility with previous configuration files
                val >> conditionsArray[i];
 				if( i == WINDOW_THICK_INDEX && conditionsArray[i] <= 1 ) conditionsArray[i] /= CM_MICRON;
 				if( i == SOURCE_SOLID_ANGLE_INDEX ) conditionsArray[i] *= FOUR_PI;
 				if( i == TUBE_CURRENT_INDEX && conditionsArray[i] >= 1 ) conditionsArray[i] /= 1000;    //  convert uA to mA
				break;
            };
		};

	};


//    int ic;
//    for( ic = 0; ic < XRF_PARAMETER_LAST; ic++ ) {
//        cout << ic << "  " << paramName[ic] << "  " << conditionsArray[ic] << endl;
//    }

    //Process spectrum

	if( ! bNumChanFound ) return -7;
	spectrum.resize( numChannels );
	int i;
	for( i=0; i<numChannels; i++ ) {
        inputFile >> spectrum[i];
		if ( !inputFile ) {
			return -8;
		};
	};

	return 0;

};


void parseEMASkeyword( const string strIn, const string delim, string &sKeyword, string &sValue) {
// strIn - in, input string
// del - in, delimiter character
// sKeyword - out, keyword up to delimiter (not inlcuding delimiter)
// sValue - out, rest of string after delimiter and following blank, if any
	int j;
	int slen;
    slen = strIn.length();
	j = strIn.find( delim );
	if( j >= 0 && j < strIn.length() ) {
		sKeyword = strIn.substr( 0, j );
		if( strIn.substr( j + 1, 1 ) == " ") j = j + 1;
		sValue = strIn.substr( j + 1, slen);
	} else {
		sKeyword = strIn.substr( 0, slen );
		sValue.erase();
	};
};


