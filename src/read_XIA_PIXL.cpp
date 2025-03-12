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
#include "read_XIA_PIXL.h"
#include "upper_trim.h"

// Dec.26, 2016
//	Begun writing, based on read_EMSA_PIXL function of Dec. 8, 2016
//  Modified April 14, 2017 to properly calculate live time from input and output count rates
//  Modified June 9, 2017
//      Improve error processing, return line number if error
//  Modified Dec. 6, 2017
//      Add input counts and output counts from count rates and times, as for PIXL DSPC histograms
//  Modified Dec. 23, 2017
//      Improve error messages, with specifics for checks before reading spectrum
//  Modified Jan. 3, 2017
//      Move tab, single quote, double quote, blank, comma, and underscore definitions to XRFconstants.h
//  Modified Nov. 6, 2019
//      Add some keyword revision because XIA changed file format (ICR: OCR: LIVE TIME: REAL TIME:)
//      Also look for kcps in ICR and OCR and multiply be 1000



	using namespace std;


//	Utility function at end of this file
int parseXIAkeyword( const string strIn, const string delim, string &sKeyword, string &sValue );

int read_XIA_PIXL ( const std::string spectrumFileName, XraySpectrum &spectrum_data_out,
				std::string &spec_acq_date, std::string &spec_title,
				std::string &spec_sample, std::string &spec_unitID ) {


//		Open input file

	ifstream inputFile(spectrumFileName.c_str(), ios::in);
	if ( !inputFile ) {
		return -999999;
	};
	int line_number = 0;

//		Read line from file, interpret data, and place in arguments

    string strRead, sKeyword;
	int numChannels;
    bool bNumChanFound;
    bNumChanFound = false;

	float icr = 0;
	float ocr = 0;
	float real_time = 0;
	float live_time = 0;

	int keyword_length = 1;

	while ( !( !inputFile ) && keyword_length > 0 ) {

        //Process XIA keywords and check file format
		getline( inputFile, strRead );
		line_number++;
		if ( !inputFile ) break;
		string sValue;
		//  Get rid of trailing CR if file is Windows line endings on Linux or Mac
		if( (int) (strRead.data()[strRead.length()-1]) == 13 ) strRead.erase( strRead.length()-1,1);
		keyword_length = parseXIAkeyword(strRead, "=", sKeyword, sValue);
        sKeyword = upper_trim(sKeyword);
//        cout << "|" << sKeyword << "|    " << sValue.length() << "  |" << sValue << "|" << endl;
		if (sKeyword == "FILE VERSION") {
            if ( sValue.length()<=0 || sValue != "2" ) {
                return -999999;
            };
        } else if (sKeyword == "MCA DATA") {
            if ( sValue.length()<=0 ) {
                return -999999;
            };
            if ( sValue.find( "ProSpect") < 0 || sValue.find( "ProSpect") >= sValue.length() ) {
                return -999999;
            };

//	Handle text keywords with returned strings

        } else if (sKeyword == "TITLE") {
        	 spec_title = sValue;
        } else if (sKeyword == "SAMPLE") {
            spec_sample = sValue;
        } else if (sKeyword == "CURRENT DATE") {
            spec_acq_date = sValue;
        } else if (sKeyword == "USER NAME") {
            spec_unitID = sValue;

//	No energy calibrations in this file format yet

         } else if (sKeyword == "NUMBER MCA BINS") {
            if(sValue.length()>0) {
                istringstream val( sValue );
                val >> numChannels;
                bNumChanFound = true;
            }
       } else if (sKeyword == "INPUT COUNT RATE" || sKeyword == "ICR:") {
            if( sValue.length() > 0 ) {
                istringstream val( sValue );
                float temp;
                val >> temp;
                if( !val ) return -line_number;
                icr = temp;
                if( sValue.find( "kcps" ) ) icr = icr * 1000;
//                cout << "ICR " << sKeyword << "    " << icr << endl;
            } else {
                return -line_number;
            }
       } else if (sKeyword == "OUTPUT COUNT RATE" || sKeyword == "OCR:") {
            if( sValue.length() > 0 ) {
                istringstream val( sValue );
                float temp;
                val >> temp;
                if( !val ) return -line_number;
                ocr = temp;
                if( sValue.find( "kcps" ) ) ocr = ocr * 1000;
//                cout << "OCR " << sKeyword << "    " << ocr << endl;
            } else {
                return -line_number;
            }
       } else if (sKeyword == "REALTIME" || sKeyword == "REAL TIME:") {
            if( sValue.length() > 0 ) {
                istringstream val( sValue );
                float temp;
                val >> temp;
                if( !val ) return -line_number;
                real_time = temp;
//                cout << "RT " << sKeyword << "    " << real_time << endl;
            } else {
                return -line_number;
            }
       } else if (sKeyword == "LIVETIME" || sKeyword == "LIVE TIME:") {
            if( sValue.length() > 0 ) {
                istringstream val( sValue );
                float temp;
                val >> temp;
                if( !val ) return -line_number;
                live_time = temp;
//                cout << "LT " << sKeyword << "    " << live_time << endl;
            } else {
                return -line_number;
            }

//  End of keywords, start reading spectrum
        } else if ( keyword_length <= 0 ) {
            break;
        } else {
            //  Ignore keywords not processed
        }
	};  //  End while

//		Process spectrum

	if( ! bNumChanFound ) {
        cout << "*** Error - number of channels not found. ***" << endl;
        return -line_number;
    }
	if( numChannels <= 0 ) return 0; // Done if no spectrum data
	int nc_check = 0;
	inputFile >> nc_check;
    line_number++;
	if( nc_check != numChannels )  {
        cout << "*** Error - number of channels preceding channel data does not match keyword value. ***" << endl;
        return -line_number;
    }
    vector <float> spec_temp( numChannels );
	int i;
	for( i=0; i<numChannels; i++ ) {
		if ( !inputFile ) return -line_number;
        float temp;
        inputFile >> temp;
        if( !inputFile ) return -line_number;
        spec_temp[i] = temp;
	};
	XraySpectrum spec_data_temp( spec_temp, 0, 0 ); //  No energy calibration
//	cout << "LT calc " << real_time * ocr / icr << "     " << live_time * icr << "     " << real_time * ocr << endl;
	spec_data_temp.live_time( real_time * ocr / icr );  //  ProSpect Software Software Revision: 1.0.24
	spec_data_temp.real_time( real_time );
	spec_data_temp.header_info_change().triggers = live_time * icr;
	spec_data_temp.header_info_change().events = real_time * ocr;
    spectrum_data_out = spec_data_temp;

	return 0;

};


int parseXIAkeyword( const string strIn, const string delim, string &sKeyword, string &sValue ) {
// strIn - in, input string
// del - in, delimiter character
// sKeyword - out, keyword up to delimiter (not inlcuding delimiter)
// sValue - out, rest of string after delimiter and following blank, if any
	const int slen = strIn.length();
	int j;
	j = strIn.find( delim );
	if( j >= 0 && j < slen ) {
		sKeyword = strIn.substr( 0, j );
		if( slen > j + 1 && strIn.substr( j+1, 1 ) == BLANK_CHARACTER ) j++;
		sValue = strIn.substr( j + 1, slen );
	} else {
        sKeyword.clear();
		sValue = strIn.substr( 0, slen );
	};
	return sKeyword.length();
};
