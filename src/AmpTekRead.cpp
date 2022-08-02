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
#include <string>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include "AmpTekRead.h"
#include "Fit.h"

//	Modified Oct. 15, 2013
//		to avoid substr error if DESCRIPTION line is missing the last three characters " - "
//	Modified Oct. 24, 2013
//		refine TAG check to only look at 9 characters (in case TAG is live_data_1, etc.)
//  Modified Aug. 10, 2016
//      Check status of input file and string streams while reading
//          to avoid infinite loops and catch invalid characters
// Modified Dec. 8, 2016
//  Remove training CR if found (file with Windows line endings read on Linux or Mac)
// Modified Feb. 20, 2017
//  Return zero energy start and unity energy per channel if no calibration (originally returned nan)

	using namespace std;

int amptek_read ( string inputFileName, AmpTekSpec &sp ) {

//		open input file
	ifstream inputFile(inputFileName.c_str(), ios::in);
	if ( !inputFile ) {
		return -1;
	};

//		read file, interpret data, and place in structure

//		file identification string
	string inputLine;
	getline( inputFile, inputLine );
    //  Get rid of trailing CR if file is Windows line endings on Linux or Mac
    if( (int) (inputLine.data()[inputLine.length()-1]) == 13 ) inputLine.erase( inputLine.length()-1,1);
	if ( inputLine != "<<PMCA SPECTRUM>>" ) {
		cout << "Can't interpret file, first line is: ";
		cout << inputLine << endl;
		return -2;
	};
	sp.fileID = inputLine;

//		read spectrum parameters

	bool parameters = true;
	while ( parameters ) {
//			read a full line at a time
		getline( inputFile, inputLine );
        if ( !inputFile ) {
            return -3;
        };
        //  Get rid of trailing CR if file is Windows line endings on Linux or Mac
        if( (int) (inputLine.data()[inputLine.length()-1]) == 13 ) inputLine.erase( inputLine.length()-1,1);
//			check for start of calibration info or data
		if ( inputLine == "<<CALIBRATION>>" || inputLine == "<<DATA>>" || inputLine == "<<END>>" ) {
			parameters = false;
			break;
		};
//			set up string stream to interpret flags and number data
		istringstream inputSS(inputLine);
		string token;
		string minusChar;
		inputSS >> token;
		bool tokenOK = false;
		if ( token == "TAG" ) {
			tokenOK = true;
//				swallow dash
			inputSS >> minusChar;
			inputSS >> sp.dataType;
			if ( sp.dataType.substr( 0, 9 ) != "live_data" ) cout << "*** Warning - unkown data type: " << sp.dataType << endl;
		};
		if ( token == "DESCRIPTION" ) {
			tokenOK = true;
//				put entire rest of line (if any) into description
			if( inputLine.length() > 14 ) {		//	Oct. 15, 2013 sometimes inputLIne is only 11 chars, " - " is missing
				sp.description = inputLine.substr( 14 );
			}
		};
		if ( token == "GAIN" ) {
			tokenOK = true;
//				swallow dash
			inputSS >> minusChar;
			inputSS >> sp.gain;
		};
		if ( token == "THRESHOLD" ) {
			tokenOK = true;
//				swallow dash
			inputSS >> minusChar;
			inputSS >> sp.threshold;
		};
		if ( token == "LIVE_MODE" ) {
			tokenOK = true;
//				swallow dash
			inputSS >> minusChar;
			inputSS >> sp.live_mode;
		};
		if ( token == "PRESET_TIME" ) {
			tokenOK = true;
//				swallow dash
			inputSS >> minusChar;
			inputSS >> sp.preset_time;
		};
		if ( token == "LIVE_TIME" ) {
			tokenOK = true;
//				swallow dash
			inputSS >> minusChar;
			inputSS >> sp.live_time;
		};
		if ( token == "REAL_TIME" ) {
			tokenOK = true;
//				swallow dash
			inputSS >> minusChar;
			inputSS >> sp.real_time;
		};
		if ( token == "START_TIME" ) {
			tokenOK = true;
//				swallow dash
			inputSS >> minusChar;
			inputSS >> sp.start_time;
		};
		if ( token == "SERIAL_NUMBER" ) {
			tokenOK = true;
//				swallow dash
			inputSS >> minusChar;
			inputSS >> sp.serial_number;
		};
		if ( ! tokenOK ) {
			cout << "*** Warning - Unrecognized token: " << token <<endl;
		};
        if ( !inputSS ) {
            return -3;
        };
	};

//		read calibration information
	if ( inputLine == "<<CALIBRATION>>" || inputLine == "<<END>>" ) {
		bool calibration = true;
		while ( calibration ) {
//				read a full line at a time
			getline( inputFile, inputLine );
            if ( !inputFile ) {
                return -3;
            };
            //  Get rid of trailing CR if file is Windows line endings on Linux or Mac
            if( (int) (inputLine.data()[inputLine.length()-1]) == 13 ) inputLine.erase( inputLine.length()-1,1);
			if ( inputLine == "<<DATA>>" ) {
				calibration = false;
				break;
			};
//				set up string stream to interpret flags and number data
			istringstream inputSS(inputLine);
			string token;
			string minusChar;
			inputSS >> token;
			bool tokenOK = false;
			if ( token == "LABEL" ) {
				tokenOK = true;
//					swallow dash
				inputSS >> minusChar;
				inputSS >> sp.cal_label;
                if ( !inputSS ) {
                    return -3;
                };
				if ( sp.cal_label != "Channel" ) cout << "*** Warning - unknown calibration label: " << sp.cal_label << endl;
			};
			if ( ! tokenOK ) {
//					assume its a channel - energy pair
				istringstream inputNumbers(inputLine);
				float ch;
				inputNumbers >> ch;
				float en;
				inputNumbers >> en;
                if ( !inputNumbers ) {
                    return -3;
                };
				sp.cal_channel.push_back( ch );
				sp.cal_energy.push_back( en );
				tokenOK = true;
			};
		};
	};

//		read spectral data as channels of counts

	if ( inputLine == "<<DATA>>" ) {
		bool data = true;
		while ( data ) {
//				read a full line at a time
			getline( inputFile, inputLine );
            if ( !inputFile ) {
                return -3;
            };
            //  Get rid of trailing CR if file is Windows line endings on Linux or Mac
            if( (int) (inputLine.data()[inputLine.length()-1]) == 13 ) inputLine.erase( inputLine.length()-1,1);
			if ( inputLine == "<<END>>" ) {
				data = false;
				break;
			};
			istringstream inputNumbers(inputLine);
			float counts;
			inputNumbers >> counts;
            if ( !inputNumbers ) {
                return -3;
            };
			sp.spectrum.push_back( counts );
		};
	};

//		calculate slope and intercept for calibration (it information is present)
    if( sp.cal_channel.size() > 2 ) {
        //void fit( vector <float> x, vector <float> y, vector <float> sig, float &a,
        //	float &b, float &siga, float &sigb, float &chi2, float &q)
        float a, b, siga, sigb, chi2, q;
        vector <float> sig(0);
        fit ( sp.cal_channel, sp.cal_energy, sig, a, b, siga, sigb, chi2, q );
    //		convert from keV to eV
        sp.ev_ch = 1000 * b;
        sp.ev_start = 1000 * a;
	}

	return 0;

};
