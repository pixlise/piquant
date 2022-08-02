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
#include <string>
#include "read_EMSA_PIXL.h"
#include "XRFconditions.h"
#include "parse_records.h"
#include "upper_trim.h"

//	Reads X-ray Fluorescence spectrum files in EMSA/MAS format (ISO 22029 2012)
//	This version is explicitly for the Planetary Instrument for X-ray Lithochemistry (PIXL)
//	It has many user-defined keywords (following the standard) and is set up for two parallel detectors

// March 9, 2016
//	Begun writing, based on borehole_read function of Dec. 11, 2015
// Modified June 1, 2016
//  Fix a few things to read example files from David Flannery (new PIXL breadboard format)
//  Handle more than one word (separated by blanks) in OWNER field
//  Return if no channels (NPOINTS zero)  [error if missing, it is a required keyword in the ISO standard]
// Modified Dec. 8, 2016
//  Remove training CR if found (file with Windows line endings read on Linux or Mac)
//  Return spectrum objects even if no channels to capture energy calibration info from configuration file
// Modified Jan. 18, 2017
//  Use vector for conditionsArray
// Modified Feb. 1, 2017
//  Change SOLID_ANGLE_INDEX to GEOMETRY_INDEX for better compatibility with new EMSA format
// Modified March 24, 2017
//  Several changes to make fully compatible with EMSA format configuration files
//      Change X-ray tube emission current to microAmps (from milliAmps)
//      Change optic keyword to handle numbers (esp zero for no optic and 3 for PIXL default optic) as well as a file name
//              (as of this date a file name is returned but not handled in the main code)
//  Modified May 14, 2017
//      Add keyword for eMin in conditions, to read from configuration file
//  Modified June 9, 2017
//      Improve error processing, return line number if error
//  Modified June 26, 2017
//      Fix handling of optic file keyword for number or file name, including use of quotes
//  Modified Sept. 14, 2017
//      Read location information and store in XraySpectrum object
//      Move title info to XraySpectrum object
//  Modified Sept. 29, 2017
//      Return units for each parameter read in (and check carefully to be sure this is right)
//  Modified Jan. 3, 2017
//      Move tab, single quote, double quote, blank, comma, and underscore definitions to XRFconstants.h
//      Change parseCommas function to use parse_records and use comma instead of blank as delimiter for values
//  Modified Jan. 17, 2018
//      Change XUNITS and YUNITS in MSA read to be case-insensitive
//      Process eV or keV XUNITS
//  Modified Jan. 31, 2018
//      Process new ##EVENTS and ##TRIGGERS keywords and make XIA livetime calculation
//      Add checks for less than or equal to zero for NPOINTS (0 OK for config files), XPERCHAN, LIVETIME, TRIGGERS, and EVENTS
//  Modified Mar. 7, 2018
//      Fix termination of quoted strings to search for double quotes instead of backslash (entered by mistake)
//  Modified July 25, 2018
//      Read real time and put in XraySpectrum objects
//      Add write_EMSA_PIXL to write spectrum as an msa file
//  Modified Sep. 19, 2018
//      Allow zero values for LIVETIME, REALTIME, TRIGGERS, and EVENTS (for testing MSA files without collecting a spectrum)
//  Modified May 6, 2019
//      Write measured spectrum for bulk sum (not calculated spectrum, pass flag for bulk sum)
//  Modified May 13, 2019
//      In XraySpectrum, all non-spectrum information put in separate structure
//      Date, time, etc. now in that structure, not passed as argument
//      Get rid of spectrum_hold struct, enter in XraySpectrum object when read in
//      Fix XIA live time calculation to include overflows and underflows
//      Parse more keywords for aux info and header info from MSA files
//  Modified Feb. 5, 2020
//      X position value put into Y location in spectrum, fixed
//  Modified May 22, 2020 to put conditions vector and optic file name in struct
//                          add file name for X-ray tube spectrum input from external calculation
//                          add optic_file_name and tube_file_name to get_EMSA_keyword for error messages
//  Modified Dec. 4, 2020   fix minor bug in removing CR from end of lines
//  Modified June 9, 2021   Write geometry factor if it is non-zero

	using namespace std;

//	Utility functions at end of this file
void parseEMSAkeyword( const string strIn, const string delim, string &sKeyword, vector <string> &sValue );
int parse_EMSA_description( const int index, const string &s );

int read_EMSA_PIXL ( const std::string spectrumFileName, XRFconditionsInput &conditionsStructEMSA,
				std::vector <XraySpectrum> &spectrum_vector ) {

//  Initialize output conditions structure
    conditionsStructEMSA.conditionsVector.resize( XRF_PARAMETER_LAST, 0 );

//		Open input file

	ifstream inputFile(spectrumFileName.c_str(), ios::in);
	if ( !inputFile ) {
		return -1;
	};

//		Read line from file, interpret data, and place in arguments

    string strRead, sKeyword;
	int numChannels = 0;
    bool bNumChanFound;
    bNumChanFound = false;
	int numColumns = 0;
	int line_number = 0;
    bool kev_units = false;
	spectrum_vector.clear();
	//  Need to hold these so they can be entered together into XraySpectrum object
	vector <float> ev_ch;
	vector <float> ev_start;
	vector <vector <float>> spectrum_hold_vec;
	Spec_Aux_Info spec_info_hold;
	bool livetime_XIA = false;
	int triggers_line = -1;

	while ( !( !inputFile ) ) {

        //Process required EMSA keywords and check file format
		getline( inputFile, strRead );
		line_number++;
		if ( !inputFile ) {
            break;
		}
		vector <string> sValue;
		//  Get rid of trailing CR if file is Windows line endings on Linux or Mac
		if( strRead.length() > 0 && (int) (strRead.data()[strRead.length()-1]) == 13 ) strRead.erase( strRead.length()-1,1);
		parseEMSAkeyword(strRead, ":", sKeyword, sValue);
        sKeyword = upper_trim(sKeyword);
//			work-around for early typographical error in spectrum save
		if( sKeyword == "SOLID_ANLGE" ) sKeyword = "SOLID_ANGLE";
//		cout << strRead << "   keyword '" << sKeyword << "'  " << (sValue.size()>0?sValue[0]:" ") << endl;
		if (sKeyword == "#FORMAT") {
            if ( sValue.size()<=3 || sValue[0] != "EMSA/MAS" || sValue[1] != "spectral" || sValue[2] != "data" || sValue[3] != "file" ) {
                return -line_number;
            };
        } else if (sKeyword == "#VERSION") {
            if ( sValue.size()<=1 || sValue[0] != "TC202v2.0" || sValue[1] != "PIXL" ) {
                return -line_number;
            };
//	Need to check for PIXL, return -6                                *****
       } else if (sKeyword == "#SIGNALTYPE") {
            if ( sValue.size()<=0 || upper_trim(sValue[0]) != "XRF" ) {
                return -line_number;
            };
       } else if (sKeyword == "#DATATYPE") {
            if ( sValue.size()<=0 || ( upper_trim(sValue[0]) != "Y" && upper_trim(sValue[0]) != "YY" ) ) {
                return -line_number;
            };
       } else if (sKeyword == "#COMMENT") {
        	 if( sValue.size()>0) {
                int k;
                string temp_str;
                for(k=0; k<sValue.size(); k++ ) temp_str += sValue[k] + BLANK_CHARACTER;
                spec_info_hold.comments.push_back( temp_str );
            }

//	Handle text keywords with returned strings

        } else if (sKeyword == "#TITLE") {
        	 if( sValue.size()>0) {
                int k;
                string temp_str;
                for(k=0; k<sValue.size(); k++ ) temp_str += sValue[k] + BLANK_CHARACTER;
                spec_info_hold.titles.push_back( temp_str );
            }
        } else if (sKeyword == "#DATE") {
            if(sValue.size()>0) spec_info_hold.date = sValue[0];
        } else if (sKeyword == "#TIME") {
            if(sValue.size()>0) spec_info_hold.time = sValue[0];
        } else if (sKeyword == "#OWNER") {
            spec_info_hold.owner.clear();
            if( sValue.size()>0) {	// Only keep first title
                int k;
                for(k=0; k<sValue.size(); k++ ) spec_info_hold.owner += sValue[k] + BLANK_CHARACTER;
            }
        } else if (sKeyword == "#NPOINTS") {
            if(sValue.size()>0) {
                istringstream val( sValue[0] );
                val >> numChannels;
                if( !val || numChannels < 0 ) return -line_number;
                bNumChanFound = true;
            }
        } else if (sKeyword == "#NCOLUMNS") {
            if(sValue.size()>0) {
                istringstream val( sValue[0] );
                val >> numColumns;
                if( numColumns<1 || numColumns>2 ) {
                    return -line_number;
                } else {
                    spectrum_vector.resize(numColumns);
                    ev_ch.resize(numColumns);
                    ev_start.resize(numColumns);
                    spectrum_hold_vec.resize(numColumns);
                }
            }
        } else if (sKeyword == "#XUNITS") {
            if ( sValue.size()<=0 ) return -line_number;
            if( upper_trim(sValue[0]) == "EV" ) {
                kev_units = false;
            } else if( upper_trim(sValue[0]) == "KEV" ) {
                kev_units = true;
            } else {
                return -line_number;
            };
            // Don't overwrite label from XLABEL keyword, only use if no label set so far
            //if( spec_x_label.length() == 0 ) spec_x_label = sValue[0];
        } else if (sKeyword == "#YUNITS") {
            if (  sValue.size()<=0 || upper_trim(sValue[0]) != "COUNTS") {
                return -line_number;
            };
            // Don't overwrite label from YLABEL keyword, only use if no label set so far
            //if( spec_y_label.length() == 0 ) spec_y_label = sValue[0];
        } else if (sKeyword == "#XLABEL") {
            //if(sValue.size()>0) spec_x_label = sValue[0];
            continue;
        } else if (sKeyword == "#YLABEL") {
            //if(sValue.size()>0) spec_y_label = sValue[0];
            continue;
        } else if (sKeyword == "##OPTICFILE") {
            if(sValue.size()>0) {
                istringstream val( sValue[0] );
                int temp;
                val>> temp;
                if( !val ) {
                    //  If optic file is not a number, set optic type in conditions array
                    conditionsStructEMSA.optic_file_name = sValue[0];
                    conditionsStructEMSA.conditionsVector[TEST_OPTIC_TYPE_INDEX] = 4;
                } else {
                    //  If optic file is a number, set optic type in conditions array
                    conditionsStructEMSA.conditionsVector[TEST_OPTIC_TYPE_INDEX] = temp;
                }
            }
        } else if (sKeyword == "##ANODE") {
            if(sValue.size()>0) {
                istringstream val( sValue[0] );
                int temp;
                val>> temp;
                if( !val ) {
                    //  If anode file is not a number, put the entire string in conditions structure for later parsing
                    string anode_str = sValue[0];
                    if( sValue.size() > 1 ) {
                        unsigned int jj;
                        for( jj=1; jj<sValue.size(); jj++ ) anode_str += "," + sValue[jj];
                    }
                    conditionsStructEMSA.anode_element_list = anode_str;
                    conditionsStructEMSA.conditionsVector[ANODE_Z_INDEX] = 0;
                } else {
                    //  If anode file is a number, set anode Z in conditions array
                    conditionsStructEMSA.conditionsVector[ANODE_Z_INDEX] = temp;
                }
            }
        } else if (sKeyword == "##TUBEFILE") {
            if(sValue.size()>0) {
                    conditionsStructEMSA.tube_file_name = sValue[0];
            }

//	Handle values that are returned as arguments (not in conditions array)

//	Need to handle multiple live times and energy calibrations                                *****
        } else if (sKeyword == "#XPERCHAN") {
            if(sValue.size()>=numColumns) {
                int k;
                for( k=0; k<numColumns; k++) {
                    istringstream val( sValue[k] );
                    float temp;
                    val >> temp;
                    if( !val || temp <= 0 ) return -line_number;
                    ev_ch[k] = temp;
                }
            } else {
                return -line_number;
            }
        } else if (sKeyword == "#OFFSET") {
            if(sValue.size()>=numColumns) {
                int k;
                for( k=0; k<numColumns; k++) {
                    istringstream val( sValue[k] );
                    float temp;
                    val >> temp;
                    if( !val ) return -line_number;
                    ev_start[k] = temp;
                }
            } else {
                return -line_number;
            }
        } else if (sKeyword == "#LIVETIME") {
            if(sValue.size()>=numColumns) {
                int k;
                for( k=0; k<numColumns; k++) {
                    istringstream val( sValue[k] );
                    float temp;
                    val >> temp;
                    if( !val || temp < 0 ) return -line_number;
                    spectrum_vector[k].live_time( temp );
                }
            } else {
                return -line_number;
            }
        } else if (sKeyword == "#REALTIME") {
            if(sValue.size()>=numColumns) {
                int k;
                for( k=0; k<numColumns; k++) {
                    istringstream val( sValue[k] );
                    float temp;
                    val >> temp;
                    if( !val || temp < 0 ) return -line_number;
                    spectrum_vector[k].real_time( temp );
                }
            } else {
                return -line_number;
            }
        } else if (sKeyword == "##TRIGGERS") {
            triggers_line = line_number;
            livetime_XIA = true;
            if(sValue.size()>=numColumns) {
                int k;
                for( k=0; k<numColumns; k++) {
                    istringstream val( sValue[k] );
                    float temp;
                    val >> temp;
                    if( !val || temp < 0 ) return -line_number;
                    spectrum_vector[k].header_info_change().triggers = temp;
                }
            } else {
                return -line_number;
            }
        } else if (sKeyword == "##EVENTS") {
            if(sValue.size()>=numColumns) {
                int k;
                for( k=0; k<numColumns; k++) {
                    istringstream val( sValue[k] );
                    float temp;
                    val >> temp;
                    if( !val || temp < 0 ) return -line_number;
                    spectrum_vector[k].header_info_change().events = temp;
                }
            } else {
                return -line_number;
            }
        } else if (sKeyword == "##OVERFLOWS") {
            if(sValue.size()>=numColumns) {
                int k;
                for( k=0; k<numColumns; k++) {
                    istringstream val( sValue[k] );
                    float temp;
                    val >> temp;
                    if( !val || temp < 0 ) return -line_number;
                    spectrum_vector[k].header_info_change().overflows = temp;
                }
            } else {
                return -line_number;
            }
        } else if (sKeyword == "##UNDERFLOWS") {
            if(sValue.size()>=numColumns) {
                int k;
                for( k=0; k<numColumns; k++) {
                    istringstream val( sValue[k] );
                    float temp;
                    val >> temp;
                    if( !val || temp < 0 ) return -line_number;
                    spectrum_vector[k].header_info_change().underflows = temp;
                }
            } else {
                return -line_number;
            }
        } else if (sKeyword == "##BASE_EVENTS") {
            if(sValue.size()>=numColumns) {
                int k;
                for( k=0; k<numColumns; k++) {
                    istringstream val( sValue[k] );
                    float temp;
                    val >> temp;
                    if( !val || temp < 0 ) return -line_number;
                    spectrum_vector[k].header_info_change().baseline_samples = temp;
                }
            } else {
                return -line_number;
            }
        } else if (sKeyword == "##RESETS") {
            if(sValue.size()>=numColumns) {
                int k;
                for( k=0; k<numColumns; k++) {
                    istringstream val( sValue[k] );
                    float temp;
                    val >> temp;
                    if( !val || temp < 0 ) return -line_number;
                    spectrum_vector[k].header_info_change().preamp_resets = temp;
                }
            } else {
                return -line_number;
            }
        } else if (sKeyword == "##OVER_ADCMAX") {
            if(sValue.size()>=numColumns) {
                int k;
                for( k=0; k<numColumns; k++) {
                    istringstream val( sValue[k] );
                    float temp;
                    val >> temp;
                    if( !val || temp < 0 ) return -line_number;
                    spectrum_vector[k].header_info_change().saturates = temp;
                }
            } else {
                return -line_number;
            }
        } else if (sKeyword == "#XPOSITION") {
            if(sValue.size()>0) {
                istringstream val( sValue[0] );
                float temp;
                val >> temp;
                spec_info_hold.x = temp;
            }
        } else if (sKeyword == "#YPOSITION") {
            if(sValue.size()>0) {
                istringstream val( sValue[0] );
                float temp;
                val >> temp;
                spec_info_hold.y = temp;
            }
        } else if (sKeyword == "#ZPOSITION") {
            if(sValue.size()>0) {
                istringstream val( sValue[0] );
                float temp;
                val >> temp;
                spec_info_hold.z = temp;
            }
        } else if (sKeyword == "##IPOSITION") {
            if(sValue.size()>0) {
                istringstream val( sValue[0] );
                float temp;
                val >> temp;
                spec_info_hold.i = temp;
            }
        } else if (sKeyword == "##JPOSITION") {
            if(sValue.size()>0) {
                istringstream val( sValue[0] );
                float temp;
                val >> temp;
                spec_info_hold.j = temp;
            }
        } else if (sKeyword == "##RTT") {
            if(sValue.size()>0) {
                istringstream val( sValue[0] );
                float temp;
                val >> temp;
                spec_info_hold.rtt = temp;
            }
        } else if (sKeyword == "##PMC") {
            if(sValue.size()>0) {
                istringstream val( sValue[0] );
                float temp;
                val >> temp;
                spec_info_hold.pmc = temp;
            }
        } else if (sKeyword == "##DETECTOR_ID") {
            if(sValue.size()>0) {
                spec_info_hold.det_ID = sValue[0];
            }
        } else if (sKeyword == "#SPECTRUM") {
            break;
        } else {

//		Process values returned in conditions array

            int i;
            for( i=0; i<XRF_PARAMETER_LAST; i++ ) {
                if( i == TEST_OPTIC_TYPE_INDEX ) continue;  //  handled above since it might be a file name
                if( i == ANODE_Z_INDEX ) continue;  //  handled above since it might be an element list
                string upper_test = upper_trim( get_EMSA_keyword(i) );
                if (sKeyword == upper_test) {
                    if(sValue.size()>0) {
                        istringstream val( sValue[0] );
                        val >> conditionsStructEMSA.conditionsVector[i];  //  This will automatically handle if it is a number
                        if( ! val ) {
                            conditionsStructEMSA.conditionsVector[i] = parse_EMSA_description( i, sValue[0] );
                            //  These are always non-negative integers since they encode items from a list
                            if( conditionsStructEMSA.conditionsVector[i] < 0 ) return -line_number;
                        }
                        if( i == TUBE_CURRENT_INDEX ) conditionsStructEMSA.conditionsVector[i] /= 1000;   //  convert from microAmps to milliAmps
                    }
                    break;
                };
            };
        };

	};  //  End while !( !inputFile )

//    int ic;
//    for( ic = 0; ic < XRF_PARAMETER_LAST; ic++ ) {
//        cout << ic << "  " << get_EMSA_keyword(ic) << "  " << conditionsStructEMSA.conditionsVector[ic] << endl;
//    }

//		Process spectrum

	if( ! bNumChanFound && numColumns > 0 ) return -999999;
	if( numColumns <= 0 ) return 0; // Done if no columns of data
	spectrum_vector.resize(numColumns);
	int k;
	for( k=0; k<numColumns; k++ ) {
        if( kev_units ) {
            //  Convert from keV to eV for energy calibration to work properly
            spectrum_vector[k].calibration( ev_start[k] * 1000,
               ev_ch[k] * 1000 );
        } else {
            spectrum_vector[k].calibration( ev_start[k],
               ev_ch[k] );
        }
        //  Put all of the auxiliary info into the XraySpectrum objects
        spectrum_vector[k].aux_info_replace( spec_info_hold );
        if( livetime_XIA ) {
            //  See writeup in file "JPL-XIA_PIXL_FPGA_Specification_v2.06.pdf", page 9
            spectrum_vector[k].header_info_change().live_time_DSPC = spectrum_vector[k].live_time();
            const Spec_Header_Info &header = spectrum_vector[k].header_info();
            if( header.triggers > 0 ) {
                float total_counts_in = header.events + header.overflows + header.underflows;
                spectrum_vector[k].live_time( header.live_time_DSPC * total_counts_in / header.triggers );
            } else if( header.live_time_DSPC != 0 ) return -triggers_line;
        }
        spectrum_hold_vec[k].resize( numChannels );
	}
	if( numChannels <= 0 ) return 0; // Done if no spectrum data
	int i;
	for( i=0; i<numChannels; i++ ) {
		string data_line;
		vector <string> sValue;
		getline( inputFile, data_line );
		line_number++;
		//cout << "data_line " << data_line << endl;
		if ( !inputFile ) {
            return -line_number;
		}
//        parseCommas( data_line, sValue );
        string delim( COMMA_CHARACTER );
        delim += BLANK_CHARACTER;
        parse_records( delim, data_line, sValue );
		if ( sValue.size() <= 0 || sValue[0] == "#ENDOFDATA" ) {
            return -line_number;
		}
        if( sValue.size() < numColumns ) {
            return -line_number;
        }
        int k_prime = -1;
        for( k=0; k<numColumns; k++ ) {
            k_prime++;
            //  Skip empty records (caused by commas, tabs, and blanks together in data)
            if( sValue[k_prime].length() <= 0 && sValue.size() > k_prime + 1 ) k_prime++;
            istringstream val( sValue[k_prime] );
            float temp;
            val >> temp;
            if( !val ) return -line_number;
            spectrum_hold_vec[k][i] = temp;
		}
	};

	for( k=0; k<numColumns; k++ ) {
        spectrum_vector[k].meas( spectrum_hold_vec[k] );
	}

	return 0;
}

//  Only writes spectrum information, not configuration information
int write_EMSA_PIXL ( const XraySpectrum &spectrum,
                const std::string spectrumFileName,
                const bool meas_only ) {
    //  Write spectrum information to MSA file for additional processing
    //  Mostly intended for calculated spectra but works OK for measured spectra
    //  Will include calculated spectrum if present, otherwise includes measured spectrum
    //  Doesn't include fit or component information
    //  See if there is a measured spectrum to include in plot
    bool measured_present = true;
    int is;
    float sum = 0;
    for( is=0; is<spectrum.numberOfChannels(); is++ )
            sum += spectrum.meas()[is];
    if( sum <= 0 ) measured_present = false;
    bool calc_present = spectrum.calc().size() >= spectrum.numberOfChannels();
    if( ( !calc_present && !measured_present ) || ( meas_only && !measured_present ) ) {
        cout << "Error - spectrum does not have measured or calculated data to write. " << endl;
        return -2;
    }
    //  Open file and prepare for writing
    ofstream spectrumFile ( spectrumFileName.c_str() );
    if( !spectrumFile ) {
        cout << "Error opening plot file " << spectrumFileName << endl;
        return -1;
    }
    spectrumFile << "#FORMAT      : EMSA/MAS spectral data file" << endl;
    spectrumFile << "#VERSION     : TC202v2.0 PIXL" << endl;
    int it;
    for( it=0; it<spectrum.aux_info().titles.size(); it++ ) {
        spectrumFile << "#TITLE       : ";
        spectrumFile << spectrum.aux_info().titles[it];
        spectrumFile << endl;
    }
    spectrumFile << "#DATE        : ";
    spectrumFile << spectrum.aux_info().date;
    spectrumFile << endl;
    spectrumFile << "#TIME        : ";
    spectrumFile << spectrum.aux_info().time;
    spectrumFile << endl;
    spectrumFile << "#NPOINTS     : ";
    spectrumFile << spectrum.numberOfChannels();
    spectrumFile << endl;
    spectrumFile << "#NCOLUMNS    : 1" << endl;
    spectrumFile << "#XUNITS      :  eV" << endl;
    spectrumFile << "#YUNITS      :  COUNTS" << endl;
    spectrumFile << "#DATATYPE    :  Y" << endl;
    spectrumFile << "#XPERCHAN    : ";
    spectrumFile << spectrum.calibration().energyPerChannel();
    spectrumFile << endl;
    spectrumFile << "#OFFSET      : ";
    spectrumFile << spectrum.calibration().energyStart();
    spectrumFile << endl;
    spectrumFile << "#SIGNALTYPE  :  XRF" << endl;
//    spectrumFile << "#XPOSITION   : ";
//    spectrumFile << spectrum.x();
//    spectrumFile << endl;
//    spectrumFile << "#YPOSITION   : ";
//    spectrumFile << spectrum.y();
//    spectrumFile << endl;
//    spectrumFile << "#ZPOSITION   : ";
//    spectrumFile << spectrum.z();
//    spectrumFile << endl;
    spectrumFile << "#LIVETIME    : ";
    spectrumFile << spectrum.live_time();
    spectrumFile << endl;
    spectrumFile << "#REALTIME    : ";
    spectrumFile << spectrum.real_time();
    spectrumFile << endl;
    if( spectrum.geometry() != 0 ) {
        spectrumFile << "##GEOMETRY   : ";
        spectrumFile << spectrum.geometry();
        spectrumFile << endl;
    }
    spectrumFile << "#SPECTRUM    :" << endl;
    for( is=0; is<spectrum.numberOfChannels(); is++ ) {
        if( calc_present && (!meas_only) ) spectrumFile << spectrum.calc()[is] << endl;
        else if( measured_present ) spectrumFile << spectrum.meas()[is] << endl;
    }
    spectrumFile << "#ENDOFDATA    : " << endl;

    return 0;
}


const string get_EMSA_keyword( const int index ) {
// TODO: get rid of statics!!!
    //		set up keywords for conditions array using EMSA keywords and user-defined keywords
    static vector <string> paramName( XRF_PARAMETER_LAST );
    paramName[ANODE_Z_INDEX] = "##ANODE";
    paramName[KV_INDEX] = "#BEAMKV";
    paramName[TUBE_INC_ANGLE_INDEX] = "##TUBEINCANG";
    paramName[TUBE_TAKEOFF_ANGLE_INDEX] = "##TUBETAKEOF";
    paramName[TUBE_BE_WINDOW_INDEX] = "##TUBEWINDOW";
    paramName[TUBE_CURRENT_INDEX] = "#EMISSION";
    paramName[FILTER_Z_INDEX] = "##FILTERZ";
    paramName[FILTER_THICK_INDEX] = "##FILTERTH";
    paramName[EXCIT_ANGLE_INDEX] = "##INCANGLE";
    paramName[EMERG_ANGLE_INDEX] = "#ELEVANGLE";
    paramName[AZIMUTH_ANGLE_INDEX] = "#AZIMANGLE";
    paramName[XTILT_ANGLE_INDEX] = "#XTILTSTGE";
    paramName[YTILT_ANGLE_INDEX] = "#YTILTSTGE";
    paramName[X_POSITION_INDEX] = "#XPOSITION";
    paramName[Y_POSITION_INDEX] = "#YPOSITION";
    paramName[Z_POSITION_INDEX] = "#ZPOSITION";
    paramName[SOURCE_SOLID_ANGLE_INDEX] = "##INCSR";
    paramName[DET_SOLID_ANGLE_INDEX] = "#SOLIDANGLE";
    paramName[GEOMETRY_INDEX] = "##GEOMETRY";
    paramName[PATH_TYPE_INDEX] = "##ATMOSPHERE";
    paramName[INC_PATH_LENGTH_INDEX] = "##PATHINCLEN";
    paramName[EMERG_PATH_LENGTH_INDEX] = "##PATHEMGLEN";
    paramName[WINDOW_TYPE_INDEX] = "##WINDOWTYPE";
    paramName[WINDOW_THICK_INDEX] = "##WINDOWTH";
    paramName[DETECTOR_TYPE_INDEX] = "#EDSDET";
    paramName[DET_RESOLUTION_INDEX] = "##DETRES";
    paramName[DET_BE_WINDOW_INDEX] = "#TBEWIND";
    paramName[DET_ACTIVE_THICK_INDEX] = "#TACTLYR";
    paramName[TEST_OPTIC_TYPE_INDEX] = "##OPTICFILE";
    paramName[MINIMUM_ENERGY_INDEX] = "##MINIMUM_EN"; //   So that default values can be set if needed (kept for compatibility)
    paramName[ENERGY_CORRECTION_SLOPE_INDEX] = "##DL_SLOPE"; //  Modified July 31, 2018 to add linear energy calibration correction
    paramName[ENERGY_CORRECTION_OFFSET_INDEX] = "##DL_OFFSET";
    string Optic_file_name("Optic file name");
    string Xray_tube_file_name("X-ray tube file name");

//    cout << "get_EMSA_keyword   " << index << "  " << XRF_PARAMETER_FIRST << "  " << XRF_PARAMETER_LAST << endl;
    if( XRF_PARAMETER_FIRST < index && index < XRF_PARAMETER_LAST ) return paramName[index];
    else if( index == XRF_PARAMETER_OPTIC_FILE ) return Optic_file_name;
    else if( index == XRF_PARAMETER_TUBE_FILE ) return Xray_tube_file_name;
    string dummy( "bad index" );
    return dummy;
};

const string get_EMSA_units( const int index, const int value ) {
// TODO: get rid of statics!!!
    //		return units for each conditions parameter read using an EMSA keyword (or user keyword)
    static vector <string> units_msa( XRF_PARAMETER_LAST );
    units_msa[ANODE_Z_INDEX] = "(Z)";
    units_msa[KV_INDEX] = "kV";
    units_msa[TUBE_INC_ANGLE_INDEX] = "deg";
    units_msa[TUBE_TAKEOFF_ANGLE_INDEX] = "deg";
    units_msa[TUBE_BE_WINDOW_INDEX] = "mm";
    units_msa[TUBE_CURRENT_INDEX] = "mA";
    units_msa[FILTER_Z_INDEX] = "(Z)";
    units_msa[FILTER_THICK_INDEX] = "micron";
    units_msa[EXCIT_ANGLE_INDEX] = "deg";
    units_msa[EMERG_ANGLE_INDEX] = "deg";
    units_msa[AZIMUTH_ANGLE_INDEX] = "deg";
    units_msa[XTILT_ANGLE_INDEX] = "deg";
    units_msa[YTILT_ANGLE_INDEX] = "deg";
    units_msa[X_POSITION_INDEX] = "mm";
    units_msa[Y_POSITION_INDEX] = "mm";
    units_msa[Z_POSITION_INDEX] = "mm";
    units_msa[SOURCE_SOLID_ANGLE_INDEX] = "sr";
    units_msa[DET_SOLID_ANGLE_INDEX] = "sr";
    units_msa[GEOMETRY_INDEX] = "";
    units_msa[PATH_TYPE_INDEX] = "";
    units_msa[INC_PATH_LENGTH_INDEX] = "cm";
    units_msa[EMERG_PATH_LENGTH_INDEX] = "cm";
    units_msa[WINDOW_TYPE_INDEX] = "";
    units_msa[WINDOW_THICK_INDEX] = "micron";
    units_msa[DETECTOR_TYPE_INDEX] = "";
    units_msa[DET_RESOLUTION_INDEX] = "eV";
    units_msa[DET_BE_WINDOW_INDEX] = "micron";
    units_msa[DET_ACTIVE_THICK_INDEX] = "mm";
    units_msa[TEST_OPTIC_TYPE_INDEX] = "";
    units_msa[MINIMUM_ENERGY_INDEX] = "eV"; //   So that default values can be set if needed (kept for compatibility)
    units_msa[ENERGY_CORRECTION_SLOPE_INDEX] = "eV/keV"; //  Modified July 31, 2018 to add linear energy calibration correction
    units_msa[ENERGY_CORRECTION_OFFSET_INDEX] = "eV";

    static string bad_value( "bad" );
    static string none_value( "none" );

//    cout << "get_EMSA_units   " << index << "  " << XRF_PARAMETER_FIRST << "  " << XRF_PARAMETER_LAST << endl;
    if( XRF_PARAMETER_FIRST < index && index < XRF_PARAMETER_LAST ) {
        if( value != -999 ) {
            //  Attempt to return more useful text for conditions that are enumerated choices
            if( index == TEST_OPTIC_TYPE_INDEX ) {
                if( value == 0 ) return none_value;
                static vector <string> optic_types( 6 );
                optic_types[1] = "none";
                optic_types[2] = "boxcar";
                optic_types[3] = "oldBB";
                optic_types[4] = "file";
                optic_types[5] = "newBB";
                if( value > 0 && value < optic_types.size() ) return optic_types[value];
                return bad_value;
            } else if( index == PATH_TYPE_INDEX ) {
                static vector <string> atm_types( 5 );
                atm_types[0] = "vac";
                atm_types[1] = "He";
                atm_types[2] = "Mars";
                atm_types[3] = "HeCO2";
                atm_types[4] = "air";
                XrayAtmosphere atm = (XrayAtmosphere) value;
                if( atm == VACUUM ) return atm_types[0];
                else if( atm == HELIUM ) return atm_types[1];
                else if( atm == MARS ) return atm_types[2];
                else if( atm == HE_MARS ) return atm_types[3];
                else if( atm == AIR || atm == EARTH ) return atm_types[4];
                else return bad_value;
            } else if( index == WINDOW_TYPE_INDEX ) {
                XrayWindowMaterials win = (XrayWindowMaterials) value;
                if( win == NO_WINDOW ) return none_value;
                static vector <string> win_types( 9 );
                win_types[1] = "B4C";
                win_types[2] = "Plas";
                win_types[3] = "CFRP";
                win_types[4] = "Zr";
                win_types[5] = "Al";
                win_types[6] = "Nylon";
                win_types[7] = "Nyl+Zr";
                win_types[8] = "Al2O3";
                if( win == B4C ) return win_types[1];
                else if( win == PLASTIC ) return win_types[2];
                else if( win == CFRP ) return win_types[3];
                else if( win == ZR ) return win_types[4];
                else if( win == AL ) return win_types[5];
                else if( win == NYLON ) return win_types[6];
                else if( win == NYLON_ZR ) return win_types[7];
                else if( win == AL2O3 ) return win_types[8];
                else return bad_value;
            } else if( index == DETECTOR_TYPE_INDEX ) {
                DetectorType det = (DetectorType) value;
                if( det == NO_DETECTOR ) return none_value;
                static vector <string> det_types( 5 );
                det_types[1] = "SiPIN";
                det_types[2] = "SDD";
                det_types[3] = "CdTe";
                det_types[4] = "HP-Ge";
                if( det == SI_PIN ) return det_types[1];
                else if( det == SI_SDD ) return det_types[2];
                else if( det == CD_TE ) return det_types[3];
                else if( det == HP_GE ) return det_types[4];
               else return bad_value;
            } else {
                return units_msa[index];
            }
        } else {
            return units_msa[index];
        }
    }
    static string dummy( "bad index" );
    return dummy;
};

void parseEMSAkeyword( const string strIn, const string delim, string &sKeyword, vector <string> &sValue ) {
// strIn - in, input string
// del - in, delimiter character
// sKeyword - out, keyword up to delimiter (not including delimiter)
// sValue - out, rest of string after delimiter and following blank, if any
//          entries in the sValue vector are separated by blanks in the input string
	const int slen = strIn.length();
	int j;
	j = strIn.find( delim );
	if( j >= 0 && j < slen ) {
		sKeyword = strIn.substr( 0, j );
        string delim( COMMA_CHARACTER );
        delim += BLANK_CHARACTER;
		parse_records( delim, strIn.substr( j + 1, slen ), sValue );
	} else {
		sKeyword = strIn.substr( 0, slen );
	};
};


void parseCommas( const string strIn, vector <string> &sValue ) {
// strIn - in, input string
// sValue - out, vector of substrings
//          entries in the sValue vector are separated by blanks in the input string
	const int slen = strIn.length();
	sValue.clear();
	if( slen <= 0 ) return;
	int j;
	j = -1;
	string tab( TAB_CHARACTER );
    // Separate values between blanks and move into output vector
    while( j+1 < slen ) {
        while( j+1<slen && strIn.substr( j + 1, 1 ) == BLANK_CHARACTER ) j++; //  skip any leading blanks
        while( j+1<slen && strIn.substr( j + 1, 1 ) == tab ) j++; //  skip any tabs
        if( strIn.substr( j + 1, 1 ) == DOUBLE_QUOTE_CHARACTER || strIn.substr( j + 1, 1 ) == SINGLE_QUOTE_CHARACTER ) {
            //  Keep part in quotes together
            j++;    //  Move past delimiter
            string temp_str;
            while( j+1<slen && strIn.substr( j + 1, 1 ) != DOUBLE_QUOTE_CHARACTER && strIn.substr( j + 1, 1 ) != SINGLE_QUOTE_CHARACTER ) {
                temp_str.append( strIn.substr( j + 1, 1 ) );
                j++;
            }
            sValue.push_back( temp_str );
            j++;    //  Move past delimiter
        } else {
            //  Separate at blanks
            int sep = strIn.find(BLANK_CHARACTER, j + 1);
            if( sep >= 0 && sep < slen ) {
                sValue.push_back( strIn.substr( j + 1, sep-j-1 ) );
                j = sep;
            } else {
                sValue.push_back( strIn.substr( j + 1, slen) );    //  To end of string
                j = slen;
            }
        }
    }
};

int parse_EMSA_description( const int index, const string &s ) {
//  Convert words given in EMSA format under certain keywords into numeric values
//  These numeric values must match the values in fpSetupConditions.cpp, function condQuant
//  So they are stored in header file XRFconditions.h

    string u = upper_trim( s );
    switch( index ) {
        case PATH_TYPE_INDEX:
            //XrayAirPath.h   enum XrayAtmosphere { Vacuum, Helium, Mars, HeMars, Air, Earth };
            if( u.substr(0,3) == "VAC" ) return VACUUM;
            else if( u.substr(0,2) == "HE" ) return HELIUM;
            else if( u == "MARS" ) return MARS;
            else if( u == "HE_MARS" ) return HE_MARS;
            else if( u == "AIR" ) return AIR;
            else if( u == "EARTH" ) return EARTH;
            else return -100 - index;
         case WINDOW_TYPE_INDEX:
            // XrayWindow.h  enum XrayWindowMaterials { B4C, Plastic, Brass, Zr, Al, Nylon, NylonZr, Al2O3 };
            if( u == "NONE" ) return NO_WINDOW;
            if( u == "B4C" ) return B4C;
            else if( u.substr(0,4) == "PLAS" ) return PLASTIC;
            else if( u == "CFRP" ) return CFRP;
            else if( u == "ZR" ) return ZR;
            else if( u == "AL" ) return AL;
            else if( u == "NYLON" ) return NYLON;
            else if( u == "NYLONZR" ) return NYLON_ZR;
            else if( u == "AL2O3" ) return AL2O3;
            else return -100 - index;
        case DETECTOR_TYPE_INDEX:
            //XrayDetector.h   enum DetectorType { Si = 0, SiPIN, SiSDD, CdTe, HPGe };
            if( u == "SIBEW" ) return SI_PIN;
            else if( u == "SDBEW" ) return SI_SDD;
            else if( u == "SIBEW" ) return SI_PIN;
            else if( u == "CDBEW" ) return CD_TE;
            else if( u == "GEBEW" ) return HP_GE;
            else return -100 - index;
        default:
            return 0;
    }
};

