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
#include "Element.h"
#include "XrayEdge.h"
#include "upper_trim.h"
#include "XRFconstants.h"
#include "parse_records.h"
#include "XRFutilities.h"
#include "read_tube_spectrum.h"

//  This function reads an X-ray tube spectrum calculated by PENELOPE (or other external program)

//  Written May 19, 2020



using namespace std;

int read_tube_spectrum( const std::string &tubeSpectrumFileName, std::vector <XrayLines> &tube_lines_out, float &tube_kV,
                            std::vector <float> &brem_energy, std::vector <float> &brem_spec, std::string &title ) {

    //  Initialize output arguments
    tube_kV = 0;
    tube_lines_out.clear();
    brem_energy.clear();
    brem_spec.clear();

    //  Avoid too many warnings or errors is file is wrong format
    int max_warnings = 10;
    int n_warnings = 0;
    int max_errors = 10;
    int n_errors = 0;

//		open file of calculated tube spectrum information
	ifstream tubeInputFile( tubeSpectrumFileName.c_str(), ios::in );
	if ( !tubeInputFile ) {
		return -1;
	};
	bool error = false;
    int line_number = 0;

    //  Read title line (first line, can be anything)
    getline( tubeInputFile, title );
    line_number++;

//		read and process lines in file
	while ( ! ( !tubeInputFile ) ) {
        bool entry_error = false;
        string input_str;
        getline( tubeInputFile, input_str );
        line_number++;
		//  Get rid of trailing CR if file is Windows line endings on Linux or Mac
		if( input_str.length()>0 && (int) (input_str.data()[input_str.length()-1]) == 13 ) input_str.erase( input_str.length()-1,1);
//		cout << "read_tube_spectrum   line # " << line_number << "  " << input_str.size() << " : " << input_str << ":  ! " << !tubeInputFile << endl;
		if ( !tubeInputFile ) break;    //  Stop at end of file
        if( input_str.size() <= 0 ) continue;	//  Skip empty line
        //  Parse line into comma separated fields
        vector <string> records;
        int result = parse_records( COMMA_CHARACTER, input_str, records );
        if( result < 0 ) {
            if( n_errors <= max_errors ) cout << "*** Error parsing comma separated entries on line " << line_number << ". ***" << endl;
            error = true;
            n_errors++;
            continue;
        }
        if( records.size() <= 0 ) continue;  //  No entries, skip this line

        //  Check keyword
        string keyword = upper_trim( records[0] );
//        cout << "read_tube_spectrum   Keyword " << keyword << "    number of records " << records.size() << endl;
		if( keyword == "COMMENT" ) continue;
        //	Characteristic emission lines from X-ray tube anode
        else if( keyword == "LINES" ) {
            //  Check units to verify file format
            if( records[4] != "ph/sec/sr/mA" ) {
                cout << "*** Error: invalid units or format for LINES on line " << line_number << ". ***" << endl;
                error = true;
                n_errors++;
                continue;
            }
            //  Read number of emission line entries
            vector <float> float_values;
            result = convert_to_float( records, 1, 1, float_values );
            if( result < 0 ) {
                entry_error = true;
                continue;
            }
            int n_lines = float_values[0];
            if( n_lines < 0 ) entry_error = true;
            //  Set up Element object for X-ray tube anode
            Element tube_anode;
            bool good_element = check_element_input( records[2], tube_anode );
            if( ! good_element ) {
                cout << "*** Error: invalid element symbol " << records[2] << " on line " << line_number << ". ***" << endl;
                error = true;
                n_errors++;
            }
            //  Read tube acceleration voltage used in calculation
            if( records[3].length() > 0 ){
                result = convert_to_float( records, 3, 3, float_values );
                tube_kV = float_values[0];
                if( tube_kV < 0 ) entry_error = true;
            }
            //  Set up list of emission lines for tube anode
            vector <EdgeIndex> edgeList;
            XrayEdge::numberOfEdges ( edgeList, tube_anode, tube_kV * 1000 );
            unsigned int edgeIndex;
            for ( edgeIndex=0; edgeIndex<edgeList.size(); edgeIndex++ ) {
                XrayEdge newEdge ( tube_anode, edgeList[edgeIndex] );
                XrayLines newLine ( newEdge );
                int lineIndex;
                for ( lineIndex=0; lineIndex<newLine.numberOfLines(); lineIndex++ ) {
                    newLine.factor( lineIndex, 0 );
                }
                tube_lines_out.push_back( newLine );
            }
            //  Read and process characteristic emission lines from anode
//            cout << "TUBE LINES  n " << n_lines << "   kV " << tube_kV << "   anode " << tube_anode.symbol() << endl;
            int i_line = 0;
            while( i_line<n_lines ) {
                getline( tubeInputFile, input_str );
                line_number++;
                if ( !tubeInputFile ) break;    //  Stop at end of file
                if( input_str.size() <= 0 ) continue;	//  Skip empty line
                int result = parse_records( COMMA_CHARACTER, input_str, records );
                if( result < 0 ) {
                    if( n_errors <= max_errors ) cout << "*** Error parsing comma separated entries on line " << line_number << ". ***" << endl;
                    error = true;
                    n_errors++;
                    continue;
                }
                if( records.size() <= 0 ) continue;  //  No entries, skip this line
                i_line++;   //  An actual emission line entry, error or no, so count it
                //  Check element against anode
                Element tube_anode_check;
                good_element = check_element_input( records[0], tube_anode_check );
                if( ! good_element ) {
                    if( n_warnings <= max_warnings ) cout << "Warning: element on line " << line_number << " is invalid or does not match anode - skipped." << endl;
                    n_warnings++;
                }
                //  Find the input emission line in the list for this anode
                bool found = false;
                int lineIndex = 0;
                for ( edgeIndex=0; edgeIndex<tube_lines_out.size(); edgeIndex++ ) {
                    for ( lineIndex=0; lineIndex<tube_lines_out[edgeIndex].numberOfLines(); lineIndex++ ) {
                        string symbol_IUPAC = tube_lines_out[edgeIndex].symbolIUPAC( lineIndex );
                        string symbol_Siegbahn = tube_lines_out[edgeIndex].symbolSiegbahn( lineIndex );
                        //  Convert to symbol format in input file
                        unsigned int p = symbol_IUPAC.find( "," );
                        if( p < symbol_IUPAC.length() ) symbol_IUPAC.replace( p, 1, "_" );
                        p = symbol_Siegbahn.find( "," );
                        if( p < symbol_Siegbahn.length() ) symbol_Siegbahn.replace( p, 1, "_" );
                        if( symbol_IUPAC == records[1] || symbol_Siegbahn == records[2] ) {
                            found = true;
                            break;
                        }
                    }
                    if( found ) break;
                }
                if( found ) {
                    //  Read the calculated intensity for this emission line
                    result = convert_to_float( records, 3, 3, float_values );
                    if( result < 0 ) {
                        entry_error = true;
                    } else {
                        //  Set line factor to get proper intensity in ph/sec/sr/mA
                        float relative_intensity = tube_lines_out[edgeIndex].relative( lineIndex );
                        float intensity_factor = float_values[0] / relative_intensity;
                        tube_lines_out[edgeIndex].factor( lineIndex, intensity_factor );
                    }
                }   //  Some very weak lines included in calculations, ignore
            }
        //  Continuum emission (Bramsstrahlung) spectrum
        } else if( keyword == "CONTINUUM" ) {
//            cout << "TUBE CONTINUUM  rec " << records[1] << "  " << records[2] << "  " << records[3] << endl;
             //  Check units to verify file format
            if( records[2] != "eV" || records[3] != "ph/sec/keV/sr/mA" ) {
                cout << "*** Error: invalid units or format for LINES on line " << line_number << ". ***" << endl;
                error = true;
                n_errors++;
                continue;
            }
            //  Read number of continuum entries
            vector <float> float_values;
            result = convert_to_float( records, 1, 1, float_values );
            int n_brem = 0;
            if( result < 0 ) {
                entry_error = true;
            } else {
                n_brem = float_values[0];
                if( n_brem < 0 ) entry_error = true;
            }
//            cout << "TUBE CONTINUUM  n " << n_brem << endl;
            //  Read and process continuum entries
            int i = 0;
            while( i<n_brem ) {
                getline( tubeInputFile, input_str );
                line_number++;
                if ( !tubeInputFile ) break;    //  Stop at end of file
                if( input_str.size() <= 0 ) continue;	//  Skip empty line
                int result = parse_records( COMMA_CHARACTER, input_str, records );
                if( result < 0 ) {
                    if( n_errors <= max_errors ) cout << "*** Error parsing comma separated entries on line " << line_number << ". ***" << endl;
                    error = true;
                    n_errors++;
                    continue;
                }
                if( records.size() <= 0 ) continue;  //  No entries, skip this line
                result = convert_to_float( records, 0, 1, float_values );
                i++;
                brem_energy.push_back( float_values[0] );
                brem_spec.push_back( float_values[1] );
            }
        } else if( n_warnings <= max_warnings ) {
            cout << "Warning - unrecognized keyword on line " << line_number << "." << endl;
            n_warnings++;
        }   //  Check keywords
        //  Check for error among entries and write error message
        if( entry_error && n_errors <= max_errors ) {
            cout << "*** Invalid value on line " << line_number << ". ***" << endl;
            error = true;
            n_errors++;
        }
        if( n_warnings >= max_warnings || n_errors >= max_errors ) break;
	}   //		read and process lines

	tubeInputFile.close();

	if( error ) return -1000 - n_errors;
	return 0;

};
