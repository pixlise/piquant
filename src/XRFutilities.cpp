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

#include <sstream>
#include <time.h>
#include <ctype.h>  //  For isdigit
#include "upper_trim.h"
#include "XRFconstants.h"
#include "XRFcontrols.h"
#include "XRFutilities.h"

//	Modified April 7, 2019	To match Beam Geometry Tool
//							Add suffix functions
//							Prevent check_file_extension crash for very short file names
//  Modified May 6, 2018
//      compile problem in convert_to_float, and some new functions were left out of header
//  Modified Dec. 16, 2020
//      oxideFormulaString moved to XrayMaterial class

using namespace std;

bool extract_path( const string &full_path_in, string &path_out,
        string &filename_out) {
    //  Separate path and file name from full path
    //  Path separator characters are defined in XRFconstants.h, currently \ for Windows and / for Unix (and Mac)
    int slash_pos = full_path_in.rfind( PATH_SEPARATOR_WINDOWS );  //    Try backslash first as it is less commonly used for anything else
    if ( slash_pos <= 0 || slash_pos >= full_path_in.length()-1 ) slash_pos = full_path_in.rfind( PATH_SEPARATOR_UNIX );
    if ( slash_pos <= 0 || slash_pos >= full_path_in.length()-1 ) slash_pos = -1;
    if( slash_pos >= 0 ) {
        //  Leave separator in path so it can be used with another file name
        path_out = full_path_in.substr( 0, slash_pos+1 );
        filename_out = full_path_in.substr( slash_pos+1, 99999 );
        return true;
    } else {
        path_out.clear();
        filename_out = full_path_in;
        return false;
    }
    return false;
};


bool check_file_extension( const string &file_name_in , const string &file_extension ) {
    //  Check for a file name with the extension specified, case-insensitive
    //  FILE_EXTENSION_CHARS is defined in #include XRFcontrols.h
    if( file_name_in.length() < FILE_EXTENSION_CHARS ) return false;
    string spectrum_upper_ext = upper_trim( file_name_in.substr( file_name_in.length() - FILE_EXTENSION_CHARS, FILE_EXTENSION_CHARS ) );
    string ext_upper = upper_trim( "." + file_extension );
    if( 0 <= spectrum_upper_ext.rfind( ext_upper ) && spectrum_upper_ext.rfind( ext_upper ) < FILE_EXTENSION_CHARS )
        return true;
    return false;
};



//  Get the local clock time and convert to ASCII (human-readable ISO format) to put in the header
//  This function isolates the system and library calls for this purpose, which may change
string datetime()
{
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];

    time (&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer,80,"%Y-%m-%d %H:%M:%S",timeinfo);
    return std::string(buffer);
}


//  Utility function to strip off suffix from number
string strip_suffix( const string &str_in ) {
    string str_out;
    unsigned int i;
    for( i=0; i<str_in.length(); i++ ) {
        string number_test( str_in.substr( i, 1 ) );
        if( isdigit( number_test.c_str()[0] ) || number_test == "." || number_test == "," || number_test == "-" || number_test == "+" || number_test == "e" ) str_out += str_in.substr( i, 1 );
        else break;
    }
    return str_out;
}


//  Utility function to convert text records to float numbers
int convert_to_float( const vector <string> &records_in, const int start, const int stop, vector <float> &float_out ) {
	float_out.clear();
	if( start < 0 || start >= records_in.size() ) return -100;
	if( stop < 0 || stop >= records_in.size() ) return -100;
	int k;
	for( k=start; k<=stop; k++ ) {
		istringstream value( strip_suffix( records_in[k] ) );
		float temp = UNLIKELY_VALUE;	//	Unlikely value
		value >> temp;
		if ( ! value || temp == UNLIKELY_VALUE ) return -k;
		else float_out.push_back( temp );
	}
	return 0;
}

bool check_element_input( const std::string symbol_Z_in, Element &element_out ) {
    //  Set up Element object from input string
    bool error = false;
    //  Check for valid symbol
    if( Element::check_symbol( symbol_Z_in ) ) {
        Element temp( symbol_Z_in );
        element_out = temp;
    } else {
        //  If this is not a valid element symbol, check for an atomic number
        istringstream temp_instr( symbol_Z_in );
        int z_test = 0;
        temp_instr >> z_test;
//        cout << "check_element_input " << symbol_Z_in << "  " << !temp_instr << "  " << z_test << "  " << Element::check_Z( z_test ) << endl;
        if( ! temp_instr || ! Element::check_Z( z_test ) ) {
            error = true;
        } else {
            Element temp( z_test );
            element_out = temp;
        }
    }
    return !error;
}
