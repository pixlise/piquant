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
//  map_spectrum_file_increment.cpp
//  PIQUANT
//
//  Created by W. T. Elam on 7/18/2017.
//  Copyright (c) 2017 APL/UW. All rights reserved.
//

#include <iostream>
#include <sstream>
#include <math.h>
#include "XRFconstants.h"
#include "map_spectrum_file_increment.h"

//  Written July 18, 2017
//  Modified Sept. 27, 2017
//      Also process Seq type increments
//  Modified Jan. 3, 2017
//      Move tab, single quote, double quote, blank, comma, and underscore definitions to XRFconstants.h
//  Modified Jan. 30, 2017
//      Search for "e pt" for first EM scan, mudstone, yesterday (comment this out when finished)
//  Modified Mar. 2, 2018
//      Return sequence number for new file name
	using namespace std;

int map_spectrum_file_increment( std::string &spec_file_incr, int &sequence_number ) {

//  Increments the spectrum file name to the next spectrum composing the map
//  Trys several approaches to see if one fits the spectrum file name and yields a valid number
    const int no_number_found = - int( pow( PI, 13 ) ); //  A very unlikely number
    int n = no_number_found;
    int seq_no_start = -1;
    int seq_no_len = -1;

    int trial;
    for( trial=1; trial<3; trial++ ) {
        //  Look for the right pattern and isolate the characters that might be a number
        switch( trial ) {
            case 1: {   //  brackets to limit scope of variables declared in this case
                //  Look for a sequence number between the last underscore and the dot (PIXL breadboard Labview program)
                int underscore_pos = spec_file_incr.rfind( UNDERSCORE_CHARACTER );
                if( underscore_pos < 0 && underscore_pos >= spec_file_incr.length() ) continue;   //  Skip to next trial
                seq_no_start = underscore_pos + 1;
                int dot_pos = spec_file_incr.rfind( "." );
                if( dot_pos < 0 || dot_pos >= spec_file_incr.length() ) continue;   //  Skip to next trial
                seq_no_len = dot_pos - seq_no_start;
                break;  }   //  And avoid error "jump to case label"
            case 2: {    //  brackets to limit scope of variables declared in this case
                //  Look for a sequence number between the characters "Seq" and an underscore
                int seq_pos = spec_file_incr.rfind( "Seq" );
                if( seq_pos < 0 && seq_pos >= spec_file_incr.length() ) continue;   //  Skip to next trial
                seq_no_start = seq_pos + 3;
                int underscore_pos = spec_file_incr.substr( seq_no_start ).find( UNDERSCORE_CHARACTER );
            //    cout << "  _  " << seq_pos << "  " << seq_no_pos << "  " << underscore_pos << "  " << seq_no_len << endl;
                if( underscore_pos < 0 && underscore_pos >= spec_file_incr.substr( seq_no_start ).length() ) continue;   //  Skip to next trial
                seq_no_len = underscore_pos;
            //    cout << "  " << seq_pos << "  " << seq_no_pos << "  " << underscore_pos << "  " << seq_no_len << endl;
                break;  }   //  And avoid error "jump to case label"
/*            case 3: {   //      Search for "e pt" for first EM scan, mudstone
                //  Look for a sequence number between the characters "e pt" and an underscore
                int underscore_pos = spec_file_incr.rfind( "e pt" );
                cout << "incr " << spec_file_incr << "  " << underscore_pos << endl;
                if( underscore_pos < 0 && underscore_pos >= spec_file_incr.length() ) continue;   //  Skip to next trial
                seq_no_start = underscore_pos + 4;
                int dot_pos = spec_file_incr.rfind( UNDERSCORE_CHARACTER );
                cout << "incr " << seq_no_start << "  " << dot_pos << endl;
                if( dot_pos < 0 || dot_pos >= spec_file_incr.length() ) continue;   //  Skip to next trial
                seq_no_len = dot_pos - seq_no_start;
                break;  }   //  And avoid error "jump to case label"
*/
        }
        string ns = spec_file_incr.substr( seq_no_start, seq_no_len );
        stringstream ns_in( ns );
        ns_in >> n;
        if( !ns_in ) {
            n = no_number_found;
            continue;   //  Invalid number, skip to next trial
        }
    }
    //  See if a valid number was found
    if( n != no_number_found ) {
        sequence_number = n + 1;
        //  Increment it and put in file name
        ostringstream ns_incr_os;
        ns_incr_os << sequence_number;
        if( ! (!ns_incr_os) ) {
//           cout << "incr " << spec_file_incr << "  " << sequence_number << endl;
            spec_file_incr = spec_file_incr.substr( 0, seq_no_start ) + ns_incr_os.str() + spec_file_incr.substr( seq_no_start + seq_no_len );
        } else {
            return -2;
        }
    } else {
        return -1;
    }
    return 0;
}
