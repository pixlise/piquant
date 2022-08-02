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
//  parse_records.cpp
//  PIQUANT
//
//  Created by W. T. Elam on 1/14/2017.
//  Copyright (c) 2017 APL/UW. All rights reserved.
//

#include <iostream>
#include "parse_records.h"
#include "XRFconstants.h"

//  Written Jan. 15, 2017
//      Parse a string of values separated by any of the given delimiters
//      Return a vector of strings with each value separated
//      Handle quoted strings using single or double quotes
//      Handle end of string during quotes (terminate quotes and return value)
//      Absorb leading and trailing blanks
//      Treat tab as a single blank unless tab is a delimiter
//  Modified June 9, 2017
//      Fix bug in finding tabs in delimiter list
//  Modified Jan. 3, 2017
//      Move tab, single quote, double quote, blank, comma, and underscore definitions to XRFconstants.h

using namespace std;

int parse_records( const std::string delimiters_in, const std::string str_in,
        std::vector <std::string> &records_out) {

    const int slen = str_in.length();
	records_out.clear();
	if( slen <= 0 ) return 0;
	//  If delimiters not given use comma and space
	string delimiters( delimiters_in );
	if( delimiters.length() <= 0 ) {
        delimiters += COMMA_CHARACTER;
        delimiters += BLANK_CHARACTER;
    }
	//  Process quotes and tabs if they are not already delimiters
	string single_quote( SINGLE_QUOTE_CHARACTER );
	string double_quote( DOUBLE_QUOTE_CHARACTER );
	string tab( TAB_CHARACTER );
	string quotes;
	if( delimiters.find( single_quote ) < 0
            || delimiters.find( single_quote ) > delimiters.length() )
            quotes += single_quote;
	if( delimiters.find( double_quote ) < 0
            || delimiters.find( double_quote ) > delimiters.length() )
            quotes += double_quote;
    //  Don't treat tab as a blank if it is a delimiter
    int l;
    for( l=0; l<delimiters.length(); l++ ) {
        if( delimiters.substr( l, 1 ) == tab ) tab.clear();
    }
    string blank( BLANK_CHARACTER );
    bool blank_delimiter = false;
    int b = delimiters.find( blank );
    if( b >= 0 && b < delimiters.length() ) blank_delimiter = true;
	int j = 0;
	string temp;
    // Separate values between delimiters and move into output vector
    while( j < slen ) { //   record loop
        bool extra_record = false;  //  In case there is a delimiter at the end of the string
        //  skip any leading blanks or tabs at beginning of record
        bool skip = true;
        while( j<slen && skip ) {
            skip = ( str_in.substr( j, 1 ) == blank ) || ( str_in.substr( j, 1 ) == tab );
            if( skip ) j++;
        }
        //  Is this record quoted?   (avoid out_of_range if this is the end of the string)
        int quoted = (j<slen) ? quotes.find( str_in.substr( j, 1 ) ) : -1;
        //  If so, scan for the same quote mark while moving characters to record
        if( quoted >= 0 && quoted < quotes.length() ) {
            j++;    //  Skip over quote mark
            while( j<slen && str_in.substr( j, 1 ) != quotes.substr( quoted, 1 ) ) {
                temp += str_in.substr( j, 1 );
                j++;
            }
            j++;    //  Skip over final quote mark
            skip = !blank_delimiter;   //  Skip trailing blanks or tabs unless blank is a delimiter
            while( j<slen && skip ) {
                skip = ( str_in.substr( j, 1 ) == blank ) || ( str_in.substr( j, 1 ) == tab );
                if( skip ) j++;
            }
            //  If the next character is not a delimiter (or the end), report an error
            int d = (j<slen) ? delimiters.find( str_in.substr( j, 1 ) ) : 0;
            if( d < 0 || d >= delimiters.length() ) {
                return -j;
            } else {
                j++;    //  Skip over delimiter
                if( j == slen && str_in.substr( j-1, 1 ) != blank && str_in.substr( j-1, 1 ) != tab )
                        extra_record = true;    //  last character was a non-blank delimiter
            }
        } else {
        //  If not quoted, scan for the next delimiter while moving characters to record
            while( j<slen ) {
                int d = delimiters.find( str_in.substr( j, 1 ) );
                if( blank_delimiter && (str_in.substr( j, 1 ) == tab ) ) d = 0;
                if( d >= 0 && d < slen ) {
                    j++;    //  Skip over delimiter
                    if( j == slen && str_in.substr( j-1, 1 ) != blank && str_in.substr( j-1, 1 ) != tab )
                        extra_record = true;    //  last character was a non-blank delimiter
                    break;
                } else {
                    temp += str_in.substr( j, 1 );
                    j++;
                }
            }
        }
        records_out.push_back( temp );
        temp.clear();
        //  Put an empty record at the end if the last character was a non-blank delimiter
        if( extra_record ) records_out.push_back( temp );
    }   //   record loop


    return 0;
}
