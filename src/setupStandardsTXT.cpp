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
#include "XRFutilities.h"
#include "setupStandardsTXT.h"

//  Oct. 25, 2013
//      Added minimum amount of an element to be included in standard composition
//          as calling argument, with default of zero
//  May 10, 2015
//      Changes from MTXRF codes - return comments from input file
//  Modified July 7, 2016
//      Change input file name from hard coded to argument
//  Modified July 6, 2016 to write to file rather than std out
//  Modified Nov. 9, 2016 to use ostream instead of ofstream (so it can be used with cout if desired)
//  Modified march 22, 2017 to add path from standards input file to all data files (assumes data is in same folder)
//  Modified Oct. 27, 2017
//      Standards list struct moved to "XRFstandards.h" with single spectrum file name (not vector)
//      Change this routine to use the element list entries from parse_element_list (included in above struct)
//      Use parse_element_string to interpret element symbols (including qualifiers)
//      Bring into line with new standards input routine that reads a CSV file
//  Modified Mar. 7, 2018
//      Use utility function to extract file path for spectrum file
//      Only add standards file path if no path included in spectrum file name read from standards file
//      Use utility to check file extensions
//  Modified July 31, 2018
//      Use spectrum file name as name of standard since no names in TXT standards file
//  Modified Nov. 24, 2020
//      Put element percent into element list as given percent

using namespace std;

int setupStandardsTXT( const std::string &standardsInputFileName, std::ostream &termOutFile,
                      std::vector <StandardInformation> &standards_out, const float miniumum_amount ) {

    standards_out.clear();

//		open file of standard compositions and spectrum file names
	string StdCalEntryName = standardsInputFileName;
	termOutFile << "Reading standard compositions from file " << StdCalEntryName << endl;
    if( miniumum_amount != 0 ) termOutFile << "Ignoring elements less than " << 10000 * miniumum_amount << " ppm" << endl;
	ifstream StdCalEntryFile(StdCalEntryName.c_str(), ios::in);
	if ( !StdCalEntryFile ) {
		termOutFile << "Cannot open standards file " << StdCalEntryName << endl;
		return -1;
	};
	bool error = false;
	vector <string> comment_list;
	vector <Element> element_list;
	vector <float> given;
    int line_number = 0;
    bool element_list_not_read = true;
//		read entries in standard compositions file
	while ( ! ( !StdCalEntryFile ) ) {
        bool entry_error = false;
//			read spectrum file name
        string input_str;
        getline( StdCalEntryFile, input_str );
        line_number++;
		//  Get rid of trailing CR if file is Windows line endings on Linux or Mac
		if( (int) (input_str.data()[input_str.length()-1]) == 13 ) input_str.erase( input_str.length()-1,1);
//		termOutFile << "setupStandardsTXT   in " << line_number << "  " << input_str.size() << " : " << input_str << ":  ! " << !StdCalEntryFile << endl;
		if ( !StdCalEntryFile ) break;
        if( input_str.size() <= 0 ) continue;
//			skip this line if comment
		if ( input_str.substr(0,2) == COMMENT_STRING ) {
			comment_list.push_back( input_str );
			continue;
		};
		if ( element_list_not_read ) {
            //  The element list must be the first non-comment entry
            istringstream el_input(input_str);
            int ne_in;
            el_input >> ne_in;
            if( ne_in <= 0 ) {
                termOutFile << "*** Element list has zero or negative number of entries. ***" << endl;
                error = true;
            } else {
                string skip;
                getline( StdCalEntryFile, skip );
                line_number++;
                element_list_not_read = false;
                continue;
            }
		};
		StandardInformation standard_entry;
        //  Get path from input file name (argument 1)
        string spectrumPathName;
        string standardFileOnly;
        bool path_found = extract_path( standardsInputFileName, spectrumPathName, standardFileOnly );
        //  Check for path in input spectrum file name and add path from standards file if none
        string path_dummy;
        string spectrumFileOnly;
        path_found = extract_path( input_str, path_dummy, spectrumFileOnly );
        if( path_found ) {
            standard_entry.spectrumFileName = input_str;
        } else {
            standard_entry.spectrumFileName = spectrumPathName + input_str;
        }
        //  Use spectrum file name as name of standard since no names in TXT standards file
        standard_entry.names.push_back( spectrumFileOnly );
//			read number of elements and percents
		int ne;
		StdCalEntryFile >> ne;
        line_number++;
//			skip this line and this entry if no elements (not a valid standard)
		if ( ne <= 0 ) {
			string skip;
			getline( StdCalEntryFile, skip );
			continue;
		};
//			read list of elements and percents for this sample
        int i;
		for ( i=0; i<ne; i++ ) {
			string e;
			StdCalEntryFile >> e;
			if( ! StdCalEntryFile ) {
                termOutFile << "Read error on standards file, line number " << line_number << endl;
				entry_error = true;
                //  Skip the rest of this line, see if anything can be salvaged
                string skip;
                getline( StdCalEntryFile, skip );
                break;
			}
			Element el;
            ElementListEntry el_entry;
            bool element_error = parse_element_string( e, el_entry );
			if( ! element_error ) {
				el = el_entry.element;
			} else {
				termOutFile << "Invalid element symbol or qualifier on line " << line_number << ",  " << e << endl;
				entry_error = true;
			};
			float p;
			StdCalEntryFile >> p;
			if ( p < 0 || p > 100 ) {
				termOutFile << "Invalid percent on line " << line_number << ",  Element " << el.symbol() << " " << p << endl;
				entry_error = true;
			};
			if( p > 0 && p >= miniumum_amount ) {
				element_list.push_back( el );
				given.push_back( p / 100 );
				el_entry.percent = p;
				el_entry.given = p;
				standard_entry.element_list.push_back( el_entry );
			}
		};
		if ( ! entry_error ) {
            XrayMaterial temp_mat( element_list, given );
            standard_entry.mat = temp_mat;
            standard_entry.comments = comment_list;
            standards_out.push_back( standard_entry );
            element_list.clear();
            given.clear();
            comment_list.clear();
            standard_entry.element_list.clear();
		} else {
            error = true;
		}
	};
	StdCalEntryFile.close();

	if( error ) {
		return -2;
	} else {
		return 0;
	};
};
