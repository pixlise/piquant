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
//  read_spectrum_file.cpp
//  PIQUANT
//
//  Created by W. T. Elam on 1/18/2017.
//  Copyright (c) 2017 APL/UW. All rights reserved.
//

#include <iostream>
#include <fstream>
#include "read_spectrum_file.h"
#include "AmpTekRead.h"
#include "borehole_read.h"
#include "read_EMSA_PIXL.h"
#include "read_XIA_PIXL.h"
#include "upper_trim.h"
#include "XRFcontrols.h"
#include "XRFutilities.h"

//  Written Jan. 16, 2017
//      Parse list of element symbols and qualifiers for PIQUANT Subprocess
//      Qualifiers are separated from the element symbol by an underscore
//      Write helpful information to cout if errors and return a boolean (true if any errors)
//  Modified April 12, 2017 to handle upper and lower case file extensions
//  Modified April 14, 2017 to write file name here and write live time and energy calibration from file
//  Modified June 9, 2017
//      Write error messages with line numbers if appropriate
//  Modified July 18, 2017
//      Change "offset" to "eV start" in read OK message at end, for consistency
//  Modified Nov. 1, 2017 (preparation for implementing CSV standards file)
//      Change all checks for file extension (.txt, .msa, etc.) to check only last 4 characters in name, and be case insensitive
//  Modified Jan. 3, 2018
//      Fix energy start display when spectrum read OK
//  Modified Jan. 24, 2018
//      Add total counts to information written about each spectrum read
//  Modified Jan. 26, 2018
//      Fixed total counts output (was after new line)
//  Modified Mar. 2, 2018
//      Store spectrum file name in XraySpectrum objects
//  Modified Mar. 7, 2018
//      Use utility function to extract file name without path for spectrum file
//      Use utility to check file extensions
//  Modified Sep. 19, 2018
//      Write message about linear energy correction (REMOVED, already in main program)
//  Modified May 22, 2020 to put conditions vector and optic file name in struct
//                          add file name for X-ray tube spectrum input from external calculation
//  Modified Jan 13, 2021 to use XRFconditionsInput struct and to resize conditionsVector in this struct before using it in borehole_read

	using namespace std;

void print_spectrum_summary(const std::vector <XraySpectrum> &spectra, std::ostream &termOutFile);

const int read_spectrum_file( std::ostream &termOutFile, const std::string &spectrumPathName,
        std::vector <XraySpectrum> &spectra, XRFconditionsInput &condStruct_spec ) {

//  Determines the spectrum file type and calls appropriate routine to read it

    termOutFile << "Reading spectrum from file " << spectrumPathName << endl;

    //  Extract the file name from the full path for storing with the spectrum
    string spectrumFileName;
    string spectrumPathOnly;
    /*bool path_found =*/ extract_path( spectrumPathName, spectrumPathOnly, spectrumFileName );

    int result = 0;
    bool error = false;
    spectra.clear();
    vector <float> spectrum;
    float ev_start = 0, ev_ch = 0, live_time = 0, x = 0, y = 0, z = 0;
    string spectrum_upper_ext = upper_trim( spectrumPathName.substr( spectrumPathName.length() - FILE_EXTENSION_CHARS, FILE_EXTENSION_CHARS ) );
    if ( check_file_extension( spectrumPathName, "MCA" ) || check_file_extension( spectrumPathName, "MCS" ) ) {
        //  Decide if this is an AmpTek file or a ProSpect file (Ketek or XIA)
        ifstream mca_file_check( spectrumPathName.c_str() );
        string check_format;
        getline( mca_file_check, check_format );
        if( ! mca_file_check ) {
            termOutFile << "Can't open mca/mcs spectrum file " << spectrumPathName << endl;
            error = true;
        }
        mca_file_check.close();
        if( check_format.substr( 0, 17 ) == "<<PMCA SPECTRUM>>" ) {
            //			read spectrum from Amptek .mca file using function
            AmpTekSpec specData;
            result = amptek_read( spectrumPathName, specData );
            if ( result != 0 ) {
                termOutFile << "Can't read AmpTek spectrum file, result = " << result << "  for file " << spectrumPathName << endl;
                error = true;
            };
            XraySpectrum temp( specData.spectrum, specData.ev_start, specData.ev_ch );
            temp.live_time( specData.live_time );
            spectra.push_back( temp );
        } else if( check_format.substr( 0, 12 ) == "File Version" ) {
            //  Read spectrum file in format from ProSpect, from Ketek or XIA digital pulse processor
            XraySpectrum spec_data;
            string spec_acq_date;
            string spec_title;
            string spec_sample;
            string spec_unitID;
            result =  read_XIA_PIXL ( spectrumPathName, spec_data, spec_acq_date, spec_title,
                                        spec_sample, spec_unitID );
            if ( result != 0 ) {
                termOutFile << "Can't read Ketek/XIA spectrum file, result = " << result << "  for file name " << spectrumPathName << endl;
                if( result == -999999 ) {
                    termOutFile << "Can't open file or invalid file format." << endl;
                } else {
                    termOutFile << "Error on line " << -result << "." << endl;
                }
                error = true;
            };
            spectra.push_back( spec_data );
        } else {
            termOutFile << "*** Spectrum file first line not recognized." << endl;
            termOutFile << "It should be <<PMCA SPECTRUM>> or File Version but is " << check_format.substr( 0, 17 ) << endl;
        }
    } else if ( check_file_extension( spectrumPathName, "MSA" ) ) {
        //      open and read the ISO 22029 2012 EMSA format file
        vector <XraySpectrum> spectrum_vec;
        string acq_date, acq_time, x_label, y_label, unit;
        result = read_EMSA_PIXL( spectrumPathName, condStruct_spec, spectra );
        if ( result != 0 ) {
            termOutFile << "Can't read msa configuration file, result = " << result << "  for file name " << spectrumPathName << endl;
            if( result == -999999 ) {
                termOutFile << "Invalid file format or missing required keyword." << endl;
            } else {
                termOutFile << "Error on line number = " << -result << "." << endl;
            }
            error = true;
        };
    } else if ( check_file_extension( spectrumPathName, "XSP" ) ) {
         // Open and read the borehole XRF file format (older version of EMSA format)
        vector <string> spectrum_titles;
        result = borehole_read ( spectrumPathName, condStruct_spec.conditionsVector, spectrum, ev_start, ev_ch, live_time, spectrum_titles, x, y, z );
        if ( result != 0 ) {
            termOutFile << "Can't read XSP spectrum file, result = " << result << "  for file name " << spectrumPathName << endl;
            error = true;
        };
        XraySpectrum temp( spectrum, ev_start, ev_ch );
        temp.live_time( live_time );
        temp.aux_info_change().x = x;
        temp.aux_info_change().y = y;
        temp.aux_info_change().z = z;
        int m;
        for( m=0; m<spectrum_titles.size(); m++ ) temp.aux_info_change().titles.push_back( spectrum_titles[m] );
        spectra.push_back( temp );
    } else {
        termOutFile << "Can't read spectrum file, unrecognized file type, for file name " << spectrumPathName << endl;
        error = true;
    }

    if( ! error ) {
        //  Store file name
        int is;
        for( is=0; is<spectra.size(); is++ ) spectra[is].file_name( spectrumFileName );
        print_spectrum_summary(spectra, termOutFile);
    }

    return result;
}

void print_spectrum_summary(const std::vector <XraySpectrum> &spectra, std::ostream &termOutFile)
{
    termOutFile << "Spectrum read OK, " << spectra.size() << (spectra.size() != 1 ? " detectors" : " detector") << endl;

    for(int is=0; is<spectra.size(); is++ ) {
        termOutFile << "Detector " << is;
        termOutFile.precision(2);
        termOutFile << "  live time " << spectra[is].live_time();
        termOutFile << "    energy calibration ";
        termOutFile.precision(1);
        termOutFile << "  eV start = " << spectra[is].calibration().energyStart();
        termOutFile.precision(4);
        termOutFile << "  eV/ch = " << spectra[is].calibration().energyPerChannel();
        termOutFile.precision(0);
        termOutFile << "    total counts = " << spectra[is].total_counts();
        termOutFile << endl;
    }
}
