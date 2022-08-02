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

#include <iostream>
#include <fstream>
#include <sstream>
#include "quantWriteMap.h"
#include "quantWriteResults.h"
#include "XrayMaterial.h"
#include "upper_trim.h"
#include "XRFconstants.h"
#include "XRFcontrols.h"
#include "XRFutilities.h"


using namespace std;

//  Appends one line to a map file

//  Written July 18, 2017
//  Modified Sept. 14, 2017 to generate fake data and to write map CSV files
//      (CSV files are presently read in Python and converted to PNG images)
//  Modified Nov. 8, 2017
//      Remove fake data and move to write_EDRhistogram_data function
//  Modified Mar. 2, 2018
//      Process -q option to specify outputs to map file
//          P=percent (default), I=intensity, E=error, K=ECF used, L=fit coefficient, Z=atomic number instead of symbol,
//          T=total counts, X=reduced chi squared, C=energy cal, R=det. res, N=number of iterations,
//          F=file name, S=sequence number
//      Load quantification results into element list in quantWriteResults and use here
//          to get spectrum fit info for each element and to put elements in same order as input list
//  Modified Mar. 7, 2018
//      Make quant map outputs case insensitive
//  Modified May 13, 2019
//      In XraySpectrum, all non-spectrum information put in separate structure
//      Make quant map outputs case sensitive and add options for all spectrum info (UC is used in quantitative calcs, LC is aux info)
//      Put output info in the order that they appear in quant map outputs
//  Modified May 14, 2019
//`     Change default quant map outputs to be entered in this file, not done here when options are blank
//      Finished radical rewrite to enable writing all of the possible information, including aux info, using upper and lower case letters
//      Correct all occurrences of quant map options to quant map outputs (to be consistent with argument list)
//      Add live time (V) and real time (M), and all aux and header info using lower case letters
//  Modified May 16, 2019
//      Add counts in 1-7.25 keV region to quant map outputs, move from quantWriteResults to XraySpectrum
//      Remove negative intensities from map and log output (leave coefficients as-is for now)
//  Modified May 21, 2019
//      Implement Evaluate action, add G and H output options (given and relative error vs given, from element list)
//  Modified May 21, 2019       (after Version 2.40)
//      Fix bug - no comma between energy start and energy per channel
//  Modified Nov. 4, 2019
//      Change err column to use total error entry (to include ECF and eventually certificate uncertainty)
//  Modified Dec. 16, 2019
//      Change map output options: seq# changed to Q, S added for element sum (to compare to 100%)
//      Element sum added to arguments so the calculation in quantWriteResults can be used
//  Modified Nov. 24, 2020
//      Write matrix effect factor, output option "W"
//  Modified Feb. 24, 2021  Add "U" option for spectrum aux info title (also used for standard names in calibrate and evaluate)


int quantWriteMapHeader(std::ostream &map_out_stream,
    const string &title,
    const string &quant_map_outputs,
    const vector <ElementListEntry> &element_list,
    const bool oxidesOutput)
{
    map_out_stream.setf( ios::fixed, ios::floatfield );

    // TIMTIME: Spectrum may not have a title in future, we may be reading from binary file or CSV file
    //  Add a title using the title from the input file or the file name if no title
    //if( spectrum.aux_info().titles.size() > 0 ) map_out_stream << spectrum.aux_info().titles[0] << endl;
    //else map_out_stream << "First spectrum file " << spectrum.file_name() << endl;

    // Instead, just printing what is passed to us
    map_out_stream << title << endl;

    unsigned int qmo_len_case = quant_map_outputs.length();
    //  See if the user wants the atomic number instead of the symbol in the headers
    bool atomic_number = quant_map_outputs.find("Z") < qmo_len_case;
    vector <string> header_labels;
    unsigned int opt_index;
    unsigned int ie;

    //  Generate strings of Element or oxide symbols to use in output headers
    for ( ie=0; ie<element_list.size(); ie++ ) {
        if( element_list[ie].qualifier == IGNORE
                || element_list[ie].qualifier == EXCLUDE
                || element_list[ie].qualifier == MATRIX ) continue;
        stringstream ofs_Z;
        ofs_Z << element_list[ie].element.Z();
        if( !oxidesOutput ) {
            //  Label element
            if( atomic_number ) {
                header_labels.push_back( ofs_Z.str() );
            } else {
                header_labels.push_back( element_list[ie].element.symbol() );
            }
        } else {
            // Label oxide
            if( atomic_number ) {
                header_labels.push_back( ofs_Z.str() + XrayMaterial::formula_string( element_list[ie].element, element_list[ie].stoichiometry, true ) );
            } else {
                header_labels.push_back( XrayMaterial::formula_string( element_list[ie].element, element_list[ie].stoichiometry ) );
            }
        }
    };

    //  Write headers for all of the information requested in the map file
    for( opt_index=0; opt_index<qmo_len_case; opt_index++ ) {
        //  Get next letter in output selections
        string single_opt = quant_map_outputs.substr( opt_index, 1 );
        if( opt_index != 0 ) map_out_stream << ", ";  //  Don't put this for the first entry on each line

        //  Individual element information (grouped by category with that category included together for all elements)
        if( single_opt == "P" ) {
            //  Write headers for percents (default if no entries in quant map outputs)
            for ( ie=0; ie<header_labels.size(); ie++ ) map_out_stream << (ie==0?"":", ") << header_labels[ie] << "_%";
        }
        else if( single_opt == "I" ) {
            //  Write headers for intensities
            for ( ie=0; ie<header_labels.size(); ie++ ) map_out_stream << (ie==0?"":", ") << header_labels[ie] << "_int";
        }
        else if( single_opt == "E" ) {
            //  Write headers for fit relative errors
            for ( ie=0; ie<header_labels.size(); ie++ ) map_out_stream << (ie==0?"":", ") << header_labels[ie] << "_err";
        }
        else if( single_opt == "L" ) {
            //  Write headers for fit coefficients
            for ( ie=0; ie<header_labels.size(); ie++ ) map_out_stream << (ie==0?"":", ") << header_labels[ie] << "_coeff";
        //  Diagnostic information
        }
        else if( single_opt == "K" ) {
            //  Write headers for Element Calibration Factors used for quantification
            for ( ie=0; ie<header_labels.size(); ie++ ) map_out_stream << (ie==0?"":", ") << header_labels[ie] << "_ECF";
        }
        else if( single_opt == "G" ) {
            //  Write headers for Element Calibration Factors used for quantification
            for ( ie=0; ie<header_labels.size(); ie++ ) map_out_stream << (ie==0?"":", ") << header_labels[ie] << "_Given";
        }
        else if( single_opt == "H" ) {
            //  Write headers for Element Calibration Factors used for quantification
            for ( ie=0; ie<header_labels.size(); ie++ ) map_out_stream << (ie==0?"":", ") << header_labels[ie] << "_errG";
        }
        else if( single_opt == "W" ) {
            //  Write headers for matrix effect factor found during fundamental parameters calculation
            for ( ie=0; ie<header_labels.size(); ie++ ) map_out_stream << (ie==0?"":", ") << header_labels[ie] << "_M";
        }
        else if( single_opt == "T" ) {
            //  Write header for total counts
            map_out_stream << "total_counts";
        }
        else if( single_opt == "X" ) {
            //  Write header for reduced chi squared
            map_out_stream << "chisq";
        }
        else if( single_opt == "C" ) {
            //  Write header for energy calibration
            map_out_stream << "eVstart, eV/ch";
        }
        else if( single_opt == "R" ) {
            //  Write header for detector resolution
            map_out_stream << "res";
        }
        else if( single_opt == "N" ) {
            //  Write header for number of iterations
            map_out_stream << "iter";
        }
        else if( single_opt == "F" ) {
            //  Write header for file name
            map_out_stream << "filename";
        }
        else if( single_opt == "S" ) {
            //  Write header for sequence number in file name
            map_out_stream << "sum_%";
        }
        else if( single_opt == "Q" ) {
            //  Write header for sequence number in file name
            map_out_stream << "seq#";
        }
        else if( single_opt == "V" ) {
            //  Write header for sequence number in file name
            map_out_stream << "livetime";
        }
        else if( single_opt == "M" ) {
            //  Write header for sequence number in file name
            map_out_stream << "realtime";
        }
        else if( single_opt == "7" ) {
            //  Write header for sequence number in file name
            map_out_stream << "region_counts";
        }
        //  Auxiliary information
        else if( single_opt == "x" ) map_out_stream << "X";
        else if( single_opt == "y" ) map_out_stream << "Y";
        else if( single_opt == "z" ) map_out_stream << "Z";
        else if( single_opt == "i" ) map_out_stream << "I";
        else if( single_opt == "j" ) map_out_stream << "J";
        else if( single_opt == "s" ) map_out_stream << "SCLK";
        else if( single_opt == "r" ) map_out_stream << "RTT";
        else if( single_opt == "d" ) map_out_stream << "DPC";
        else if( single_opt == "p" ) map_out_stream << "PMC";
        else if( single_opt == "e" ) map_out_stream << "Events";
        else if( single_opt == "t" ) map_out_stream << "Triggers";
        else if( single_opt == "o" ) map_out_stream << "Overflows";
        else if( single_opt == "u" ) map_out_stream << "Underflows";
        else if( single_opt == "b" ) map_out_stream << "baseline_samples";
        else if( single_opt == "a" ) map_out_stream << "Resets";
        else if( single_opt == "s" ) map_out_stream << "Saturates";
        else if( single_opt == "l" ) map_out_stream << "Fast_livetime";
        else if( single_opt == "n" ) map_out_stream << "USN";
        else if( single_opt == "U" ) map_out_stream << "Title";
        else {
            cout << "*** Invalid quant map output selection: " << single_opt << "   ****" << endl;
            return -10;
        }
    }   //  Get next letter in output selections

    map_out_stream << endl;
    return 0;
}


void quantWriteMapRow(std::ostream &map_out_stream,
    const string &quant_map_outputs,
    const vector <ElementListEntry> &element_list,

    const XrayDetector &detector,
    const XraySpectrum &spectrum,
    float element_sum)
{
    unsigned int qmo_len_case = quant_map_outputs.length();
    vector <int> header_indices;
    unsigned int opt_index;
    unsigned int ie;

    map_out_stream.setf( ios::fixed, ios::floatfield );

    //  Generate indices of elements to match headers
    for ( ie=0; ie<element_list.size(); ie++ ) {
        if( element_list[ie].qualifier == IGNORE
                || element_list[ie].qualifier == EXCLUDE
                || element_list[ie].qualifier == MATRIX ) continue;
        header_indices.push_back( ie );
    };

    //  Write the line of information to match the headers

    for( opt_index=0; opt_index<qmo_len_case; opt_index++ ) {
        if( opt_index != 0 ) map_out_stream << ", ";  //  Don't put this for the first entry on each line
       //  Get next letter in output selections
        string single_opt = quant_map_outputs.substr( opt_index, 1 );

        if( single_opt == "P" ) {
            //  Write Element percents (default if no entries in quant map outputs)
            map_out_stream.precision(4);
            for ( ie=0; ie<header_indices.size(); ie++ ) map_out_stream << (ie==0?"":", ") << element_list[ header_indices[ie] ].percent;
        }
        else if( single_opt == "I" ) {
            //  Write Element intensities
            map_out_stream.precision(1);
            for ( ie=0; ie<header_indices.size(); ie++ ) {
                if( element_list[ header_indices[ie] ].intensity >= 0 ) map_out_stream << (ie==0?"":", ") << element_list[ header_indices[ie] ].intensity;
                else map_out_stream << (ie==0?"":", ") << 0.0;
            }
        }
        else if( single_opt == "E" ) {
            //  Change err column to use total error entry (to include ECF and eventually certificate uncertainty)
            //  Write Element fit relative errors
            map_out_stream.precision(4);
            for ( ie=0; ie<header_indices.size(); ie++ ) map_out_stream << (ie==0?"":", ") << element_list[ header_indices[ie] ].total_err;
        }
        else if( single_opt == "L" ) {
            //  Write Element fit coefficients
            map_out_stream.precision(4);
            for ( ie=0; ie<header_indices.size(); ie++ ) map_out_stream << (ie==0?"":", ") << element_list[ header_indices[ie] ].coefficient;
        //  Diagnostic information
        }
        else if( single_opt == "K" ) {
            map_out_stream.precision(3);
            //  Write Element Calibration Factors used for quantification
            for ( ie=0; ie<header_indices.size(); ie++ ) map_out_stream << (ie==0?"":", ") << element_list[ header_indices[ie] ].ecf;
        }
        else if( single_opt == "G" ) {
            map_out_stream.precision(4);
            //  Write Element Given Percents (used for Evaluate)
            for ( ie=0; ie<header_indices.size(); ie++ ) map_out_stream << (ie==0?"":", ") << element_list[ header_indices[ie] ].given;
        }
        else if( single_opt == "H" ) {
            map_out_stream.precision(1);
            //  Write Relative Error vs Given (used for Evaluate)
            for ( ie=0; ie<header_indices.size(); ie++ ) map_out_stream << (ie==0?"":", ") << element_list[ header_indices[ie] ].rel_err_given;
        }
        else if( single_opt == "W" ) {
            map_out_stream.precision(3);
            //  Write Relative Error vs Given (used for Evaluate)
            for ( ie=0; ie<header_indices.size(); ie++ ) map_out_stream << (ie==0?"":", ") << element_list[ header_indices[ie] ].matrix;
        }
        else if( single_opt == "T" ) {
            //  Write total counts
            map_out_stream.precision(0);
            map_out_stream << spectrum.total_counts();
        }
        else if( single_opt == "X" ) {
            //  Write reduced chi squared
            map_out_stream.precision(2);
            map_out_stream << spectrum.chisq();
        }
        else if( single_opt == "C" ) {
            //  Write energy calibration
            map_out_stream.precision(1);
            map_out_stream << spectrum.calibration().energyStart();
            map_out_stream.precision(4);
            map_out_stream << ", ";
            map_out_stream << spectrum.calibration().energyPerChannel();
        }
        else if( single_opt == "R" ) {
            //  Write detector resolution
            map_out_stream.precision(0);
            map_out_stream << detector.resolution();
        }
        else if( single_opt == "N" ) {
            //  Write number of iterations
            map_out_stream << spectrum.iterations();
        }
        else if( single_opt == "F" ) {
            //  Write file name
            map_out_stream << spectrum.file_name();
        }
        else if( single_opt == "S" ) {
            //  Write header for sequence number in file name
            map_out_stream.precision(2);
            map_out_stream << element_sum;
        }
        else if( single_opt == "Q" ) {
            //  Write header for sequence number in file name
            map_out_stream << spectrum.seq_number();
        }
        else if( single_opt == "V" ) {
            //  Write live time
            map_out_stream.precision(2);
            map_out_stream << spectrum.live_time();
        }
        else if( single_opt == "M" ) {
            //  Write real time
            map_out_stream.precision(2);
            map_out_stream << spectrum.real_time();
        }
        else if( single_opt == "7" ) {
            //  Write real time
            map_out_stream.precision(0);
            map_out_stream << spectrum.region_counts();
        }
        //  Auxiliary information
        else if( single_opt == "x" ) map_out_stream << spectrum.aux_info().x;
        else if( single_opt == "y" ) map_out_stream << spectrum.aux_info().y;
        else if( single_opt == "z" ) map_out_stream << spectrum.aux_info().z;
        else if( single_opt == "i" ) map_out_stream << spectrum.aux_info().i;
        else if( single_opt == "j" ) map_out_stream << spectrum.aux_info().j;
        else if( single_opt == "s" ) map_out_stream << spectrum.aux_info().sclk;
        else if( single_opt == "r" ) map_out_stream << spectrum.aux_info().rtt;
        else if( single_opt == "d" ) map_out_stream << spectrum.aux_info().dpc;
        else if( single_opt == "p" ) map_out_stream << spectrum.aux_info().pmc;
        else if( single_opt == "e" ) map_out_stream << spectrum.header_info().events;
        else if( single_opt == "t" ) map_out_stream << spectrum.header_info().triggers;
        else if( single_opt == "o" ) map_out_stream << spectrum.header_info().overflows;
        else if( single_opt == "u" ) map_out_stream << spectrum.header_info().underflows;
        else if( single_opt == "b" ) map_out_stream << spectrum.header_info().baseline_samples;
        else if( single_opt == "a" ) map_out_stream << spectrum.header_info().preamp_resets;
        else if( single_opt == "s" ) map_out_stream << spectrum.header_info().saturates;
        else if( single_opt == "l" ) map_out_stream << spectrum.header_info().live_time_DSPC;
        else if( single_opt == "n" ) map_out_stream << spectrum.aux_info().usn;
        else if( single_opt == "U" ) {
            if( spectrum.aux_info().titles.size() > 0 ) map_out_stream << spectrum.aux_info().titles[0];
            else map_out_stream << " ";
        }
    }   //  Get next letter in output selections

    //  End of line
    map_out_stream << endl;
}
