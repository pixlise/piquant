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
#include <iomanip>
#include <string>
#include <sstream>
#include "XRFconstants.h"
#include "XRFcontrols.h"
#include "parse_records.h"
#include "histogram_from_SDD_data.h"

//  Converts the output of the PIXL SEND_ADD_DATA command to an X-ray histogram
//  The binary data is broken up into a series of 16-bit integers
//  Each 16-bit integer is represented as a decimal number and stored in a CSV file
//  Each line of the CSV file is the data from one SEND_SDD_DATA command
//  It usually has the data from two detectors
//  Binary data format is from "JPL-XIA_PIXL_FPGA_Specification_v2.06.pdf"

//  Written Nov. 8, 2017
//  Completed testing using simulated iFSW data on Dec. 6, 2017
//      (Still need to test with data from actual hardware)
//  Modified Dec. 10, 2017
//      Move maximum number of error messages definition to XRFcontrols.h
//  Modified Dec. 15, 2017
//      Add check for minimum number of channels in histogram
//  Modified Jan. 3, 2017
//      Move tab, single quote, double quote, blank, comma, and underscore definitions to XRFconstants.h
// Modified May 13, 2019
//     In XraySpectrum, all non-spectrum information put in separate structure
//     Fix XIA live time calculation to include overflows and underflows, also use real time (not DSPC live time)

using namespace std;

//  Setup description of binary data locations
#define SDD_DATA_OFFSET 1       //  Start of the histogram data from the DPP FPGA as described in above document
#define SDD_DATA_INCREMENT 2    //  Each unity increment moves to the next 16-bit integer in the input row
                                //  Set to 2 if each entry has an associated 16-bit address
#define SDD_TIME_UNITS 500e-9       //  500 nanoseconds (in seconds)
#define SDD_DATA_SHIFT16 65536       //  Two to the sixteenth power
#define SDD_DATA_HISTOGRAMS_PER_LINE 2       //  Two histograms on each line (will only process one if only one is present)
//  Tag word 0xAA55 precedes statistics data
#define SDD_TAGWORD1_LENGTH 1
#define SDD_DPPSTATUS_LENGTH 1
#define SDD_RUNSTATUS_LENGTH 1
#define SDD_TAGWORD1_VALUE 0xAA55
//  Measured real time (while GATE=0), in 500 ns units, 48-bits (3 words, low word first)
#define SDD_REALTIME_LENGTH 3
//  Measured trigger live time (time under threshold*), in 500 ns units, 48-bits (3 words, low word first)
#define SDD_LIVETIME_LENGTH 3
//  Total number of events in the spectrum 32-bits (2 words, low word first)
#define SDD_EVTSINRUN_LENGTH 2
//  Total number of triggers (threshold crossings*) detected 32-bits (2 words, low word first)
#define SDD_TRIGGERS_LENGTH 2
//  Total number of overflows detected 32-bits (2 words, low word first)
#define SDD_OVERFLOWS_LENGTH 2
//  Total number of underflows detected 32-bits (2 words, low word first)
#define SDD_UNDERFLOWS_LENGTH 2
//  Total number of baseline samples acquired 32-bits (2 words, low word first)
#define SDD_BASEEVENTS_LENGTH 2
//  Total number of preamplifier resets detected (ADC excursions below ADCMIN) 32-bits (2 words, low word first)
#define SDD_PRERESETS_LENGTH 2
//  Total number of ADC excursions above ADCMAX 32-bits (2 words, low word first)
#define SDD_SATURATES_LENGTH 2
//  Seven (7) reserved locations in SRAM
#define SDD_RESERVED_LENGTH 7
//  Number of bins in the spectrum
#define SDD_MCALIMHI_LENGTH 1
//  Tag word 0x55AA precedes spectrum data
#define SDD_TAGWORD2_LENGTH 1
#define SDD_TAGWORD2_VALUE 0x55AA
//  First 32-bit bin value in histogram
#define SDD_BINWORD_LENGTH 2

//  Helper function to reduce code size and complexity (defined at end of file)
bool parse_one_sdd_entry( int &index, const int length, const vector <string> &records,
            const int line_number, double &value_out );
//  Note that this function moves the first argument forward to the next record

int histogram_from_SDD_data( const std::string sdd_file_name,
                const std::string fake_edr_file_name, std::vector <XraySpectrum> &spec_vec_out ) {

    spec_vec_out.clear();

//		open file of 16-bit integers is CSV format from PIXL SEND_SDD_DATA command
	cout << "Reading PIXL SDD data (in CSV format) from file " << sdd_file_name << endl;
	ifstream sdd_data_file(sdd_file_name.c_str(), ios::in);
	if ( !sdd_data_file ) {
		cout << "Cannot open SDD Data CSV file " << sdd_file_name << endl;
		return -1;
	};

//  See if we want to write a fake EDR file for EM display by Fang Zhong's interface
    int line_number = 0;
    int error_number = 0;
//		read and process lines in standard compositions file
	while ( ! ( !sdd_data_file ) ) {
        if( error_number > MAX_ERROR_MESSAGES ) {
        cout << "*** Processing terminated after too many errors. ***" << endl;
        return -1;
        }
        string input_str;
        getline( sdd_data_file, input_str );
        line_number++;
        //  Get rid of trailing CR if file is Windows line endings on Linux or Mac
        if( (int) (input_str.data()[input_str.length()-1]) == 13 ) input_str.erase( input_str.length()-1,1);
//		cout << "histogram_from_SDD_data   in " << line_number << "  " << input_str.size() << " : " << input_str << ":  ! " << !sdd_data_file << endl;
        //  Stop at end of file
        if ( !sdd_data_file ) break;
        //  Skip empty line
        if( input_str.size() <= 0 ) continue;
        //  Parse line into comma separated fields
        vector <string> records;
        int result = parse_records( COMMA_CHARACTER, input_str, records );
        if( result < 0 ) {
            cout << "*** Error parsing comma separated entries on line " << line_number << ". ***" << endl;
            error_number++;
            continue;
        }
        cout << "Records " << records.size() << endl;
        if( records.size() <= 0 ) continue;  //  No entries, skip this line

        //  ***************************************************************************
        //  Interpret input integers and convert to internal format for X-ray histogram
        //      (Histogram is X-ray spectrum prior to energy calibration)
        //  ***************************************************************************

        //  Attempt to process more than one histogram on each line
        int sdd_data_position = SDD_DATA_OFFSET;
        int i_hist_line;
        for( i_hist_line=0; i_hist_line<SDD_DATA_HISTOGRAMS_PER_LINE;i_hist_line++ ) {
            XraySpectrum temp_spec;

//            int ir;
//            cout << std::hex;
//            for( ir=0; ir<64; ir++ ) cout << ir << "  " << ir*2+sdd_data_position << "  " << records[ir*2+sdd_data_position] << endl;
//            cout << std::dec;

            //  Tag word 0xAA55 precedes statistics data
            double value;
            bool entry_error = parse_one_sdd_entry( sdd_data_position, SDD_TAGWORD1_LENGTH, records,
                line_number, value );
            if( entry_error || int(value) != SDD_TAGWORD1_VALUE ) {
                cout << "*** Error - tag word preceding statistics is missing or has incorrect value: ";
                cout << int( value ) << " on line " << line_number << " should be " << SDD_TAGWORD1_VALUE<< ". ***" << endl;
                if( i_hist_line > 0 ) cout << "The problem may be an incorrect number of channels in the previous histogram on this line." << endl;
                error_number++;
                continue;
            }

            //  Skip RUN STATUS and DPP STATUS registers
            sdd_data_position += SDD_DPPSTATUS_LENGTH * SDD_DATA_INCREMENT;
            sdd_data_position += SDD_RUNSTATUS_LENGTH * SDD_DATA_INCREMENT;

            //  Measured real time (while GATE=0), in 500 ns units, 48-bits (3 words, low word first)
            entry_error = parse_one_sdd_entry( sdd_data_position, SDD_REALTIME_LENGTH, records,
                line_number, value );
            if( entry_error ) error_number++;
            value *= SDD_TIME_UNITS;
            temp_spec.real_time( float( value ) );

            //  Measured trigger live time (time under threshold), in 500 ns units, 48-bits (3 words, low word first)
            entry_error = parse_one_sdd_entry( sdd_data_position, SDD_LIVETIME_LENGTH, records,
                line_number, value );
            if( entry_error ) error_number++;
            value *= SDD_TIME_UNITS;
            temp_spec.header_info_change().live_time_DSPC = float( value ); //  This is not the actual live time, see below

            //  Total number of events in the spectrum 32-bits (2 words, low word first)
            entry_error = parse_one_sdd_entry( sdd_data_position, SDD_EVTSINRUN_LENGTH, records,
                line_number, value );
            if( entry_error ) error_number++;
            temp_spec.header_info_change().events = int( value );

            //  Total number of triggers (threshold crossings) detected 32-bits (2 words, low word first)
            entry_error = parse_one_sdd_entry( sdd_data_position, SDD_TRIGGERS_LENGTH, records,
                line_number, value );
            if( entry_error ) error_number++;
            temp_spec.header_info_change().triggers = int( value );

            //  Total number of overflows detected 32-bits (2 words, low word first)
            entry_error = parse_one_sdd_entry( sdd_data_position, SDD_OVERFLOWS_LENGTH, records,
                line_number, value );
            if( entry_error ) error_number++;
            temp_spec.header_info_change().overflows = int( value );

            //  Total number of underflows detected 32-bits (2 words, low word first)
            entry_error = parse_one_sdd_entry( sdd_data_position, SDD_UNDERFLOWS_LENGTH, records,
                line_number, value );
            if( entry_error ) error_number++;
            temp_spec.header_info_change().underflows = int( value );

            //  Total number of baseline samples acquired 32-bits (2 words, low word first)
            entry_error = parse_one_sdd_entry( sdd_data_position, SDD_BASEEVENTS_LENGTH, records,
                line_number, value );
            if( entry_error ) error_number++;
            temp_spec.header_info_change().baseline_samples = int( value );

            //  Total number of preamplifier resets detected (ADC excursions below ADCMIN) 32-bits (2 words, low word first)
            entry_error = parse_one_sdd_entry( sdd_data_position, SDD_PRERESETS_LENGTH, records,
                line_number, value );
            if( entry_error ) error_number++;
            temp_spec.header_info_change().preamp_resets = int( value );

            //  Total number of ADC excursions above ADCMAX 32-bits (2 words, low word first)
            entry_error = parse_one_sdd_entry( sdd_data_position, SDD_SATURATES_LENGTH, records,
                line_number, value );
            if( entry_error ) error_number++;
            temp_spec.header_info_change().saturates = int( value );

            //  Skip reserved locations between statistics and histogram
            sdd_data_position += SDD_RESERVED_LENGTH * SDD_DATA_INCREMENT;

            //  Number of bins in the spectrum
            entry_error = parse_one_sdd_entry( sdd_data_position, SDD_MCALIMHI_LENGTH, records,
                line_number, value );
            if( entry_error ) error_number++;
            //  Simulated data from Rboert Denise (Dec. 13, 2017) has MCALIMHI as the index of the high limit
            //      of the MCA channels (and thus 4095).  But it might actually be the number of MCA Channels (and thus 4096)
            //  If so, the second histogram will get an incorrect tag word message (see above for TAGWORD1)
            int mca_limit_high = int( value ) + 1;

            //  Tag word 0x55AA precedes spectrum data
            entry_error = parse_one_sdd_entry( sdd_data_position, SDD_TAGWORD2_LENGTH, records,
                line_number, value );
            if( entry_error || int( value ) != SDD_TAGWORD2_VALUE ) {
                cout << "*** Error tag word preceding channel data is missing or has incorrect value: ";
                cout << int( value ) << " on line " << line_number << " should be " << SDD_TAGWORD2_VALUE << ". ***" << endl;
                error_number++;
                continue;
            }

            //  Check the number of channels in the histogram
            const int minimum_channels = 2;
            if( mca_limit_high < minimum_channels ) {
                cout << "*** Error - not enough channels (" << mca_limit_high << ") in histogram " << i_hist_line+1;
                cout << " on line " << line_number << ", should be at least " << minimum_channels << ". ***" << endl;
                error_number++;
            }
            //  Don't continue processing all of the bin values if there is already a problem
            //  Avoid too many error messages if wrong format file is opened
            if( error_number > MAX_ERROR_MESSAGES ) continue;
            //  Now read all of the bin values in the histogram
            vector <float> meas_histogram( mca_limit_high, 0 );
            int i_bin;
            for( i_bin=0; i_bin<mca_limit_high; i_bin++ ) {
                entry_error = parse_one_sdd_entry( sdd_data_position, SDD_BINWORD_LENGTH, records,
                    line_number, value );
                if( entry_error ) {
                    cout << "*** Error reading histogram " << i_hist_line+1 << ", channel " <<  i_bin << " on line " << line_number << ". ***" << endl;
                    error_number++;
                    if( entry_error && error_number > MAX_ERROR_MESSAGES ) break;
                }
                meas_histogram[i_bin] = float( value );

            }
            if( entry_error ) continue;
            temp_spec.meas( meas_histogram );
            //  Calculate actual live time from real time, input count rate, and output count rate
            float icr = temp_spec.header_info().triggers / temp_spec.header_info().live_time_DSPC;
            float ocr = ( temp_spec.header_info().events + temp_spec.header_info().overflows + temp_spec.header_info().underflows ) / temp_spec.real_time();
            temp_spec.live_time( temp_spec.real_time() * ocr/ icr );
            spec_vec_out.push_back( temp_spec );
        }   //  for( i_hist_line=0; i_hist_line<2;i_hist_line++ )
    }   //  while ( ! ( !sdd_data_file ) )

    if( error_number > 0 ) return -3;
    return 0;

};

//  Helper function to reduce code size and complexity
bool parse_one_sdd_entry( int &index, const int length, const vector <string> &records,
            const int line_number, double &value_out ) {
    if( records.size() <= index + (length-1) * SDD_DATA_INCREMENT ) {
        cout << "*** Error - unexpected end of line while reading line " << line_number << ". ***" << endl;
        return true;
    }
    value_out = 0;
    int word_count;
    bool error_found = false;
    for( word_count=0; word_count<length; word_count++ ) {
        istringstream record_stream( records[index] );
        int value16 = 0;
        record_stream >> value16;
        if ( ! record_stream || value16 < 0 || value16 >= SDD_DATA_SHIFT16 ) {
            cout << "Missing or invalid value on line " << line_number << ", entry number " << index << ", " << records[index] << endl;
            error_found = true;
        } else {
            value_out += value16 * pow( SDD_DATA_SHIFT16, word_count );
        }
        index += SDD_DATA_INCREMENT;
    }
    return error_found;
};
