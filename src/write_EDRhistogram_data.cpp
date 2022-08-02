// Copyright (c) 2018-2022 California Institute of Technology (‚ÄúCaltech‚Äù) and
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
#include "XRFconstants.h"
#include "XRFcontrols.h"
#include "write_EDRhistogram_data.h"



using namespace std;

//  Appends one line to a map file

//  Written Nov. 8, 2017    (adapted from code in quantWriteMap.cpp)


int write_EDRhistogram_data( const int histogram_count,
            const XraySpectrum &spectrum1, const XraySpectrum &spectrum2,
            const std::string  &edr_file_name ) {

    //  Write element composition information and some diagnostic information
    //  Just writes CSV files for this first version

    //  Open the map file for append
    if( edr_file_name.length() <= 0 ) return -1;
    ofstream edr_out_stream;
    edr_out_stream.open( edr_file_name.c_str(), ios::out | ios::app );
    if( !edr_out_stream ) return -2;
	edr_out_stream.setf( ios::fixed, ios::floatfield );

    int is;

    //  Clear the file and add a header if this is the first histogram
    if( histogram_count <= 0 ) {

        //  Close and re-open the file for output (not append) to overwrite if file exists
        edr_out_stream.close();
        edr_out_stream.open( edr_file_name.c_str(), ios::out );
        edr_out_stream.setf( ios::fixed, ios::floatfield );

        //  Add a title and header if this is the first histogram
        edr_out_stream << "Histogram EDR file (EM version, not final flight version)." << endl;

        //  Fake data for EM software development   for Fang Zhong   Sept. 14, 2017
        //  Format for Histogram EDR from "EDR_Organization_proposal_Aug302017.pptx"
        //      Typical histogram EDR:  All histograms from multiple DPO with one RTT hashkey; each histogram in a row as below
        //      Row 1: 					          detector 1                         detector 2
        //      Row 2: SCLK, RTT, PIXL DP Category, USN, PIXL motion counter, data1, data2, data3, ÖÖÖÖ..

        edr_out_stream << "SCLK";   //  Spacecraft clock
        edr_out_stream << ", " << "RTT";    //  Round Trip Tracking token (32-bits)
        edr_out_stream << ", " << "PDPC";   //  PIXL Data Product Category
        edr_out_stream << ", " << "USN";    //  PIXL Universal Sequence Number
        edr_out_stream << ", " << "PMC";    //  PIXL Motion Counter
        edr_out_stream << ", " << "DPPSTATUS";  //  Detector 1
        edr_out_stream << ", " << "REALTIME_1";
        edr_out_stream << ", " << "LIVETIME_1";
        edr_out_stream << ", " << "EVTSINRUN_1";
        edr_out_stream << ", " << "TRIGGERS_1";
        edr_out_stream << ", " << "OVERFLOWS_1";
        edr_out_stream << ", " << "UNDERFLOW_1";
        edr_out_stream << ", " << "BASEEVENTS_1";
        edr_out_stream << ", " << "PRERESETS_1";
        edr_out_stream << ", " << "SATURATES_1";
        edr_out_stream << ", " << "MCALIM_1"; //  # channels in histogram
        for( is=0; is<spectrum1.meas().size(); is++ ) edr_out_stream << ", " << "H1_CH" << is;
        edr_out_stream << ", " << "REALTIME_2";  //  Detector 2
        edr_out_stream << ", " << "LIVETIME_2";
        edr_out_stream << ", " << "EVTSINRUN_2";
        edr_out_stream << ", " << "TRIGGERS_2";
        edr_out_stream << ", " << "OVERFLOWS_2";
        edr_out_stream << ", " << "UNDERFLOW_2";
        edr_out_stream << ", " << "BASEEVENTS_2";
        edr_out_stream << ", " << "PRERESETS_2";
        edr_out_stream << ", " << "SATURATES_2";
        edr_out_stream << ", " << "MCALIM_2";
        for( is=0; is<spectrum2.meas().size(); is++ ) edr_out_stream << ", " << "H2_CH" << is;

        edr_out_stream << endl;

    }


    //  Write the line of information to match the header

    //  Temporary code to write fake data for EM software development
    //  CSV file with data product header info, DPP register info, and histogram channels
    //  Intended to be a close approximation to a PIXL EDR

    //  PILX DP header info (see above for EDR format example)
    //  Get the local clock time and convert to ASCII (human-readable ISO format) to put in the header
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];
    time (&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buffer,80,"%H:%M:%S",timeinfo);
    std::string time_str(buffer);
    int pixl_dp_category = 18;  //  Histogram Nominal
    const float ns_per_sec = 2000000;

    edr_out_stream << time_str;
    edr_out_stream << ", " << "0x" << hex << spectrum1.aux_info().rtt << dec;
    edr_out_stream << ", " << pixl_dp_category;
    edr_out_stream << ", " << spectrum1.aux_info().usn;
    edr_out_stream << ", " << spectrum1.aux_info().pmc;

    //  First detector (use actual measured spectrum)
    edr_out_stream << ", " << "0x0007";  //  DPPSTATUS
    edr_out_stream << ", " << spectrum1.real_time() * ns_per_sec;       //  REALTIME 48-bit Real Time in 500 ns units
    edr_out_stream << ", " << spectrum1.live_time() * ns_per_sec;       //  LIVETIME    48-bit Live Time in 500 ns units
    edr_out_stream << ", " << spectrum1.header_info().events;                 //  EVTSINRUN
    edr_out_stream << ", " << spectrum1.header_info().triggers;                 //  TRIGGERS
    edr_out_stream << ", " << spectrum1.header_info().overflows;                    //  OVERFLOWS
    edr_out_stream << ", " << spectrum1.header_info().underflows;                   //  UNDERFLOW
    edr_out_stream << ", " << spectrum1.header_info().baseline_samples;             //  BASEEVENTS
    edr_out_stream << ", " << spectrum1.header_info().preamp_resets;                //  PRERESETS
    edr_out_stream << ", " << spectrum1.header_info().saturates;                    //  SATURATES
    edr_out_stream << ", " << spectrum1.meas().size();                  //  MCALIMHI (number of channels in histogram)
    for( is=0; is<spectrum1.meas().size(); is++ ) edr_out_stream << ", " << spectrum1.meas()[is];    //  HISTOGRAM


    //  Second detector histogram
    edr_out_stream << ", " << spectrum2.live_time() * ns_per_sec;       //  REALTIME 48-bit Real Time in 500 ns units
    edr_out_stream << ", " << spectrum2.live_time() * ns_per_sec;       //  LIVETIME    48-bit Live Time in 500 ns units
    edr_out_stream << ", " << spectrum2.header_info().events;                 //  EVTSINRUN
    edr_out_stream << ", " << spectrum2.header_info().triggers;                 //  TRIGGERS
    edr_out_stream << ", " << spectrum2.header_info().overflows;                    //  OVERFLOWS
    edr_out_stream << ", " << spectrum2.header_info().underflows;                   //  UNDERFLOW
    edr_out_stream << ", " << spectrum2.header_info().baseline_samples;             //  BASEEVENTS
    edr_out_stream << ", " << spectrum2.header_info().preamp_resets;                //  PRERESETS
    edr_out_stream << ", " << spectrum2.header_info().saturates;                    //  SATURATES
    edr_out_stream << ", " << spectrum2.meas().size();                  //  MCALIMHI (number of channels in histogram)
    for( is=0; is<spectrum2.meas().size(); is++ ) edr_out_stream << ", " << spectrum2.meas()[is];    //  HISTOGRAM

    //  End of line and close map file
    edr_out_stream << endl;
    edr_out_stream.close();

    return 0;

};

