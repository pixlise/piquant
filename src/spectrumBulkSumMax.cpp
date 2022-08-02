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
#include "spectrumBulkSumMax.h"
#include "setup_spectrum_parameters.h"
#include "read_spectrum_file.h"
#include "quantCombineSpectra.h"


int spectrumBulkSumMax(const std::string &map_spec_file,

    const XRFconditionsInput &condStruct_config,

    const ARGUMENT_LIST &arguments,
    const bool oxidesOutput,

    const XraySpectrum &configSpectrum,
    int n_map_spectra,

    int sequence_number,
    vector <float> &bulk_sum,
    vector <float> &max_value,
    float &sum_live_time,
    float &sum_geometry,
    int &geometry_count,

    XraySpectrum &singleSpectrum,
    bool &error
)
{
    // Init a place to write outputs to
    std::ostringstream termOutFile;

    // Spectrum we've read
    vector <XraySpectrum> spectrum_vec;

    int result = 0;
    error = false;

    XRFconditionsInput condStruct_map;
    result = read_spectrum_file( termOutFile, map_spec_file,
            spectrum_vec, condStruct_map );
    if ( result != 0 ) {
        termOutFile << "read_spectrum_file failed, result = " << result << "   file " << endl;
        error = true;
        return -1;
    };
    if( condStruct_map.conditionsVector[GEOMETRY_INDEX] != 0 ) {
        sum_geometry += condStruct_map.conditionsVector[GEOMETRY_INDEX];
        geometry_count++;
    }
    //  Set up energy calibration, background parameters, and measurement conditions
    setup_spectrum_parameters( arguments, configSpectrum.calibration(), spectrum_vec,
            condStruct_config, condStruct_map, cout );
    if( spectrum_vec.size() <= 0 ) {
        termOutFile << "No spectra in file " << map_spec_file << endl;
        error = true;
        return -1;
    } else {
        //  Combine the spectrum information from several detectors (or the selected detector) into the variable where they will be used
        //      NB: quantCombineSpectra modifies the spectra in the input list to match them to a single energy axis
        //          for proper plotting
        result = quantCombineSpectra( spectrum_vec, singleSpectrum, arguments.detector_select );
        if( result < 0 ) {
            error = true;
            return -1;
        };
    }
    singleSpectrum.seq_number( sequence_number );

    //  Just read the spectra and calculate the sum spectrum and the maximum value spectrum
    if( n_map_spectra <= 0 ) {
        //  Initialize the sum and max channel storage
        bulk_sum.resize( singleSpectrum.meas().size() );
        max_value.resize( singleSpectrum.meas().size() );
        int is;
        for( is=0; is<bulk_sum.size(); is++ ) {
            bulk_sum[is] = 0;
            max_value[is] = 0;
        }
    } else {
        if( singleSpectrum.meas().size() != bulk_sum.size() ) {
            termOutFile << "Spectrum in file " << map_spec_file << " is not the same size as previous spectra." << endl;
            error = true;
            return -1;
        }
    }
    //  Eventually want to use quantCombineSpectra to do this and get energy alignment correct
    //  Will have to sort out live time and solid angle sums in quantCombineSpectra, as well as
    //      a way to specify that the new detector histograms get added to the existing sum
    int is;
    for( is=0; is<singleSpectrum.meas().size(); is++ ) {
        bulk_sum[is] += singleSpectrum.meas()[is];
        if( singleSpectrum.meas()[is] > max_value[is] ) max_value[is] = singleSpectrum.meas()[is];
    }
    sum_live_time += singleSpectrum.live_time();

    cout << termOutFile.str();
    return 0;
}
