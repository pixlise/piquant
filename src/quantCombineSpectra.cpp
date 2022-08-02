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

#include "quantCombineSpectra.h"
#include "rebin.h"
#include "XRFconstants.h"
#include "XRFcontrols.h"

//  Started Sept. 30, 2017
//      Combine two (or more) detectors (simple channel-by-channel sum in this version)
//      Also implements detector selection in argument list
//  Modified Dec. 15, 2017
//      Add rebin of spectra before combining using individual energy calibrations
//      Put rebinned spectra into list so they all plot correctly against a single energy axis
//  Modified Jan. 15, 2018
//      Fix logic error if number of channels in list spectrum is not equal to number in sum
//  Modified Jan. 26, 2018
//      Use average live time when combining multiple-detector spectra
//      Added total counts for summed detectors to terminal output

//      NB: This function modifies the spectra in the input list to match them to a single energy axis
//          for proper plotting (so that peak alignment can be checked visually)

using namespace std;

int quantCombineSpectra( std::vector <XraySpectrum> &spectrum_list_in,
                XraySpectrum &combinedSpectrum_out, const int detector_selection ) {

    if( spectrum_list_in.size() <= 0 ) {
        return -1;
    }
    int det_sel_verified = detector_selection;
    if( det_sel_verified >= 0 && det_sel_verified >= spectrum_list_in.size() ) {
        cout << "*** Error - invalid detector selection: " << det_sel_verified;
        cout << " (only " << spectrum_list_in.size() << " detectors found)" << endl;
        return -2;
    }

    if( det_sel_verified >= 0 ) {
        combinedSpectrum_out = spectrum_list_in[det_sel_verified];
        cout << "Detector " << det_sel_verified << " selected." << endl;
        return 0;
    }

    if( det_sel_verified < 0 && spectrum_list_in.size() == 1 ) {
        combinedSpectrum_out = spectrum_list_in[0];
        return 0;
    }

    //  Choose a spectrum in the list for the basis of the combined spectrum
    int basis_spec_index = -1;
    int isv;
    for( isv=0; isv<spectrum_list_in.size(); isv++ ) {
        if( ! spectrum_list_in[isv].calibration().good() ) continue;
        if( spectrum_list_in[isv].numberOfChannels() < 2 ) continue;
        basis_spec_index = isv;
        break;
    }
    if( basis_spec_index < 0 ) {
        //  In case we are only plotting, put something in output spectrum to provide non-zero # channels
        combinedSpectrum_out = spectrum_list_in[0];
        cout << "*** Could not combine spectra, all spectra in list are missing energy calibration or do not have enough channels. ***" << endl;
        return -3;
    }

    //  Initialize the combined spectrum
    combinedSpectrum_out = spectrum_list_in[basis_spec_index];
    int ns = combinedSpectrum_out.numberOfChannels();
    //  ASet up the vector to hold the summed spectrum
    vector <float> new_spec( ns );
    //  Set up vector of energy bins for combined spectrum
    vector <float> new_energy( ns );
    int is;
    for( is=0; is<ns; is++ ) {
        new_spec[is] = combinedSpectrum_out.meas()[is];
        new_energy[is] = combinedSpectrum_out.energy( is );
    }
    //  Add up the live times and real times
    float live_time_sum = combinedSpectrum_out.live_time();
    float real_time_sum = combinedSpectrum_out.real_time();

    //  Loop over all spectra in the list(except the basis spectrum) and add them into the output spectrum
    const int list_errors_return = -4;
    int list_errors = 0;
    for( isv=0; isv<spectrum_list_in.size(); isv++ ) {
        if( isv == basis_spec_index) continue;  //  Skip the basis spectrum
        if( ! spectrum_list_in[isv].calibration().good() ) {
                cout << " ** Spectrum #" << isv << " could not be summed, it does not have an energy calibration." << endl;
                list_errors++;
                if( list_errors > MAX_ERROR_MESSAGES ) return list_errors_return;
        }
        if( spectrum_list_in[isv].meas().size() < ns / 10 + 2 ) {   //  10x expansion is too much!
                cout << " ** Spectrum #" << isv << " could not be summed, it does not have enough channels." << endl;
                list_errors++;
                if( list_errors > MAX_ERROR_MESSAGES ) return list_errors_return;
        }
        //  Vector to hold this spectrum after rebin
        vector <float> new_list_spec( ns );
        if( ! ( spectrum_list_in[isv].calibration() == combinedSpectrum_out.calibration() )
                    || spectrum_list_in[isv].meas().size() != ns ) {
            //  Re-bin the spectrum using its energy calibration to match the combined energy calibration
            //  Set up vectors for energy bins of spectrum to be added in to combination
            vector <float> energy_list( spectrum_list_in[isv].meas().size() );
            for( is=0; is<ns; is++ ) energy_list[is] = spectrum_list_in[isv].energy( is );
            int result = rebin( energy_list, spectrum_list_in[isv].meas(), new_energy, new_list_spec );
            if( result < 0 ) {  //  This should never happen after above checks
                cout << " ** Spectrum #" << isv << " could not be summed, it could not be changed to a common energy scale." << endl;
                list_errors++;
                if( list_errors > MAX_ERROR_MESSAGES ) return list_errors_return;
            }
            //  Replace the energy calibration and measured spectrum in the list object for proper plotting
            spectrum_list_in[isv].calibration( combinedSpectrum_out.calibration() );
            spectrum_list_in[isv].meas( new_list_spec );
        } else {
            for( is=0; is<ns; is++ ) new_list_spec[is] = spectrum_list_in[isv].meas()[is];
        }
        //  Add the spectra in the input list together
        for( is=0; is<ns; is++ ) new_spec[is] += new_list_spec[is];
        live_time_sum += spectrum_list_in[isv].live_time();
        real_time_sum += spectrum_list_in[isv].real_time();
    }

    if( list_errors > 0 ) {
        return list_errors_return;
    }

    combinedSpectrum_out.meas( new_spec );
    combinedSpectrum_out.live_time( live_time_sum );
    combinedSpectrum_out.real_time( real_time_sum );

    if( spectrum_list_in.size() > 2 ) cout << "All " << spectrum_list_in.size();
    else if( spectrum_list_in.size() == 2 ) cout << "Both";
    if( spectrum_list_in.size() > 1 ) {
        cout << " detectors summed (after matching channels using individual energy calibrations)";
        cout << ", total counts = " << combinedSpectrum_out.total_counts();
        cout << endl;
    }

    return 0;

};
