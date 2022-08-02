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

#include "setup_spectrum_parameters.h"

//  Modified May 10, 2021    Add parameters for -bh and -bx background options, eliminate adj_calc_bkg (from -a option)
//                          Move setup of energy calibration, detector resolution, and Compton convolution here

//  Consolidate code to copy input conditions structure
void copy_conditions_struct( const XRFconditionsInput &condStruct_in, XRFconditionsInput &condStruct_out ) {
	condStruct_out.conditionsVector.resize( XRF_PARAMETER_LAST, 0 );
     //  Use configuration file conditions for any missing in the spectrum file
    int icv;
    for( icv=0; icv<condStruct_in.conditionsVector.size(); icv++ )
            if( condStruct_out.conditionsVector[icv] == 0 )
                condStruct_out.conditionsVector[icv] = condStruct_in.conditionsVector[icv];
    //  Transfer optic file name and X-ray tube external input file name if not already set
    if( condStruct_out.optic_file_name.length() <= 0 ) condStruct_out.optic_file_name = condStruct_in.optic_file_name;
    if( condStruct_out.tube_file_name.length() <= 0 ) condStruct_out.tube_file_name = condStruct_in.tube_file_name;
    return;
};

//  Consolidate code to choose calibration and other spectrum parameters
void setup_spectrum_parameters( const ARGUMENT_LIST &arguments, const XrayEnergyCal &config_cal,
        vector <XraySpectrum> & spectrum_vec_out,
        const XRFconditionsInput &condStruct_config, XRFconditionsInput &condStruct_out, std::ostream &logout ) {

    if(condStruct_out.conditionsVector.size() < XRF_PARAMETER_LAST) {
        cout << "setup_spectrum_parameters: conditionsVector size " << condStruct_out.conditionsVector.size() << " less than expected " << XRF_PARAMETER_LAST << endl;
        return;
    }

    int iv;
    for( iv=0; iv<spectrum_vec_out.size(); iv++ ) {
        //     Setup adjustments to energy calibration and detector resolution in fits (-f and -g options)
        spectrum_vec_out[iv].adjust_energy( arguments.fit_adjust_energy );
        spectrum_vec_out[iv].adjust_width( arguments.fit_adjust_width );
        //  Setup convolution of Compton components (now brute force, so very expensive in compute time)
        spectrum_vec_out[iv].convolve_Compton( arguments.convolve_Compton );
        //  Put in background control parameters from argument list (will be zero size if none)
        spectrum_vec_out[iv].put_bkg_parameters( arguments.bkg_args );
        if( arguments.bkg_args.size() > 0 ) {
            //  Write a message to document the background parameters
            logout << "Background arguments:";
            int ib;
            for( ib=0; ib<arguments.bkg_args.size(); ib++ ) logout << "  " << arguments.bkg_args[ib];
            logout << endl;
        }
        vector <float> bkg_test;
        spectrum_vec_out[iv].get_bkg_parameters( bkg_test );
        //  Put in background control parameters from argument list for high-energy background
        spectrum_vec_out[iv].put_bh_parameters( arguments.bh_args );
        if( arguments.bh_args.size() > 0 ) {
            //  Write a message to document the background parameters
            logout << "High-energy background arguments:";
            int ib;
            for( ib=0; ib<arguments.bh_args.size(); ib++ ) logout << "  " << arguments.bh_args[ib];
            logout << endl;
        }
        //  Put in background control parameters from argument list for crossover background
        spectrum_vec_out[iv].put_bx_parameters( arguments.bx_args );
        if( arguments.bx_args.size() > 0 ) {
            //  Write a message to document the background parameters
            logout << "Background crossover arguments:";
            int ib;
            logout.precision(2);
            for( ib=0; ib<arguments.bx_args.size(); ib++ ) logout << "  " << arguments.bx_args[ib];
            logout << endl;
        }
        // If an energy calibration was given in the argument list, it overrides
        if( arguments.eV_ch > 0 ) {
            spectrum_vec_out[iv].calibration( arguments.eV_start, arguments.eV_ch );
            logout << "Using energy calibration from option argument  ";
            logout.precision(1);
            logout << "  eV start = " << spectrum_vec_out[iv].calibration().energyStart();
            logout.precision(4);
            logout << "  eV/ch = " << spectrum_vec_out[iv].calibration().energyPerChannel();
            logout << endl;
        }
        //  If a good energy calibration is not present, try the one from configuration file
        if( ! spectrum_vec_out[iv].calibration().good() ) {
            spectrum_vec_out[iv].calibration( config_cal );
            logout << "Using energy calibration from configuration file  ";
            logout.precision(1);
            logout << "  eV start = " << spectrum_vec_out[iv].calibration().energyStart();
            logout.precision(4);
            logout << "  eV/ch = " << spectrum_vec_out[iv].calibration().energyPerChannel();
            logout << endl;
        }
        //  Use configuration file conditions (if any) for any missing in the spectrum file
        if( condStruct_config.conditionsVector.size() >= condStruct_out.conditionsVector.size() ) {
            int icv;
            for( icv=0; icv<condStruct_out.conditionsVector.size(); icv++ )
                    if( condStruct_out.conditionsVector[icv] == 0 )
                        condStruct_out.conditionsVector[icv] = condStruct_config.conditionsVector[icv];
        }
        //  Implement linear energy calibration correction from Chris
        if( condStruct_out.conditionsVector[ENERGY_CORRECTION_OFFSET_INDEX] != 0 || condStruct_out.conditionsVector[ENERGY_CORRECTION_SLOPE_INDEX] != 0 ) {
            spectrum_vec_out[iv].calibration_change().linearCorrection( condStruct_out.conditionsVector[ENERGY_CORRECTION_OFFSET_INDEX], condStruct_out.conditionsVector[ENERGY_CORRECTION_SLOPE_INDEX] );
            logout << "Applying linear energy correction at low energy: ";
            logout << "   offset " << spectrum_vec_out[iv].calibration().linearCorrectionOffset();
            logout << "   slope " << spectrum_vec_out[iv].calibration().linearCorrectionSlope();
            logout << "  (eV per keV, stops at ";
            logout << ( - spectrum_vec_out[iv].calibration().linearCorrectionOffset() / spectrum_vec_out[iv].calibration().linearCorrectionSlope() ) * 1000;
            logout << " ).";
            logout << endl;
        }
    }   //  Loop over spectrum_vec_out
    if( condStruct_out.optic_file_name.length() <= 0 ) condStruct_out.optic_file_name = condStruct_config.optic_file_name;
    if( condStruct_out.tube_file_name.length() <= 0 ) condStruct_out.tube_file_name = condStruct_config.tube_file_name;
    if( condStruct_out.anode_element_list.length() <= 0 ) condStruct_out.anode_element_list = condStruct_config.anode_element_list;

    return;
}
