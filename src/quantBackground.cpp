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

#include "quantBackground.h"
#include "quantComponents.h"
#include "XRFconstants.h"
#include "XRFcontrols.h"
#include "snip.h"
#include "split_component.h"
#include "scale_under_peaks.h"


using namespace std;

//		calculate the background in an X-ray spectrum
//      using the SNIP nonlinear digital filter

//  Written Feb. 20, 2017
//      From code in PIQUANT Version 1 main programs
//  Modified May 14, 2017
//      Add use of option arguments to control where background removal starts, width, & # of iterations
//              (-b,start_ch,width_ch,iterations)
//      Use an average of SNIP with and without LSQ (one is too large and one is too small, try average)
//  Modified July 27, 2018
//      Change background to 2-zone SNIP
//  Modified Aug. 1, 2018
//      Change end channel checks to avoid accepting small non-zero channel number when source kV is zero
//  Modified July 2, 2019
//      Allow for background to be adjusted manually instead of in least-squares fit, using -b option parameter
//      Ignore background parameters with value <= 0
//  Modified July 3, 2019
//      Background is now always put into a component, (but is excluded and component not used if bkg not fit)
//      This function now makes all of the background decisions based on the option parameters
//  Modified Oct. 21, 2020
//      Add multiple background components (split for independent fitting)
//  Modified Dec. 4, 2020   Use SNIP background to fit anomalous background at low energies
//  Modified Dec. 28, 2020  Modify SNIP parameters to get better low energy background, single-region, no fit
//  Modified Jan. 5, 2021   Use BKG_SNIP_CROSSOVER to set number of bkg components if calculated bkg is used
//  Modified Feb. 1, 2021   Revert to all SNIP bkg, with fit at low energies and not at high energies (crossover in arguments)
//  Modified Apr. 6, 2021   Add CONTINUUM component type and sort out how to handle background and Compton escape (remove Det shelf component)
//  Modified Apr. 11, 2021  Fix quantBckground to handle plot cmd properly (just use SNIP, don't insert components and update calc)
//  Modified Apr. 28, 2021  Remove option to have -b,0 perform SNIP bkg with default parameters (must now specify parameters to use SNIP)
//                          Add back crossover, but using SNIP at high energies and adj calc bkg at low energies (via option -a)
//                          Change default to this new crossover, since it works much better for trace elements
//  Modified May 10, 2021   Add -bh and -bx background options, eliminate -a option, implement more controls and simplify code using helper functions
//                          Use scale_under_peaks function when scale factor is negative
//  Modified June 10, 2021  Add argument to use single SNIP bkg for plots and not use update_calc to get bkg in spectrum (to avoid residual, etc.)
//                          Include final decision on background defaults

//  Helper functions to consolidate common code for low and high energy backgrounds (code at end)
void perform_SNIP( const vector <float> &bkg_params, const  XraySpectrum &spectrum, vector <float> &bkg_out, const int end_chan_in = 0 );
void set_up_params( const vector <float> &default_params, const vector <float> bk_args, const XRFconditions &conditions_in,
                        const XraySpectrum &spectrum, vector <float> &bk_params_out );

int quantBackground( const XRFconditions &conditions_in, XraySpectrum &spectrum, bool plot ) {

    //const vector <float> bkg_SNIP_defaults = { 0, 10, 40, 910, 2800, 14 };
    const vector <float> bkg_SNIP_defaults = { 0, 12, 60, 910, 2800, 16 };  //  from Chris Heirwegh, April 28, 2021
//    vector <float> bkg_defaults = { -1,-5 };  //  Version 3.1.2 default bkg, May 10, 2021
//    vector <float> bh_defaults = { 0, 12, 60, 910, 2800, 16, 1 };  //  Version 3.1.2 default bkg, May 10, 2021
//    vector <float> bx_defaults = { 7000, 500 };  //  Version 3.1.2 default bkg, May 10, 2021
//    vector <float> bx_def_for_bh = { 7000, 500 };  //  In case bx defaults above are changed and -bh only is specified

    //  Final decision for surface operations   June 10, 2021   by Tim Elam and Chris Heirwegh, from PIXL Elemental Calibration data set
    vector <float> bkg_defaults = { -1,-5 };
    vector <float> bh_defaults = { 0,10,60,910,1260,6,1 };
    vector <float> bx_defaults = { 7150,150 };
    vector <float> bx_def_for_bh = { 7150,150 };  //  In case bx defaults above are changed and -bh only is specified

    //  For compatibility with quantCalculate and older PIQUANT versions, use single bkg if only -b option is present
    bool single_bkg = false;    //   One bkg for the entire spectrum instead of crossover background
    vector <float> bkg_single_SNIP( spectrum.numberOfChannels(), 0 );

    //  For plots only, use a single SNIP background and don't use update_calc to get background into spectrum
    if( plot ) {
        bkg_defaults = bkg_SNIP_defaults;
        bh_defaults.clear();
        bx_defaults.clear();
        bx_def_for_bh.clear();
        single_bkg = true;
    }
    //  Set up the background subtraction using the option parameters (if any)

    //  Get any background arguments stored with the spectrum
    vector <float> bkg_args;
    spectrum.get_bkg_parameters( bkg_args );
    vector <float> bh_args;
    spectrum.get_bh_parameters( bh_args );
    vector <float> bx_args;
    spectrum.get_bx_parameters( bx_args );
    //  Adjust the end channel to the present spectrum
    const int endCh = spectrum.channel( conditions_in.source.kV() * 1000 );
    if( bkg_args.size() > 0 && bh_args.size() == 0 && bx_args.size() == 0 ) single_bkg = true;
    vector <SpectrumComponent> bkg_components;
    vector <XrayLines> empty_lines;

    //  Set up parameters for low-energy background (or full-spectrum background if no crossover)
    vector <float> bkg_params(7);
    set_up_params( bkg_defaults, bkg_args, conditions_in, spectrum, bkg_params );
    //float bkg_multiplier = 1;
    if( bkg_params[0] < 0 ) {     //  Calculate a continuum background and use it
        //  Add a CONTINUUM component
        int result = makeComponents( CONTINUUM, empty_lines, bkg_components );
        if( result < 0 ) {
            return -100 + result;
        }
        unsigned int index = bkg_components.size() - 1;
        bkg_components[index].spectrum.resize( spectrum.numberOfChannels() );
        //  Set up scaling for calculated background
        if( bkg_params[1] > 0 ) {   //  Fixed amplitude scaling via option
            bkg_components[index].fit = false;
            bkg_components[index].coefficient = bkg_params[1];
        } else if( bkg_params[1] == 0 ) {   //  Amplitude scaled via least-squares fit to spectrum
            bkg_components[index].fit = true;
        } else if( bkg_params[1] < 0 ) {   //  Amplitude scaled via scale-under-peaks algorithm (done after calculation, in quantCalculate)
            bkg_components[index].fit = false;
            bkg_components[index].scale_under = -bkg_params[1];  //  Sigma multiplier for scale-under-peaks algorithm in quantCalculate
        }
    } else {    //  Use SNIP background
        //  Add a SNIP component
        int result = makeComponents( SNIP_BKG, empty_lines, bkg_components );
        if( result < 0 ) {
            return -100 + result;
        }
        unsigned int index = bkg_components.size() - 1;
        //  Compute SNIP background
        perform_SNIP( bkg_params, spectrum, bkg_single_SNIP, endCh );
        bkg_components[index].spectrum = bkg_single_SNIP;
        //  Set up scaling for SNIP background
        if( bkg_params[6] > 0 ) {   //  Fixed amplitude scaling via option
            bkg_components[index].fit = false;
            bkg_components[index].coefficient = bkg_params[6];
        } else if( bkg_params[6] == 0 ) {   //  Amplitude scaled via least-squares fit to spectrum
            bkg_components[index].fit = true;
        } else if( bkg_params[6] < 0 ) {   //  Amplitude scaled via scale-under-peaks algorithm (done after calculation, in quantCalculate)
            bkg_components[index].fit = false;
            bkg_components[index].coefficient = scale_under_peaks( bkg_components[index].spectrum, spectrum.meas(), spectrum.sigma(), fabs(bkg_params[6]) );
        }
    }

    //  Process the -bx arguments and crossover information
    vector <float> bx_params(2);
    set_up_params( bx_defaults, bx_args, conditions_in, spectrum, bx_params );
    //  Force crossover if -bh option but no -bx option and crossover is not default
    if( !single_bkg && bh_args.size() > 0 && bx_args.size() == 0 && bx_defaults.size() == 0 ) bx_params = bx_def_for_bh;

    //  Set up parameters for high-energy background
    vector <float> bh_params(7);
    if( !single_bkg ) {
        set_up_params( bh_defaults, bh_args, conditions_in, spectrum, bh_params );
        //bkg_multiplier = 1;
        if( bh_params[0] < 0 ) {     //  Calculate a continuum background and use it
            //  Add a CONTINUUM component
            int result = makeComponents( CONTINUUM, empty_lines, bkg_components );
            if( result < 0 ) {
                return -100 + result;
            }
            unsigned int index = bkg_components.size() - 1;
            bkg_components[index].spectrum.resize( spectrum.numberOfChannels() );
            //  Set up scaling for calculated background
            if( bh_params[1] > 0 ) {   //  Fixed amplitude scaling via option
                bkg_components[index].fit = false;
                bkg_components[index].coefficient = bh_params[1];
            } else if( bh_params[1] == 0 ) {   //  Amplitude scaled via least-squares fit to spectrum
                bkg_components[index].fit = true;
            } else if( bh_params[1] < 0 ) {   //  Amplitude scaled via scale-under-peaks algorithm (done after calculation, in quantCalculate)
                bkg_components[index].fit = false;
                bkg_components[index].scale_under = -bh_params[1];  //  Sigma multiplier for scale-under-peaks algorithm in quantCalculate
            }
        } else {    //  Use SNIP background
            //  Add a SNIP component
            int result = makeComponents( SNIP_BKG, empty_lines, bkg_components );
            if( result < 0 ) {
                return -100 + result;
            }
            unsigned int index = bkg_components.size() - 1;
            //  Compute SNIP background
            perform_SNIP( bh_params, spectrum, bkg_components[index].spectrum, endCh );
            //  Set up scaling for SNIP background
            if( bh_params[6] > 0 ) {   //  Fixed amplitude scaling via option
                bkg_components[index].fit = false;
                bkg_components[index].coefficient = bh_params[6];
            } else if( bh_params[6] == 0 ) {   //  Amplitude scaled via least-squares fit to spectrum
                bkg_components[index].fit = true;
            } else if( bh_params[6] < 0 ) {   //  Amplitude scaled via scale-under-peaks algorithm (done after calculation, in quantCalculate)
                bkg_components[index].fit = false;
                bkg_components[index].coefficient = scale_under_peaks( bkg_components[index].spectrum, spectrum.meas(), spectrum.sigma(), fabs(bh_params[6]) );
            }
        }
    }

    //  Set up the crossover components
    if( !single_bkg && bx_params[0] > 0 ) {
        //  Crossover background, one component for low energies and one for high energies
        vector <float> bkg_split_energies = { bx_params[0] - bx_params[1], bx_params[0] + bx_params[1] };
        spectrum.put_bkg_split( bkg_split_energies );
        unsigned int i;
        for( i=0; i<bkg_components.size(); i++ )  {
            bkg_components[i].bkg_index = i;
            unsigned int k;
            for( k=0; k<bkg_components[i].spectrum.size(); k++ ) {
                float e = spectrum.energy( k );
                float split = split_weight( e, bkg_split_energies, bkg_components[i].bkg_index );
                //  Split up background
                bkg_components[i].spectrum[k] = bkg_components[i].spectrum[k] * split;
            }
        }
    }

    //  Add the components into the spectrum
    unsigned int i;
    for( i=0; i<bkg_components.size(); i++ )  {
        bkg_components[i].enabled = true;
        if( bx_args.size() > 0 ) bkg_components[i].plot = true; //  Plot all components even though it looks crazy with splits
        else bkg_components[i].plot = false;    //  Only plot overall background, not individual components
        spectrum.add_component( bkg_components[i] );
        //        cout << "quantBkg        " << i << "  " << componentDescription( bkg_components[i] );
        //        cout << "      fit " << bkg_components[i].fit;
        //        cout << "      plot " << bkg_components[i].plot;
        //        cout << "      n " << spectrum.numberOfComponents();
        //        cout << "        " << i << "  " << SpectrumComponent_toString( bkg_components[i] );
        //        cout << endl;
    }

    //  Put the full background into the spectrum
    if( plot ) spectrum.bkg( bkg_single_SNIP );
    else spectrum.update_calc();
	return 0;

};


//  Helper function to consolidate common code for low and high energy backgrounds (code at end)

void set_up_params( const vector <float> &default_params, const vector <float> bk_args, const XRFconditions &conditions_in,
                        const XraySpectrum &spectrum, vector <float> &bk_params_out ) {
    //  Set up defaults
    unsigned int ib;
    unsigned int nb = default_params.size();
    if( bk_params_out.size() < nb ) nb = bk_params_out.size();
    if( nb > 0 ) for( ib=0; ib<nb; ib++ ) bk_params_out[ib] = default_params[ib];
    //  Replace the default parameters with any arguments present
    nb = bk_args.size();
    if( bk_params_out.size() < nb ) nb = bk_params_out.size();
    if( nb > 0 ) for( ib=0; ib<nb; ib++ ) bk_params_out[ib] = bk_args[ib];
    //  Handle the situation where option is -b,0,s where s is the scaling parameter => default SNIP using s
    if( bk_args.size() == 2 && bk_args[0] == 0 && bk_params_out.size() > 6 ) {
        bk_params_out[1] = default_params[1];
        bk_params_out[6] = bk_args[1];
    }
    //  If bkg is SNIP, adjust the parameters to the present spectrum if any are zero
    //  Start at channel corresponding to minimum energy
    if( bk_params_out.size() > 3 && bk_params_out[0] == 0 && spectrum.calibration().good() ) bk_params_out[0] = spectrum.channel( conditions_in.eMin );
    if( bk_params_out.size() > 3 && bk_params_out[0] >= 0 && bk_params_out[1] == 0 )
        //  Get filter width from detector (possibly using default value for resolution)
        bk_params_out[1] = conditions_in.detector.resolution() / spectrum.calibration().energyPerChannel() + 1;
    return;
};


void perform_SNIP( const vector <float> &bkg_params, const  XraySpectrum &spectrum, vector <float> &bkg_out, const int end_chan_in ) {
    //  Uses parameters to perform SNIP background calculation, with checks and defaults for zero parameters
    //  If an input argument was given in the options, use it preferentially
    int startCh = 0;
    if( bkg_params.size() > 0 && bkg_params[0] > 0 ) {
        startCh = bkg_params[0];
    }
    //  Find first non-zero channel in spectrum (avoid initial zero channels)
    if( startCh <= 0 ) {
        int i;
        for( i=2; i<spectrum.numberOfChannels(); i++ ) {
            if( spectrum.meas()[i] > 0 ) {
                startCh = i + 2;
                break;
            };
        };
    };
    //  Avoid out-of-range errors
    if( startCh < 1 ) startCh = 1;
    if( startCh > spectrum.numberOfChannels()-1 ) startCh = spectrum.numberOfChannels()-1;
    int endCh = end_chan_in;
    //  Avoid accepting small non-zero channel number when source kV is zero
    if( endCh <= startCh ) endCh = spectrum.numberOfChannels()-10;    //  Avoid possible extra information in last few channels
    if( endCh <= startCh ) endCh = startCh + 10;
    if( endCh > spectrum.numberOfChannels()-1 ) endCh = spectrum.numberOfChannels()-1;
    //  Get filter width
    int width_chan = 0;
    if( bkg_params.size() > 1 && bkg_params[1] > 0 ) width_chan = bkg_params[1];
    if( width_chan <= 0 ) width_chan = 12;    //  Last resort for uncalibrated spectrum, about 125 eV at 10 eV/channel
    int iterations = 0;
    //  Use value from option arguments if any
    if( bkg_params.size() > 2 && bkg_params[2] > 0 ) iterations = int( bkg_params[2] );
    if( iterations <= 0 ) iterations = 24;  //  Default value, used for almost everything
    //  Change background calculation to new 2-zone SNIP routine developed by Lauren O'Neil
    //  This means added a few more parameters to the arguments list to control new zone
    //  It reverts to the standard SNIP if any of the 2nd zone parameters are zero or missing
    int startCh2 = 0;
    int endCh2 = 0;
    int width2 = 0;
    if( bkg_params.size() > 5 ) {
        if( bkg_params[3] > 0 ) startCh2 = bkg_params[3];
        if( bkg_params[4] > 0 ) endCh2 = bkg_params[4];
        if( bkg_params[5] > 0 ) width2 = bkg_params[5];
    }
    bkg_out.resize( spectrum.numberOfChannels() );
    snipbg_2zone( spectrum.meas(), bkg_out, startCh, endCh, width_chan, iterations, startCh2, endCh2, width2 );
//    float sum = 0;
//    int i;
//    for( i=0; i<bkg1.size(); i++ ) sum += bkg_out[i];
//    cout << "snipbg_2zone  " << startCh << "  " << endCh << "  " << width_chan << "  " << iterations << "  " << startCh2 << "  " << endCh2 << "  " << width2 << "  " << sum << "   m " << bkg_multiplier << endl;
};
