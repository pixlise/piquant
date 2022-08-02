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

/*
 *  XraySpectrum.cpp
 *  APS1IDanalysis_SpecProcDev
 *
 *  Created by Tim Elam on 3/16/12.
 *  Copyright 2012 University of Washington. All rights reserved.
 *
 */

#include "XraySpectrum.h"
#include <math.h>
#include <sstream>
#include "toStringHelpers.h"

// re-written Feb. 7, 2017 to separate energy calibration and add components
//  Modified June 7, 2017
//      Fix error in calculation of energy per channel and offset when tilt corrections were non-zero
//      Allow background to be a fitted component (change residual to meas vs calc, was net vs calc)
//  Modified Sept. 13, 2017
//      update spectrum background using background components including coefficients
//          when components or coefficients are updated (in update_coefficients & update_component)
//  Modified Sept. 13, 2017
//      Added check for nan to components included in fit_vector function
//  Modified Sept. 30, 2017
    //  Check for very small components compared to largest and don't include in fit (to improve stability)
//  Modified July 31, 2018
//      Implement linear energy calibration correction from Chris
//  Modified Sep. 19, 2018
//      Fix bug in linear energy calibration correction (x-intercept was compared to offset, not energy_in)
//  Modified May 16, 2019
//      Add counts in 1-7.25 keV region to quant map outputs, move from quantWriteResults to XraySpectrum
//      Also get total counts and counts in 1-7.25 keV region from calculation if measured is zero
//      Allow region start and end energies to be modified
//  Modified July 2, 2019
//      Allow for background to be adjusted manually (no background component present)
//  Modified July 3, 2019
//      Change default for convolution of Compton scatter to true (for use with calculated instead of SNIP background)
//      Update background in add_component if the component is type background
//  Modified Oct. 21, 2020
//      Add multiple background components (split for independent fitting)
//      Include Compton escape component in background
//  Modified Nov. 24, 2020
//      Add matrix effect factor to update_component
//  Modified Nov. 30, 2020
//      Exclude from fit vector any components that have fit Boolean set to false
//      Also add factor to compute coefficients of non-fit components from components used for quant
//  Modified Dec. 28, 2020
//      Add adjusted coefficients (for better match when composition changed, must be set elsewhere) and function to reset coefficients to unity
//  Modified Jan. 4, 2021
//      Add residual error calculation, for contribution to uncertainty from fit residual for each component
//  Modified Feb. 26, 2021   Add detector shelf, tail, and Compton escape to background
//  Modified Mar. 9, 2021    Don't double-count any background components
//  Modified Mar. 15, 2021   Adjust coefficients for all components, not just non-fit (necessary for shelf calculation)
//  Modified Apr. 6, 2021    Add plot and bkg flags in components (for separate control to cover lots of cases, which are still changing)
//  Modified Apr. 11, 2021   Add function to change coefficient using index (used for bkg adjustment factor)
//  Modified May 10, 2021    Add storage for -bh and -bx background options, eliminate adj_calc_bkg (from -a option)   [changes only in header]
//                           Comment out bkg_SNIP for now, may need it in the future
//  Modified June 9, 2021   Fix bug in energy per channel calculation that was disturbing convolution normalization  (header change only)
//                          Add geometry factor so it can be written to bulk sum MSA files  (header change only)


using namespace std;

// ********************************
//  ***** XrayEnergyCal class *****
// ********************************

XrayEnergyCal::XrayEnergyCal() {
	energyStart_save = 0;
	energyPerChannel_save = 0;
	quad_save = 0;
	offset_save = 0;
	tilt_save = 0;
}

XrayEnergyCal::XrayEnergyCal( const float energyStart_in, const float energyPerChannel_in, const float quad_cal_in ) {
	energyStart_save = energyStart_in;
	energyPerChannel_save = 0;
	if( energyPerChannel_in > 0 )energyPerChannel_save  = energyPerChannel_in;
	quad_save = quad_cal_in;
	offset_save = 0;
	tilt_save = 0;
}


const float XrayEnergyCal::energy_calc( const float channel_in, const bool corrected ) const {
    if( energyPerChannel_save > 0 ) {
        float temp_energy = energyStart_save + channel_in * energyPerChannel_save
        + channel_in * channel_in * quad_save;
        if( corrected ) temp_energy += channel_in * tilt_save + offset_save;
        //  Implement linear energy calibration correction from Chris
        temp_energy += linear_correction( temp_energy );
        return temp_energy;
    } else {
        return channel_in;
    }
};

const float XrayEnergyCal::channel_calc( const float energy_in, const bool corrected ) const {
    float energy_calc = energy_in;
    float offset = energyStart_save;
    float eV_ch = energyPerChannel_save;
    if( corrected ) {
        offset += offset_save;
        eV_ch += tilt_save;
    }
    //  Implement linear energy calibration correction from Chris
    //  Assume correction is small so it doesn't have to be inverted
    energy_calc -= linear_correction( energy_calc );
    if( eV_ch > 0 ) {
        if( quad_save < 1e-8f * energyPerChannel_save ) {
            float ch = ( energy_calc - offset ) / eV_ch;
            return ch;
        } else {
            double b = eV_ch;
            double c = offset - energy_calc;
            double b2_4ac = b*b - 4 * quad_save * c;
            if( b2_4ac < 0 ) return 0;
            double num = - b + sqrt( b2_4ac );
            if( num <= 0 ) return 0;
            double guess_plus = num / ( 2 * quad_save );
            return guess_plus;
        }
    } else {
        return 0;
    }
};

const float XrayEnergyCal::linear_correction( const float energy_in ) const {
    //  Implement linear energy calibration correction from Chris
    //  eV_shift = slope * (keV_energy) + offset
    if( energyCorrectionSlope_save == 0 ) return 0;
    float x_intercept = ( - energyCorrectionOffset_save / energyCorrectionSlope_save ) * 1000;
    if( energy_in > x_intercept ) return 0;
    float correction = energyCorrectionSlope_save * ( energy_in / 1000 ) + energyCorrectionOffset_save;
    return correction;
};

std::string XrayEnergyCal::toString() const
{
    ostringstream os;
    os << "XrayEnergyCal:" << endl;
    os << "  energyStart_save" << energyStart_save << endl;
    os << "  energyPerChannel_save" << energyPerChannel_save << endl;
    os << "  quad_save" << quad_save << endl;
    os << "  offset_save" << offset_save << endl;
    os << "  tilt_save" << tilt_save << endl;
    os << "  energyCorrectionOffset_save" << energyCorrectionOffset_save << endl;
    os << "  energyCorrectionSlope_save" << energyCorrectionSlope_save << endl;

    return os.str();
}



// *******************************
//  ***** XraySpectrum class *****
// *******************************

XraySpectrum::XraySpectrum() {
	live_time_save = 0;;
	real_time_save = 0;
};

XraySpectrum::XraySpectrum( const vector <float> &counts_in, const float energyStart_in,
						   const float energyPerChannel_in, const float quad_cal_in ) {
	XrayEnergyCal new_calibration( energyStart_in, energyPerChannel_in, quad_cal_in );
	if( new_calibration.good() ) spectrum_calibration = new_calibration;
	live_time_save = 0;;
	real_time_save = 0;
	if( counts_in.size() <= 0 ) return;
	setup_measured( counts_in );
};

XraySpectrum::XraySpectrum( const vector <float> &counts_in ) {
	live_time_save = 0;;
	real_time_save = 0;
	if( counts_in.size() <= 0 ) return;
	setup_measured( counts_in );
};

XraySpectrum::XraySpectrum( const int nc_in, const float counts_in[] ) {
	live_time_save = 0;;
	real_time_save = 0;
	if( nc_in <= 0 ) return;
	vector <float> temp_meas( nc_in );
	int i;
	for( i=0; i<nc_in; i++ ) temp_meas[i] = counts_in[i];
	setup_measured( temp_meas );
};

void XraySpectrum::meas( std::vector <float> &counts_in ) {
	setup_measured( counts_in );
};

void XraySpectrum::bkg( const std::vector <float> &background_in, const float multiplier ) {
	move_spectrum( background_in, background, multiplier );
	vector <float> temp_net( measured_data.size() );
	int i;
	for( i=0; i<measured_data.size(); i++ ) {
        temp_net[i] = measured_data[i] - background[i];
	}
	move_spectrum( temp_net, measured_net );
};

void XraySpectrum::bkg( float background_in[] ) {
	vector <float> temp_bkg( measured_data.size() );
	int i;
	for( i=0; i<temp_bkg.size(); i++ ) temp_bkg[i] = background_in[i];
	bkg( temp_bkg );
};

void XraySpectrum::calc( std::vector <float> &calculation_in ) {
	move_spectrum( calculation_in, calculation );
	//  New residual (if measured data available)
    float new_chisq = 0;
    unsigned int i;
	if( measured_data.size() > 0 ) {
        vector <float> temp_res( measured_data.size() );
        for( i=0; i<measured_data.size(); i++ ) {
            temp_res[i] = measured_data[i] - calculation[i];
            //  Calculate reduced chi squared for new fit
            new_chisq += temp_res[i] * temp_res[i] / ( measured_sigma[i] * measured_sigma[i] );
        }
        move_spectrum( temp_res, residual_calc );
        //  Find number of included components in fit (number of degrees of freedom in fit)
        int nce = 0;
        int ic;
        for( ic=0; ic<components.size(); ic++ ) {
            if( components[ic].included ) nce++;
        }
        residual_chisq = new_chisq / ( measured_data.size() - nce );
    }
    //  If measured spectrum total counts and region counts are zero, use calculated values
    if( total_counts_save <= 0 && calculation.size() > 0 ) {
        total_counts_save = 0;
        for( i=0; i<calculation.size(); i++ ) {
            total_counts_save += calculation[i];
        }
	}
    //  PIXL L5 requirements info for X-ray Subsystem
    if( region_counts_save <= 0 && calculation.size() > 0 ) {
        int range_counts_start = channel( range_counts_start_energy );
        if( range_counts_start < 0 ) range_counts_start = 0;
        if( range_counts_start >= calculation.size() ) range_counts_start = calculation.size() - 1;
        int range_counts_end = channel( range_counts_end_energy );
        if( range_counts_end < 0 ) range_counts_end = 0;
        if( range_counts_end >= calculation.size() ) range_counts_end = calculation.size() - 1;
        region_counts_save = 0;
        unsigned int is;
        for( is=range_counts_start; is<range_counts_end; is++ )
                region_counts_save += calculation[is];
    };
    //  Calculate contribution to uncertainty from fit residual for each component
    int ic;
    for( ic=0; ic<components.size(); ic++ ) {
        if( ! components[ic].included ) continue;
        if( components[ic].spectrum.size() < residual_calc.size() ) continue;
        //  Weight residual by component amplitude and normalize by amplitude squared
        //      to get relative residual for this component
        float residual_wgtd_sum = 0;
        float res_err_norm = 0;
        int is;
        for( is=0; is<residual_calc.size(); is++ ) {
            float res = fabs( residual_calc[is] );
            float sp = components[ic].spectrum[is];
            residual_wgtd_sum += res * sp;
            res_err_norm += sp * sp;
        }
        if( res_err_norm > 0 ) components[ic].residual_err = residual_wgtd_sum / res_err_norm;
    }
};

void XraySpectrum::calc( float calculation_in[] ) {
	vector <float> temp_calc( measured_data.size() );
	int i;
	for( i=0; i<temp_calc.size(); i++ ) temp_calc[i] = calculation_in[i];
	calc( temp_calc );
};

void XraySpectrum::calibration( const float energyStart_in,
                     const float energyPerChannel_in, const float quad_cal_in ) {
	XrayEnergyCal new_calibration( energyStart_in, energyPerChannel_in, quad_cal_in );
	if( new_calibration.good() ) spectrum_calibration = new_calibration;
	return;
};

void XraySpectrum::setup_measured( const vector <float> &meas_in ) {
    const int nc = meas_in.size();
	measured_data.resize( nc );
	measured_sigma.resize( nc );
	int sum = 0;
	int i;
	for( i=0; i<nc; i++ ) {
		measured_data[i] = meas_in[i];
		sum += measured_data[i];
		measured_sigma[i] = ( meas_in[i] > 0? sqrt( meas_in[i] + 2 ): sqrt(2) );
	}
	total_counts_save = sum;
    //  PIXL L5 requirements info for X-ray Subsystem
    int range_counts_start = channel( range_counts_start_energy );
    if( range_counts_start < 0 ) range_counts_start = 0;
    if( range_counts_start >= measured_data.size() ) range_counts_start = measured_data.size() - 1;
    int range_counts_end = channel( range_counts_end_energy );
    if( range_counts_end < 0 ) range_counts_end = 0;
    if( range_counts_end >= measured_data.size() ) range_counts_end = measured_data.size() - 1;
    sum = 0;
    int is;
    for( is=range_counts_start; is<range_counts_end; is++ )
            sum += measured_data[is];
    region_counts_save = sum;
};


//  Functions for handling components of calculated spectrum (also used for fits)
void XraySpectrum::add_component( const SpectrumComponent &component_in ) {
    int ic = find_component( component_in );
    //  If component exists already, replace it
    if( ic >= 0 ) {
        components[ic] = component_in;
        update_intensity( components[ic] );
     // If component does not exist already, add it
    } else {
        components.push_back( component_in );
        update_intensity( components[components.size()-1] );
    }
    return;
};

const SpectrumComponent &XraySpectrum::component( const int index_in ) const {
    if( index_in >= 0 && index_in < components.size() )
        return components[index_in];
    static SpectrumComponent temp_component;
    temp_component.type = NO_COMPONENT;
    temp_component.coefficient = 0;
    return temp_component;
};

void XraySpectrum::update_component( const SpectrumComponent &component_in ) {
    int ic = find_component( component_in );
    if( ic >= 0 ) {
        move_spectrum( component_in.spectrum, components[ic].spectrum );
        update_intensity( components[ic] );
        components[ic].matrix = component_in.matrix;
    }
    return;
};

const int XraySpectrum::index( const Element el_in ) const {
    //  Returns component used to quantify this element
    int ic = find_component( el_in );
    return ic;
};

const float XraySpectrum::coefficient( const int index_in ) const {
    //  Returns coefficient of component
    if( index_in >= 0 && index_in < components.size() ) return components[index_in].coefficient;
    return 0;
};

const float XraySpectrum::coefficient( const Element el_in ) const {
    //  Returns coefficient of component used to quantify this element
    int ic = find_component( el_in );
    if( ic >= 0 && ic < components.size() && components[ic].enabled ) return components[ic].coefficient;
    //  Return special value if no component for this element
    return COEFFICIENT_NO_COMPONENT;
};

int XraySpectrum::update_coefficient( const int index_in, const float new_coefficient ) {
    //  Changes the value of the coefficient of a component specified by index
    if( index_in >= 0 && index_in < components.size() ) {
        components[index_in].coefficient = new_coefficient;
        return 0;
    } else {
        //  Return error value if no component for this element
        return -1;
    }
};


void XraySpectrum::adjusted_coefficient( const Element el_in, const float adj_coeff_in ) {
    //  Saves the adjusted coefficient for this element
    int ic = find_component( el_in );
    if( ic >= 0 && ic < components.size() ) components[ic].adjusted_coefficient = adj_coeff_in;
    return;
};

void XraySpectrum::adjust_coefficients() {
    //  Moves adjusted coefficients into actual coefficients (to get better match to updated composition)
    int ic;
    for( ic=0; ic<components.size(); ic++ ) {
        if( components[ic].adjusted_coefficient > 0 ) {
            components[ic].coefficient = components[ic].adjusted_coefficient;
            update_intensity( components[ic] );
        }

    }
    update_non_fit_coefficients();
};

void XraySpectrum::reset_coefficients() {
    //  Changes all coefficients to unity for components included in (or affected by) fit
    int ic;
    for( ic=0; ic<components.size(); ic++ ) {
        if( components[ic].enabled && components[ic].fit ) {
            components[ic].coefficient = 1;
            update_intensity( components[ic] );
        }
    }
    update_non_fit_coefficients();
};


void XraySpectrum::disable( const Element el_in ) {
    //  Removes component used to quantify this element from the fit
    int ic = find_component( el_in );
    disable( ic );
    return;
}

void XraySpectrum::disable( const int index_in ) {
    //  Removes component with this index from the fit
    if( index_in >= 0 && index_in < components.size() ) {
        components[index_in].enabled = false;
        components[index_in].plot = false;
    }
    return;
};

void XraySpectrum::enable( const int index_in ) {
    //  Puts component with this index back into the fit (used for optic response)
    if( index_in >= 0 && index_in < components.size() ) {
        components[index_in].enabled = true;
        components[index_in].plot = true;
    }
    return;
};


const float XraySpectrum::variance( const int index_in ) const {
    //  Returns variance of coefficient
    if( index_in >= 0 && index_in < components.size() ) return components[index_in].variance;
    return 0;
};

const float XraySpectrum::residual_error( const int index_in ) const {
    //  Returns contribution to uncertainty from fit residual for this component
    if( index_in >= 0 && index_in < components.size() ) return components[index_in].residual_err;
    return 0;
};

const float XraySpectrum::intensity( const int index_in ) const {
    //  Returns total intensity of component
    if( index_in >= 0 && index_in < components.size() ) return components[index_in].intensity;
    return 0;
};

void XraySpectrum::update_background() {
    //  Sum background components to get new full background
    //  NOTE: calculation, net spectrum, and fit residual must also be updated
    //      net and residual taken care of in bkg function
    //      calculation is updated separately (this is called from update_calc)
    vector <float> temp_bkg( measured_data.size(), 0 );
    bool found_bkg = false;
    int ic;
    for( ic=0; ic<components.size(); ic++ ) {
        if( ! components[ic].bkg ) continue;
        if( ! components[ic].enabled ) continue;
        //  If the coefficient was made negative in the most recent fit, skip it
        if( components[ic].coefficient <= 0 ) continue;
        if( components[ic].spectrum.size() < measured_data.size() ) continue;
        found_bkg = true;
        int is;
        for( is=0; is<measured_data.size(); is++ ) {
            float value = components[ic].coefficient * components[ic].spectrum[is];
            temp_bkg[is] += value;
        }
    }
    if( !found_bkg ) return;
    bkg( temp_bkg );
};


void XraySpectrum::update_calc( ) {
    //  First update the background (and net spectrum)
    update_background();
    vector <float> temp_calc( measured_data.size(), 0 );
    int ic;
    for( ic=0; ic<components.size(); ic++ ) {
        if( components[ic].type == NO_COMPONENT ) continue;
        //  Don't double-count any background components
        if( components[ic].bkg ) continue;
        if( ! components[ic].enabled ) continue;
        if( components[ic].spectrum.size() < measured_data.size() ) continue;
        int is;
        for( is=0; is<measured_data.size(); is++ ) {
            float value = components[ic].coefficient * components[ic].spectrum[is];
            temp_calc[is] += value;
        }
        update_intensity( components[ic] );
    }
    //  Add the background to the calculation
    if( background.size() >= measured_data.size() ) {
        int is;
        for( is=0; is<measured_data.size(); is++ ) temp_calc[is] += background[is];
    }
    //  Put the new calculation in place (and update residual using new calculation)
	calc( temp_calc );
};

void XraySpectrum::fit_vector( std::vector <float> &componentSpectra,
                std::vector <float> &coefficients_out, std::vector <float> &centerEnergy  ) {
    //  Produces vector of all enabled component spectra for least squares fit
    //      and saves the relation between component spectra and fit coefficients
    int ns = measured_data.size();
    fit_vector_indices.resize(0);
    //  Check for very small components compared to largest and don't include
    //  (They make the fit unstable)
    float largest_int = 0;
    int ic;
    for( ic=0; ic<components.size(); ic++ ) {
        //  Keep track of which components are included in fit
        components[ic].included = false;
        //  Any changes here match similar loops in update_coefficients
        if( ! components[ic].enabled ) continue;
        if( ! components[ic].fit ) continue;
        if( components[ic].spectrum.size() < ns ) continue;
        float spec_sum = 0;
        int is;
        for( is=0; is<ns; is++ ) spec_sum += components[ic].spectrum[is];
        if( spec_sum <= 0 || isnan( spec_sum ) ) continue;   //  Avoid singular fit matrix
        if( spec_sum > largest_int ) largest_int = spec_sum;
        if( spec_sum / largest_int < 1e-10f )  continue;    //  Avoid fit instability from very small components
        //  Mark this component as included and save its index
        components[ic].included = true;
        fit_vector_indices.push_back( ic );
    }
    //  Resize the output arguments and the association vector of component indices
    componentSpectra.resize( fit_vector_indices.size()*ns );
    coefficients_out.resize( fit_vector_indices.size() );
    centerEnergy.resize( fit_vector_indices.size() );
    fit_vector_indices.resize( fit_vector_indices.size() );
    //  Load component spectra into vector
    int ic_fit;
    for( ic_fit=0; ic_fit<fit_vector_indices.size(); ic_fit++ ) {
        int ic = fit_vector_indices[ic_fit];
        coefficients_out[ic_fit] = components[ic].coefficient;
        //  Find energy of largest peak in each component to use for calculating shifts and widths
        float componentSpecMax = 0;
        centerEnergy[ic_fit] = 0;
        int is;
        for( is=0; is<ns; is++ ) {
            componentSpectra[ ic_fit*ns + is ] = components[ic].spectrum[is];
            if( components[ic].type == ELEMENT && components[ic].spectrum[is] > componentSpecMax ) {
                componentSpecMax = components[ic].spectrum[is];
                centerEnergy[ic_fit] = energy( is );
            }
        }
        fit_vector_indices[ic_fit] = ic;
//        cout << "fit_vector entry " << ic_fit << " is component number " << ic << " - " << componentDescription( components[ic] ) << endl;
    }
//        cout << "fit_vector " << fit_vector_indices.size() << "     " << coefficients_out.size() << endl;
    return;
};

int XraySpectrum::update_coefficients( const std::vector <float> &new_coefficients,
                   const std::vector <float> &new_variances ) {
    //  Update the values of the coefficients for all enabled components included in the least squares fit
    //  This list must match the spectra returned by the above function of fit spectra
    //  First find the number of enabled components and check it
    int ns = measured_data.size();
//        cout << "update_coefficients " << fit_vector_indices.size() << "     " << new_coefficients.size() << "     " <<( new_coefficients.size() != fit_vector_indices.size() ) << endl;
    if( new_coefficients.size() != fit_vector_indices.size() ) return -1;
    int ic_fit;
    for( ic_fit=0; ic_fit<new_coefficients.size(); ic_fit++ ) {
        //  Get the component index corresponding to this coefficient from the least squares fit
        int ic = fit_vector_indices[ic_fit];
        if( components[ic].spectrum.size() < ns ) return -2;
        components[ic].coefficient = new_coefficients[ic_fit];
        if( ic_fit < new_variances.size() ) components[ic].variance = new_variances[ic_fit];
        update_intensity( components[ic] );
//        cout << "update_coefficients entry " << ic_fit << " is component number " << ic << " - " << componentDescription( components[ic] ) << "   new coeff " << new_coefficients[ic_fit] << endl;
    }
    update_non_fit_coefficients();
    update_calc();
    return 0;
};

void XraySpectrum::update_non_fit_coefficients() {
    //  Update the coefficients for non-fit components
    int ic;
    for( ic=0; ic<components.size(); ic++ ) {
        if( components[ic].fit ) continue;
        float quant_coeff = coefficient( components[ic].element );
        if( quant_coeff == COEFFICIENT_NO_COMPONENT ) continue;
        //  Follow coefficient of component used to quantify this element using non-fit factor stored in component
        components[ic].coefficient = components[ic].non_fit_factor * quant_coeff;
    }
};

void XraySpectrum::clean() {
    //  Remove the component spectra to save storage space when lots of spectra processed
    int ic;
    for( ic=0; ic<components.size(); ic++ ) {
        components[ic].spectrum.clear();
        //  Force reallocation to actually free up the space
        //  (This will eventually be replaced by the shrink_to_fit member function)
        vector<float>(components[ic].spectrum).swap(components[ic].spectrum);
    }
};

void XraySpectrum::reset() {
    //  Reset everything except the measured spectrum and energy calibration
    clean();
    background.clear();
    measured_net.clear();
	calculation.clear();
	residual_calc.clear();
	residual_chisq = 0;
	components.clear();
    //  Force reallocation to actually free up the space
    //  (This will eventually be replaced by the shrink_to_fit member function)
    vector<SpectrumComponent>(components).swap(components);
};



//  Utility functions

void XraySpectrum::move_spectrum( const vector <float> &vec_in, vector <float> &vec_out, const float factor ) {
    //  Output vector is always the size of the measured data vector
    //  if it is the measured data vector, it must be sized already
    //  The factor is for coefficients that must multiply fit components
    //      (bkg is the only one this is used for right now)
    int nc = measured_data.size();
    vec_out.resize( nc );
    //  Don't exceed the size of the either vector when copying
    int nc_in = vec_in.size();
    //  If the input vector is larger, only copy part of it
    if( nc_in > nc ) nc_in = nc;
    int i;
    for( i=0; i<nc_in; i++ ) vec_out[i] = vec_in[i] * factor;
    //  Zero fill if the input vector is smaller
    if( nc_in < nc ) {
        for( i=nc_in; i<nc; i++ ) vec_out[i] = 0;
    }
};

int XraySpectrum::find_component( const Element &el_in ) const {
    int index = -1;
    int ic;
    for( ic=0; ic<components.size(); ic++ ) {
        if( components[ic].type != ELEMENT
            || ! ( components[ic].element == el_in )
            || ! components[ic].quant ) continue;
        index = ic;
        break;
    }
    return index;
};

int XraySpectrum::find_component( const SpectrumComponent &component_in ) const {
    int index = -1;
    int ic;
    for( ic=0; ic<components.size(); ic++ ) {
        if( ! matchComponent( components[ic], component_in ) ) continue;
        index = ic;
        break;
    }
    return index;
};

void XraySpectrum::update_intensity( SpectrumComponent &component_in ) {
        float sum = 0;
        int is;
        for( is=0; is<component_in.spectrum.size(); is++ ) {
            float value = component_in.coefficient * component_in.spectrum[is];
            sum += value;
        }
        component_in.intensity = sum;
        return;
};


string XraySpectrum::toString() const
{
    ostringstream os;
    os << "XraySpectrum:" << endl;
    os << "  measured_data: " << floatVecToString(measured_data) << endl;
    os << "  measured_sigma: " << floatVecToString(measured_sigma) << endl;
    os << "  background: " << floatVecToString(background) << endl;
    os << "  measured_net: " << floatVecToString(measured_net) << endl;
    os << "  calculation: " << floatVecToString(calculation) << endl;
    os << "  residual_calc: " << floatVecToString(residual_calc) << endl;
    os << "  max_value_save: " << floatVecToString(max_value_save) << endl;
    os << "  residual_chisq: " << residual_chisq << endl;
    os << "  live_time_save: " << live_time_save << endl;
    os << "  real_time_save: " << real_time_save << endl;
    os << "  total_counts_save: " << total_counts_save << endl;
    os << "  range_counts_start_energy: " << total_counts_save << endl;
    os << "  range_counts_end_energy: " << total_counts_save << endl;
    os << "  region_counts_save: " << total_counts_save << endl;

    os << "  bkg_params_save: " << floatVecToString(bkg_params_save) << endl;
    os << "  spectrum_calibration: " << "TODO" << endl;
    os << "  components: " << "sz=" << components.size() << endl;

    int c = 0;
    for(auto it = components.begin(); it != components.end(); it++)
    {
        os << "  components[" << c << "]: " << SpectrumComponent_toString(*it) << endl;
        c++;
    }

    os << "  aux_info_save: " << "TODO" << endl;
    os << "  header_info_save: " << "TODO" << endl;

    os << "  file_name_save: " << file_name_save << endl;
    os << "  seq_number_save: " << seq_number_save << endl;
    os << "  iterations_save: " << iterations_save << endl;
    os << "  adjust_energy_save: " << adjust_energy_save << endl;
    os << "  adjust_width_save: " << adjust_width_save << endl;
    os << "  convolve_Compton_save: " << convolve_Compton_save << endl;

    return os.str();
}
