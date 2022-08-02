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
 *  XraySpectrum.h
 *
 *  Created by Tim Elam on 1/18/2017 from old SpectrumData class
 *  Copyright 2017 University of Washington. All rights reserved.
 *
 */

#ifndef XraySpectrum_h
#define XraySpectrum_h

#include <vector>
#include <math.h>
#include "Element.h"
#include "quantComponents.h"    //  defines SpectrumComponents

 // re-written Feb. 7, 2017 to separate energy calibration and add components
 // Modified May 14, 2017 to include background control parameters in this object (all inline)
 // Modified June 7, 2017 to fix error in returned value for energy per channel and energy start when tilt & offset not zero
 // Modified Dec. 6, 2017 to hold info from DSPC statistics registers
 // Modified Dec. 15, 2017 to add equality boolean operator for energy calibration
 // Modified Jan. 3, 2017 to add max value (only a container for this, all inline)
 // Modified Mar. 2, 2018 to add file name, sequence number, and iterations (only a container for this, all inline)
 // Modified July 31, 2018 to implement linear energy calibration correction from Chris
 // Modified May 13, 2019 to put all non-spectrum information in separate structure for easy modification
//  Modified Dec. 4, 2020 save SNIP background to fit anomalous background at low energies
//  Modified June 9, 2021   Fix bug in energy per channel calculation that was disturbing convolution normalization  (header change only
//                          Add geometry factor so it can be written to bulk sum MSA files  (header change only)

 // The following two structures are separated according to the info that occurs once per column and once per input file

 struct Spec_Aux_Info {
    // Spectrum info not related to quantitative analysis
    string date;
    string time;
    std::vector <string> titles;
    std::vector <string> comments;
	string owner;
     // (File name and seq number still held individually since they are used in the analysis)
 	//  Location information for microXRF scans
	float x = 0;
	float y = 0;
	float z = 0;
	float i = 0;
	float j = 0;
	//  PIXL-specific information
	unsigned int sclk = 0;
	unsigned int rtt = 0;
	unsigned int usn = 0;
	unsigned int dpc = 0;
	unsigned int pmc = 0;
	string det_ID;
};

struct Spec_Header_Info {
	//  Info from DSPC (note this is the fast channel live time, not compensated for throughput)
	float live_time_DSPC = 0;
	int events = 0;
	int triggers = 0;
	int overflows = 0;
	int underflows = 0;
	int baseline_samples = 0;
	int preamp_resets = 0;
	int saturates = 0;
};

//  Put energy calibration in a separate class so it can be passed without entire spectrum object
class XrayEnergyCal {
public:
    XrayEnergyCal();
    XrayEnergyCal( const float energyStart_in, const float energyPerChannel_in, const float quad_cal_in = 0 );
	const float energy( const int channel_in ) const {
        return energy_calc( float( channel_in ), true );
    };
	bool operator==(const XrayEnergyCal& cal_in) const
        { return ( ( energyStart_save + offset_save == cal_in.energyStart() )
        && ( energyPerChannel_save + tilt_save == cal_in.energyPerChannel() ) ); };
	const float energy_uncorrected( const int channel_in ) const {
        return energy_calc( float( channel_in ), false );
    };
	const float energy( const float channel_in ) const {
        return energy_calc( channel_in, true );
    };
	const float energy_uncorrected( const float channel_in ) const {
        return energy_calc( channel_in, false );
    };
	const int channel( const float energy_in ) const {
        return int( channel_calc( energy_in, true ) + 0.5f );
    };
	const float channel_float( const float energy_in ) const {
        return channel_calc( energy_in, true );
    };
	const int channel_uncorrected( const float energy_in ) const {
        return int( channel_calc( energy_in, false ) + 0.5f );
    };
	const float channel_float_uncorrected( const float energy_in ) const {
        return channel_calc( energy_in, false );
    };
	const float energyPerChannel( const int channel_in ) const {
        return  ( 2 * channel_in * quad_save + energyPerChannel_save + tilt_save ); };
	const float energyStart() const { return energyStart_save + offset_save; };
	const float energyPerChannel() const { return energyPerChannel_save + tilt_save; };
	const float quad() const { return quad_save * ( 1 + tilt_save ); };
	const float offset() const { return offset_save; };
	const float tilt() const { return tilt_save; };
    void offset( const float offset_in ) { offset_save = offset_in; };
    void tilt( const float tilt_in ) { tilt_save = tilt_in; };
    bool good() const { return ( energyPerChannel_save > 0 && ! isnan(energyPerChannel_save) ); };
    void linearCorrection( const float lin_offset, const float lin_slope ) { energyCorrectionOffset_save = lin_offset; energyCorrectionSlope_save = lin_slope; return ; };
    const float linearCorrectionOffset( ) const { return energyCorrectionOffset_save; };
    const float linearCorrectionSlope( ) const { return energyCorrectionSlope_save; };

    std::string toString() const;

private:
	float energyStart_save = 0;
	float energyPerChannel_save = 0;
	float quad_save = 0;
	//		these are temporary corrections to the energy calibration
	//		that can be set without losing the original calibration
	float offset_save = 0;
	float tilt_save = 0;
	float energyCorrectionOffset_save = 0;
	float energyCorrectionSlope_save = 0;
    const float energy_calc( const float channel_in, const bool corrected = true ) const;
    const float channel_calc( const float energy_in, const bool corrected = true ) const;
    const float linear_correction( const float energy_in ) const;

};

class XraySpectrum {
public:
	XraySpectrum( const std::vector <float> &counts_in, const float energyStart_in,
				 const float energyPerChannel_in, const float quad_cal_in = 0 );
	XraySpectrum( const std::vector <float> &counts_in );
	XraySpectrum( const int nc_in, const float counts_in[] );
	//		must have default constructor to declare arrays
	XraySpectrum();
	//		data access functions
	const float live_time() const { return live_time_save; };
	const float real_time() const { return real_time_save; };
	const float total_counts() const { return total_counts_save; };
	const float region_counts() const { return region_counts_save; };
	const float region_start() const { return range_counts_start_energy; };
	const float region_end() const { return range_counts_end_energy; };
	int numberOfChannels() const { return measured_data.size(); };
	const std::vector <float> &meas() const { return measured_data; };
	const std::vector <float> &sigma() const { return measured_sigma; };
	const std::vector <float> &bkg() const { return background; };
	const std::vector <float> &net() const { return measured_net; };
	const std::vector <float> &calc() const { return calculation; };
	const std::vector <float> &max_value() const { return max_value_save; };
	const std::vector <float> &residual() const { return residual_calc; };
	void meas( std::vector <float> &measured_data_in );
	void bkg( const std::vector <float> &background_in, const float multiplier = 1 );
	void bkg( float background_in[] );
	void calc( std::vector <float> &calculation_in );
	void calc( float calculation_in[] );
	void max_value( std::vector <float> &max_value_in ) {
        max_value_save.resize( max_value_in.size() );
        int i;
        for( i=0; i<max_value_in.size(); i++ ) max_value_save[i] = max_value_in[i];
        return;
	};
//	const std::vector <float> &bkg_SNIP() const { return bkg_SNIP_save; };
//	void bkg_SNIP( std::vector <float> &bkg_SNIP_in ) {
//        bkg_SNIP_save.resize( bkg_SNIP_in.size() );
//        int i;
//        for( i=0; i<bkg_SNIP_in.size(); i++ ) bkg_SNIP_save[i] = bkg_SNIP_in[i];
//        return;
//	};
	void live_time( const float live_time_in ) { live_time_save = live_time_in; }
	void real_time( const float real_time_in ) { real_time_save = real_time_in; }
	const float geometry() const { return geometry_save; }
	void geometry( const float geometry_in ) { geometry_save = geometry_in; }
	void region_start( const float start_region_in ) { range_counts_start_energy = start_region_in; }
	void region_end( const float end_region_in ) { range_counts_end_energy = end_region_in; }
	//  Pass-through functions to access the energy calibration
	const XrayEnergyCal &calibration() const { return spectrum_calibration; }
	XrayEnergyCal &calibration_change() { return spectrum_calibration; };
	const float energy( const int channel_in ) const { return spectrum_calibration.energy( channel_in ); }
	const int channel( const float energy_in ) const { return spectrum_calibration.channel( energy_in ); }
    void calibration( const XrayEnergyCal &new_calibration ) { if( new_calibration.good() ) spectrum_calibration = new_calibration; }
    void calibration( const float energyStart_in, const float energyPerChannel_in, const float quad_cal_in = 0 );
    void offset( const float offset_in ) { spectrum_calibration.offset( offset_in ); }
    void tilt( const float tilt_in ) { spectrum_calibration.tilt( tilt_in ); }

    //  Functions for handling components of calculated spectrum (also used for fits)
    const int numberOfComponents() const {    return components.size(); }
    void add_component( const SpectrumComponent &component_in );
    const SpectrumComponent &component( const int index_in ) const;
    //  update_calc must be called after components are changed to update calculated spectrum and component intensities
    void update_component( const SpectrumComponent &component_in ); //  replaces calculated spectrum in the component
    const int index( const Element el_in ) const;    //  Returns index of component used to quantify this element
    const float coefficient( const int index_in ) const;    //  Returns coefficient of component
    const float coefficient( const Element el_in ) const;    //  Returns coefficient of component used to quantify this element
    int update_coefficient( const int index_in, const float new_coefficient );
    void adjusted_coefficient( const Element el_in, const float adj_coeff_in );    //  Saves the adjusted coefficient for this element
    void adjust_coefficients();             //  Moves adjusted coefficients into actual coefficients (to get better match to updated composition)
    void reset_coefficients();             //  Changes all coefficients to unity for components included in (or affected by) fit
    void disable( const Element el_in );    //  Removes component used to quantify this element from the fit
    void disable( const int index_in );    //  Removes component with this index from the fit
    void enable( const int index_in );    //  Puts component with this index back into the fit (used for optic response)
    const float variance( const int index_in ) const;    //  Returns variance of coefficient
    const float intensity( const int index_in ) const;    //  Returns total intensity of component
    const float residual_error( const int index_in ) const;    //  Returns error from fit residual for component
    void update_calc( ); //  replaces calculated spectrum by summing enabled components and updates intensity of each component
    //  produces vector of all enabled component spectra plus some other things needed for fits
    void fit_vector( std::vector <float> &componentSpectra,
                std::vector <float> &coefficients_out, std::vector <float> &centerEnergy );
    int update_coefficients( const std::vector <float> &new_coefficients,
                   const std::vector <float> &new_variances ); //  this list must match the spectra returned by the above function
    void clean();   //  Remove the component spectra to save storage space when lots of spectra processed
    void reset();   //  Reset everything except the measured spectrum and energy calibration
    const float chisq() const { return residual_chisq; };
    //  New member functions for variable number of background parameters (for 2-zone SNIP)
    void get_bkg_parameters( std::vector <float> &bkg_params_out ) const { bkg_params_out = bkg_params_save; }
    void put_bkg_parameters( const std::vector <float> &bkg_params_in ) { bkg_params_save = bkg_params_in; }
    void get_bh_parameters( std::vector <float> &bh_params_out ) const { bh_params_out = bh_params_save; }
    void put_bh_parameters( const std::vector <float> &bh_params_in ) { bh_params_save = bh_params_in; }
    void get_bx_parameters( std::vector <float> &bx_params_out ) const { bx_params_out = bx_params_save; }
    void put_bx_parameters( const std::vector <float> &bx_params_in ) { bx_params_save = bx_params_in; }
    void get_bkg_split( std::vector <float> &bkg_split_out ) const { bkg_split_out = bkg_split_save; }
    void put_bkg_split( const std::vector <float> &bkg_split_in ) { bkg_split_save = bkg_split_in; }
    //  Control behavior of spectrum fits
    const bool adjust_energy() const { return adjust_energy_save; }
    const bool adjust_width() const { return adjust_width_save; }
    void adjust_energy( const bool adj_in ) { adjust_energy_save = adj_in; }
    void adjust_width( const bool adj_in ) { adjust_width_save = adj_in; }
    const bool convolve_Compton() const { return convolve_Compton_save; }
    void convolve_Compton( const bool convolve_in ) { convolve_Compton_save = convolve_in; }

	//  File name, and sequence number for maps
    const std::string &file_name() const { return file_name_save; }
    void file_name( const std::string &file_name_in ) { file_name_save = file_name_in; }
	const int seq_number() const { return seq_number_save; }
	void seq_number( const int seq_number_in ) { seq_number_save = seq_number_in; }
	const int iterations() const { return iterations_save; }
	void iterations( const int iterations_in ) { iterations_save = iterations_in; }

	//  Access to all other auxiliary spectrum information
    const Spec_Aux_Info &aux_info() const { return aux_info_save; }
    Spec_Aux_Info &aux_info_change() { return aux_info_save; }
    void aux_info_replace( const Spec_Aux_Info &aux_in ) { aux_info_save = aux_in; }
    const Spec_Header_Info &header_info() const { return header_info_save; }
    Spec_Header_Info &header_info_change() { return header_info_save; }
    void header_info_replace( const Spec_Header_Info &header_in ) { header_info_save = header_in; }

    const std::vector <string> &std_names() const { return std_names_save; }
    void std_names( const std::vector <string> std_names_in ) { std_names_save = std_names_in; }

    std::string toString() const;

private:
	std::vector <float> measured_data;
	std::vector <float> measured_sigma;
	std::vector <float> background;
	std::vector <float> measured_net;
	std::vector <float> calculation;
	std::vector <float> residual_calc;
	std::vector <float> max_value_save;
	std::vector <float> bkg_params_save;
	std::vector <float> bh_params_save;
	std::vector <float> bx_params_save;
	std::vector <float> bkg_split_save;
//	std::vector <float> bkg_SNIP_save;
	float residual_chisq = 0;
	float live_time_save = 0;
	float real_time_save = 0;
	float geometry_save = 0;
	float total_counts_save = 0;
    float range_counts_start_energy = 1000.0f;
    float range_counts_end_energy = 7250.0f;
    float region_counts_save = 0;
	XrayEnergyCal spectrum_calibration;
	std::vector <SpectrumComponent> components;
    // (File name and seq number still held individually since they are used in the analysis)
    //  (All other info from spectrum files not related to quantification moved to separate structure)
    Spec_Aux_Info aux_info_save;
    Spec_Header_Info header_info_save;
	//  File name, and sequence number for maps
	string file_name_save;
	int seq_number_save = 0;
	int iterations_save = 0;
	bool adjust_energy_save = true;
	bool adjust_width_save = true;
	bool convolve_Compton_save = true;
    //  Saves the index of each component included in the least squares fit, in order
    //  Coefficients produced by the fit will be in the same order
 	std::vector <int> fit_vector_indices;
 	std::vector <string> std_names_save;    //  Used to avoid evaluating a standard with itself as a calibration standard

    void reload_energy();
    void setup_measured( const std::vector <float> &meas_in );
    void move_spectrum( const std::vector <float> &vec_in, std::vector <float> &vec_out, const float factor = 1 );
    int find_component( const Element &el_in ) const;
    int find_component( const SpectrumComponent &component_in ) const;
    void update_intensity( SpectrumComponent &component_in );
    void update_non_fit_coefficients();
    void update_background();   //  Used when background components are included in fit
};
#endif
