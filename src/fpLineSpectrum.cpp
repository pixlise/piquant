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
#include <vector>
#include <math.h>
#include "fpConvolve.h"
#include "XrayDetector.h"
#include "XRFcontrols.h"
#include "fpLineSpectrum.h"

using namespace std;

//	Generates calculated spectrum from Xray Lines object that has been loaded with intensity factor
//		Calculation is Gaussian with fwhm as input, integral matches line intensity
//		Calculated spectrum is counts in each channel
//	Added check for zero or negative energy at low channels     Dec. 12, 2011
//     Threshold causing problems with low concentration elements, replace with 10^-7  Oct 2012
//  Added generation for escape peaks     Oct. 31, 2013
//	Modified June 18, 2014
//		Block generation of U M lines (< 10 keV) to avoid problems with K line fits -- ** kludge warning
//  Modified May 10, 2015
//      Comment out above kludge
//  Modified Feb. 9, 2017
//      To use XrayEnergyCal for energy to channel conversions
//  Re-written Feb. 11, 2017 as fpLineSpectrum
//      to use SpectrumComponent to select XrayLines to include
//      This allows decopupling of K, L, M, and N lines to have separate fit coefficients
//      and makes future changes in which lines are included much easier
//  Modified Dec. 13, 2017
//      Add check for minimum energy to escape peaks
//  Modified July 25, 2018
//      Write out some useful information if calculated intensity is zero
//      Correct resolution of escape peaks to value at escape energy, not line energy
//  Modified Sept. 2, 2020
//      Add peak tails from incomplete charge collection
//  Modified Nov. 24, 2020
//      Move matrix effect factor from XrayLines to spectrum component
//  Modified Apr. 2, 2021  Major rearrangements to improve speed
//      Tail total intensity estimate using full energy range, not sum over individual channel intensities
//      Do a single convolution at end instead of each peak location and tail channel
//      Group lines within the detector resolution for tail and shelf calculations
//      Include detector electron loss shelf with each line (in groups)
//  Modified May 14, 2021   Move shelf factor and slope to XrayDetector and control via -T option
//  Modified May 25, 2021   Added symbol to grouped lines, for identification during debugging
//  Modified July 10, 2021  Add simple pulse pileup calculation - return peak intensity information and use average energy for grouped lines


void fpLineSpectrum( const XrayLines &lines_in, const XrayDetector detector, const float threshold_in,
				   const XrayEnergyCal &cal_in, const float eMin, std::vector <LineGroup> &pileup_list,
				   SpectrumComponent &component_out ) {
	int ns = component_out.spectrum.size();
	if ( ns <= 0 ) return;

 //		check threshold against strongest line to be sure some channels will be generated
	int j;
	float max_lines = 0;
//      Also save matrix effect factor from strongest line
    float matrix_factor = 0;
	for ( j=0; j<lines_in.numberOfLines(); j++ ) {
        if( ! checkComponent( component_out, lines_in, j ) ) continue;
        if( max_lines < lines_in.intensity(j) ) {
            max_lines = lines_in.intensity(j);
            matrix_factor = lines_in.matrix(j);
        }
    }
	if( max_lines <= 0 || isnan( max_lines ) ) {
//        cout << "Warning - all emission lines have zero (or nan) calculated intensity for component " << componentDescription( component_out ) << endl;
        return;
    }
    component_out.matrix = matrix_factor;
//    cout << "fpLS  M " << componentDescription( component_out ) << "  mf " << component_out.matrix << endl;;

    //  Consolidate the lines into a few groups by FWHM for tail and shelf calculations
    vector <LineGroup> grouped_lines;

//		calculated spectrum
	for ( j=0; j<lines_in.numberOfLines(); j++ ) {
        //  Check to see if this emission line should be included
        if( ! checkComponent( component_out, lines_in, j ) ) continue;
		float line_energy = lines_in.energy(j);
        if( line_energy < eMin ) continue;
//			turn off uranium M lines to avoid bad fits from Compton overlap
//		if (lines_in.edge().element().Z() == 92 && line_energy < 10000 ) continue;
        //  get info on escape peaks from detector
        vector <EscapeLines> escape_info;
        float non_escape_fraction = 1;
        if( escape_peaks_enable_flag ) {    //  Defined in XRFcontrols.h
            non_escape_fraction = detector.escape( line_energy, escape_info );
        }
        float int_minus_esc = lines_in.intensity(j) * non_escape_fraction;
        //  Calculate peak tail
        float tail_end_energy = detector.energy_for_C0( line_energy );
        if( tail_end_energy < eMin ) tail_end_energy = eMin;
        int peak_channel = cal_in.channel( line_energy ) + 0.5f;
        if( peak_channel < 0 || peak_channel >= ns ) continue;

        float tail_sum = 0;
        if( peak_tail_enable_flag ) {    //  Defined in XRFcontrols.h
            //  Calculate total tail intensity from lowest tail energy to peak energy
            tail_sum = lines_in.intensity(j) * detector.tail_fraction( line_energy, tail_end_energy, line_energy );
        }
        float main_peak_intensity = int_minus_esc - tail_sum;
        //  Add this peak into the line groups
        float fwhm_in = detector.resolution( line_energy );
        int found = -1;
        int ig;
        for( ig=0; ig<grouped_lines.size(); ig++ ) {
            if( fabs( grouped_lines[ig].energy - line_energy ) < fwhm_in ) {
                found = ig;
                break;
            }
        }
        if( found >= 0 && grouped_lines[found].number > 0 ) {
            //  Calculate energy via averaging with intensity weighting
            grouped_lines[found].energy = grouped_lines[found].energy * grouped_lines[found].intensity + line_energy * main_peak_intensity;
            grouped_lines[found].intensity += main_peak_intensity;
            grouped_lines[found].energy /= grouped_lines[found].intensity;
            grouped_lines[found].tail_sum += tail_sum;
            grouped_lines[found].number++;
        } else {
            LineGroup new_group;
            new_group.energy = line_energy;
            new_group.number = 1;
            new_group.intensity = int_minus_esc;
            new_group.tail_sum = tail_sum;
            new_group.symbol = lines_in.edge().element().symbol() + " " + lines_in.symbolIUPAC(j);
            grouped_lines.push_back( new_group );
        }
        //  Add main peak with Lorentzian using natural linewidth
        float line_width = lines_in.width(j);
        float peak_width = line_width * 10;  //  Arbitrary cutoff for Lorentzian, which has infinite tails
        float gamma2 = line_width*line_width / 4;       //  Half width at half maximum, the scale parameter for the Lorentzian
		int kMin = cal_in.channel( line_energy - peak_width ) - 1;
		int kMax = cal_in.channel( line_energy + peak_width ) + 1;  //  Be sure there are at least 2 channels in peak
        if ( kMin < 0 ) kMin = 0;
        if ( kMax > ns ) kMax = ns;
        //  Find Lorentzian integral empirically since we are using only a few points
        float width_integral = 0;
        int k;
		for ( k=kMin; k<=kMax; k++ ) {
			float en = cal_in.energy( k );
			if( en <= 0 ) continue;
			float diff = en - line_energy;
			width_integral += 1 / ( diff*diff + gamma2 );
        };
        if( width_integral <= 0 ) continue;
        //  Put the Lorentzian line shape into the spectrum (will be broadened by detector resolution at the end of this function)
        float norm_int = main_peak_intensity / width_integral;
		for ( k=kMin; k<=kMax; k++ ) {
			float en = cal_in.energy( k );
			if( en <= 0 ) continue;
			float diff = en - line_energy;
			component_out.spectrum[k] += norm_int / ( diff*diff + gamma2 );
		};

        if( escape_peaks_enable_flag ) {
            int i_esc;
            for( i_esc=0; i_esc<escape_info.size(); i_esc++ ) {
                float el_esc = escape_info[i_esc].energy;
                if ( el_esc < eMin) continue;
                //  Put the same Lorentzian line shape into the spectrum for the escape peak
                kMin = cal_in.channel( el_esc - peak_width ) - 1;
                if ( kMin < 0 ) kMin = 0;
                kMax = cal_in.channel( el_esc + peak_width ) + 1;
                if ( kMax > ns ) kMax = ns;
                width_integral = 0;
                for ( k=kMin; k<=kMax; k++ ) {
                    float en = cal_in.energy( k );
                    if( en <= 0 ) continue;
                    float diff = en - line_energy;
                    width_integral += 1 / ( diff*diff + gamma2 );
                };
                float norm_int = lines_in.intensity(j) * escape_info[i_esc].fraction / width_integral;
                for ( k=kMin; k<=kMax; k++ ) {
                    float en = cal_in.energy( k );
                    if( en <= 0 ) continue;
                    float diff = en - line_energy;
                    component_out.spectrum[k] += norm_int / ( diff*diff + gamma2 );
                };
            };
        };
	};

	//  Include an incomplete charge collection tail for each line (as grouped)
    if( peak_tail_enable_flag ) {    //  Defined in XRFcontrols.h
        int ig;
        for( ig=0; ig<grouped_lines.size(); ig++ ) {
//            cout << "fpLS group " << ig << "  " << grouped_lines[ig].symbol << "  " << grouped_lines[ig].energy;
//            cout << "  " << grouped_lines[ig].number << "  " << grouped_lines[ig].tail_sum << "  " << grouped_lines[ig].intensity << endl;
            //  Calculate tail for this peak from incomplete charge collection
            int i_tail;
            float line_energy = grouped_lines[ig].energy;
            float tail_end_energy = detector.energy_for_C0( line_energy );
            if( tail_end_energy < eMin ) tail_end_energy = eMin;
            int peak_channel = cal_in.channel( line_energy ) + 0.5f;
            if( peak_channel < 0 || peak_channel >= ns ) continue;
            int tail_end_channel = cal_in.channel( tail_end_energy );
            if( tail_end_channel < 0 ) tail_end_channel = 0;
            float tail_prev_en = tail_end_energy;
            for( i_tail=tail_end_channel+1; i_tail<peak_channel; i_tail++ ) {
                float tail_new_energy = cal_in.energy( i_tail );
                float tail_fraction = detector.tail_fraction( line_energy, tail_prev_en, tail_new_energy );
                tail_prev_en = tail_new_energy;
                float tail_int = grouped_lines[ig].intensity * tail_fraction;
                if( tail_int <= 0 ) continue;
                component_out.spectrum[i_tail] += tail_int;
            }
        }
    }
        //  Calculate electron escape contribution to shelf at low energies
    if( detector_shelf_enable_flag ) {    //  Defined in XRFcontrols.h
        int ig;
        for( ig=0; ig<grouped_lines.size(); ig++ ) {
            float photon_energy = grouped_lines[ig].energy;
            if( photon_energy < eMin ) continue;
            float meas_intensity = grouped_lines[ig].intensity;
            if( meas_intensity <= 0 ) continue;
            //  Find original intensity incident on detector by dividing by response at this energy
            float det_resp = detector.response( photon_energy );
            if( det_resp <= 0 ) continue;
            float incoming_int =  meas_intensity / det_resp;
            //  Get shelf contributions for this photon energy
            vector <ShelfStruct> shelf_factors;
            /* int n_shelf = */ detector.electron_shelf( photon_energy, shelf_factors );
            //  Adjustment factors for detector shelf for better quantification
            float det_shelf_factor = detector.get_shelf_factor();
            float det_shelf_slope = detector.get_shelf_slope();
            float det_shelf_slope_start = detector.get_shelf_slope_start();
            //  Loop over the possible electrons that can contribute to the shelf
            unsigned int ishelf;
            for( ishelf=0; ishelf<shelf_factors.size(); ishelf++ ) {
                //cout << "Shelf pri   " << ishelf << "   p " << shelf_factors[ishelf].probability << endl;
                //if( shelf_factors[ishelf].type != PHOTO_ACTIVE_VOLUME && shelf_factors[ishelf].type != AUGER_ACTIVE_VOLUME ) continue;
                //if( shelf_factors[ishelf].type != PHOTO_FRONT_CONTACT && shelf_factors[ishelf].type != AUGER_FRONT_CONTACT ) continue;
                float min_shelf_energy = shelf_factors[ishelf].energy_start;
                float max_shelf_energy = shelf_factors[ishelf].energy_end;
                float electron_energy = ( max_shelf_energy - min_shelf_energy );
                if( electron_energy <= 0 ) continue;
                if( min_shelf_energy < eMin ) min_shelf_energy = eMin;
                unsigned int min_shelf_channel = cal_in.channel( min_shelf_energy );
                unsigned int max_shelf_channel = cal_in.channel( max_shelf_energy );
                if( max_shelf_channel >= ns - 1 ) continue;
                if( min_shelf_channel < 0 || min_shelf_channel > ns - 1 || max_shelf_channel < 0 || max_shelf_channel > ns - 1 ) continue;
                float shelf_intensity = incoming_int * shelf_factors[ishelf].probability / electron_energy; //  Flat distribution, equal probability over electron range
                shelf_intensity *= det_shelf_factor;
                if( shelf_intensity < SHELF_THRESHOLD ) continue;   //  Defined in XRFcontrols.h
                float shelf_slope_start_loss = -det_shelf_slope_start * electron_energy;
                int ish;
                for( ish=min_shelf_channel; ish<=max_shelf_channel; ish++ ) {
                    //  Calculate shelf intensity for this channel
                    float shelf_energy = cal_in.energy( ish );
                    float loss_energy = shelf_energy - photon_energy;
                    float shelf_adjustment = 1;
                    if( loss_energy < shelf_slope_start_loss ) shelf_adjustment += ( loss_energy - shelf_slope_start_loss ) * det_shelf_slope / electron_energy;
                    if( shelf_adjustment < 0 ) continue;   //  Defined in XRFcontrols.h
                    float shelf_int_here = shelf_adjustment * shelf_intensity;
                    if( shelf_int_here < SHELF_THRESHOLD ) continue;   //  Defined in XRFcontrols.h
                    component_out.spectrum[ish] += shelf_int_here;
                }
            }
        }
    }


	//  Now convolve everything with the detector broadening
    fpConvolve( detector, cal_in, component_out.spectrum );

    //  Replace lowest intensity lines in pileup list with any stronger lines from this list
    if( PILEUP_LIST_LENGTH > 0 ) {    //  Defined in XRFcontrols.h
        int ig;
        for( ig=0; ig<grouped_lines.size(); ig++ ) {
            if( grouped_lines[ig].intensity <= 0 ) continue;
            if( pileup_list.size() < PILEUP_LIST_LENGTH ) pileup_list.push_back( grouped_lines[ig] );
            else {
                //  Find the lowest intensity line in the pileup list
                int lowest_index = -1;
                float lowest_intensity = MAXIMUM;   //  Defined in XRFconstants.h
                for ( j=0; j<pileup_list.size(); j++ ) {
                    if( pileup_list[j].intensity < lowest_intensity ) {
                        lowest_intensity = pileup_list[j].intensity;
                        lowest_index = j;
                    }
                }
                //  Replace the lowest intensity entry in the pileup list of this entry has greater intensity
                if( lowest_index >= 0 && lowest_intensity < grouped_lines[ig].intensity ) pileup_list[lowest_index] = grouped_lines[ig];
            }
        }
    }


};
