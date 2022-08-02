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

#include "quantFitSpectrum.h"
#include "Lfit.h"
#include "Fit.h"
#include "differentiate.h"
#include "XRFconstants.h"
#include "XRFcontrols.h"


using namespace std;

//		Use least squares method to fit an XRF spectrum
//          and find the coefficients for each spectrum component
//      This assumes the individual component spectra have already been calculated

//  Written Mar. 8, 2017
//      Based on stdQuant.cpp of Aug. 10, 2016
//      And specElementFit.cpp of June 27, 2016
//  Modified May 26, 2017 to leave out peak shifts and widths with zero or nan sigmas
//  Modified June 7, 2017 to fix convergence criteria (to be relative instead of absolute)
//  Modified Sept. 30, 2017
//      Check resolution change to be sure it's not too large before applying (also check for negative values)
//      Fix bug in resolution calculation for one peak - forgot to take square root
//  Modified Apr. 18, 2018
//      Increase acceptance for change to Fano factor to 40% to avoid fit problems in some noisy spectra
//  Modified May 16, 2019
//      Add capability to turn off adjustments to energy calibration and detector resolution
//  Modified July 3, 2020   Fit to net spectrum if all background components disabled (background not adjusted via fit)
//  Modified Apr. 6, 2021   Always use fit flags in components to subtract any background, never use net spectrum
//                              (So bkg can be several components, some fit and some not)
//                          This implies that everything in the background must always show up as a component
//  Modified Apr. 8, 2021   Ignore zero-energy components for done check (not an element component with a reasonable peak at a known energy)

int quantFitSpectrum( XRFconditions &conditions_in, XraySpectrum &spectrum, std::ostream &logger ) {
//		check input parameters
	if( spectrum.numberOfComponents() <= 0 ) return -801;
	if( ! spectrum.calibration().good() ) return -805;
	if( spectrum.live_time() <= 0 ) return -806;
	int nChan = spectrum.numberOfChannels();

	//  Get the individual component spectra in a single vector for the least squares fit routine
	vector <float> componentSpec;
	vector <float> coeff_save;    //  save coefficients to check changes
	vector <float> elementCenterEnergy;	//	center energy of largest peak for this element
	spectrum.fit_vector( componentSpec, coeff_save, elementCenterEnergy );
	int nc_fit = coeff_save.size();

	vector <float> coeff(nc_fit);
	vector <float> var(nc_fit);
	float chisq;

    //  Prepare spectrum with non-fit components removed
    vector <float> fit_spectrum = spectrum.meas();
    unsigned int ic;
    for( ic=0; ic<spectrum.numberOfComponents(); ic++ ) {
        if( ! spectrum.component(ic).enabled ) continue;
        if( ! spectrum.component(ic).fit ) {
            const SpectrumComponent &non_fit_component = spectrum.component(ic);
            //  Subtract non-fit component from measured spectrum
            unsigned int is;
            for( is=0; is<non_fit_component.spectrum.size(); is++ ) {
                fit_spectrum[is] -= non_fit_component.coefficient * non_fit_component.spectrum[is];
            }
        }
    }
    //  Use least squares to find coefficient values that best fit components to spectrum
    int result = 0;
    result = lfit ( fit_spectrum, spectrum.sigma(), coeff, var, chisq, componentSpec, nChan );
	if( result < 0 ) return -810 + result;
    //  **** DEBUG **** Code to write coefficients to output
//    cout << "quantFitSpectrum    fit coeff";
//    for( int ii=0; ii<coeff.size(); ii++ ) cout << "  " << ii << " " << coeff[ii];
//    cout << endl;

/*
    //  Some debug code to catch zero (or nan) intensities or coefficients that appear during fit
    //  Also numbers that are too large or too small and make fit unstable
	int j;
	for( j=0; j<nc_fit; j++ ) {
        //  Check for zero (or nan) components
        float sum = 0;
        int i;
        for( i=0; i<nChan; i++ ) sum += componentSpec[ nChan * j + i];
//        if( sum == 0 || isnan( sum ) ) {
            logger << "fit vector entry " << j;
            logger.setf( ios::scientific, ios::floatfield );
            logger << " sum " << sum << endl;
            logger.setf( ios::fixed, ios::floatfield );
//        }
//        if( isnan( coeff[j] ) ) {
            logger << "fit vector entry " << j;
            logger.setf( ios::scientific, ios::floatfield );
            logger << " is zero or not a number.   " << coeff[j] << "  " << coeff_save[j] << endl;
            logger.setf( ios::fixed, ios::floatfield );
//        }
	}
*/
	//  Update coefficients (also updates fit and calculates chi squared)
	result = spectrum.update_coefficients( coeff, var );
	if( result < 0 ) return -820 + result;

	//  See if coefficients have not changed to within specified delta (XRFcontrols.h)
	bool done = true;
    int ic_fit;
	for( ic_fit=0; ic_fit<coeff.size(); ic_fit++ ) {
        if( elementCenterEnergy[ic_fit] == 0 ) continue;    //  Skip if not an element component (with a reasonable peak at a known energy)
        if( fabs( ( coeff[ic_fit] - coeff_save[ic_fit] ) / coeff_save[ic_fit] ) > FIT_COEFF_DELTA ) {
            done = false;
//            logger << "Conv check    en " << elementCenterEnergy[ic_fit] << "  coeff " << coeff[ic_fit];
//            logger << "   diff " << ( coeff[ic_fit] - coeff_save[ic_fit] ) / coeff_save[ic_fit] << "  threshold " << FIT_COEFF_DELTA << "    done " << done << endl;
        }
	}
    if( done ) return 0;

    //  Check to see if adjustments are to be done
    if( !spectrum.adjust_energy() && !spectrum.adjust_width() ) return 1;

        //			adjust energy calibration to get good fits
        //			(necessary for accurate net intensities and quantification)

    //      Calculate peak shifts and width change of individual functions in fit
    //  storage for peak shifts, center channels, peak widths, center energy, and weights
    vector <float> peak_shift;
    vector <float> peak_channel;
    vector <float> peak_sigma;
    vector <float> peak_width;
    vector <float> peak_ref_energy;
    vector <float> peak_width_sigma;
    float sumWgt = 0;
    float sumWgtWidth = 0;
	float nominal_resolution = conditions_in.detector.resolution();
    float det_ref_energy = conditions_in.detector.fwhm_energy();
    float ev_ch = spectrum.calibration().energyPerChannel();
	for ( ic_fit=0; ic_fit<nc_fit; ic_fit++ ) {
        if( elementCenterEnergy[ic_fit] == 0 ) continue;    //  Skip if not an element component (with a reasonable peak at a known energy)
        float elementResolution = conditions_in.detector.resolution( elementCenterEnergy[ic_fit] );
		vector <float> deriv (nChan,0);
		int is;
        for ( is=0; is<nChan; is++ ) deriv[is] = coeff[ic_fit] * componentSpec[ic_fit*nChan + is] / ev_ch;
        differentiate( deriv );
		vector <float> deriv2 (nChan,0);
        for ( is=0; is<nChan; is++ ) deriv2[is] = deriv[is] / ev_ch;
        differentiate( deriv2 );
        float sumRD = 0;
        float sumDD = 0;
        float sumR2D = 0;
        float sum2D2 = 0;
        for ( is=0; is<nChan; is++ ) {
            float e = spectrum.energy( is );
            float residual = spectrum.residual()[is];
            sumRD += residual * deriv[is];
            float d2 = deriv[is] * deriv[is];
            sumDD += d2;
            //  Only include points near peak of strongest line for width calculation
            if( fabs( e - elementCenterEnergy[ic_fit] ) < elementResolution / 4 ) {
                sumR2D += residual * deriv2[is];
                sum2D2 += deriv2[is] * deriv2[is];
            }
        }
        //  Collect LSQ sums for linear fit to shifts
        float center_channel = spectrum.channel( elementCenterEnergy[ic_fit] );
        float shift = sumRD / sumDD;
        float wgt = coeff[ic_fit] * coeff[ic_fit] / var[ic_fit];
        float sigma =  sqrt( var[ic_fit] / ( coeff[ic_fit] * coeff[ic_fit] ) ); //  relative standard deviation of peak area
//        float wgt = 1;
//        logger << "en cal corr  " << elementCenterEnergy[ic_fit] << "  " << sumRD;
//        logger << "  " << sumDD << "  " << shift;
//        logger << "  " << center_channel << "  " << sigma << "  " << wgt << "  " << elementResolution / 4;

        if( sigma > 0 && !isnan(sigma) && fabs( shift ) < elementResolution / 4 ) {
            peak_shift.push_back( shift );
            peak_channel.push_back( center_channel );
            peak_sigma.push_back( sigma );
//            logger << "    included";
        }
//        logger << endl;
        sumWgt += wgt;
        //  Accumulate sums for width linear fit to sqrt(energy)  [Fano effect]
        //  Skip if energy shift is too large (width calculation is unreliable)
        if( sigma > 0 && !isnan(sigma) && fabs( shift ) < elementResolution / 3 ) {
            //float rootEn = sqrt( elementCenterEnergy[ic_fit] );
            float fwhm_increase = SQRT_EIGHT_LN_2 * sumR2D / sum2D2 / elementResolution;
//            logger << "res corr  " << rootEn << "  " << sumR2D;
//            logger << "  " << sum2D2 << "  " << elementResolution << "  " << fwhm_increase;
            if( fabs( fwhm_increase ) < elementResolution / 4 ) {
                float width = elementResolution + fwhm_increase;
                peak_width.push_back( width * width );
                peak_ref_energy.push_back( elementCenterEnergy[ic_fit] - det_ref_energy );
                peak_width_sigma.push_back( sigma );
                sumWgtWidth += wgt;
//                logger << "    included";
//                logger << "  " << width << "  " << elementCenterEnergy[ic_fit] - det_ref_energy << "  " << width * width << "  " << wgt << "  " << sigma;
          }
//            logger << endl;
        }

	}

    //  fit energy shift vs channel number of peaks to straight line with weighting
	float offset_change = 0;
	float slope_change = 0;
    float siga, sigb, chi2, q;
    if( peak_channel.size() > 1 ) {
        fit( peak_channel, peak_shift, peak_sigma, offset_change,
                slope_change, siga, sigb, chi2, q);
    } else if( peak_channel.size() > 0 ) {
        //  If there is only one channel, modify the slope for best fit and ignore offset
        offset_change = 0;
        slope_change = peak_shift[0] / peak_channel[0];
    }
    //  fit width change squared vs energy to straight line with weighting and energy offset
    //  Be sure most of the main peaks are included, not just some very small peaks
    float fwhm_res, fwhm_slope;
    if( ( sumWgtWidth / sumWgt > 0.80f ) && ( peak_ref_energy.size() > 1 ) ) {
        fit( peak_ref_energy, peak_width, peak_width_sigma, fwhm_res,
                fwhm_slope, siga, sigb, chi2, q);
    } else if( ( sumWgtWidth / sumWgt > 0.80f ) && ( peak_ref_energy.size() > 0 ) ) {
        //  If there is only one channel, modify the resolution for best fit and ignore Fano factor
        //  Get value referenced to Mn Ka
        fwhm_res = sqrt( peak_width[0] ) / conditions_in.detector.resolution( peak_ref_energy[0] + det_ref_energy ) * nominal_resolution;
        fwhm_res *= fwhm_res;
        //  Keep same value for Fano factor by calculating slope using current values
        fwhm_slope = conditions_in.detector.fano() * conditions_in.detector.energy_per_pair() * EIGHT_LN_2;      //  8 ln2, defined in XRFconstants.h
    } else {
        fwhm_res = 0;
        fwhm_slope = 0;
    }
    if( fwhm_res < 0 ) fwhm_res = 0;
    if( fwhm_slope < 0 ) fwhm_slope = 0;
    fwhm_res = sqrt(fwhm_res);
//    logger << "specFit  en cal fit  offset " << offset_change << "  slope " << 100 * slope_change / ev_ch << "%" << endl;

    //  make sure corrections are not larger than a fraction of the peak width, otherwise they are unreliable
    if( fabs( offset_change * 4 / nominal_resolution ) < 2 && fabs( slope_change / ev_ch * 100 ) < 1 ) {
//        logger << "applying energy corr    offset " << offset_change << "  slope " << 100 * slope_change / ev_ch << "%" << endl;
        if( spectrum.adjust_energy() ) {
            spectrum.calibration_change().offset( spectrum.calibration().offset() + offset_change * 0.8 );    //  reducing corrections improves convergence, since effects are nonlinear
            spectrum.calibration_change().tilt( spectrum.calibration().tilt() + slope_change * 0.8 );    //  reducing corrections improves convergence, since effects are nonlinear
        }
        //  if energy calibration is good, attempt to correct peak widths in fit
        if( fabs( offset_change * 4 / nominal_resolution ) < 0.2f && fabs( slope_change * nChan * 2 / nominal_resolution ) < 0.2f ) {
            //  Make sure resolution changes are not too large (can be caused by unstable fits to only a few peaks)
            float old_fwhm_res = conditions_in.detector.resolution();
            float old_Fano = conditions_in.detector.fano();
            float new_Fano = fwhm_slope / conditions_in.detector.energy_per_pair() / EIGHT_LN_2;
            float rel_diff_res = fabs( fwhm_res - old_fwhm_res ) / old_fwhm_res;
            float rel_diff_Fano = fabs( new_Fano - old_Fano ) / old_Fano;
//            logger << "specFit  res fit  old: r " << old_fwhm_res << " F " << old_Fano << "  new: r " << fwhm_res << " F " << new_Fano << "  diff%: r " << 100*rel_diff_res << " F " << 100*rel_diff_Fano << endl;
            if( spectrum.adjust_width() && fwhm_res != 0 && new_Fano != 0 && rel_diff_res < 0.2f && rel_diff_Fano < 0.4f ) {   //  Apr. 18, 2018  was rel_diff_Fano < 0.2f
//                logger << "applying corr    res " << fwhm_res << "  Fano " << new_Fano << endl;
                conditions_in.detector.setResolution( fwhm_res );
                conditions_in.detector.fano( new_Fano );
            }
        }
    } else {
        logger << "Energy corrections too large:  offset " << offset_change << "  eVch (%) " << slope_change / ev_ch * 100 << endl;
    }

	return 1;
}
