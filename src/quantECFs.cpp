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

//#include <iostream>
//#include <iomanip>
#include "quantECFs.h"
#include "XRFconstants.h"
#include "XRFcontrols.h"
#include "interp.h"

//  Written April 3, 2017
//      Use the information in the calibration file
//      to find the element calibration factors for the input list of elements
//  Modified Dec. 23, 2017
//      Add use of standards list with weights if available
//      Reorganize for future changes to improve ECF calculation and interpolation for missing elements
//  Modified July 25, 2018
//      Check for zero ECF from csv calibration file to prevent zero final ECF for any element
//  Modified Nov. 4, 2019
//      Calculate an uncertainty for each ECF using the standard deviation across the set of fit coefficients for the standards
//      Calculate a weighted average certificate uncertainty associated with each ECF and include it in the ecf_error entry
//      Compute weighted average of ECF fit errors and use larger of that or ECF standard deviation as ECF error
//  Modified Nov. 25, 2020
//      Re-format the output for better alignment with column headings, add commas between numbers
//  Modified Jan. 4, 2021
//      Major changes - use relative errors and standard deviations for error root-square-sum
//      Always include certificate uncertainty and larger of deviation from average or coefficient relative error
//          (coefficient relative error now includes the residual error if calibration run today or later)
//      Add interpolation vs atomic number for missing ECFs, separately for K, L, and M lines
//      Skip standards from csv file with zero or negative ECFs
//  Modified Feb. 2, 2021   Avoid evaluating a standard with itself as a calibration standard
//  Modified Mar. 12, 2021  Fix bug interpolating ECF when only one entry in interpolation list
//  Modified Apr. 5, 2021   Ignore cal file and set all ECFs to unity (errors will be calculated in quantWriteResults)
//  Modified June 9, 2021   Modify ECF interpolation vs Z so it doesn't extrapolate outside of available values (takes end value)


using namespace std;

int quantECFs( const std::vector <StandardInformation> &standards_in,
                const vector <Element> &cal_element_list, const vector <float> &cal_factors_list,
                const vector <Element> &unk_element_list, const vector <float> &unk_fractions,
                vector <float> &unk_factors_list, vector <float> &unk_ECF_rel_err_list, std::ostream &logger ) {

    //  Calculate the Element Calibration Factor (ECF) for each element in the unknown
    //  Use the standards list with weighting factors if it is populated from a csv calibration file
    //  Otherwise use just the list of elements and factors from the txt calibration file
    //  Include the unknown composition in the arguments list in case we want composition-specific ECFs in the future
    //  Set up for interpolation vs Z for missing ECFs
    vector <float> interp_Z_K;
    vector <float> interp_ECF_K;
    vector <float> interp_Uncert_K;
    //  Set up averages for each element in the unknown and an overall average across all elements
    vector <float> avgFitCoeff( unk_element_list.size(), 0 );
    vector <float> avgWeight( unk_element_list.size(), 0 );
//    vector <float> fitCoeffErrs( unk_element_list.size(), 0 );
//    vector <float> fitCoeffSDs( unk_element_list.size(), 0 );
//    vector <float> avgUncert( unk_element_list.size(), 0 );
    vector <float> totalError( unk_element_list.size(), 0 );

    //  If the calibration file was read in, use it to determine the ECFs.  If not, set them all to unity.

    bool all_unity = true;
    float overallFitCoeff = 0;
    float overallWeight = 0;
    float overallErr = 0;
    if( standards_in.size() > 0 || cal_element_list.size() > 0 ) {
        //		find ECF and calculate average fit coefficients for each element
        unsigned int ie;
        for ( ie=0; ie<unk_element_list.size(); ie++ ) {
            if( standards_in.size() <= 0 ) {
                //  Loop over all of the elements in the ECF list from txt cal file
                unsigned int il;
                for ( il=0; il<cal_element_list.size(); il++ ) {
                    if ( cal_element_list[il] == unk_element_list[ie]
                        && cal_factors_list[il] > 0 ) {
                        float weight = 1;
                        avgFitCoeff[ie] += cal_factors_list[il] * weight;
                        avgWeight[ie] += weight;
                        if( cal_factors_list[il] != 1 ) all_unity = false;
                    };
                };
            } else {
                //  Loop over all standards in list
                unsigned int is;
                for( is=0; is<standards_in.size(); is++ ) {
                    //  Loop over all of the elements in this standard
                    unsigned int i;
                    for ( i=0; i<standards_in[is].element_list.size(); i++ ) {
                        if( standards_in[is].disable ) continue;
                        if( ! ( standards_in[is].element_list[i].element == unk_element_list[ie] ) ) continue;
                        if( standards_in[is].element_list[i].ecf <= 0 ) continue;
                        float weight = standards_in[is].element_list[i].weight;
                        float trial_ecf = standards_in[is].element_list[i].ecf;
                        if( trial_ecf > 0 ) {
                            avgFitCoeff[ie] += trial_ecf * weight;
                            avgWeight[ie] += weight;
                            if( trial_ecf != 1 ) all_unity = false;
                        }
                    };
                }
            }
            if ( avgWeight[ie] > 0 ) {
                overallFitCoeff += avgFitCoeff[ie];
                overallWeight += avgWeight[ie];
                avgFitCoeff[ie] /= avgWeight[ie];
            };
            //  If a CSV standards file was used, calculate the uncertainties for the ECFs
            if( standards_in.size() > 0 ) {
                //  Loop over all standards in list
                unsigned int is;
                for( is=0; is<standards_in.size(); is++ ) {
                    if( standards_in[is].disable ) continue;
                    //  Loop over all of the elements in this standard
                    unsigned int i;
                    for ( i=0; i<standards_in[is].element_list.size(); i++ ) {
                        if( ! ( standards_in[is].element_list[i].element == unk_element_list[ie] ) ) continue;
                        float weight = standards_in[is].element_list[i].weight;
                        if( weight <= 0 ) continue;
                        float trial_ecf = standards_in[is].element_list[i].ecf;
                        if( trial_ecf > 0 ) {
                            //  Use all relative errors in root square sum calculation of total uncertainty for this ECF
                            float fitError = standards_in[is].element_list[i].ecf_sigma / 100;  //  Convert percent to fraction
                            fitError += fitError * fitError * weight;
                            //  Use larger of fit error or deviation of individual ECF from this standard to the average ECF
                            float diff = trial_ecf - avgFitCoeff[ie];
                            diff += diff * diff * weight;
                            if( diff > fitError ) fitError = diff;
                            //  Relative uncertainty in given value, from certificate, in percent
                            float givenUncertainty = standards_in[is].element_list[i].uncertainty / 100;
                            totalError[ie]  += fitError + givenUncertainty * givenUncertainty * weight;
                            //  Old error calculation, commented out, kept for reference and in case the new method proves unworkable
                            //  Convert fit error to absolute error for root square sums
    /*                        float fitError = standards_in[is].element_list[i].ecf_sigma * avgFitCoeff[ie];    //  Uncertainty in fit coefficient, in percent
                            fitCoeffErrs[ie] += fitError * fitError * weight;
                            float diff = trial_ecf - avgFitCoeff[ie];
                            fitCoeffSDs[ie] += diff * diff * weight;
                            //  Relative uncertainty in given value, from certificate, in percent
                            //  Convert to absolute uncertainty in ECF for use as part of ECF error
                            float givenUncertainty = ( standards_in[is].element_list[i].uncertainty / 100 ) * avgFitCoeff[ie];
                            avgUncert[ie] += givenUncertainty * givenUncertainty * weight;

    */                    }
                    };
                }
                if ( avgWeight[ie] > 0 ) {
                    //  Use larger of weighted average fit error or weighted standard deviation in total ECF error
    //                totalError[ie] = fitCoeffSDs[ie];
    //                if( fitCoeffErrs[ie] > totalError[ie] ) totalError[ie] = fitCoeffErrs[ie];
    //                totalError[ie] += avgUncert[ie]; //  Add in certificate uncertainty (still squared at this point)
                    overallErr += totalError[ie];
                    //  Convert to relative errors
    //                fitCoeffErrs[ie] = sqrt( fitCoeffErrs[ie] / avgWeight[ie] ) / avgFitCoeff[ie];
    //                fitCoeffSDs[ie] = sqrt( fitCoeffSDs[ie] / avgWeight[ie] ) / avgFitCoeff[ie];
    //                avgUncert[ie] = sqrt( avgUncert[ie] / avgWeight[ie] ) / avgFitCoeff[ie];
    //                totalError[ie] = sqrt( totalError[ie] / avgWeight[ie] ) / avgFitCoeff[ie];
                    totalError[ie] = sqrt( totalError[ie] / avgWeight[ie] );
                }
            }
            //  Save information for interpolation vs atomic number for missing EFs (leave out zero ECFs so interpolation works OK)
            if( avgFitCoeff[ie] > 0 ) {
                interp_Z_K.push_back( unk_element_list[ie].Z() );
                interp_ECF_K.push_back( avgFitCoeff[ie] );
                interp_Uncert_K.push_back( totalError[ie] );
            }
        };
        //		use the overall average if no fit coefficient available from standards
        if ( overallWeight > 0 ) {
            overallFitCoeff /= overallWeight;
    //        overallErr = sqrt( overallErr / overallWeight ) / overallFitCoeff;
            overallErr = sqrt( overallErr / overallWeight );
        } else {
            overallFitCoeff = 1;
        };
    /*    if ( overallWeight > 0 ) {
            for ( ie=0; ie<unk_element_list.size(); ie++ ) if ( avgWeight[ie] <= 0 ) {
                avgFitCoeff[ie] = overallFitCoeff;
                totalError[ie] = overallErr;
            }
        } else {
            for ( ie=0; ie<unk_element_list.size(); ie++ ) avgFitCoeff[ie] = 1;
        }
    */
        //  Exchange sort into Z order for interpolation function
        for ( ie=0; ie<interp_Z_K.size()-1; ie++ ) {
            int i2;
            for( i2=ie+1; i2<interp_Z_K.size(); i2++ ) {
                if( interp_Z_K[ie] > interp_Z_K[i2] ) {
                    float temp = interp_Z_K[ie];
                    interp_Z_K[ie] = interp_Z_K[i2];
                    interp_Z_K[i2] = temp;
                    temp = interp_ECF_K[ie];
                    interp_ECF_K[ie] = interp_ECF_K[i2];
                    interp_ECF_K[i2] = temp;
                    temp = interp_Uncert_K[ie];
                    interp_Uncert_K[ie] = interp_Uncert_K[i2];
                    interp_Uncert_K[i2] = temp;
                }
            }
        }
//        cout << "quantECF   Z interp  " << interp_Z_K.size() << endl;
//        for ( ie=0; ie<interp_Z_K.size(); ie++ ) cout << "  " << interp_Z_K[ie];
//        cout << endl;
//        for ( ie=0; ie<interp_Z_K.size(); ie++ ) cout << "  " << interp_ECF_K[ie];
//        cout << endl;
//        for ( ie=0; ie<interp_Z_K.size(); ie++ ) cout << "  " << interp_Uncert_K[ie];
//        cout << endl;

        //  Use linear interpolation vs Z to get missing ECFs
        for ( ie=0; ie<unk_element_list.size(); ie++ ) {
            if ( avgWeight[ie] > 0 )  continue;
            if( interp_Z_K.size() > 1 ) {
                //  Check for outside of available range
                if( unk_element_list[ie].Z() < interp_Z_K[0] ) {
                    avgFitCoeff[ie] = interp_ECF_K[0];
                    totalError[ie] = interp_Uncert_K[0];
                } else if( unk_element_list[ie].Z() > interp_Z_K[interp_Z_K.size()-1] ) {
                    avgFitCoeff[ie] = interp_ECF_K[interp_Z_K.size()-1];
                    totalError[ie] = interp_Uncert_K[interp_Z_K.size()-1];
                } else {
                    avgFitCoeff[ie] = interp( unk_element_list[ie].Z(), interp_Z_K, interp_ECF_K );
                    totalError[ie] = interp( unk_element_list[ie].Z(), interp_Z_K, interp_Uncert_K );
                }
            }
            //  Use the overall average ECF as a last resort
            if( avgFitCoeff[ie] <= 0 ) {
                avgFitCoeff[ie] = overallFitCoeff;
                totalError[ie] = overallErr;
            }
        }

    }   //  if( standards_in.size() > 0 || cal_element_list.size() > 0 )


    if( all_unity ) {
        //  Unity ECFs, errors will be calculated in quantWriteResults
        unsigned int ie;
        for ( ie=0; ie<avgFitCoeff.size(); ie++ ) avgFitCoeff[ie] = 1;

        //  Errors will actually be calculated in quantWriteResults using fit to weight percents
        logger << "All ECFs are unity and errors are from Elemental Calibration results on 30 standards." << endl;
        logger << endl;

    } else {    //  No standards read in, use unity ECF calibration (optic file adjusted to give unit value for all ECFs)

        logger << "Final element calibration factors and uncertainties for this unknown (overall ECF is ";
        logger.precision( 4 );
        logger << overallFitCoeff;
        logger.precision( 1 );
        logger << ", overall ECF relative error is " << 100 * overallErr << " %)" << endl;
        const int element_width = 12;
        logger << setw(element_width) << "Element";
        logger << setw(element_width) << ", ECF";
        logger << setw(element_width) << ", Total weight";
    //    logger << setw(element_width) << ", Fit Err";
    //    logger << setw(element_width) << ", Std Dev";
    //    logger << setw(element_width) << ", Given uncertainty";
        logger << setw(element_width) << ", ECF uncertainty (relative percent)";
        logger << endl;
        unsigned int ie;
        for ( ie=0; ie<unk_element_list.size(); ie++ ) {
                logger << setw(element_width) << unk_element_list[ie].symbol();
                logger.precision( 4 );
                logger << ",   " << setw(element_width) << avgFitCoeff[ie];
                logger.precision( 2 );
                logger << ",   " << setw(element_width) << avgWeight[ie];
                logger.precision( 1 );
    //            logger << ",   " << setw(element_width) << 100 * fitCoeffErrs[ie];
    //            logger << ",   " << setw(element_width) << 100 * fitCoeffSDs[ie];
    //            logger << ",   " << setw(element_width) << 100 * avgUncert[ie];
                logger << ",   " << setw(element_width) << 100 * totalError[ie];
                logger << endl;
        };
        logger << "    ECF uncertainty is included in total quant error and includes standard fitting errors and certificate uncertainties" << endl;
        logger << endl;

    }


    //  Copy the results to the output vector
    unk_factors_list.resize(avgFitCoeff.size());
    unk_ECF_rel_err_list.resize(totalError.size());
    int ie;
    for ( ie=0; ie<avgFitCoeff.size(); ie++ ) unk_factors_list[ie] = avgFitCoeff[ie];
    for ( ie=0; ie<totalError.size(); ie++ ) unk_ECF_rel_err_list[ie] = totalError[ie];

	return 0;
};
