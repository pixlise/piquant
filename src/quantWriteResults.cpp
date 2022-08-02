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

#include "quantWriteResults.h"
#include "Lfit.h"
#include "Fit.h"
#include "differentiate.h"
#include "XRFconstants.h"
#include "XRFcontrols.h"
#include "XRFutilities.h"
#include "spline.h"


//		Write out results for standards or unknowns

//  Written Mar. 16, 2017
//      Based on HistogramAnalysisPackage_Subprocess_PIXL_Jun2016.cpp of Jan. 2, 2017
//  Modified June 7, 2017 to use componentDescription that returns string directly
//  Modified Dec. 10, 2017
//      Use a local boolean array to keep track of which components shown with quantified elements
//  Modified Mar. 2, 2018
//      Use number of iterations stored in spectrum
//  Modified Sep. 19, 2018
//      Correct counts in range 1-7.25 keV (was total spectrum counts)
//      Write counts even if Si and Fe intensity is zero
//  Modified May 16, 2019
//      Add counts in 1-7.25 keV region to quant map outputs, move from quantWriteResults to XraySpectrum
//      Remove negative intensities from map and log output (leave coefficients as-is for now)
//  Modified May 21, 2019
//      Implement Evaluate action, add given and relative error vs given, from element list (if any given > 0)
//  Modified June 6, 2019
//      Add output of oxide oxygen percent from oxide ratio of matrix elements
//  Modified July 22, 2019
//      Change wording for oxygen from matrix element oxides
//  Modified Nov. 4, 2019
//      Add total error to element list, to include ECF standard deviation (and later certificate uncertainty)
//  Modified Dec. 16, 2019
//      Element sum added to arguments so the calculation here can be used in quantWriteMap
//  Modified Nov. 24, 2020
//      Change eval Boolean to argument for evaluation sub-command
//      Put relative error vs given into element list
//      Move matrix effect factor from spectrum component to element list
//  Modified Dec. 1, 2020
//      For non-fit components, write indicator instead of relative coefficient error
//  Modified Jan. 4, 2021
//      Fix bug in that caused element list to be populated incorrectly when component was disabled
//      Include new residual error from component into fit error for coefficient
//      Increase total error to 2-sigma (multiply error by 2)
//  Modified Jan. 4, 2021
//      Write error vs given even for disabled components if given != 0 (and add to element list for eval output)
//  Modified Apr. 5, 2021   Disable spectrum residual error estimates, use fixed uncertainties vs given from elemental calibration (mag in quantECFs)
//  Modified Apr. 28, 2021  New error values after modifications made today, significant improvements
//  Modified June 8, 2021   Change error reporting to absolute (was relative)  quantWriteResults.cpp
//  Modified June 9, 2021   Change Z range for mid-Z trace elements to eliminate Co (too much interference from Fe beta peak)
//                          Change interpolation range to better match actual error performance
//                          Change oxides and element reporting to match team wishes (e-mail from Joel 6/7/2021, 1:27 PM)
//  Modified June 11, 2021  Final errors from surface ops calibration
//  Modified June 20, 2021  Final errors from surface ops calibration with shelf bug fixed
//  Modified July 10, 2021  Add simple pulse pileup calculation and turn on Compron escape, new error values

using namespace std;

//  Based on results from PIXL Flight Model, Elemental Calibration of May 23, 2019
//  Results from file "Eval_speedup_GlassPowderScap_Mar2021.csv"
const vector <float> ERROR_GIVEN = {    0.0, 0.05, 0.5, 5.0, 100 };
const vector <float> ERROR_RELATIVE = { 298, 126, 36,   5,   5 };
const vector <float> ERROR_SPLINE( ERROR_GIVEN.size(), 0 ); //  Linear interpolation
const int TRACE_MID_Z_LO = 28;          const int TRACE_MID_Z_HI = 42;
const vector <float> TRACE_MID_Z_ERROR_GIVEN = {   0.0, 0.05, 0.5, 5.0, 100 };
const vector <float> TRACE_MID_Z_ERROR_RELATIVE = { 298,  40, 36,   5,   5 };
const vector <float> TRACE_MID_Z_ERROR_SPLINE( ERROR_GIVEN.size(), 0 ); //  Linear interpolation
const int RARE_EARTH_Z_LO = 57;          const int RARE_EARTH_Z_HI = 71;
const vector <float> RARE_EARTH_ERROR_GIVEN = {    0.0, 0.05, 0.5, 5.0, 100 };
const vector <float> RARE_EARTH_ERROR_RELATIVE = {  298,  79, 36,   5,   5 };
const vector <float> RARE_EARTH_ERROR_SPLINE( ERROR_GIVEN.size(), 0 ); //  Linear interpolation

bool display_as_pure_element( const Element el );   //  Helper function to control formatting of quant output for geologists

int quantWriteResults( const XrayMaterial &material, const XrayDetector &detector, std::vector <ElementListEntry> &element_list,
                const XraySpectrum &spectrum, const bool oxidesOutput, ostream &termOutFile, float &element_sum, bool eval ) {

    //  new energy calibration and detector resolution for best fit
    termOutFile << endl;
    termOutFile << endl;
    termOutFile.precision( 2 );
    termOutFile << "Fit results after " << spectrum.iterations() << " iterations, reduced chi sq = " << spectrum.chisq();
    termOutFile << "         live time " << spectrum.live_time() << " sec." << endl;
    termOutFile << "Final energy calibration (eV): ";
    termOutFile.precision(1);
    termOutFile << "  eV start = " << spectrum.calibration().energyStart();
    termOutFile.precision(4);
    termOutFile << "  eV/ch = " << spectrum.calibration().energyPerChannel();
    termOutFile << "  detector resolution (eV): ";
    termOutFile.precision( 0 );
    termOutFile  << detector.resolution()<< "  (at " << detector.fwhm_energy() << " eV)";
    termOutFile.precision( 3 );
    termOutFile << "  fano = " << detector.fano();
    termOutFile << endl;
    termOutFile.precision( 2 );
    termOutFile << "      Energy correction offset " << spectrum.calibration().offset();
    termOutFile << " eV   slope change " << 100 * spectrum.calibration().tilt() / spectrum.calibration().energyPerChannel() << " %" << endl;
    if( material.thickness() > 0 )
        termOutFile << "Specimen thickness=" << material.thickness() << " cm,  density=" << material.density() << " gm/cm3" << endl;
    termOutFile << endl;

    const int output_width = 8;

    //  PIXL L5 requirements info for X-ray Subsystem
    const Element si( 14 );    //   Si, Z=14
    const Element fe( 26 );     //  Fe, Z=26
    float si_int = 0;
    float fe_int = 0;
    int ic;
    for( ic=0; ic<spectrum.numberOfComponents(); ic++ ) {
        if( spectrum.component( ic ).element == si && spectrum.component( ic ).level == K )
                si_int = spectrum.component( ic ).intensity;
        if( spectrum.component( ic ).element == fe && spectrum.component( ic ).level == K )
                fe_int = spectrum.component( ic ).intensity;
    };
    termOutFile << "XRS L5 requirements info: " << endl;
    termOutFile.precision(2);
    termOutFile << "  L5-XRS-03    total count rate (" << spectrum.region_start()/1000 << " keV";
    termOutFile << " to " << spectrum.region_end()/1000 << " keV)";
    termOutFile.precision(0);
    termOutFile << " = " << spectrum.region_counts() / spectrum.live_time() << " cps (correct energy range)" << endl;
    if( si_int > 0 && fe_int > 0 ) {
        termOutFile << "  L5-XRS-10    Si intensity = " << si_int;
        termOutFile << ", Fe intensity = " << fe_int;
        termOutFile.precision( 2 );
        termOutFile << ",   Si/Fe ratio = " << si_int / fe_int << endl;
        termOutFile << endl;
    };

    //	write composition, fit coefficients, intensities, and errors from fit to standard or unknown
    termOutFile.precision(2);
    termOutFile << "Fitted elements" << endl;
    float checkSum = 0;
    const vector <Element> &elements = material.element_list();
    //  Keep track of which components are shown with quantified eleemnts
    vector <bool> shown_component_indices( spectrum.numberOfComponents(), false );
    unsigned int ie;
    //  First process elements in material that have associated quantification components
    for ( ie=0; ie<elements.size(); ie++ ) {
        //  Find the spectrum fit component that was used to quantify this element
        int ic_element = spectrum.index( elements[ie] );
        if( ic_element < 0 || spectrum.coefficient( ic_element ) == COEFFICIENT_NO_COMPONENT ) continue;    //  no component, don't list here
        shown_component_indices[ic_element] = true;
        //  Find the corresponding entry in the input element list
        int element_list_index = -1;
        unsigned int ie_list;
        for( ie_list=0; ie_list<element_list.size(); ie_list++ ) {
            if( !( elements[ie] == element_list[ie_list].element ) ) continue;
            //  Be sure this element list entry corresponds to a quantified element
            if( element_list[ie_list].qualifier == IGNORE
                    || element_list[ie_list].qualifier == EXCLUDE
                    || element_list[ie_list].qualifier == MATRIX ) continue;
            //  Be sure the quant level for this element list entry matches the spectrum fit component
            //  componentQuantLevel converts SpectrumComponent enum to corresponding ElementQuantLevel enum
            if( element_list[ie_list].quant_level != NO_QUANT_LEVEL
                && componentQuantLevel( spectrum.component( ic_element ) ) != element_list[ie_list].quant_level ) continue;
            element_list_index = ie_list;
            break;
        }
        //  Output the information for this element
        termOutFile.precision(4);
        float percent;
        if( !oxidesOutput || display_as_pure_element( elements[ie] ) ) {
            termOutFile << "   " << setw(5) << elements[ie].symbol();
            float fraction = material.fraction( elements[ie] );
            termOutFile << "  " << setw(output_width) << 100 * fraction << " %";
            checkSum += fraction;
            percent = 100 * fraction;
            //  Correct the check sum for light element fraction if this was not a pure element but was displayed as one
            if( display_as_pure_element( elements[ie] ) && element_list[element_list_index].stoichiometry.formula != PURE_ELEMENT ) {
                checkSum += material.fraction_light( elements[ie] );
            }
        } else {
            // label oxides
            termOutFile << setw(output_width) << material.formula_string( elements[ie] );
            //  fractions as oxides
            float oxideFraction = material.fraction_formula( elements[ie] );
            termOutFile.precision(2);
            termOutFile << "  " << setw(output_width) << 100 * oxideFraction << " %";
            checkSum += oxideFraction;
            percent = 100 * oxideFraction;
            //  Convert given percent from element percent to oxide percent if necessary (to correctly calculate error vs given)
            if( !element_list[element_list_index].stoichiometry.input_fractions_are_formula ) {
                element_list[element_list_index].given = 100 * XrayMaterial::calculate_fraction_formula( elements[ie], element_list[element_list_index].given/100, material.stoichiometry( elements[ie] ) );
                element_list[element_list_index].stoichiometry.input_fractions_are_formula = true;
            }
        }
        if( element_list_index >= 0 ) {
            //  Put the results of the quantification into the element list for calibration and map files
            element_list[element_list_index].percent = percent;
            element_list[element_list_index].stoichiometry = material.stoichiometry( elements[ie] );
            //  Fix the formula for this element if its fraction was for the pure element
            if( display_as_pure_element( elements[ie] ) ) element_list[element_list_index].stoichiometry.formula = PURE_ELEMENT;
            element_list[element_list_index].intensity = spectrum.intensity( ic_element );
            element_list[element_list_index].coefficient = spectrum.coefficient( ic_element );
            element_list[element_list_index].matrix = spectrum.component( ic_element ).matrix;
        }
//        float rsd = material.std_dev( elements[ie] );
//        termOutFile << "  rsd " << setw(8) << 100 * rsd;
        termOutFile << "  " << setw(output_width) << componentDescription( spectrum.component( ic_element ) );
        termOutFile.precision(1);
        if( spectrum.intensity( ic_element ) >= 0 ) termOutFile << "   int " << setw(output_width) << spectrum.intensity( ic_element );
        else termOutFile << "   int " << setw(output_width) << 0.0;
        termOutFile.precision(4);
        termOutFile << "  coeff " << setw(output_width) << spectrum.coefficient( ic_element );
        float fitError = 0;
        float ecfError = element_list[element_list_index].ecf_sigma;
        float totalError = 0;
        if( spectrum.component( ic_element ).enabled ) {
            fitError = sqrt( spectrum.variance( ic_element ) ) / spectrum.coefficient( ic_element );
//            if( spectrum.residual_error( ic_element ) > fitError ) fitError = spectrum.residual_error( ic_element );
            termOutFile.precision(1);
            termOutFile << "   re_c " << setw(output_width) << 100 * fitError << "%";
//            totalError = 2 * sqrt( fitError*fitError + ecfError*ecfError ); //  2-sigma error
//            termOutFile << "   tot_err(2s) " << setw(output_width) << 100 * totalError << "%";
            //  Estimated error from element calibration standards, error vs certificate values, expressed as a function of given percent
            float estimated_error = splint( ERROR_GIVEN, ERROR_RELATIVE, ERROR_SPLINE, 100*material.fraction_formula( elements[ie] ) ) / 100;
            if( TRACE_MID_Z_LO <= elements[ie].Z() && elements[ie].Z() <= TRACE_MID_Z_HI )
                estimated_error = splint( TRACE_MID_Z_ERROR_GIVEN, TRACE_MID_Z_ERROR_RELATIVE, TRACE_MID_Z_ERROR_SPLINE, 100*material.fraction_formula( elements[ie] ) ) / 100;
            if( RARE_EARTH_Z_LO <= elements[ie].Z() && elements[ie].Z() <= RARE_EARTH_Z_HI )
                estimated_error = splint( RARE_EARTH_ERROR_GIVEN, RARE_EARTH_ERROR_RELATIVE, RARE_EARTH_ERROR_SPLINE, 100*material.fraction_formula( elements[ie] ) ) / 100;
            totalError = sqrt( fitError*fitError + ecfError*ecfError + estimated_error*estimated_error );
//            termOutFile << "   tot_err(1s) " << setw(output_width) << 100 * totalError << "%";
            totalError *= percent / 100;  //  Convert from relative to absolute error
            termOutFile.precision(4);
            termOutFile << "   abs_err(1s) " << setw(output_width) << 100 * totalError << "%";
            termOutFile.precision(1);
//            termOutFile.precision(3);
//            if( spectrum.component( ic_element ).matrix > 0 ) termOutFile << "   M " << setw(output_width) << spectrum.component( ic_element ).matrix;
//            termOutFile << "   M " << setw(output_width) << spectrum.component( ic_element ).matrix;
            if( element_list_index >= 0 ) {
                //  Put the results of the quantification into the element list for calibration and map files
                element_list[element_list_index].intensity = spectrum.intensity( ic_element );
                element_list[element_list_index].coefficient = spectrum.coefficient( ic_element );
                element_list[element_list_index].rel_err_coeff = 100 * fitError;
                element_list[element_list_index].total_err = 100 * totalError;
            }
        } else {
            termOutFile << "   not included";
        };
        if( element_list_index >= 0 ) {
            //  Calculate relative error for quantification vs given value
            float rel_err_given = 0;
            if( element_list[element_list_index].given != 0 ) rel_err_given = 100 * ( element_list[element_list_index].percent - element_list[element_list_index].given ) / element_list[element_list_index].given;
            element_list[element_list_index].rel_err_given = rel_err_given;
            if( eval ) {
                termOutFile.precision(4);
                termOutFile << "  given " << setw(output_width) << element_list[element_list_index].given;
                termOutFile.precision(1);
                termOutFile << "  rel_err_vs_given " << setw(output_width) << rel_err_given;
            }
        }

        termOutFile.precision(1);
//        termOutFile << "   oxide " << 100*oxideFraction << "   ox " << 100*oxygen_fraction;   //  for debugging of oxide calculations
//        termOutFile << "   incl " << spectrum.component( ic_element ).included;
//        termOutFile << "   enable " << spectrum.component( ic_element ).enabled;
//        termOutFile << "   quant " << spectrum.component( ic_element ).quant;
        termOutFile << endl;
    };
    //  Now process any left-over elements without associated components or whose coefficients were zero
    termOutFile << "Matrix elements" << endl;
    for ( ie=0; ie<elements.size(); ie++ ) {
        int ic_element = spectrum.index( elements[ie] );
        if( ic_element >= 0 ) continue;    //  has peaks in spectrum, not a matrix element
        float fraction = material.fraction_formula( elements[ie] );
        checkSum += fraction;
        //  For oxides output, the contribution from carbon and oxygen is already summed in the fractions
        if( elements[ie].atomicNumber() == 6 && oxidesOutput ) checkSum -= material.added_carbon();
        if( elements[ie].atomicNumber() == 8 && oxidesOutput ) checkSum -= material.added_oxygen();
        termOutFile << "   " << setw(2) << material.formula_string( elements[ie], material.stoichiometry( elements[ie] ) );
        termOutFile.precision(4);
        termOutFile << "  " << setw(output_width) << 100 * fraction << " %";
        if( elements[ie].atomicNumber() == 6 ) {
            termOutFile.precision(2);
            termOutFile << "     (" << 100 * material.added_carbon() << " % from stoichiometry)";
        } else if( elements[ie].atomicNumber() == 8 ) {
            termOutFile.precision(2);
            termOutFile << "     (" << 100 * material.added_oxygen() << " % from stoichiometry)";
        } else if( material.stoichiometry( elements[ie] ).formula != PURE_ELEMENT ) {
            termOutFile << "     (" << 100 * material.fraction( elements[ie] ) << " %";
            termOutFile << " " <<  elements[ie].symbol() << ")";
        }
        termOutFile << endl;
    };
    termOutFile << endl;
    termOutFile.precision(2);
    termOutFile << "    Element sum " << checkSum * 100 << " %" << endl;
    //      Write out Compton and Rayleigh info
    termOutFile << endl;
//    termOutFile << "Compton  ";
//    termOutFile << net_Compton;
//    termOutFile << "        Rayleigh  ";
//    termOutFile << net_Rayleigh;
//    termOutFile << endl;

    //  Now list any components that were not used for quantification
    termOutFile << "Other spectrum components" << endl;
    for( ic=0; ic<spectrum.numberOfComponents(); ic++ ) {
        if( shown_component_indices[ic] ) continue;  //  Already listed in first group above
        termOutFile << "   " << ic;
        termOutFile << "  " << setw(output_width) << componentDescription( spectrum.component( ic ) );
        termOutFile.precision(1);
        termOutFile << "  int " << setw(output_width) << spectrum.component( ic ).intensity;
        termOutFile.precision(4);
        termOutFile << "  coeff " << setw(output_width) << spectrum.component( ic ).coefficient;
        if( spectrum.component( ic ).enabled ) {
            if( spectrum.component( ic ).fit ) {
                float fitError = sqrt( spectrum.variance( ic ) ) / spectrum.coefficient( ic );
                termOutFile.precision(2);
                termOutFile << "   re_c " << setw(output_width) << 100 * fitError << "%";
            } else {
                termOutFile << "   not fit";
                if( spectrum.component( ic ).non_fit_factor != 0 ) termOutFile << " (tracks quant component)";
            }
        } else {
            termOutFile << "   not included";
        };
//        termOutFile << "   incl " << spectrum.component( ic ).included;
//        termOutFile << "   enable " << spectrum.component( ic ).enabled;
//        termOutFile << "   quant " << spectrum.component( ic ).quant;
        termOutFile << endl;
    };

    termOutFile << endl;
    termOutFile << endl;
    element_sum = checkSum * 100;   //  Convert to percents for map and evaluate
	return 0;

};


bool display_as_pure_element( const Element el ) {
	const int max_Z_c = 101;
	//  Change oxides vs element reporting to match team wishes (e-mail from Joel 6/7/2021, 1:27 PM)
	//  Na2O, MgO, Al2O3, SiO2, P2O5, SO3, Cl, K2O, CaO, TiO2, Cr2O3, MnO, FeO-T
	//  11    12   13     14    15    16   17  19   20   22    24     25   26
	//  Everything else report as elements
	//  Note that these values are inverted
	//      (0=display as element, 1=display as oxide)
	const int oxides[max_Z_c] = { 0,	//	E-mail from Joel Hurowitz, 12/14/2020, 8:55 AM
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,	//	1-10    Zero means element does not form a carbonate
        1,  1,  1,  1,  1,  1,  0,  0,  1,  1,	//	11-20
        0,  1,  0,  1,  1,  1,  0,  0,  0,  0,	//	21-30
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,	//	31-40
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,	//	41-50
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,	//	51-60
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,	//	61-70
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,	//	71-80
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,	//	81-90
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0 	//	91-100
    };
    int temp_Z = el.Z();
    if( temp_Z > 0 && temp_Z < max_Z_c ) {
        return (oxides[temp_Z] == 0 );
    } else {
        return true;
    }
}
