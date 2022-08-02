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

#include "XrayMaterial.h"
#include <math.h>
#include <sstream>
#include "XRFconstants.h"
#include "toStringHelpers.h"


//  This class is intended to allow calculation of X-ray properties
//      of a material with arbitrary composition
//  Plan is for it to eventually replace all objects in the X-ray beam
//      such as windows, paths, filters, etc.
//      so they are not limited to a single element
//  Skeleton class written Dec. 22, 2016 for dust error calculations
//  Modified Dec. 26, 2016 to add default constructor and fix density bug
//  Modified March 29, 2017 to include add-element function (used in quantUnknown)
//      Also handle oxides individually using oxide ratio, remove global oxides flag
//  Modified Nov. 1, 2017 to add uncertainty in composition
//  Modified July 18, 2018 to add oxygen to element list when oxide ratios change from all zeros
//      Also add flag to specify element fractions as input, not oxide fractions (still handles oxide ratios)
//      Keep separate track of the oxygen added from oxides and return it via a function
//      Add function to return fraction input (without any modification)
//  Modified July 25, 2018
//      Prevent nan if sum of fractions is zero
//  Modified July 22, 2019
//      Change light element oxide ratios to zero (below Z=10)
//  Modified July 22, 2019
//      Add functions to return average Z, average A, and average Z/A
//      Add photoelectric cross section
//  Modified Nov. 24, 2019
//      Change some of the utility functions for oxide calculations to static public functions
//  Modified Nov. 24, 2020
//      Added individual element photoelectric cross sections (for detector shelf calculations)
//  Modified Nov. 14, 2020
//      Add capability for more than just oxygen to be included by stoichiometry (just carbonates for now)
//  Modified Jan. 7, 2020
//      Fix bug with formula_fraction
//  Modified June 9, 2021   Change oxides and element reporting to match team wishes (e-mail from Joel 6/7/2021, 1:27 PM)
//  Modified June 27, 2021  Add command line option to normalize element sum to 100% (or any value)
//  Modified July 9, 2021   Add command line option to change Fe oxide ratio (-Fe)


using namespace std;

//  Initialize static class member
float XrayMaterial::default_modified_iron_oxide_ratio = -1;


XrayMaterial::XrayMaterial() {
}

XrayMaterial::XrayMaterial(const Element &element_in, const bool oxides_in, const bool oxides_frac_flag ) {
    element_list_input.push_back( element_in );
    fractions_input.push_back( 1 );
    LightElements temp;
    if( oxides_in ) {
        temp.formula = OXIDE;
        temp.formula_ratio = default_oxide_ratio( element_in );
        temp.input_fractions_are_formula = oxides_frac_flag;
    }
    formula_info.push_back( temp );
    uncertainties.push_back( 0 );
    populate_element_list();
    calculate_element_fractions();
};

XrayMaterial::XrayMaterial(const std::vector <Element> &element_list_in, const std::vector <float> &fractions_in, const bool oxides_in, const bool oxides_frac_flag) {
    int ie;
    for( ie=0; ie<element_list_in.size(); ie++ ) {
        element_list_input.push_back( element_list_in[ie] );
        fractions_input.push_back( fractions_in[ie] );
        LightElements temp;
        if( oxides_in ) {
            temp.formula = OXIDE;
            temp.formula_ratio = default_oxide_ratio( element_list_in[ie] );
            temp.input_fractions_are_formula = oxides_frac_flag;
        }
        formula_info.push_back( temp );
        uncertainties.push_back( 0 );
    }
    populate_element_list();
    calculate_element_fractions();
}

XrayMaterial::XrayMaterial(const std::vector <Element> &element_list_in, const std::vector <float> &fractions_in, const std::vector  <LightElements> formula_info_in ) {
    LightElements default_formula_info;
    int ie;
    for( ie=0; ie<element_list_in.size(); ie++ ) {
        element_list_input.push_back( element_list_in[ie] );
        ie<fractions_in.size()?fractions_input.push_back( fractions_in[ie] ) : fractions_input.push_back( 0 );
        ie<formula_info_in.size()?formula_info.push_back( formula_info_in[ie] ) : formula_info.push_back( default_formula_info );
        uncertainties.push_back( 0 );
    }
    populate_element_list();
    calculate_element_fractions();
}

XrayMaterial::XrayMaterial(const std::vector <Element> &element_list_in, const std::vector <float> &fractions_in,
                const std::vector  <LightElements> formula_info_in, const std::vector  <float> uncertainties_in ) {
    LightElements default_formula_info;
    int ie;
    for( ie=0; ie<element_list_in.size(); ie++ ) {
        element_list_input.push_back( element_list_in[ie] );
        ie<fractions_in.size()?fractions_input.push_back( fractions_in[ie] ) : fractions_input.push_back( 0 );
        ie<formula_info_in.size()?formula_info.push_back( formula_info_in[ie] ) : formula_info.push_back( default_formula_info );
        ie<uncertainties_in.size()?uncertainties.push_back( uncertainties_in[ie] ) : uncertainties.push_back( 0 );
    }
    populate_element_list();
    calculate_element_fractions();
}

XrayMaterial::XrayMaterial(const int n_el, const int element_Z_in[], const float fractions_in[], const bool oxides_in, const bool oxides_frac_flag ) {
    int ie;
    for( ie=0; ie<n_el; ie++ ) {
        if( element_Z_in[ie] < 1 || element_Z_in[ie] > Element::maxZ() ) continue;
        Element temp_el( element_Z_in[ie] );
        element_list_input.push_back( temp_el );
        fractions_input.push_back( fractions_in[ie] );
        LightElements temp;
        if( oxides_in ) {
            temp.formula = OXIDE;
            temp.formula_ratio = default_oxide_ratio( temp_el );
            temp.input_fractions_are_formula = oxides_frac_flag;
        }
        formula_info.push_back( temp );
        uncertainties.push_back( 0 );
    }
    populate_element_list();
    calculate_element_fractions();
}


//      X-ray properties functions

float XrayMaterial::transmission( const float energy_in, const float csc ) const {
//      Calculate X-ray transmission at a given energy
	if ( m_thickness <= 0 || elements.size() <= 0 ) return 1;
	if ( energy_in <= 0 ) return 0;
	float mu_x = cross_section( energy_in );
	mu_x *= m_thickness * csc;
	float resp = 0;
	if( mu_x < EXP_FLOAT_TEST ) resp = exp ( -mu_x );
	return resp;
}

float XrayMaterial::absorption( const float energy_in, const float csc) const {
//      Calculate X-ray absorption at a given energy
    return 1 - transmission( energy_in, csc );
}

float XrayMaterial::cross_section( const float energy_in ) const {
//      Calculate X-ray absorption cross-section at a given energy (cm2/gm)
	if ( energy_in <= 0 ) return 0;
	float sigma = 0;
	float sum = 0;
	int i;
	//  Normalize composition so that absorption does not vary inappropriately
	for ( i=0; i<fractions.size(); i++ ) sum += fractions[i];
	for ( i=0; i<elements.size(); i++ ) {
		sigma += fractions[i] * absorption_tables[i].total( energy_in );
	};
	if( sum == 0 ) return 0;
	return sigma / sum;
};

float XrayMaterial::cross_section( const Element el, const float energy_in ) const {
    //  Returns the absorption cross-section for a single element at the given energy
	if ( energy_in <= 0 ) return 0;
	float sigma = 0;
    int ie = find_element( el, elements );
    if( ie >= 0 ) sigma = absorption_tables[ie].total( energy_in );
    return sigma;
};

float XrayMaterial::photo( const float energy_in ) const {
//      Calculate X-ray photoelectric cross-section at a given energy (cm2/gm)
	if ( energy_in <= 0 ) return 0;
	float sigma = 0;
	float sum = 0;
	int i;
	//  Normalize composition so that absorption does not vary inappropriately
	for ( i=0; i<fractions.size(); i++ ) sum += fractions[i];
	for ( i=0; i<elements.size(); i++ ) {
		sigma += fractions[i] * absorption_tables[i].photo( energy_in );
	};
	if( sum == 0 ) return 0;
	return sigma / sum;
};

float XrayMaterial::photo( const Element el, const float energy_in ) const {
//      Calculate X-ray photoelectric cross-section at a given energy (cm2/gm)
	if ( energy_in <= 0 ) return 0;
	float sigma = 0;
    int ie = find_element( el, elements );
    if( ie >= 0 ) sigma = absorption_tables[ie].photo( energy_in );
	return sigma;
};

const XrayXsectTable &XrayMaterial::cross_section_table( const Element el ) const {
    //  Returns the absorption cross-section table for a single element
    //  XrayXsectTable no_table;    Moved to private data to avoid passing a deleted object
    int ie = find_element( el, elements );
    if( ie < 0 ) {
    //  Return the default table
        return no_table;
    };
    return absorption_tables[ie];
};

float XrayMaterial::incoherent( const float energy_in, const float theta_in ) const {
//      Calculate X-ray incoherent scatter cross-section at a given energy and angle (cm2/gm)
	if ( energy_in <= 0 ) return 0;
	float sigma = 0;
	float sum = 0;
	int i;
	//  Normalize composition so that absorption does not vary inappropriately
	for ( i=0; i<fractions.size(); i++ ) sum += fractions[i];
	for ( i=0; i<elements.size(); i++ ) {
		sigma += fractions[i] * scatter_tables[i].incoherent( energy_in, theta_in );
	};
	if( sum == 0 ) return 0;
	return sigma / sum;
};

float XrayMaterial::incoherent( const float energy_in, const float theta_in, const float scattered_energy_in ) const {
//      Calculate X-ray incoherent scatter cross-section at a given energy and angle
//      This is the doubly-differential cross-section vs solid angle and energy  (cm2/gm/eV)
	if ( energy_in <= 0 ) return 0;
	float sigma = 0;
	float sum = 0;
	int i;
	//  Normalize composition so that absorption does not vary inappropriately
	for ( i=0; i<fractions.size(); i++ ) sum += fractions[i];
	for ( i=0; i<elements.size(); i++ ) {
		sigma += fractions[i] * scatter_tables[i].incoherent( energy_in, theta_in, scattered_energy_in );
	};
	if( sum == 0 ) return 0;
	return sigma / sum;
};


float XrayMaterial::coherent( const float energy_in, const float theta_in ) const {
//      Calculate X-ray coherent scatter cross-section at a given energy and angle (cm2/gm)
	if ( energy_in <= 0 ) return 0;
	float sigma = 0;
	float sum = 0;
	int i;
	//  Normalize composition so that absorption does not vary inappropriately
	for ( i=0; i<fractions.size(); i++ ) sum += fractions[i];
	for ( i=0; i<elements.size(); i++ ) {
		sigma += fractions[i] * scatter_tables[i].coherent( energy_in, theta_in );
	};
	if( sum == 0 ) return 0;
	return sigma / sum;
};


//      Data retrieval functions - get info from element_list (may be longer than input list because of oxygen, matrix, etc.)

float XrayMaterial::fraction( const Element el ) const {
    //  Returns the actual fraction of an element (not the oxide fraction)
    float fraction_out = 0;
    int ie = find_element( el, elements );
    if( ie >= 0 ) fraction_out = fractions[ie];
    return fraction_out;
}

float XrayMaterial::fraction( const int z_in ) const {
    if( ! Element::check_Z( z_in ) ) return 0;
    Element temp_el( z_in );
    return fraction( temp_el );
};


float XrayMaterial::oxide_ratio( const Element el ) const {
    //  Returns the oxide ratio associated with an element
    float oxide_ratio_out = 0;
    int ie = find_element( el, elements );
    if( ie >= 0 && ie < formula_info.size() ) oxide_ratio_out = formula_info[ie].formula_ratio;
    return oxide_ratio_out;
}

const LightElements &XrayMaterial::stoichiometry( const Element el ) const {
    //  Returns light element info included with this analyte element via stoichiometry
    int ie = find_element( el, elements );
    if( ie >= 0 && ie < formula_info.size() ) return formula_info[ie];
    return dummy_light_elements;
};

float XrayMaterial::uncertainty( const Element el ) const {
    //  Returns the uncertainty associated with an element
    //  This is just carried along, it is not used in this class
    float uncertainty_out = 0;
    int ie = find_element( el, element_list_input );
    if( ie >= 0 && ie < uncertainties.size() ) uncertainty_out = uncertainties[ie];
    return uncertainty_out;
}

float XrayMaterial::fraction_formula( const Element el ) const {
    float f = 0;
    int ie = find_element( el, elements );
    //  Get formula info from input element list
    LightElements formula_temp; //  With default values, pure element
    int ie_f = find_element( el, element_list_input );
    if( ie_f >= 0 && ie_f < formula_info.size() ) formula_temp = formula_info[ie_f];
    if( ie >= 0 )
        f = calculate_fraction_formula( elements[ie], fractions[ie], formula_temp );
    return f;
}

float XrayMaterial::fraction_input( const Element el ) const {
    //  Returns the input fraction of an element (without modification)
    float fraction_out = 0;
    int ie = find_element( el, element_list_input );
    if( ie >= 0 ) fraction_out = fractions_input[ie];
    return fraction_out;
};

float XrayMaterial::fraction_oxygen( const Element el ) const {
    //  Returns the oxygen fraction associated with an element
    float fraction_out = 0;
    int ie = find_element( el, element_list_input );
    if( ie >= 0 ) {
        fraction_out = calculate_fraction_oxygen( el, fractions_input[ie], formula_info[ie] );
    }
    return fraction_out;
};

float XrayMaterial::fraction_carbon( const Element el ) const {
    //  Returns the oxygen fraction associated with an element
    float fraction_out = 0;
    int ie = find_element( el, element_list_input );
    if( ie >= 0 ) {
        fraction_out = calculate_fraction_carbon( el, fractions_input[ie], formula_info[ie] );
    }
    return fraction_out;
};

float XrayMaterial::fraction_light( const Element el ) const {
    //  Returns the totaal light element fraction associated with an element
    float fraction_out = 0;
    int ie = find_element( el, element_list_input );
    if( ie >= 0 ) {
        fraction_out = calculate_fraction_light( el, fractions_input[ie], formula_info[ie] );
    }
    return fraction_out;
};




//      Data change functions - put data into element_list_input

void XrayMaterial::add_element( const Element element_in, const float fraction_in, const bool oxide_in ) {
    LightElements temp;
    if( oxide_in ) {
        temp.formula = OXIDE;
        temp.formula_ratio = default_oxide_ratio( element_in );
    }
    int ie = find_element( element_in, element_list_input );
    if( ie < 0 ) {
        element_list_input.push_back( element_in );
        fractions_input.push_back( fraction_in );
        formula_info.push_back( temp );
        uncertainties.push_back( 0 );
    } else {
        fractions_input[ie] = fraction_in;
        formula_info[ie] = temp;
        uncertainties[ie] = 0;
    }
    populate_element_list();
    calculate_element_fractions();
    return;
};

void XrayMaterial::add_element( const Element element_in, const float fraction_in, const LightElements formula_in ) {
    int ie = find_element( element_in, element_list_input );
    if( ie < 0 ) {
        element_list_input.push_back( element_in );
        fractions_input.push_back( fraction_in );
        formula_info.push_back( formula_in );
        uncertainties.push_back( 0 );
    } else {
        fractions_input[ie] = fraction_in;
        formula_info[ie] = formula_in;
        uncertainties[ie] = 0;
    }
    populate_element_list();
    calculate_element_fractions();
    return;
};

void XrayMaterial::fraction( const Element el, const float val ) {
    //  Sets the fraction of an element
    //  Note that whether it is setting the element or oxide fraction
    //     depends on the value of the oxide ratio for this element
    //  Unless the INPUT_FRACTIONS_OXIDES flag is false, the it is always element fraction
    //  Ignore invalid values, fraction must be between zero and one (inclusive)
    if( val < 0 ) return;
    int ie = find_element( el, element_list_input );
    if( ie >= 0 && ie < fractions_input.size() ) {
        fractions_input[ie] = val;
        calculate_element_fractions();
    }
    return;
}

void XrayMaterial::normalize( const float normalize_in ) {
    //  Normalize input fractions to input value (ignore if input is zero)
    if( normalize_in <= 0 ) return;
    //  Update fractions to insure formulas are included properly
    calculate_element_fractions();
    //  Sum over all elements to get normalization correct including additions from stoichiometry
    float sum = 0;
    unsigned int ie;
    for( ie=0; ie<fractions.size(); ie++ ) sum += fractions[ie];
    for( ie=0; ie<element_list_input.size(); ie++ ) fractions_input[ie] *= normalize_in / sum;
    calculate_element_fractions();
    return;
};

void XrayMaterial::uncertainty( const Element el, const float val ) {
    //  Sets the uncertainty of an element
    //  This is just carried along, it is not used in this class
    if( val < 0 ) return;
    int ie = find_element( el, element_list_input );
    if( ie >= 0 && ie < uncertainties.size() ) {
        uncertainties[ie] = val;
    }
    return;
}

void XrayMaterial::oxide_ratio( const Element el, const float val ) {
    //  Sets the oxide ratio of an element
    //  Negative value implies use default oxide ratio
    float temp_ratio = val;
    if( temp_ratio < 0 ) temp_ratio = default_oxide_ratio( el );
    int ie = find_element( el, element_list_input );
    if( ie >= 0 && ie < formula_info.size() ) {
        formula_info[ie].formula = OXIDE;
        formula_info[ie].formula_ratio = temp_ratio;
        populate_element_list();
        calculate_element_fractions();
    }
    return;
}

void XrayMaterial::stoichiometry( const Element el, const LightElements formula_in ) {
    //  Sets the oxide ratio of an element
    //  Negative value implies use default oxide ratio
    float temp_ratio = formula_in.formula_ratio;
    if( temp_ratio < 0 ) temp_ratio = default_formula_ratio( el, formula_in );
    int ie = find_element( el, element_list_input );
    if( ie >= 0 && ie < formula_info.size() ) {
        formula_info[ie] = formula_in;
        formula_info[ie].formula_ratio = temp_ratio;
        populate_element_list();
        calculate_element_fractions();
    }
    return;
}

void XrayMaterial::convert_to_oxides() {
    //  Sets all elements that don't already have a formula to OXIDE and oxide ratios to the default values
    //  Converts all input fractions to formula fractions if not already
    int ie;
    for( ie=0; ie<element_list_input.size(); ie++ ) {
        if( formula_info[ie].formula != PURE_ELEMENT ) continue;
        formula_info[ie].formula = OXIDE;
        formula_info[ie].formula_ratio = default_oxide_ratio( element_list_input[ie] );
        formula_info[ie].input_fractions_are_formula = true;
        fractions_input[ie] = calculate_fraction_formula( element_list_input[ie], fractions_input[ie], formula_info[ie] );
    }
    populate_element_list();
    calculate_element_fractions();
    return;
};




//      Private functions

int XrayMaterial::find_element( const Element el_in, const vector <Element> &e_list ) const {
    int ie;
    for( ie=0; ie<e_list.size(); ie++ ) {
        if( e_list[ie] == el_in ) {
            return ie;
        }
    }
    return -1;
}


void XrayMaterial::calculate_element_fractions() {
//  Calculates the actual element fractions given the input element list and the oxide ratios
//  Don't change element list if it is already set up (so tables do not have to be re-loaded)
    int ie;
    for( ie=0; ie<fractions.size(); ie++ ) fractions[ie] = 0;
	Element oxygen( 8 );
	float oxygen_fraction = 0;
	Element carbon( 6 );
	float carbon_fraction = 0;
    for( ie=0; ie<element_list_input.size(); ie++ ) {
        fractions[ie] = fractions_input[ie];
        if( element_list_input[ie] == oxygen || element_list_input[ie] == carbon ) continue;
        //  If this element has a light element formula associated with it, re-calculate element fractions
        if( formula_info[ie].formula != PURE_ELEMENT ) {
            if( formula_info[ie].input_fractions_are_formula ) fractions[ie] = calculate_fraction_element( element_list_input[ie], fractions_input[ie], formula_info[ie] );
            //  Note that we have to use the element fraction here, not the formula fraction
            oxygen_fraction += calculate_fraction_oxygen( element_list_input[ie], fractions[ie], formula_info[ie] );
            carbon_fraction += calculate_fraction_carbon( element_list_input[ie], fractions[ie], formula_info[ie] );
        }
    }
    //  Increase oxygen and carbon amounts from direct input (may or may not be zero) by sums from formulas
    for( ie=0; ie<elements.size(); ie++ ) {
        if( elements[ie] == oxygen ) fractions[ie] += oxygen_fraction;
        if( elements[ie] == carbon ) fractions[ie] += carbon_fraction;
    }
    oxygen_added = oxygen_fraction;
    carbon_added = carbon_fraction;

    //update density and mass thickness
    //  Calculate theoretical density if density has not been given
    if( ! fixed_density ) mass_density = calculate_theoretical_density();
    m_thickness = thickness_in * mass_density;
}

void XrayMaterial::populate_element_list() {
//  Populates a new element list with element fractions (not oxide fractions) and adds oxygen as appropriate
    const int n = element_list_input.size();
    elements.clear();
    elements.resize( n );
    fractions.resize( n, 0 );
    absorption_tables.clear();
    scatter_tables.clear();
	Element oxygen( 8 );
	Element carbon( 6 );
	int oxygen_index = -1;
	int carbon_index = -1;
	bool add_oxygen = false;
	bool add_carbon = false;
    int ie;
    for( ie=0; ie<element_list_input.size(); ie++ ) {
        elements[ie] = element_list_input[ie];
         //  If this is oxygen, note its index
        if( element_list_input[ie] == oxygen ) oxygen_index = ie;
         //  If this is carbon, also note its index
        if( element_list_input[ie] == carbon ) carbon_index = ie;
        //  Check to see if we need to add oxygen or carbon to the list
        switch( formula_info[ie].formula ) {
            case PURE_ELEMENT:  break;
            case OXIDE:  add_oxygen = true; break;
            case CARBONATE:  add_oxygen = true; add_carbon = true; break;
        }
    }
    //  Add oxygen and carbon to the element list if needed and not already present
        //  Fractions will be set in calculate_element_fractions
    if( add_oxygen && oxygen_index < 0 ) {
        elements.push_back( oxygen );
        fractions.push_back( 0 );
    }
    if( add_carbon && carbon_index < 0 ) {
        elements.push_back( carbon );
        fractions.push_back( 0 );
    }
//		load absorption tables
    absorption_tables.clear();
    scatter_tables.clear();
	for ( ie=0; ie<elements.size(); ie++ ) {
		XrayXsectTable abs( elements[ie] );
		absorption_tables.push_back( abs );
        ScatterXsectTable scat( elements[ie] );
        scatter_tables.push_back( scat );
	};
}


float XrayMaterial::calculate_fraction_element( const Element el, const float formula_fraction, const LightElements formula_info_in ) {
    float formula_weight = calculate_formula_weight( el, formula_info_in );
    return formula_fraction * el.atomicWeight() / formula_weight;
}

float XrayMaterial::calculate_fraction_formula( const Element el, const float element_fraction, const LightElements formula_info_in ) {
    //  Calculates the fraction of the total formula (including the analyte element) given the pure element fraction
    float formula_weight = calculate_formula_weight( el, formula_info_in );
    return element_fraction * formula_weight / el.atomicWeight();
};

float XrayMaterial::calculate_fraction_oxygen( const Element el, const float element_fraction, const LightElements formula_info_in ) {
    Element oxygen( 8 );
    float oxf = calculate_atomic_ratio( el, oxygen, formula_info_in );
    return element_fraction * oxf * oxygen.atomicWeight() / el.atomicWeight();
}

float XrayMaterial::calculate_fraction_carbon( const Element el, const float element_fraction, const LightElements formula_info_in ) {
    Element carbon( 6 );
    float oxf = calculate_atomic_ratio( el, carbon, formula_info_in );
    return element_fraction * oxf * carbon.atomicWeight() / el.atomicWeight();
}

float XrayMaterial::calculate_fraction_light( const Element el, const float element_fraction, const LightElements formula_info_in ) {
    //  Calculates the total fraction of the light elements from the formula (not including the analyte element) given the pure element fraction
    float formula_weight = calculate_formula_weight( el, formula_info_in );
    return element_fraction * ( formula_weight - el.atomicWeight() ) / el.atomicWeight();
};

float XrayMaterial::calculate_formula_weight( const Element el, const LightElements formula_info_in ) {
    //  Calculates formula weight for use in calculating fractions of elements that appear in the formulas
    Element oxygen( 8 );
    Element carbon( 6 );
    float formula_weight = el.atomicWeight();
    formula_weight += calculate_atomic_ratio( el, oxygen, formula_info_in ) * oxygen.atomicWeight();
    formula_weight += calculate_atomic_ratio( el, carbon, formula_info_in ) * carbon.atomicWeight();
    return formula_weight;
};

float XrayMaterial::calculate_atomic_ratio( const Element el, const Element formula_el, const LightElements formula_info_in ) {
    //  Calculates the atomic ratio of a selected element to the analyte element given the formula info
    Element oxygen( 8 );
    Element carbon( 6 );
    Element iron( 26 );
    float formula_ratio = formula_info_in.formula_ratio;
    float atomic_ratio = 0;
    switch( formula_info_in.formula ) {
        case PURE_ELEMENT:  break;
        case OXIDE:  if( formula_el == oxygen ) atomic_ratio = formula_ratio; break;
        case CARBONATE: //  CO3 with -2 oxidation state => formulas per element atom same as oxygen
            if( el == iron ) formula_ratio = 1;   //  Ignore oxide ratio for iron and make it Fe2+ if carbonate
            if( formula_el == oxygen ) atomic_ratio = 3 * formula_ratio;
            if( formula_el == carbon ) atomic_ratio = formula_ratio;
            break;
        default:   break;
    }
    return atomic_ratio;
};



float XrayMaterial::default_oxide_ratio( const Element el ) {
	const int max_Z_ox = 101;
	const int oxidationState[max_Z_ox] = { 0,	//	Sargent Welch periodic table 1962
//        1,  0,  1,  2,  3,  4,  3, -2, -1,  0,	//	1-10
        0,  0,  0,  0,  0,  0,  0, 0, 0,  0,	//	1-10    No associated oxide since these will always be matrix elements
        1,  2,  3,  4,  5,  6, -1,  0,  1,  2,	//	11-20
//        3,  4,  5,  3,  2,  3,  2,  2,  2,  2,	//	21-30	//	Fe2O3
        3,  4,  5,  3,  2,  2,  2,  2,  2,  2,	//	21-30	//	FeO
        3,  4,  3,  4, -1,  0,  1,  2,  3,  4,	//	31-40
        5,  6,  7,  4,  3,  2,  1,  2,  3,  4,	//	41-50
        3,  4, -1,  0,  1,  2,  3,  3,  4,  3,	//	51-60
        3,  3,  3,  3,  3,  3,  3,  3,  3,  3,	//	61-70
        3,  4,  5,  6,  7,  4,  4,  4,  3,  2,	//	71-80
        1,  2,  3,  2,  0,  0,  1,  2,  3,  4,	//	81-90
        5,  6,  5,  4,  3,  3,  3,  3,  0,  0,	//	91-100
    };
    int temp_Z = el.Z();
    //  Allow default value for iron to be modified by command line option
    if( temp_Z == 26 && default_modified_iron_oxide_ratio >= 0 ) return default_modified_iron_oxide_ratio;
    if( temp_Z > 0 && temp_Z < max_Z_ox && oxidationState[ temp_Z ] > 0 ) {
        return float( oxidationState[ temp_Z ] ) / 2; //  assumes oxygen is -2;
    } else {
        return 0;
    }
}
float XrayMaterial::default_carbonate_ratio( const Element el ) {
	const int max_Z_c = 101;
	const float carbonate_atomic_ratio[max_Z_c] = { 0,	//	E-mail from Joel Hurowitz, 12/14/2020, 8:55 AM
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,	//	1-10    Zero means element does not form a carbonate
        0,  1,  0,  0,  0,  0,  0,  0,  0,  1,	//	11-20   Mg and Ca
        0,  0,  0,  0,  1,  1,  0,  0,  0,  0,	//	21-30   Mn and Fe
        0,  0,  0,  0,  0,  0,  0,  1,  0,  0,	//	31-40   Sr
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,	//	41-50
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,	//	51-60
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,	//	61-70
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,	//	71-80
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,	//	81-90
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0 	//	91-100
    };
    int temp_Z = el.Z();
    if( temp_Z > 0 && temp_Z < max_Z_c ) {
        return carbonate_atomic_ratio[ temp_Z ];
    } else {
        return 0;
    }
}

float XrayMaterial::calculate_theoretical_density() const {
//  Calculate density using theoretical solid density of each element
//  Since these are weight fractions, we know the mass and are actually averaging the volume
//      of each element, so the sum is over inverse density
	float d = 0;
	int ie;
	for ( ie=0; ie<elements.size(); ie++ ) {
        d += fractions[ie] / elements[ie].density();
    }
    if( d <= 0 ) return 0;
	else return 1 / d;
}

float XrayMaterial::avgZ() const {
	float sum = 0;
	float z = 0;
	int i;
	//  Normalize composition to be sure average is accurate
	for ( i=0; i<fractions.size(); i++ ) sum += fractions[i];
	for ( i=0; i<elements.size(); i++ ) {
		z += fractions[i] * elements[i].Z();
	};
	if( sum == 0 ) return 0;
	return z / sum;
};

float XrayMaterial::avgA() const {
	float sum = 0;
	float a = 0;
	int i;
	//  Normalize composition to be sure average is accurate
	for ( i=0; i<fractions.size(); i++ ) sum += fractions[i];
	for ( i=0; i<elements.size(); i++ ) {
		a += fractions[i] * elements[i].atomicWeight();
	};
	if( sum == 0 ) return 0;
	return a / sum;
};

float XrayMaterial::avgZoverA() const {
	float sum = 0;
	float z_a = 0;
	int i;
	//  Normalize composition to be sure average is accurate
	for ( i=0; i<fractions.size(); i++ ) sum += fractions[i];
	for ( i=0; i<elements.size(); i++ ) {
		z_a += fractions[i] * elements[i].Z() / elements[i].atomicWeight();
	};
	if( sum == 0 ) return 0;
	return z_a / sum;
};

float XrayMaterial::avgAoverZ() const {
	float sum = 0;
	float z_a = 0;
	int i;
	//  Normalize composition to be sure average is accurate
	for ( i=0; i<fractions.size(); i++ ) sum += fractions[i];
	for ( i=0; i<elements.size(); i++ ) {
		z_a += fractions[i] * elements[i].atomicWeight() / elements[i].Z();
	};
	if( sum == 0 ) return 0;
	return z_a / sum;
};


string XrayMaterial::formula_string( const Element formula_element_in ) const {
    //  Returns a string that represents the formula unit included for this analyte
    //      element via stoichiometry (version for instantiated objects)
    int ie = find_element( formula_element_in, element_list_input );
    if( ie >= 0 && ie < formula_info.size() ) return formula_string( formula_element_in, formula_info[ie] );
    else return "";
};

float XrayMaterial::default_formula_ratio( const Element el, const LightElements formula_info_in ) {
    //  Finds the default formula ratio for oxides or carbonates
    float temp_ratio = 0;
    if( formula_info_in.formula == OXIDE ) temp_ratio = default_oxide_ratio( el );
    if( formula_info_in.formula == CARBONATE ) temp_ratio = default_carbonate_ratio( el );
    return temp_ratio;
};

string XrayMaterial::formula_string( const Element formula_element_in, const LightElements formula_info_in, const bool suffix_only ) {
    //  Returns a string that represents the formula unit included for this analyte
    //      element via stoichiometry (static version, all info in arguments)
    string output_str;
    if( !suffix_only ) output_str += formula_element_in.symbol();
    int oxide_ratio_int = 2 * formula_info_in.formula_ratio + 0.01f;
    if( formula_info_in.formula == PURE_ELEMENT ) return output_str;
    else if( formula_info_in.formula == OXIDE ) {
        switch( oxide_ratio_int ) {
            case 0:             break;
            case 1: output_str +="2O";      break;
            case 2: output_str += "O";      break;
            case 3: output_str += "2O3";    break;
            case 4: output_str += "O2";     break;
            case 5: output_str += "2O5";    break;
            case 6: output_str += "O3";     break;
            case 8: output_str += "O4";     break;
            default:    output_str += "O_Err";
                break;
            }
    } else if( formula_info_in.formula == CARBONATE ) {
        switch( oxide_ratio_int ) {
            case 0:             break;
            case 1: output_str +="2CO3";      break;
            case 2: output_str += "CO3";      break;
            case 3: output_str += "2(CO3)3";    break;
            case 4: output_str += "(CO3)2";     break;
            case 5: output_str += "2(CO3)5";    break;
            case 6: output_str += "(CO3)3";     break;
            default:    output_str += "C_Err";
                break;
            }
    } else output_str += "Undef";
    //  Indicate that Fe is total amount
    if( formula_element_in.Z() == 26 ) output_str += "-T";
	return output_str;
};

string XrayMaterial::toString() const
{
    ostringstream os;
    os << "XrayMaterial:" << endl;
    os << "Inputs=[" << endl;
    unsigned int ii;
    for( ii=0; ii<formula_info.size(); ii++ ) {
        os << formula_string( element_list_input[ii] ) << "  ";
        os << (formula_info[ii].input_fractions_are_formula?"formula %":"element %");
        os << ":  " << 100*fraction_input( element_list_input[ii] ) << " " << endl;
    }
    os << "]" << endl;
    os << "  uncertainties=[" << floatVecToString(uncertainties) << endl;
    os << "  fixed_density=" << fixed_density << endl;
    os << "  mass_density=" << mass_density << endl;
    os << "  thickness_in=" << thickness_in << endl;
    os << "  oxygen_added=" << oxygen_added << endl;
    os << "  carbon_added=" << carbon_added << endl;
    os << "  elements=[" << elementVecToString(elements) << "]" << endl;
    os << "  fractions=[" << floatVecToString(fractions) << endl;
    os << "  m_thickness=" << m_thickness << endl;
    os << "  absorption_tables=" << "TODO" << endl;
    os << "  scatter_tables=" << "TODO" << endl;
    return os.str();
}
