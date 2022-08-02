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

#ifndef XRAYMATERIAL_H
#define XRAYMATERIAL_H

#include <vector>
#include "Element.h"
#include "XrayXsectTable.h"
#include "ScatterXsectTable.h"

//  This class is intended to allow calculation of X-ray properties
//      of a material with arbitrary composition

//  More versatile struct and enum for stoichiometrically included light elements added Dec. 14, 2020
//  The populate_element_list and calculate_element_fractions functions have to be modified when new entries are added
enum LightElementFormula { PURE_ELEMENT=0, OXIDE, CARBONATE };
struct LightElements {
    //  Formula for light element inclusion
    LightElementFormula formula = PURE_ELEMENT;
    //  Atomic ratio of light element part of formula to analyte element
    //  Ratio is element oxidation state divided by oxygen oxidation state (2)
    //  For Na2O, ratio is 0.5, for FeO ratio is 1, for Fe2O3 ratio is 1.5, and for Fe3O4 ratio is 1.3333
    //  For CaCO3 ratio is 1
    float formula_ratio = 0;
    //  Fractions input are for entire formula, not just the element fraction
    bool input_fractions_are_formula = false;
    //  Add another float here to include waters of hydration, someday
};


class XrayMaterial
{
    public:
//      Constructors
        XrayMaterial();
        XrayMaterial(const Element &element_in, const bool oxides_in = false, const bool oxides_frac_flag = false);
        XrayMaterial(const std::vector <Element> &element_list_in, const std::vector <float> &fractions_in, const bool oxides_in = false, const bool oxides_frac_flag = false );
        XrayMaterial(const std::vector <Element> &element_list_in, const std::vector <float> &fractions_in, const std::vector  <LightElements> formula_info_in );
        XrayMaterial(const std::vector <Element> &element_list_in, const std::vector <float> &fractions_in,
                const std::vector  <LightElements> formula_info_in, const std::vector  <float> uncertainties_in );
        XrayMaterial(const int n_el, const int element_Z_in[], const float fractions_in[], const bool oxides_in = false, const bool oxides_frac_flag = false );

//      X-ray properties functions
        float transmission( const float energy_in, const float csc = 1 ) const;   //  csc is for slant path
        float absorption( const float energy_in, const float csc = 1 ) const;
        float photo( const float energy_in ) const;
        float photo( const Element el, const float energy_in ) const;
        float cross_section( const float energy_in )const;
        float cross_section( const Element el, const float energy_in ) const;
        const XrayXsectTable &cross_section_table( const Element el ) const;
        float incoherent( const float energy_in, const float theta_in ) const;
        float incoherent( const float energy_in, const float theta_in, const float scattered_energy_in ) const;
        float coherent( const float energy_in, const float theta_in ) const;

 //     Data access functions
        const int number_of_original_elements() const { return element_list_input.size(); };
        const int number_of_elements() const { return elements.size(); };
        const vector <Element> &original_element_list() const { return element_list_input; };
        const vector <Element> &element_list() const { return elements; };
        void add_element( const Element element_in, const float fraction_in, const bool oxide_in = false );
        void add_element( const Element element_in, const float fraction_in, const LightElements formula_in );
        float fraction( const Element el )const;
        float fraction( const int z_in ) const;
        void fraction( const Element el, const float val );
        void normalize( const float normalize_in );
        const vector <float> &fraction_list() const { return fractions; };
        float oxide_ratio( const Element el )const;
        void oxide_ratio( const Element el, const float val );
        const LightElements &stoichiometry( const Element el ) const;
        void stoichiometry( const Element el, const LightElements formula_in );
        string formula_string( const Element formula_element_in ) const;
        float uncertainty( const Element el ) const;
        void uncertainty( const Element el, const float val );
        float fraction_input( const Element el ) const;
        float fraction_formula( const Element el ) const;
        float fraction_oxygen( const Element el ) const;
        float fraction_carbon( const Element el ) const;
        float fraction_light( const Element el ) const;
        float added_oxygen( ) const { return oxygen_added; };
        float added_carbon( ) const { return carbon_added; };
        void convert_to_oxides();   //  Converts all input fractions to oxides using the default oxide ratios
        //  Set flag to force re-calculation of mass thickness if density or thickness changed
        float density() const { return mass_density; }    //  density in cm2/gm
        void density( const float val ) { if( val > 0 ) { fixed_density = true; mass_density = val; calculate_element_fractions(); }; }
        float thickness() const { return thickness_in; }  //  thickness in cm
        float mass_thickness() const { return m_thickness; }
        void thickness( const float val ) { if( val > 0 ) { thickness_in = val; calculate_element_fractions(); };}
        float avgZ() const;
        float avgA() const;
        float avgZoverA() const;
        float avgAoverZ() const;

        static float default_oxide_ratio( const Element el );
        static float default_carbonate_ratio( const Element el );
        static float default_formula_ratio( const Element el, const LightElements formula_info_in );
        static float calculate_fraction_element( const Element el, const float formula_fraction, const LightElements formula_info_in );
        static float calculate_fraction_formula( const Element el, const float element_fraction, const LightElements formula_info_in );
        //  Note that the following functions take the element fraction, not the formula fraction
        static float calculate_fraction_oxygen( const Element el, const float element_fraction, const LightElements formula_info_in );
        static float calculate_fraction_carbon( const Element el, const float element_fraction, const LightElements formula_info_in );
        static float calculate_fraction_light( const Element el, const float element_fraction, const LightElements formula_info_in );
        static float calculate_formula_weight( const Element el, const LightElements formula_info_in );
        static float calculate_atomic_ratio( const Element el, const Element formula_el, const LightElements formula_info_in );
        static string formula_string( const Element formula_element_in, const LightElements formula_info_in, const bool suffix_only = false );
        static float default_iron_oxide_ratio() { return default_modified_iron_oxide_ratio; };
        static void default_iron_oxide_ratio( const float value_in ) { default_modified_iron_oxide_ratio = value_in; };

//      debug function (usually commented out)
//        void debug_print();

        string toString() const;

    private:

    //  data
        std::vector <Element> element_list_input;
        std::vector <float> fractions_input;
        //  Replace oxide ratio with more versatile light element info    Dec. 14, 2020
        //std::vector <float> oxide_ratios;
        std::vector <LightElements> formula_info;
        std::vector <float> uncertainties;
        bool fixed_density = false;
        float mass_density = 0;
        float thickness_in = 0;
        float oxygen_added = 0;
        float carbon_added = 0;
        std::vector <Element> elements;
        std::vector <float> fractions;
        float m_thickness = 0;
        static float default_modified_iron_oxide_ratio;
        //  absorption cross-section tables (in cm2/gm)
        vector <XrayXsectTable> absorption_tables;
        vector <ScatterXsectTable> scatter_tables;
        LightElements dummy_light_elements; //  Default values for return if element not found in stoichiometry function
        XrayXsectTable no_table;    //  Dummy table for default return in cross_section_table for a missing element

    //  private functions
        int find_element( const Element el_in, const std::vector <Element> &e_list ) const;
        void calculate_element_fractions();
        void populate_element_list();
        float calculate_theoretical_density() const;
};

#endif // XRAYMATERIAL_H
