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

//
//  parse_arguments.cpp
//  PIQUANT
//
//  Created by W. T. Elam on 1/14/2017.
//  Copyright (c) 2017 APL/UW. All rights reserved.
//

#include <iostream>
#include <sstream>
#include "parse_element_list.h"
#include "parse_records.h"
#include "upper_trim.h"
#include "XRFutilities.h"
#include "XRFconstants.h"

//  Written Jan. 16, 2017
//      Parse list of element symbols and qualifiers for PIQUANT Subprocess
//      Qualifiers are separated from the element symbol by an underscore
//      Write helpful information to cout if errors and return a boolean (true if any errors)
//  Modified Sept. 30, 2017
//      Fix minor bug in error message, it was writing the wrong string
//  Modified Oct. 27, 2017
//      Separate parser for individual element string with qualifiers
//      Add information for holding standards description to element list entries
//  Modified Dec. 6, 2017
//      Add function to append or replace entry in element list, keeping selected info
//  Modified Dec. 6, 2017
//      Add MATRIX qualifier (to match standards input file parser)
//  Modified Jan. 3, 2017
//      Ignore trailing spaces and tabs for element symbols
//      Move tab, single quote, double quote, blank, comma, and underscore definitions to XRFconstants.h
//  Modified Jan. 26, 2018
//      Add net peak intensity to element list and calibration file output
//  Modified May 21, 2019
//      Implement Evaluate action, add given and relative error vs given to element list
//  Modified July 3, 2019
//      Allow matrix element fraction or percent to be entered in element list (as C_M=23.7%)
//  Modified Nov. 4, 2019
//      Add total error to element list, to include ECF standard deviation (and later certificate uncertainty)
//  Modified Dec. 16, 2019
//      Add element qualifier OUTPUT ("O") to force the element to be included in the evaluate list (with zeros if not in any standard in this run)
//  Modified Aug. 24, 2019
//      Allow element followed by an equal sign to define a material without X-ray qualifiers
//  Modified Nov. 24, 2020
//      Added Nov. 24, 2020  to hold matrix effect factor from FP calculation
//  Modified Dec. 15, 2020
//      Added check for "CO3" in element list to indicate some elements to be included as carbonates instead of oxides

using namespace std;


const bool parse_element_list( const std::string &element_list_in, std::vector <ElementListEntry> &elements_out, bool &carbonates, bool oxides ) {

    bool error = false; //  Set this to true and keep going if errors, to check all of input string
    elements_out.clear();
    vector <string> records;
    string delimiters( COMMA_CHARACTER );
    delimiters += BLANK_CHARACTER;
    int result = parse_records( delimiters, element_list_in, records );
    if( result < 0 ) {
        cout << "Error separating element list into entries at position " << -result << endl;
        error = true;
    }
    int ie;
    for( ie=0; ie<records.size(); ie++ ) {
//        cout << "record    \"" << records[ie] << "\"" << endl;
        //  Check for exact match to "CO3", if so set flag and skip parsing this entry
        if( records[ie] == "CO3" ) {
            carbonates = true;
            continue;
        }
        ElementListEntry temp_entry;
        error = parse_element_string( records[ie], temp_entry );
        if( carbonates ) {
            temp_entry.stoichiometry.formula = CARBONATE;
            temp_entry.stoichiometry.formula_ratio = XrayMaterial::default_formula_ratio( temp_entry.element,temp_entry.stoichiometry );
        }
        //  See if we should make this entry an oxide
        if( temp_entry.stoichiometry.formula_ratio == 0 && oxides ) {
            temp_entry.stoichiometry.formula = OXIDE;
            temp_entry.stoichiometry.formula_ratio = XrayMaterial::default_formula_ratio( temp_entry.element,temp_entry.stoichiometry );
        }
        elements_out.push_back( temp_entry );
    }

    return error;

}


const bool parse_element_string( const std::string &element_string_in, ElementListEntry &element_entry_out ) {

    bool error = false;
    const string underscore( UNDERSCORE_CHARACTER );
    const string equal_str( EQUAL_CHARACTER );

    string temp_symbol;
    string level_qual_str;
    int u = element_string_in.find( underscore );
    if( u < 0 || u >= element_string_in.length() ) {
        //  No underscore
        //  Check for element followed by equal sign then percent (or fraction etc.)
        int eq = element_string_in.find( equal_str );
        if( eq >= 0 && eq < element_string_in.length() ) {
            temp_symbol = element_string_in.substr( 0, eq );
            if( eq+1 < element_string_in.length() )
                    level_qual_str = upper_trim( element_string_in.substr( eq ) );
        } else {
            temp_symbol = element_string_in.substr(0,2);
            if( temp_symbol.length() == 2 && temp_symbol.substr(1,1) == BLANK_CHARACTER ) temp_symbol.erase( 1, 1 );
            if( temp_symbol.length() == 2 && temp_symbol.substr(1,1) == TAB_CHARACTER ) temp_symbol.erase( 1, 1 );
        }
    } else {
        temp_symbol = element_string_in.substr( 0, u );
        //  Some level or qualifier symbol follows the element, so save it for parsing later
        if( u+1 < element_string_in.length() ) level_qual_str = upper_trim( element_string_in.substr( u+1 ) );
    }
    if( Element::check_symbol( temp_symbol ) ) {
        Element temp( temp_symbol );
        element_entry_out.element = temp;
    } else {
        //  If this is not a valid element symbol, check for an atomic number
        istringstream temp_instr( temp_symbol );
        int z_test = 0;
        temp_instr >> z_test;
        if( ! temp_instr || ! Element::check_Z( z_test ) ) {
            cout << "Invalid element symbol " << temp_symbol << endl;
            error = true;
        } else {
            Element temp( z_test );
            element_entry_out.element = temp;
        }
    }
    //  Take apart and parse the level symbol and qualifier symbol (could be either or both)
//    cout << "Element list  s " << temp_symbol << "  " << element_entry_out.element.symbol() << "  ";
    if( level_qual_str.length() > 0 ) {
        //  Check for matrix element with fraction or percent value after equal sign
        string percent_str;
        if( level_qual_str.length() >= 2 && level_qual_str.substr(0,2) == "M=" ) {
            element_entry_out.qualifier = MATRIX;
            percent_str = level_qual_str.substr( 2 );
        } else if( level_qual_str.length() >= 1 && level_qual_str.substr(0,1) == equal_str ) {
            element_entry_out.qualifier = MATRIX;
            percent_str = level_qual_str.substr( 1 );
        }
        if( percent_str.length() > 0 ) {
            //  Parse the fraction or percent and put into element list entry
			enum ElementEntryNormalization { PPM = -1, FRACTION, PERCENT };
			ElementEntryNormalization composition_format = PERCENT;
            float percent_in = 0;
            if( percent_str.length() > 0 ) {
                istringstream value( strip_suffix(percent_str) );
                value >> percent_in;
                if( 0 <= percent_str.rfind( "%" ) && percent_str.rfind( "%" ) < percent_str.length() ) {
                    composition_format = PERCENT;
                } else if( 0 <= percent_str.rfind( "f" ) && percent_str.rfind( "f" ) < percent_str.length() ) {
                    composition_format = FRACTION;
                } else if( 0 <= percent_str.rfind( "F" ) && percent_str.rfind( "F" ) < percent_str.length() ) {
                    composition_format = FRACTION;
                } else if( 0 <= percent_str.rfind( "p" ) && percent_str.rfind( "p" ) < percent_str.length() ) {
                    composition_format = PPM;
                } else if( 0 <= percent_str.rfind( "P" ) && percent_str.rfind( "P" ) < percent_str.length() ) {
                    composition_format = PPM;
                }
                if( composition_format == PPM ) percent_in *= PPM_PERCENT;
                else if( composition_format == FRACTION ) percent_in *= 100;
                if ( ! value || percent_in <= 0 || percent_in > 100 ) {
                    error = true;
                    cout << "Invalid matrix percent in element list,  Element " << element_entry_out.element.symbol() << ", " << percent_str << endl;
                } else {
                    element_entry_out.percent = percent_in;
                }
            }
        } else {
            string temp_level_str = level_qual_str.substr(0,1);
            string temp_qual_str;
            //  If there are two symbols, treat the second one as a qualifier
            if( level_qual_str.length() > 1 ) temp_qual_str = level_qual_str.substr(1,1);
            if( temp_level_str == "K" ) element_entry_out.quant_level = K_LEVEL;
            else if( temp_level_str == "L" ) element_entry_out.quant_level = L_LEVEL;
            else if( temp_level_str == "M" ) element_entry_out.quant_level = M_LEVEL;
            else if( temp_level_str == "N" ) element_entry_out.quant_level = N_LEVEL;
            //  Level symbol not found, it may be a qualifier symbol (or an error)
            else temp_qual_str = temp_level_str;
            //  Parse qualifier symbol (if any)
            if( temp_qual_str.length() > 0 ) {
                if( temp_qual_str == "I" ) element_entry_out.qualifier = IGNORE;
                else if( temp_qual_str == "F" ) element_entry_out.qualifier = FORCE;
                else if( temp_qual_str == "X" ) element_entry_out.qualifier = EXCLUDE;
                else if( temp_qual_str == "M" ) element_entry_out.qualifier = MATRIX;
                else if( temp_qual_str == "O" ) element_entry_out.qualifier = OUTPUT;
                else {
                    cout << "Invalid quantification lines or qualifier " << temp_qual_str;
                    cout << " for element " << temp_symbol << endl;
                    error = true;
                }
            }
        }
//        cout << "     lv " << temp_level_str << "   " << element_entry_out.quant_level;
//        cout << "     q " << temp_qual_str << "   " << element_entry_out.qualifier;
    }
//    cout << endl;

    return error;

};


void add_element_list_entry( const ElementListEntry &element_entry, std::vector <ElementListEntry> &element_list_out,
            const bool ignore_qualifier ) {

//      Append or replace entry in element list, keeping selected info

    //  See if there is already an entry that matches this one (element and line qualifier)
    int element_entry_index = -1;
    int iel;
    for( iel=0; iel<element_list_out.size(); iel++ ) {
        if( ! ( element_list_out[iel].element == element_entry.element ) ) continue;
        if( ! ( element_list_out[iel].quant_level == element_entry.quant_level ) ) continue;
        //  Ignore qualifier match if a new non-default qualifier is given
        if( ( ! ignore_qualifier ) && ! ( element_list_out[iel].qualifier == element_entry.qualifier ) ) continue;
        if( element_list_out[iel].type != ELEMENT && element_list_out[iel].type != element_entry.type ) continue;
        element_entry_index = iel;
        break;
    }
    if( element_entry_index < 0 ) {
        //  Add new entry to element list
        element_list_out.push_back( element_entry );
    } else {
        //  replace the existing element list entry with the new one
        ElementListEntry default_entry;
        ElementListEntry new_entry = element_entry;
        ElementListEntry old_entry = element_list_out[element_entry_index];
        //  Transfer selected info to the new entry if it is not specified (check for change from defaults in new entry before transferring)
        if( new_entry.percent == default_entry.percent && old_entry.percent != default_entry.percent ) new_entry.percent = old_entry.percent;
        if( new_entry.stoichiometry.formula == default_entry.stoichiometry.formula && old_entry.stoichiometry.formula != default_entry.stoichiometry.formula ) new_entry.stoichiometry.formula = old_entry.stoichiometry.formula;
        if( new_entry.stoichiometry.formula_ratio == default_entry.stoichiometry.formula_ratio && old_entry.stoichiometry.formula_ratio != default_entry.stoichiometry.formula_ratio ) new_entry.stoichiometry.formula_ratio = old_entry.stoichiometry.formula_ratio;
        if( new_entry.stoichiometry.input_fractions_are_formula == default_entry.stoichiometry.input_fractions_are_formula && old_entry.stoichiometry.input_fractions_are_formula != default_entry.stoichiometry.input_fractions_are_formula ) new_entry.stoichiometry.input_fractions_are_formula = old_entry.stoichiometry.input_fractions_are_formula;
        if( new_entry.uncertainty == default_entry.uncertainty && old_entry.uncertainty != default_entry.uncertainty ) new_entry.uncertainty = old_entry.uncertainty;
        if( new_entry.weight == default_entry.weight && old_entry.weight != default_entry.weight ) new_entry.weight = old_entry.weight;
        if( new_entry.ecf == default_entry.ecf && old_entry.ecf != default_entry.ecf ) new_entry.ecf = old_entry.ecf;
        if( new_entry.ecf_sigma == default_entry.ecf_sigma && old_entry.ecf_sigma != default_entry.ecf_sigma ) new_entry.ecf_sigma = old_entry.ecf_sigma;
        if( new_entry.intensity == default_entry.intensity && old_entry.intensity != default_entry.intensity ) new_entry.intensity = old_entry.intensity;
        if( new_entry.coefficient == default_entry.coefficient && old_entry.coefficient != default_entry.coefficient ) new_entry.coefficient = old_entry.coefficient;
        if( new_entry.rel_err_coeff == default_entry.rel_err_coeff && old_entry.rel_err_coeff != default_entry.rel_err_coeff ) new_entry.rel_err_coeff = old_entry.rel_err_coeff;
        if( new_entry.given == default_entry.given && old_entry.given != default_entry.given ) new_entry.given = old_entry.given;
        if( new_entry.rel_err_given == default_entry.rel_err_given && old_entry.rel_err_given != default_entry.rel_err_given ) new_entry.rel_err_given = old_entry.rel_err_given;
        if( new_entry.total_err == default_entry.total_err && old_entry.total_err != default_entry.total_err ) new_entry.total_err = old_entry.total_err;
        if( new_entry.matrix == default_entry.matrix && old_entry.matrix != default_entry.matrix ) new_entry.matrix = old_entry.matrix;
        element_list_out[element_entry_index] = new_entry;
    }

    return;

};
