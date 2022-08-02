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
#include <sstream>
#include "upper_trim.h"
#include "fpMain.h"
#include "quantComponents.h"
#include "XRFcontrols.h"
#include "toStringHelpers.h"

// Feb.5, 2017
//	Begun writing, new function
//	Sets up the list of components that sum to make the full spectrum
//  Sets and checks the components that are used to quantify each element
//  Maps XrayLines objects to each component,
//      matching associated elements and physical processes that generate the signal
//  Also generates and parses the text descriptions of the components
//      for output to and input from calibration files
//  Do all this here so it's in one place if it is changed to improve quantification
//  At present the lines are grouped by edge level (K, L, M, or N) and all lines
//      with that principal quantum number are fit together (with one fit coefficient)
// Modified May 26, 2017
//  Change componentDescription to return a static string
// Modified Dec. 9, 2017
//  Add componentQuantLevel function to convert EdgeLevel of component to ElementQuantLevel in element list
//  Clean up selection of quant components (was setting quant flag for scatter components)
//  Modified Jan. 3, 2017
//      Move tab, single quote, double quote, blank, comma, and underscore definitions to XRFconstants.h
//  Modified July 2, 2019
//      Allow for background to be adjusted manually instead of in least-squares fit, using -b option parameter
//  Modified July 22, 2019
//      Treat any elements that have no components in the spectrum as matrix elements, and print warning
//  Modified Dec. 11, 2019
//      Make quantComponents and quantDefaults skip IGNORE elements (they have already been taken care of)
//  Modified Oct. 21, 2020
//      Add index for multiple background components (split for independent fitting)
//      Add separate component for Compton escape  (fixed bugs Oct. 26, 2020)
//  Modified Nov. 9, 2020
//      Fixed bug in quantDefaults where it was missing match for the component with index zero
//  Modified Nov. 24, 2020
//      Add matrix effect factor
//  Modified Nov. 30, 2020
//      Add fit Boolean so that L (or M) components could be fixed relative to K (or L) components
//      Also add factor to relate coefficients of non-fit components to components used for quant
//  Modified Feb. 26, 2021  Add DETECTOR_SHELF and DETECTOR_CE spectrum component types
//  Modified Apr. 6, 2021   Add CONTINUUM component type and sort out how to handle background and Compton escape (remove Det shelf component)
//  Modified Apr. 28, 2021  Add La and Lb1 component types to checkComponent (for extra Rh L lines, especially Lb1 fit separately)
//  Modified May 10, 2021   Added scale_under to use scale-under-peaks algorithm for calculated background (every time it is calculated in quantCalculate)
//                          Modify matchComponent to handle multiple SNIP_BKG components
//  Modified July 10, 2021  Add pulse pileup component


using namespace std;

int setupComponents( const std::vector <XrayLines> sourceLines, const std::vector <XrayLines> pureLines,
                std::vector <SpectrumComponent> &components_out ) {

    bool error = false;
    int result = 0;
    //  Set up components for the calculated spectrum
    //  Put in components for fluorescence from list of specimen emission lines
    result = makeComponents( ELEMENT, pureLines, components_out );
    if( result < 0 ) {
        cout << "Element makeComponents failed, result is " << result << endl;
        error = true;
    }
    //  Put in components for Rayleigh scatter from list of source lines
    result = makeComponents( RAYLEIGH, sourceLines, components_out );
    if( result < 0 ) {
        cout << "Rayleigh makeComponents failed, result is " << result << endl;
        error = true;
    }
    //  Put in components for Compton scatter from list of source lines
    result = makeComponents( COMPTON, sourceLines, components_out );
    if( result < 0 ) {
        cout << "Compton makeComponents failed, result is " << result << endl;
        error = true;
    }
    //  Put in a component for pulse pileup if indicated in XRFcontrols.h
    if( PILEUP_LIST_LENGTH > 0 ) {
        result = makeComponents( PILEUP, sourceLines, components_out );
        if( result < 0 ) {
            cout << "Pulse pileup makeComponents failed, result is " << result << endl;
            error = true;
        }
    }

    if( error ) return -1;
    return 0;

};


int makeComponents( const SpectrumComponentType type_in, const std::vector <XrayLines> &lines_in,
                std::vector <SpectrumComponent> &components_out, const int n_bkg ) {
    SpectrumComponent temp_component;
    //  Create a special component for the continuum background if input emission lines are empty
    if( lines_in.size() == 0 && type_in == PRIMARY_CONTINUUM ) {
        temp_component.type = type_in;
        //  Use defaults for everything else
        components_out.push_back( temp_component );
        return 0;
    };
    //  Create component for the continuum background if input emission lines are empty
    if( lines_in.size() == 0 && type_in == CONTINUUM && n_bkg > 0 ) {
        temp_component.type = type_in;
        temp_component.bkg = true;
        //  Include desired number of background components (default 1)
        int ic;
        for( ic=0; ic<n_bkg; ic++ ) {
            //  Set index sequentially
            temp_component.bkg = true;
            temp_component.bkg_index = ic;
            temp_component.plot = false;
            //  Use defaults for everything else
            components_out.push_back( temp_component );
        }
        return 0;
    };
    //  Create component for the SNIP background if input emission lines are empty
    if( lines_in.size() == 0 && type_in == SNIP_BKG && n_bkg > 0 ) {
        temp_component.type = type_in;
        temp_component.bkg = true;
        //  Include desired number of background components (default 1)
        int ic;
        for( ic=0; ic<n_bkg; ic++ ) {
            //  Set index sequentially
            temp_component.bkg = true;
            temp_component.bkg_index = ic;
            temp_component.plot = false;
            //  Use defaults for everything else
            components_out.push_back( temp_component );
        }
        return 0;
    };
    //  Create a separate component for Compton escape
    if( lines_in.size() == 0 && type_in == DETECTOR_CE ) {
        temp_component.type = type_in;
        temp_component.bkg = true;
        temp_component.fit = false;
        //  Use defaults for everything else
        components_out.push_back( temp_component );
        return 0;
    };
    //  Create a separate component for pulse pileup
    if( lines_in.size() == 0 && type_in == PILEUP ) {
        temp_component.type = type_in;
        temp_component.bkg = false;
        temp_component.fit = false;
        //  Use defaults for everything else
        components_out.push_back( temp_component );
        return 0;
    };
    //  Make sure there is a component to include every X-ray emission line
    //  A component may include more that one XrayLines object
    //      (or, maybe in the future, part of an object, like alpha or beta lines)
    //  The lines will be from the source (for scatter) or the specimen (for fluorescence)
    int il;
    for( il=0; il<lines_in.size(); il++ ) {
        //  See if the component is already listed
        bool found = false;
        int ic;
        for( ic=0; ic<components_out.size(); ic++ ) {
            if( components_out[ic].type == type_in
                && checkComponent( components_out[ic], lines_in[il], il ) ) {
                found = true;
                break;
            }
        }
        if( found ) continue;
        //  Make a new component and add it to the output list
        EdgeLevel level = lines_in[il].edge().level();
        Element el = lines_in[il].edge().element();
        temp_component.type = type_in;
        temp_component.element = el;
        temp_component.level = level;
        components_out.push_back( temp_component );
    }
    return 0;
};


//  Checks match between element and line-to-component map (assumes input XrayLines object matches component type)
bool checkComponent( const SpectrumComponent &component_in, const XrayLines &lines_in, const int line_index_in ) {
    //  Line index is currently ignored, included here in case included lines someday depends in index (like alpha/beta)
    if( ( component_in.type == ELEMENT || component_in.type == COMPTON || component_in.type == RAYLEIGH
            || component_in.type == PRIMARY_LINES || component_in.type == La || component_in.type == Lb1 )
            && component_in.element == lines_in.edge().element()
            && component_in.level == lines_in.edge().level() ) return true;
    return false;
};


bool matchComponent( const SpectrumComponent &component_1, const SpectrumComponent &component_2 ) {
    bool identical = true;
    if( component_1.type != component_2.type ) identical = false;
    if( ! ( component_1.element == component_2.element ) ) identical = false;
    if( component_1.level != component_2.level ) identical = false;
    if( ( component_1.type == CONTINUUM || component_1.type == SNIP_BKG )
            && ( component_1.bkg_index != component_2.bkg_index ) ) identical = false;
    return identical;
};


//  Set flags of components that will be used to quantify their associated elements
//  Mark any components that are to be excluded based on the element list inputs
int quantComponents( const std::vector <ElementListEntry> element_list_in,
				std::vector <SpectrumComponent> &components_out ) {
    int ie;
    for( ie=0; ie<element_list_in.size(); ie++ ) {
        if( element_list_in[ie].qualifier == IGNORE ) continue;
        int ic;
        for( ic=0; ic<components_out.size(); ic++ ) {
            //  Skip components that are not element emission lines
            if( components_out[ic].type != ELEMENT ) continue;
            //  Find out if this component will be used to quantify its associated element
            Element el = components_out[ic].element;
            //  Check for element match
            if( ! ( el == element_list_in[ie].element ) ) continue;
            //  See if all components for this element are to be excluded
            if( element_list_in[ie].quant_level == NO_QUANT_LEVEL && element_list_in[ie].qualifier == EXCLUDE ) {
                components_out[ic].enabled = false;
                continue;
            }
            //  Check the element emission line for a match
            ElementQuantLevel quant_level_check = componentQuantLevel( components_out[ic] );
            if( element_list_in[ie].quant_level != quant_level_check ) continue;
            //  Be sure this component is not excluded or ignored
            if( element_list_in[ie].qualifier == EXCLUDE ) {
                components_out[ic].enabled = false;
            } else {
                //  OK, set it as the component used for quantification of its element
                if( element_list_in[ie].qualifier != IGNORE ) {
                    components_out[ic].quant = true;
                    break;
                }
            }
        }
    }
    return 0;
};


//  Choose default components to quantify any elements that do not already have an associated component
//  Zero the quantity of any elements that are not matrix and have no component, and print warning
int quantDefaults( std::vector <ElementListEntry> &element_list_in,
				std::vector <SpectrumComponent> &components_out ) {
    int ie;
    for( ie=0; ie<element_list_in.size(); ie++ ) {
        if( element_list_in[ie].qualifier == IGNORE ) continue;
        Element element_in = element_list_in[ie].element;
        bool quant_found = false;
        int k_index = -1;
        int l_index = -1;
        int m_index = -1;
        int n_index = -1;
        int ic;
        for( ic=0; ic<components_out.size(); ic++ ) {
            //  Check for element match
            if( ! ( element_in == components_out[ic].element ) ) continue;
            if( components_out[ic].quant ) {
                quant_found = true;
                break;
            } else {
                switch( components_out[ic].level ) {
                    case K: k_index = ic;   break;
                    case L: l_index = ic;   break;
                    case M: m_index = ic;   break;
                    case N: n_index = ic;   break;
                    default:                break;
                }
            }
        }
        if( ! quant_found ) {
            //  Quantify with the highest level that has an available component
            //  Components won't be created if they are not excited under the conditions
            if( k_index >= 0 ) components_out[k_index].quant = true;
            else if( l_index >= 0 ) components_out[l_index].quant = true;
            else if( m_index >= 0 ) components_out[m_index].quant = true;
            else if( n_index >= 0 ) components_out[n_index].quant = true;
            else {      //  No component to quantify this element
                //  If it is already a matrix element or excluded, all is OK
                if( element_list_in[ie].qualifier == MATRIX || element_list_in[ie].qualifier == EXCLUDE ) continue;
                //  Mark it as a matrix element and print a warning
                element_list_in[ie].qualifier = MATRIX;
                cout << "*** Warning - there are no emission lines in the spectrum for ";
                cout << element_list_in[ie].element.symbol();
                cout << " (it will be treated as a matrix element)." << endl;
            }
        }
    }
    return 0;
};


std::string componentDescription( const SpectrumComponent &component_in ) {
    string s;
    s += component_in.element.symbol();
    s += UNDERSCORE_CHARACTER;
    switch( component_in.level ) {
        case NO_EDGE:   break;
        case K:   s += "K";    break;  //  EdgeLevel Kl avoids conflict with EdgeIndex K
        case L:    s += "L";    break;
        case M:    s += "M";    break;
        case N:    s += "N";    break;
        default:                break;
    };
    switch( component_in.type ) {
        case NO_COMPONENT:      s = "none";    break;
        case ELEMENT:           break;  //  this is the default case if no additional description
        case COMPTON:           s += "_inc";    break;   //  Compton scatter = incoherent scatter
        case RAYLEIGH:          s += "_coh";    break;   //  Rayleigh scatter = incoherent scatter
        case CONTINUUM:        s = "calc bkg";    break;    //  Get rid of element symbol and edge level, do not apply to this component
        case SNIP_BKG:        s = "SNIP bkg";    break;    //  Get rid of element symbol and edge level, do not apply to this component
        case PRIMARY_LINES:     s += "_pri";    break;    //  Primary spectrum from anode element lines
        case PRIMARY_CONTINUUM: s = "continuum";    break;    //  Get rid of element symbol and edge level, do not apply to this component
        case La:                s += "_coh_La";    break;   //  For debugging extra intensity in tube scatter peaks from L lines
        case Lb1:               s += "_coh_Lb1";    break;   //  For debugging extra intensity in tube scatter peaks from L lines
        case DETECTOR_CE:       s = "DetCE";    break;   //  For debugging extra intensity from Compton escape in detector
        case OPTIC_TRANS:       s = "Optic";    break;   //  For plotting optic response during primary spectrum calculation
        case PILEUP:            s = "Pileup";    break;   //  For pulse pileup calculation
    };
    if( component_in.type == CONTINUUM ) {
        //  Add index for multiple background components
        ostringstream temp_ostr;
        temp_ostr << component_in.bkg_index;
        s += temp_ostr.str();
    }
    return s;
};


int parse_component( const std::string component_string_in,
				SpectrumComponent &component_out ) {
    const string underscore( UNDERSCORE_CHARACTER );
    string temp_symbol = upper_trim( component_string_in );
    const int l = component_string_in.length();
    //  Check for a background component (which has no element symbol)
    if( temp_symbol.size() == 3 && temp_symbol == "BKG" ) {
        component_out.type = CONTINUUM;
        component_out.level = NO_EDGE;
        return 0;
    }
    component_out.type = ELEMENT;   //  type defaults to ELEMENT
    component_out.level = K;   //  EdgeLevel defaults to K
    string temp_level_str;
    string temp_type_str;
    //  Find element symbol, instantiate the Element object, and put it in the component
    int u = component_string_in.find( underscore );
    if( u < 0 || u >= l ) {
        temp_symbol = component_string_in;
    } else {
        temp_symbol = component_string_in.substr( 0, u );
        if( u+1 < l ) temp_level_str = upper_trim( component_string_in.substr( u+1, 1 ) );
        if( u+2 < l ) temp_type_str = upper_trim( component_string_in.substr( u+2, l ) );
    }
    if( Element::check_symbol( temp_symbol ) ) {
        Element temp_el( temp_symbol );
        component_out.element = temp_el;
    } else {
        //  If this is not a valid symbol, check for an atomic number
        istringstream temp_instr( temp_symbol );
        int z_test = 0;
        temp_instr >> z_test;
        if( ! temp_instr || ! Element::check_Z( z_test ) ) {
            cout << "Invalid element symbol " << temp_symbol << endl;
            return -1;
        } else {
            Element temp_el( z_test );
            component_out.element = temp_el;
        }
    }
    if( temp_level_str.length() > 0 ) {
        if( temp_level_str == "K" ) component_out.level = K;
        else if( temp_level_str == "L" ) component_out.level = L;
        else if( temp_level_str == "M" ) component_out.level = M;
        else if( temp_level_str == "N" ) component_out.level = N;
    }
    int k = temp_type_str.length();
    if( k > 0 ) {
        int f = temp_type_str.find( "INC" );
        if( f > 0 && f < k ) component_out.type = COMPTON;
        f = temp_type_str.find( "COH" );
        if( f > 0 && f < k ) component_out.type = RAYLEIGH;
    }
    return 0;

};


const ElementQuantLevel componentQuantLevel( const SpectrumComponent &component_in ) {
    ElementQuantLevel level_out = NO_QUANT_LEVEL;
    if( component_in.level == K ) level_out = K_LEVEL;
    if( component_in.level == L ) level_out = L_LEVEL;
    if( component_in.level == M ) level_out = M_LEVEL;
    if( component_in.level == N ) level_out = N_LEVEL;
    return level_out;
};

std::string SpectrumComponent_toString(const SpectrumComponent &comp)
{
    ostringstream os;
    os << "SpectrumComponent:" << endl;

    os << "  type: ";
    switch(comp.type)
    {
    case NO_COMPONENT:
        os << "NO_COMPONENT";
        break;
    case ELEMENT:
        os << "ELEMENT";
        break;
    case COMPTON:
        os << "COMPTON";
        break;
    case RAYLEIGH:
        os << "RAYLEIGH";
        break;
    case CONTINUUM:
        os << "CONTINUUM";
        break;
    case SNIP_BKG:
        os << "SNIP_BKG";
        break;
    case PRIMARY_LINES:
        os << "PRIMARY_LINES";
        break;
    case PRIMARY_CONTINUUM:
        os << "PRIMARY_CONTINUUM";
        break;
    case La:
        os << "La";
        break;
    case Lb1:
        os << "Lb1";
        break;
    case DETECTOR_CE:
        os << "Compton Escape";
        break;
    case OPTIC_TRANS:
        os << "Optic Transmission";
        break;
    case PILEUP:
        os << "Pulse pileup";
        break;
    }
    os << endl;

    os << "  element: " << comp.element.toString() << endl;

    os << "  level: ";
    switch(comp.level)
    {
    case NO_EDGE:
        os << "NO_EDGE";
        break;
    case K:
        os << "K";
        break;
    case L:
        os << "L";
        break;
    case M:
        os << "M";
        break;
    case N:
        os << "N";
        break;
    case O:
        os << "O";
        break;
    case P:
        os << "P";
        break;
    case Q:
        os << "Q";
        break;
    }
    os << endl;

    os << "quant: " << comp.quant << endl;

    os << "spectrum: " << floatVecToString(comp.spectrum) << endl;

    os << "coefficient: " << comp.coefficient << endl;
    os << "variance: " << comp.variance << endl;
    os << "intensity: " << comp.intensity << endl;
    os << "residual_err: " << comp.residual_err << endl;
    os << "enabled: " << comp.enabled << endl;
    os << "ignore: " << comp.ignore << endl;
    os << "bkg: " << comp.bkg << endl;
    os << "bkg_index: " << comp.bkg_index << endl;
    os << "fit: " << comp.fit << endl;
    os << "plot: " << comp.plot << endl;
    os << "non-fit factor: " << comp.non_fit_factor << endl;
    os << "matrix effect factor: " << comp.matrix << endl;
    os << "included: " << comp.included << endl;

    return os.str();
}
