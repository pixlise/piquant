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

#ifndef parse_element_list_h
#define parse_element_list_h

#include <string>
#include <vector>
#include "Element.h"
#include "XrayMaterial.h"

//  Modified Oct. 27, 2017
//      Add information for holding standards description to element list entries
//  Modified Dec. 6, 2017
//      Add element calibration factors
//      Add function to append or replace entry in element list, keeping selected info
//  Modified Dec. 8, 2017
//      Add MATRIX qualifier and SpectrumComponentType
//  Modified Dec. 16, 2019
//      Add element qualifier OUTPUT ("O") to force the element to be included in the evaluate list (with zeros if not in any standard in this run)
//  Modified Dec. 15, 2020
//      Change oxide ratios from float to LightElements struct from XrayMaterial.h
//  Modified Feb. 26, 2021  Add DETECTOR_SHELF and DETECTOR_CE spectrum component types

//  This is here to prevent circular header references
enum SpectrumComponentType { NO_COMPONENT = -1, ELEMENT, COMPTON, RAYLEIGH, CONTINUUM, SNIP_BKG, PRIMARY_LINES, PRIMARY_CONTINUUM, La, Lb1, DETECTOR_CE, OPTIC_TRANS, PILEUP };

enum ElementQuantLevel { NO_QUANT_LEVEL = -1, K_LEVEL, L_LEVEL, M_LEVEL, N_LEVEL };
enum ElementQualifiers { NO_QUALIFIER = 0, IGNORE, FORCE, EXCLUDE, MATRIX, OUTPUT };
struct ElementListEntry {
    Element element;
    ElementQuantLevel quant_level = NO_QUANT_LEVEL;
    ElementQualifiers qualifier = NO_QUALIFIER;
    //  The following entries just serves as a convenient place to put this information
    SpectrumComponentType type = NO_COMPONENT;
    float percent = -1; //  No percent entered for this element
//    float oxide_ratio = 0; //  No oxide ratio entered for this element
    LightElements stoichiometry; //  Replaces oxide ratio so carbonates can be specified
    float uncertainty = 0;  //  Relative error of given element percent (expressed as percent)
    float weight = 1;   //  Default is equal weights
    float ecf = -1;   //  Element calibration factor for this element and line, -1 => none available
    float ecf_sigma = 0;   //  Uncertainty in ecf (expressed as relative percent)
    float intensity = 0;   //  Net peak intensity, added Jan. 26, 2018
    float coefficient = -1;   //  Spectrum fit coefficient, added Mar. 2, 2018
    float rel_err_coeff = 0;   //  Standard deviation (sqrt of variance) of spectrum fit coefficient (expressed as relative percent), added Mar. 2, 2018
    float total_err = 0;   //  To include ECF standard deviation in error output (and later certificate uncertainty), added Nov. 4, 2019
    float given = 0;    //  For use during evaluate, to pass given element percent to output (added May 21, 2019)
    float rel_err_given = 0;    //  For use during evaluate (added May 21, 2019)
    float matrix = 0;
    //  If any more entries are added, function add_element_list_entry must be modified to include them in replacement test
};

const bool parse_element_list( const std::string &element_list_in, std::vector <ElementListEntry> &elements_out, bool &carbonates, bool oxides = true );

const bool parse_element_string( const std::string &element_string_in, ElementListEntry &element_entry_out );

void add_element_list_entry( const ElementListEntry &element_entry, std::vector <ElementListEntry> &element_list_out,
            const bool ignore_qualifier = true );

//  Returns true if any errors, false if none

#endif

