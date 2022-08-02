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

#ifndef quantComponents_h
#define quantComponents_h

#include <vector>
#include <string>
#include "Element.h"
#include "XrayLines.h"  //  Defines EdgeLevel
#include "XrayMaterial.h"
#include "XRFconditions.h"
#include "parse_element_list.h" //  Defines ElementQualifiers


//	Sets up the list of components that sum to make the full spectrum
//  Sets and checks the components that are used to quantify each element
//  Maps XrayLines objects to each component,
//      matching associated elements and physical processes that generate the signal
//  Also generates and parses the text descriptions of the components
//      for output to and input from calibration files
//  Do all this here so it's in one place if it is changed to improve quantification

//  Maybe convert this into a class someday, with == and > operators to replace matchComponent and checkComponent


#define COEFFICIENT_NO_COMPONENT -9999  //  Special value for coefficient return to indicate that there is no coefficient

struct SpectrumComponent {
    SpectrumComponentType type = NO_COMPONENT;  //  Defined in parse_element_list.h (to avoid circular headers)
    //  Information for mapping associated elements to components
    Element element;
    //  Information for mapping lines to components
    EdgeLevel level = NO_EDGE;
    bool quant = false;     //  if true this is the component used to quantify its associated element
//      Notes for possible future changes
//  A combined component will be specified in the element qualifiers, so when to set this flag is problematic right now
//  If it is set after the components are put into the vector using makeComponents as it is now, then the extra ones will have to be removed
//  A component will include more that one XrayLines object
//  This is all done by the generation in makeComponents and the tests in checkComponent
//  Note that any new flags or other variables have to be initialized by assigning default values here
//      and added to the matchComponent function for checking matching components
//  The code in checkComponent has to be changed to properly handle checking for matches
//  For example, set this flag if we DON'T want to decouple K, L, M or N lines
//    bool combined;
//      (or, maybe in the future, part of an XrayLines object, like alpha or beta lines)
//    bool alpha;
//    bool beta;
    //  storage for the computed contribution of this component to the full spectrum
    std::vector <float> spectrum;
    float coefficient = 1;  //  coefficient from fit to measured spectrum
    float variance = 0;  //  variance of coefficient
    float intensity = 0;  //  Integrated intensity of this component
    float residual_err = 0;  //  Contribution to uncertainty from fit residual for this component
    bool enabled = true;
    bool ignore = false;
    float matrix = 0;   //  To hold matrix effect factor from FP calculation so it can be written to map file if desired
    bool bkg = false;    //  Include this component in the spectrum background
    int bkg_index = 0;  //  For splitting the background function into several components for independent fitting
    bool plot = true;    //  Include this component in the plot (separated from enabled for debugging and other special cases, set with disable in XraySpectrum)
    bool fit = true;    //  Added so that L (or M) components could be fixed relative to K (or L) components
    float scale_under = 0;    //  Added to use scale-under-peaks algorithm for calculated background (every time it is calculated in quantCalculate)
    float non_fit_factor = 0;   //  Used to set coefficient of non-fit components, ratio to coefficient of a fit component
    float adjusted_coefficient = -1;   //  Used to get coefficients to better match updated concentration for next calculation (used in detector shelf calculation)
    bool included = true;   //  Only to be used in XraySpectrum to form fit vector, will be set there as needed

};

std::string SpectrumComponent_toString(const SpectrumComponent &comp);


int setupComponents(  const std::vector <XrayLines> sourceLines, const std::vector <XrayLines> pureLines,
                std::vector <SpectrumComponent> &components_out );

int makeComponents( const SpectrumComponentType type_in, const std::vector <XrayLines> &lines_in,
                std::vector <SpectrumComponent> &components_out, const int n_bkg = 1 );

int quantComponents( const std::vector <ElementListEntry> element_list_in,
				std::vector <SpectrumComponent> &components_out );

int quantDefaults( std::vector <ElementListEntry> &element_list_in,
				std::vector <SpectrumComponent> &components_out );

bool checkComponent( const SpectrumComponent &component_in, const XrayLines &lines_in, const int line_index_in );

bool matchComponent( const SpectrumComponent &component_1, const SpectrumComponent &component_2 );

std::string componentDescription( const SpectrumComponent &component_in );

int parse_component( const std::string component_string_in,
				SpectrumComponent &component_out );

const ElementQuantLevel componentQuantLevel( const SpectrumComponent &component_in );

#endif
