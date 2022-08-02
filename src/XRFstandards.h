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

#ifndef XRFSTANDARDS_H_INCLUDED
#define XRFSTANDARDS_H_INCLUDED

#include "XrayMaterial.h"
#include "XraySpectrum.h"
#include "parse_element_list.h"

struct StandardInformation {
    std::vector <std::string> names;
	XrayMaterial mat;
    std::string  spectrumFileName;
    std::vector <std::string> comments;
  //  This will only be used for the first entry, to capture comments before the first standard
    std::vector <std::string> preceding_comments;
    XraySpectrum spectrum;
    std::vector <ElementListEntry> element_list;
    //  Added July 31, 2018 to indicate that user weights were read in from a CSV standards file
    //  Used to reject ECFs from standards with too low concentration or too large uncertainty
    //      when standards were read in from an old TXT standards input file
    bool user_weights = false;
    bool carbonates = false;
    bool input_fractions_are_formula = false;
    bool disable = false;       //  Used to avoid evaluating a standard with itself as a calibration standard
};


#endif // XRFSTANDARDS_H_INCLUDED
