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

#ifndef parse_arguments_h
#define parse_arguments_h

#include <string>
#include <vector>
#include "XRFcontrols.h"


enum PIQUANT_SUBCOMMAND {
    ENERGY_CAL = 1,
    PLOT,
    PRIMARY,
    CALCULATE,
    CALIBRATE,
    QUANTIFY,
    EVALUATE,
    MAP,
    COMPARE,
    FIT_ONE_STANDARD,
    BULK_SUM_MAX,
    EM_SDD_DATA,
    PRINT_VERSION,
    OPTIC_RESPONSE
};

struct ARGUMENT_LIST {
    //  They should always appear in the argument list in this order
    //  Not all arguments will be present for all subcommands
    std::string configuration_file;
    std::string standards_file;
    std::string calibration_file;
    std::string spectrum_file;
    std::string element_list;
    std::string plot_file;
    std::string map_file;
    std::string terminal_text_file;
    std::string invalid_arguments;
    std::string quant_map_outputs;  //  Added mar 2, 2018
    float eV_start = 0; //  Initializers added May 14, 2017
    float eV_ch = 0;
    //  Background arguments Added May 14, 2017
    //  Converted to vector of float July 27, 2018
    std::vector <float> bkg_args;
    //  Add -bh and -bx background options May 10, 2021
    std::vector <float> bh_args;
    std::vector <float> bx_args;
    int detector_select = -1;
    int max_map_arg = -1;
    bool fit_adjust_energy = true;
    bool fit_adjust_width = true;
    bool convolve_Compton = true;
    int map_threads = 1;
    bool standard_selected = false;
    int standard_selection = 0;
    std::string standard_name;
    std::string cal_eval_file;
    bool carbonates = false;
    float min_wgt_eval = MINIMUM_WEIGHT_EVALUATE;
    std::vector <float> detector_shelf_parameters;
    float normalization = 0;
    float iron_oxide_ratio = -1;
};

int parse_arguments( const int argc, const char * argv[],
            PIQUANT_SUBCOMMAND &cmd, ARGUMENT_LIST &arguments );

#endif
