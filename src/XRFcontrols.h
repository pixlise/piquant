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

#ifndef XRFcontrols_h
#define XRFcontrols_h



// settings to control operation of program   Mar. 8, 2017 W. T. Elam  APL/UW

//      flags to control operations (with default values)

//  Minimum amount for standards, elements with amounts below this value will be left out of the standard composition
//#define MINIMUM_AMOUNT 0.001f  //  Value is in percent, 0.001 % or 10 ppm
#define MINIMUM_AMOUNT 0.00f  //  Turn this off since it is controlled via the standards input file and the element list

//  Old values for text calibration files and no element control list
//#define calibration_minimum_fraction 0.01f
//#define calibration_maximum_rsd 0.02f
//#define calibration_maximum_rsd 0.1f
#define calibration_minimum_fraction 0.000009f
#define calibration_maximum_rsd 1e6f    //  Turn this off since it is controlled via the element list
#define calibration_minimum_z  11       //  Leave out everything below Na from text calibration file element list
//#define fp_enable_flag true   //  This is not used in PIQUANT Version 2, fp is always on
#define escape_peaks_enable_flag true
#define peak_tail_enable_flag true
#define detector_shelf_enable_flag true //  See notes in quantCalculate.cpp
#define Compton_escape_enable_flag true
#define composition_normalization_value 0
//#define minimum_quant_Z 10
//#define minimum_quant_Z 19
#define ratio_Compton_flag false

#define MAX_ITERATIONS 40

#define FIT_COEFF_DELTA 0.001f  //  0.1% relative

#define NEGLIGIBLE_FRACTION 1e-8f;

#define MINIMUM_ITERATIONS 3    //  Last one adjusts fit removing negative components

//#define SEC_FLUOR_THRESHOLD 0.001f  //  only calculates sec fluor if exciting element is above this fraction (1 => no sec fluor)
//#define SEC_FLUOR_THRESHOLD 1  //  Turn off sec fluor Jan. 7, 2021 (until quant accuracy gets better)
#define SEC_FLUOR_THRESHOLD 0.05f  //  Try turning back on to see if anything improves

#define FILE_EXTENSION_CHARS 5  //  Characters checked at end of file name to determine extension (.txt, .msa, etc.)

#define MAX_ERROR_MESSAGES 10    //  Avoid too many error messages if wrong file format is opened

#define MINIMUM_WEIGHT_EVALUATE 0.15f    //  Minimum weight for standard element to be included in evaluate
//      can be changed by option input using -w

#define COEFF_RATIO_L_K 1.0f            //  Used to set coefficient of non-fit components, ratio to coefficient of a fit component
//#define COEFF_RATIO_L_K 2.7f            //  Used to set coefficient of non-fit components, ratio to coefficient of a fit component
#define COEFF_RATIO_M_L 1.0f

#define SHELF_THRESHOLD 0   //  Smallest per-channel shelf counts that will be included calculated spectrum
#define SHELF_THRESHOLD_FACTOR 0.00f   //  Determines when detector shelf calculation terminates, when new contributions are less than this max value

//#define BKG_SNIP_CROSSOVER 2500,3000    //  Determines where SNIP background at low energy turns into calculated background
#define BKG_SNIP_CROSSOVER 0

//  Parameters to adjust the shape of the calculated continuum scatter background by applying a linear ramp
#define BKG_RAMP_CONSTANT 0.9    //  Ramp value is this value at zero energy
#define BKG_RAMP_SLOPE 0.004f    //  Ramp starts at zero and increased by this slope every eV in energy

//  Length of list of peaks for pulse pileup calculation (zero => no pileup included)
//      Note that the time taken for the pileup calculation depends on the square of this length
#define PILEUP_LIST_LENGTH  8

#endif
