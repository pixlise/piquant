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

//  From Log_OpticResp_Teflon_Rework_FM_Thermal_Dec2020.txt     Dec. 29, 2020  3:30pm
//  Post-rework Thermal Test     Teflon thickness 2 mm    PIQUANT fits to Teflon spectra
//  This response is from Yellow Piece temperature 12.5C (Elemental Calibration was at 10C)
//  Temperature dependence was small, less than about 10% at highest energies from 12.5C to 30C

//  Recalculated Apr. 28, 2021 with Be window thickness of 157 microns, for unity ECFs
//      Started with fit to Teflon spectrum above, then increased zero energy value to 2.2x 4 keV value, to get Na thru Cl ECFs to unity
//      30 keV value made much smaller and derivative zeroed to avoid problems at highest energies (fit to SNIP bkg in Teflon spectrum gets this better)

const int N_FM_OpticResp = 12;
const float X_FM_OpticResp[N_FM_OpticResp] = {  0,  4000,  6000,  8000,  10000,  12000,  14000,  16000,  18000,  20000,  25000,  30000 };
//const float Y_FM_OpticResp[N_FM_OpticResp] = {  15.7810,  6.8613,  9.3393,  6.0398,  4.0440,  2.1134,  1.2380,  0.7581,  0.5089,  0.3082,  0.0117,  0.0117 };
//const float D_FM_OpticResp[N_FM_OpticResp] = {  -1.3604e-06,  2.7208e-06,  -3.1971e-06,  1.4014e-06,  -4.5320e-07,  5.0906e-07,  -1.7288e-10,  8.4921e-08,  6.5435e-09,  -3.8315e-08, 0, 0 };
//      Tweaked by hand to get Sr thru Zr ECFs to unity and for better fit to Teflon bkg above 25 keV
//const float Y_FM_OpticResp[N_FM_OpticResp] = {  15.7810,      6.8613,      9.3393,       6.0398,      4.0440,       2.1134,      1.2380,       0.6,         0.377,       0.3,    0.04,   0.0117 };
//      Tweaked 6 keV value to get Ca closer and better bkg under Cr
const float Y_FM_OpticResp[N_FM_OpticResp] = {  15.7810,      6.8613,      7.8,       6.0398,      4.0440,       2.1134,      1.2380,       0.6,         0.377,       0.3,    0.04,   0.0117 };
const float D_FM_OpticResp[N_FM_OpticResp] = {  -1.3604e-06,  2.7208e-06,  -3.1971e-06,  1.4014e-06,  -4.5320e-07,  5.0906e-07,  -1.7288e-10,  8.4921e-08,  6.5435e-09,  0,      0,      0 };



//--------------------------------------------------------------------------------------------------
//  The optic response below this line (PIXL_FM_OPTIC_OLD, OpticFile #7) was incorrectly calculated with wrong Be window thickness for flight X-ray tube

//  Matched by hand and eye above 20 keV to Teflon spectrum (30 keV value taken from 25 keV, derivatives zeroed)
//const int N_FM_OpticResp = 12;
//const float X_FM_OpticResp[N_FM_OpticResp] = {  0,            4000,        6000,         8000,        10000,        12000,       14000,       16000,       18000,       20000,        25000,       30000 };
//const float Y_FM_OpticResp[N_FM_OpticResp] = {  3.9338,       7.8675,      10.0934,      6.4421,      4.2282,       2.2448,      1.3170,      0.8082,      0.5341,      0.3270,       0.0105,      0.0105 };
//const float D_FM_OpticResp[N_FM_OpticResp] = {  -3.0708e-07,  6.1416e-07,  -2.6823e-06,  1.2992e-06,  -3.5853e-07,  4.8075e-07,  1.8820e-08,  7.2626e-08,  4.2546e-08,  0,            0,           0 };

//  Adjusted optic response to get unity fit coefficients for elements in 5-minute elemental calibration spectra (pure compounds)
//      Below 3 keV, multiply by 1.85, from 3 keV to 12 keV multiply by 0.89, above 12 keV ~linearly decreasing to 0.76 at 18 keV   (factor = 0.89 - 0.13 * ( E - 12 ) / 6 )
const int N_FM_OpticResp_7 = 14;
const float X_FM_OpticResp_7[N_FM_OpticResp_7] = {  0,            2999,        3000,         4000,        6000,         8000,        10000,        12000,       14000,       16000,       18000,       20000,        25000,       30000 };
//const float Y_FM_OpticResp[N_FM_OpticResp_7] = {  7.29,        14.56,        5.68,         7.00,        8.90,        5.73,        3.76,         2.00,        1.13,        0.649,       0.406,       0.234,        0.0105,      0.0105 };
//  April 2, 2021    manual trim of values at 6 & 8 keV to remove some structure, didn't affect quantification of FM Elemental Calibration (all standards checked)
const float Y_FM_OpticResp_7[N_FM_OpticResp_7] = {  7.29,        14.56,        5.68,         7.00,        8.00,        6.00,        3.76,         2.00,        1.13,        0.649,       0.406,       0.234,        0.0105,      0.0105 };
const float D_FM_OpticResp_7[N_FM_OpticResp_7] = {  0,           0,            -3.0708e-07,  6.1416e-07,  -2.6823e-06,  1.2992e-06,  -3.5853e-07,  4.8075e-07,  1.8820e-08,  7.2626e-08,  4.2546e-08,  0,            0,           0 };
