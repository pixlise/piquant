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

#ifndef XRFconstants_h
#define XRFconstants_h

//	various constants for use in x-ray fluorescence calculations

#define MINIMUM 1.0e-30f
#define MAXIMUM 1.0e+30f
//	From Physics Today Buyer's Guide, August 2001
//	(CODATA recommended values 1998)
#define PI 3.141592653589793f	//	CRC Handbook 1996 page A-1   corrected Aug. 14, 2020, was 3.141592653598793f
#define AVOGADRO 6.02214199e23f
#define HC 12.39841856f	//	keV Angstrom (from h=4.135 667 27 e-15 eV sec; c=299 792 458 meter/sec)
#define FOURPIINV 0.079577472f
#define RADDEG 0.017453293f
#define DEGRAD 57.2957795f
#define ELECTRON_CHARGE 1.602176462E-19f
#define RE2 7.940787e-26f	//	classical electron radius squared [cm2]
#define ME 510998.902f	//	elecron rest mass in eV
#define ALPHA_INV 137.03599976f	//	inverse of fine structure constant
#define EIGHT_LN_2 5.545177f//  ratio of fwhm to sigma in Gaussian is sqrt( 8 ln2 ) ~= 2.35482, usually used when squared
#define SQRT_EIGHT_LN_2 2.35482f//  ratio of fwhm to sigma in Gaussian is sqrt( 8 ln2 ) ~= 2.35482, usually used when squared
#define FOUR_PI 12.566370614359172954f

#define CM_MM 0.1f
#define MM_CM 10.0f
#define NM_CM 0.0000001f     //  nanometers to centimeters
#define CM_MICRON 0.0001f
#define MICRON_CM 10000.0f
#define MICRON_MM 1000.0f
#define GAS_MOLE_VOLUME 22413.996f	//	cm3 = milliLiters
#define PPM_PERCENT 0.0001f //  Convert parts-per-million to percent

#define RESOLUTION_REFERENCE_ENERGY 5984	//	Mn K alpha, usual energy for quoting detector resolution

//	sigma = fwhm * 1 / [ 2 * sqrt( ln(2) ]
#define FWHM_SIGMA 1.66510922f		//	this is [ 2 * sqrt( ln(2) ], not 1 / [ 2 * sqrt( ln(2) ]
#define SIGMA_FWHM 0.6005612f
#define GAUSSIAN_INTEGRAL 1.064463936f  //   Gaussian integral is sqrt(PI/4ln2)*fwhm
#define SIGMA_FWHM 0.6005612f

//		added July 8, 2011  for thin film specimens
#define EXP_FLOAT_TEST	70	//	~= ln(10^-30)   1-exp(-70) = zero for floats
#define THIN_SEC_FLUOR_TEST	1	//	used to turn off secondary fluorescence for very thin films

//  Character definitions   added Jan. 3, 2018
#define TAB_CHARACTER    "\t"
#define SINGLE_QUOTE_CHARACTER    "'"
#define DOUBLE_QUOTE_CHARACTER    "\""
#define UNDERSCORE_CHARACTER    "_"
#define BLANK_CHARACTER    " "
#define SLASH_CHARACTER    "/"
#define BACKSLASH_CHARACTER    "\\"
#define COMMA_CHARACTER    ","
#define EQUAL_CHARACTER    "="
//  These are operating system dependent
#define PATH_SEPARATOR_WINDOWS    "\\"
#define PATH_SEPARATOR_UNIX    "/"
//  This string is used to denote a comment at the beginning of a line in most files
#define COMMENT_STRING "//"

//  This float value is returned from some functions when value is not valid
#define UNLIKELY_VALUE 9.2693e+30f

#endif
