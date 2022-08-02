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

#include <string>
#include <math.h>
#include "Element.h"
#include "XrayEdge.h"
#include "XrayLines.h"
#include "interp.h"
#include <sstream>

using namespace std;

//	Elam Ravel Sieber database class

//	the following three constants are used only to translate IUPAC to Siegbahn notation
//	they cannot be used to obtain line indices (which vary among the different edges)

//  Modified July 20, 2018
//      Relative intensities for K lines changed to use values from Scofield
//      Had to change how list of possibilities to check was created
//      Old way: call XrayEdge::number (checked edge energies for entry)
//      New way: call XrayEdge::numberOccupied (checks electron occupancy)
//      This also is consistent with strongly typed enum EdgeIndex


#define MAXLINES 27
const string XrayLines::LINE_NAMES_Siegbahn[MAXLINES] = {
	"Ka1",	"Ka2",	"Ka3",	"Kb1",	"Kb2",		"Kb3",	"Kb4",		"Kb5",
	"Lb3",		"Lb4",		"Lg2",		"Lg3",
	"Lb1",		"Ln",		"Lg1",		"Lg6",
	"La1",		"Lb2,15",	"La2",		"Lb5",		"Lb6",		"Ll",
	"Ma",		"Mb",		"Mg",		"Mz",		"M2-N4"
};
const string XrayLines::LINE_NAMES_IUPAC[MAXLINES] = {
	"K-L3",	"K-L2",	"K-L1",	"K-M3",	"K-N2,3",	"K-M2",	"K-N4,5",	"K-M4,5",
	"L1-M3",	"L1-M2",	"L1-N2",	"L1-N3",
	"L2-M4",	"L2-M1",	"L2-N4",	"L2-O4",
	"L3-M5",	"L3-N4,5",	"L3-M4",	"L3-O4,5",	"L3-N1",	"L3-M1",
	"M5-N6,7",	"M4-N6",	"M3-N5",	"M4,5-N2,3",	"M2-N4"
};

//	this is for returning a blank string when no name was available
const string XrayLines::blank("");


XrayLines::XrayLines ( ) {
//		default constructor just allocates variables, for vector allocation
	lineCount = 0;
//		let default constructors take care of everything else
};

XrayLines::XrayLines ( const float energy ) {
//		constructor for single-entry line
//		single entry at given energy, relative intensity unity, blank name, H K edges
	Element fakeH( 1 );
	XrayEdge fakeEdge ( fakeH, K1 );
	edgeVacant = fakeEdge;
	lineCount = 0;
    XrayLinesInfo temp;
	temp.edgeOccupied = fakeEdge;
	temp.lineIUPAC = " " ;
	temp.lineEnergy = energy;
	temp.lineIntensity = 1;
	temp.lineFactor = 1;
    lineList.push_back( temp );
	lineCount = 1;
};

//      modified constructor to produce subset of available lines for special fits   Aug 20, 2012
XrayLines::XrayLines ( const XrayEdge& newEdge, const float separation,
                      const float energy_low_limit, const float energy_high_limit ) {
//		constructor for list of x-ray emission lines from a vacancy
//			in the energy level associated with an absorption edge
//		keeps energies and relative intensities in memory to save search and
//			interpolation time during repetitive calculations
	edgeVacant = newEdge;
	lineCount = 0;
//		find the number of edges with energy below this edge
	vector<EdgeIndex>  possibilities;
    newEdge.numberOccupied( possibilities, newEdge.element() );
    //  Fix up a couple of differences between occupancy table and Scofield configuration
    //  Scofield assumes that Al and Si have K-M3 (Kb1) lines even though the occupancy table only has M1 and M2 electrons
    if( newEdge.element().Z() == 13 || newEdge.element().Z() == 14 ) {
        possibilities.push_back( M3 );
    }
    float sum = 0;
    int i;
    for ( i=0; i<possibilities.size(); i++ ) {
//				check possible edges for transitions
        XrayEdge testEdge ( edgeVacant.element(), possibilities[i] );
//				see if this line has any intensity (also get IUPAC symbol)
        string testSymbol;
//				intSymbol returns the intensity of a line and its IUPAC symbol,
//					given the two edges between which the transition occurs
        float testInt = intSymbol ( testSymbol, edgeVacant, testEdge );
//        cout << testEdge.symbol() << "  " << testInt << "  " << testSymbol << endl;
        if ( testInt > 0.0 ) {
//				create an entry for an emission line, transition from occupied edge level
//					to vacant edge level
            XrayLinesInfo temp;
            temp.edgeOccupied = testEdge;
            temp.lineIUPAC = testSymbol;
            temp.lineEnergy = edgeVacant.energy()-testEdge.energy();
            temp.lineIntensity = testInt;
            temp.lineFactor = 1.0;
            if( energy_low_limit <= temp.lineEnergy
                   && temp.lineEnergy <= energy_high_limit ) {
                sum += temp.lineIntensity;
                lineList.push_back( temp );
            };
        };
    };
    lineCount = lineList.size();
//			normalize relative intensities
    for ( i=0; i<lineCount; i++ ) {
        lineList[i].lineIntensity /= sum;
    };
    if ( separation > 0 ) merge_peaks( separation, lineList );
    lineCount = lineList.size();
};

XrayLines::XrayLines ( const XrayLines& newLines, const float separation ) {
	edgeVacant = newLines.edge();
	lineCount = newLines.numberOfLines();
	lineList.resize(lineCount);
	int i;
	for ( i=0; i<lineCount; i++ ) {
		lineList[i].edgeOccupied = newLines.edgeSource(i);
		lineList[i].lineIUPAC = newLines.symbolIUPAC(i);
		lineList[i].lineEnergy = newLines.energy(i);
		lineList[i].lineIntensity = newLines.relative(i);
		lineList[i].lineFactor = newLines.factor(i);
	};
    if ( separation > 0 ) merge_peaks( separation, lineList );
    lineCount = lineList.size();
};

const float XrayLines::energy (  const string& symbol ) const {
//		return the energy of the line whose IUPAC or Siegbahn designation
//			matches the given symbol (check IUPAC first)
	int i;
	i = (*this).index(symbol);
	if ( i >= 0 ) return (*this).energy(i);
	return 0.0;
};

const string& XrayLines::symbolSiegbahn ( const int index ) const {
	int i;
	string test = (*this).symbolIUPAC(index);
	for ( i=0; i<MAXLINES; i++ ) {
		if ( test == LINE_NAMES_IUPAC[i] ) return LINE_NAMES_Siegbahn[i];
	};
	return blank;
};


const int XrayLines::index ( const string& symbol ) const {
//		finds the line whose IUPAC or Siegbahn designation matches the given symbol
	int i;
// check IUPAC symbols first
	for ( i=0; i<lineCount; i++ ) {
		if ( symbol == (*this).symbolIUPAC(i) ) return i;
	};
// now check Siegbahn symbols
	for ( i=0; i<lineCount; i++ ) {
		if ( symbol == (*this).symbolSiegbahn(i) ) return i;
	};
	return -1;
};

//          Line width from width of energy levels
const float XrayLines::width ( const int index ) const {
    if( index < 0 || index >= lineList.size() ) return 0.00001f;
    float w1 = edgeVacant.width();
    float w2 = lineList[index].edgeOccupied.width();
    float wid = sqrt( w1*w1 + w2*w2 );
    return wid;
};


float XrayLines::intSymbol ( std::string& symbol, const XrayEdge& upper, const XrayEdge& lower ) {
//	Finds the relative intensity and IUPAC symbol

//		The indices for each edge are defined in class XrayEdge in an enum

//  July 20, 2018
//  Relative intensities for K lines from  James H. Scofield, "Exchange corrections of K x-ray emission rates",
//      Phys. Rev. A 9 (3), March 1974, 1041-49.  (Table V on page 1074)
//  The data for Cr (Z=24) and Cu (Z=29) have been modified to fit in line with the trends in Z by interpolating between the adjacent elements
//      See the article by Iain Campbell for the Group 4 report of the Fundamental parameters initiative
//      Published in IRPS Bulletin (Newsletter of the International Radiation Physics Society) Vol 24 No 1 pp17-30.
//      This recommendation is on page 21, first full paragraph in the left column.
//      It is attributed to Schönfeld, E. and Janβen, H. Physikalische-Technische Bundesanstalt Report PTB-Ra-37 (1995).
//          and Schönfeld, E. and Janβen, H. Nucl. Instrum. Meth. A369 (1996) 527.
//	Use linear interpolation between values from table
const int NUMBER_Z_K_L2 = 50;
const float DATA_Z_K_L2[50] = { 10,  13,    14,    15,    16,    17,    18,    19,    20,    22,    23,    24,    25,    26,    28,    29,    30,    32,    33,    34,    35,    36,    37,    38,    40,    42,    47,    50,    51,    54,    56,    60,    63,    64,    65,    68,    70,    72,    73,    74,    78,    79,    80,    81,    82,    85,    90,    92,    96,    98 };
const float DATA_K_L2[50] = { 0.5028,0.5033,0.5037,0.5048,0.5053,0.5056,0.5049,0.5055,0.5061,0.5076,0.5083,0.5091,0.5099,0.5107,0.5124,0.5133,0.5142,0.5149,0.5153,0.5158,0.5181,0.5186,0.5195,0.5205,0.5225,0.5247,0.5305,0.5343,0.5356,0.5398,0.5428,0.5491,0.5542,0.5559,0.5577,0.5634,0.5673,0.5714,0.5736,0.5757,0.585, 0.5874,0.5899,0.5924,0.595,0.6033,0.6182,0.6247,0.6387,0.6462 };
//  Note that Scofield assumes that Al and Si have K-M3 (Kb1) lines even though the occupancy table only has M1 and M2 electrons
const int NUMBER_Z_K_M = 49;
const float DATA_Z_K_M[NUMBER_Z_K_L2] = {       13,    14,    15,    16,    17,    18,    19,    20,    22,    23,    24,    25,    26,    28,    29,    30,    32,    33,    34,    35,    36,    37,    38,    40,    42,    47,    50,    51,    54,    56,    60,    63,    64,    65,    68,    70,    72,    73,    74,    78,    79,    80,    81,    82,    85,    90,    92,    96,    98 };
//const float DATA_Kb_OVER_Ka[NUMBER_Z_K_M] = { 0.0134,0.0294,0.0472,0.0659,0.0862,0.1088,0.1211,0.1315,0.1355,0.1367,0.1337,0.1385,0.1391,0.1401,0.1379,0.141,0.1504,0.156,0.1624,0.1683,0.1727,0.178,0.1831,0.1913,0.1981,0.213,0.223,0.2266,0.2368,0.2433,0.2504,0.2549,0.257,0.2575,0.2612,0.2634,0.2666,0.2682,0.2698,0.2758,0.2772,0.2788,0.2804,0.2821,0.2873,0.2952,0.2975,0.3019,0.3037 };
const float DATA_Kb3_OVER_Kb1[NUMBER_Z_K_M] = { 0.5057,0.5052,0.5048,0.5047,0.5041,0.5041,0.5042,0.5043,0.5054,0.506, 0.507, 0.5073,0.5079,0.5093,0.5105,0.5108,0.5105,0.5113,0.5116,0.5116,0.5111,0.5113,0.5115,0.512, 0.5125,0.5138,0.5148,0.5151,0.5157,0.516, 0.5167,0.517, 0.5171,0.5171,0.517, 0.5175,0.5176,0.5176,0.5176,0.5173,0.5172,0.517, 0.5167,0.5165,0.5158,0.5134,0.5122,0.509, 0.507 };
const float DATA_KbM_OVER_Ka1[NUMBER_Z_K_M] = { 0.0201,0.0443,0.071, 0.0992,0.1298,0.1638,0.1824,0.1982,0.2043,0.2063,0.2077,0.2092,0.2102,0.2119,0.2127,0.2135,0.2229,0.2277,0.2331,0.2372,0.2381,0.2423,0.2463,0.2543,0.2617,0.2775,0.2857,0.2882,0.2951,0.2997,0.3086,0.3147,0.3166,0.3185,0.324, 0.3274,0.3307,0.3323,0.3338,0.3399,0.3414,0.343, 0.3444,0.3459,0.3503,0.3577,0.3606,0.3665,0.3695 };
const int NUMBER_Z_K_N = 33;
const float DATA_Z_K_N[NUMBER_Z_K_N] = {        32,    33,    34,    35,    36,    37,    38,    40,    42,    47,    50,    51,    54,    56,    60,    63,    64,    65,    68,    70,    72,    73,    74,    78,    79,    80,    81,    82,    85,    90,    92,    96,    98 };
const float DATA_KbN_OVER_Ka1[NUMBER_Z_K_N] = { 0.0049,0.0086,0.0131,0.0183,0.024 ,0.0281,0.032, 0.037, 0.0403,0.0484,0.0564,0.0597,0.0695,0.0756,0.0792,0.0813,0.0832,0.0826,0.0843,0.0853,0.0883,0.0898,0.0913,0.0972,0.0987,0.1004,0.1023,0.1043,0.1105,0.1205,0.1233,0.129, 0.1315 };

//  *************  Relative intensities for K lines changed to use values from Scofield, July 20, 2018 ********************
//	Relative emission rates, fits from Kaleidagraph, low-Z extrapolations by hand and eye
//	data from Salem, Panossian, and Krause, Atomic Data and Nuclear Data Tables Vol. 14 No.2 August 1974, pp92-109.
//	M shell data is from T. P. Schreiber and A. M. Wims, X-ray Spectrometry Vol. 11, No. 2, 1982, pp42-45.

//	Arrays contain:				(except for K-N32 transitions, which contain Z & value pairs)
//		value for Z=0
//		slope vs Z for extrapolation below Zmin of polynomial fit
//		Zmin for polynomial fit
//		Zmax for polynomial fit
//		5 polynomial coefficients, for Z**0 to Z**4

//	K relative emission rates (relative to K-L3 (Ka1) as 1.0)

//const float DATA_K_L2[9] = { 0.502f, 0.0f, 22.0f, 100.0f, 0.47659f, 0.0013468f, -1.1051e-05f, 1.4555e-07f, 0.0f };	//	Ka2

const float DATA_K_L1[9] = { 0.0f, 1.83e-5f, 60.0f, 92.0f, -0.013296f, 0.00062767f, -1e-05f, 5.4395e-08f, 0.0f };	//	Ka3

//const float DATA_K_M3[9] = { 0.107f, 1.33e-3f, 64.0f, 100.0f, -1.0379f, 0.041531f, -0.00045837f, 1.7078e-06f, 0.0f };	//	Kb1

//	K-N32, Z=36 to 100	//	Kb2
//	Use actual values, linearly extrapolate to zero at Z=30
//	(structure shows up in matrix-corrected Scofield calculations)
//const float DATA_Z_K_N32[33] = {
//	36.0f, 38.0f, 40.0f, 42.0f, 44.0f, 46.0f, 48.0f, 50.0f, 52.0f, 54.0f, 56.0f, 58.0f, 60.0f, 62.0f, 64.0f, 66.0f, 68.0f, 70.0f, 72.0f, 74.0f,
//	76.0f, 78.0f, 80.0f, 82.0f, 84.0f, 86.0f, 88.0f, 90.0f, 92.0f, 94.0f, 96.0f, 98.0f, 100.0f
//};
//const float DATA_K_N32[33] = {
//	0.019f, 0.030f, 0.037f, 0.041f, 0.045f, 0.048f, 0.053f, 0.055f, 0.058f, 0.064f,
//	0.070f, 0.076f, 0.083f, 0.086f, 0.089f, 0.089f, 0.088f, 0.087f, 0.085f, 0.086f,
//	0.087f, 0.091f, 0.096f, 0.102f, 0.108f, 0.113f, 0.117f, 0.120f, 0.123f, 0.125f,
//	0.128f, 0.132f, 0.135f
//};

const float DATA_K_N45[9] = { 0.0f, 1.32e-5f, 64.0f, 92.0f, 0.0073753f, -0.00018264f, 6.6587e-07f, 9.3221e-09f, 0.0f };	//	Kb4

const float DATA_K_M45[9] = { 0.0f, 4.72e-5f, 64.0f,  92.0f, -0.0098129f, 0.00020066f, 0.0f, 0.0f, 0.0f };	//	Kb5

//	This is the ratio of K-M2 (Kb3) to K-M3 (Kb1)
//const float DATA_K_M2_M3[9] = { 0.518f, 0.0f, 48.0f, 100.0f, 0.56232f, -0.0014291f, 1.1312e-05f, 0.0f, 0.0f };	//	Kb3/Kb1


//	L1 relative emission rates (relative to L1-M3 (Lb3) as 100.0)

const float DATA_L1_M2[9] = { 70.6f, 0.0f, 42.0f, 96.0f, 201.57f, -4.8189f, 0.039388f, 2.3756e-05f, 0.0f };	//	Lb4

const float DATA_L1_N2[9] = { 9.6f, 0.15f, 64.0f,  96.0f, -76.165f, 5.2791f, -0.096091f, 0.00057673f, 0.0f };	//	 Lg2

const float DATA_L1_N3[9] = { 8.31f, 0.281f, 32.0f,  96.0f, 16.793f, -0.11096f, 0.0042381f, 0.0f, 0.0f };		//	Lg3


//	L2 relative emission rates (relative to L2-M4 (Lb1) as 100.0)

const float DATA_L2_M1[9] = { 13.3f, -0.203f, 28.0f, 96.0f, 18.917f, -0.56982f, 0.006189f, -2.087e-05f, 0.0f };	//	Ln

const float DATA_L2_N4[9] = { -23.6f, 0.672f, 40.0f,  96.0f, -175.28f, 9.5323f, -0.18014f, 0.0015239f, -4.7415e-06f };	//	Lg1

const float DATA_L2_O4[9] = { 0.0f, 0.0f, 74.0f, 96.0f, -157.55f, 4.3458f, -0.037263f, 0.00010052f, 0.0f };	//	Lg6


//	L3 relative emission rates (relative to L3-M5 (La1) as 100.0)

//	Note that the L3_N45 fit is split into two Z ranges
const float DATA_L3_N45_1[9] = { 0.0f, 0.0f, 40.0f, 70.0f, -259.32f, 11.946f, -0.16561f, 0.00074045f, 0.0f };	//	Lb2,15 (part 1)
const float DATA_L3_N45_2[9] = { 0.0f, 0.0f, 70.0f, 96.0f, -3165.7f, 147.79f, -2.5697f, 0.01986f, -5.7501e-05f };	//	Lb2,15 (part 2)

const float DATA_L3_M4[9] = { 11.0f, 0.0f, 40.0f, 96.0f, 11.052f, 0.0014163f, 0.0f, 0.0f, 0.0f };	//	La2

const float DATA_L3_O45[9] = { 0.0f, 0.0f, 72.0f, 96.0f, -44.376f, 0.88467f, -0.0037163f, 0.0f, 0.0f };	//	Lb5

const float DATA_L3_N1[9] = { -0.985f, 0.031f, 60.0f, 96.0f, -1.0706f, 0.032066f, 0.0f, 0.0f, 0.0f };	//	Lb6

// This is L3-M1 to L3-M45 ratio times 100 (Ll/La)*100, in two Z ranges
const float DATA_Ll_La_1[9] = { 11.0f, 0.0f, 26.0f, 40.0f, 826.43f, -92.536f, 3.93f, -0.074403f, 0.00052853f };	//	Ll (part 1)
const float DATA_Ll_La_2[9] = { 0.0f, 0.0f, 40.0f, 96.0f, 14.145f, -0.47213f, 0.0070206f, -4.1231e-05f, 1.126e-07f };	//	Ll (part 2)

//		calculate relative emisson rates for lines using data in above arrays
//		Salem, Panossian, and Krause, Atomic Data and Nuclear Data Tables Vol. 14 No.2 August 1974, pp92-109.
//			(use helper routine to evaluate polynomials stored in some arrays)

//		select the line based on the indices of the two levels
	float z = float( upper.element().Z() );
	switch ( upper.index() ) {
		case K:
			switch (lower.index()) {
				case L1:		symbol = "K-L1";	return linePolyCalc ( z,  DATA_K_L1 );	// K-L1   Ka3   continue to use Salem value, no value from Scofield
//				case L2:		symbol = "K-L2";	return linePolyCalc ( z,  DATA_K_L2 );	// K-L2   Ka2 //  Old Salem value, not used
				case L2:		symbol = "K-L2";	return interp ( z, DATA_Z_K_L2, DATA_K_L2, NUMBER_Z_K_L2 );	// K-L2   Ka2   Scofield value
				case L3:		symbol = "K-L3";	return 1.0;								// K-L3   Ka1
//				case M3:		symbol = "K-M3";	return linePolyCalc ( z,  DATA_K_M3 );	// K-M3   Kb1 //  Old Salem value, not used
                                //  Kb'1 is sum of K-M transitions, so solve for Kb1
				case M3:	{	symbol = "K-M3";                                            // K-M3   Kb1
                                float int_Kb5 = linePolyCalc ( z,  DATA_K_M45 );	// (K-M4,5 Kb5   use Salem value)
                                float int_K_M = interp ( z, DATA_Z_K_M, DATA_KbM_OVER_Ka1, NUMBER_Z_K_M );
                                float ratio_Kb3_over_Kb1 = interp ( z, DATA_Z_K_M, DATA_Kb3_OVER_Kb1, NUMBER_Z_K_M );
                                return ( int_K_M - int_Kb5 ) / ( 1 + ratio_Kb3_over_Kb1 );
                            };
				case M2:	{	symbol = "K-M2";											// K-M2   Kb3   Scofield value
//										Scofield table gives the ratio of K-M2 (Kb3) to K-M3 (Kb1), so calculate Kb1 as above then use ratio
                                float int_Kb5 = linePolyCalc ( z,  DATA_K_M45 );	// (K-M4,5 Kb5   use Salem value)
                                float int_K_M = interp ( z, DATA_Z_K_M, DATA_KbM_OVER_Ka1, NUMBER_Z_K_M );
                                float ratio_Kb3_over_Kb1 = interp ( z, DATA_Z_K_M, DATA_Kb3_OVER_Kb1, NUMBER_Z_K_M );
                                float int_Kb1 = ( int_K_M - int_Kb5 ) / ( 1 + ratio_Kb3_over_Kb1 );
                                return int_Kb1 * ratio_Kb3_over_Kb1;
                            };
//										data is the ratio of K-M2 (Kb3) to K-M3 (Kb1)
//                                return linePolyCalc ( z,  DATA_K_M2_M3 ) * linePolyCalc ( z,  DATA_K_M3 );    //  Old Salem value, not used
				case M4:		symbol = "K-M4,5";	return linePolyCalc ( z,  DATA_K_M45 );	// K-M4,5 Kb5   continue to use Salem value, no value from Scofield
				case N2:	{    symbol = "K-N2,3";											// K-N2,3 Kb2
                                //  Kb'2 is sum of K-N transitions, so solve for Kb2
                                float int_Kb4 = linePolyCalc ( z,  DATA_K_N45 );	// K-N4,5 Kb4   (use Salem value, no value from Scofield)
                                float int_K_N = interp ( z, DATA_Z_K_N, DATA_KbN_OVER_Ka1, NUMBER_Z_K_N );
                                return int_K_N - int_Kb4;   //  Subtract Kb4 to get Kb2 from total of K-N transitions
                            };
//						use linear interpolation between values stored in array (too much structure for poly fit)
//					return interp ( z, DATA_Z_K_N32, DATA_K_N32, 33 ) ; //  Old Salem value, not used
				case N4:	symbol = "K-N4,5";	return linePolyCalc ( z,  DATA_K_N45 );	// K-N4,5 Kb4   (use Salem value, no value from Scofield)
				default: return 0.0;
			};
		case L1:
			switch (lower.index()) {
				case M3:		symbol = "L1-M3";	return 100.0;							// L1-M3   Lb3
				case M2:		symbol = "L1-M2";	return linePolyCalc ( z,  DATA_L1_M2 );	// L1-M2   Lb4
				case N2:	symbol = "L1-N2";	return linePolyCalc ( z,  DATA_L1_N2 );	// L1-N2   Lg2
				case N3:	symbol = "L1-N3";	return linePolyCalc ( z,  DATA_L1_N3 );	// L1-N3   Lg3
				default:	return 0.0;
			};
		case L2:
			switch (lower.index()) {
				case M4:		symbol = "L2-M4";	return 100.0;							// L2-M4   Lb1
				case M1:		symbol = "L2-M1";	return linePolyCalc ( z,  DATA_L2_M1 );	// L2-M1   Ln
				case N4:	symbol = "L2-N4";	return linePolyCalc ( z,  DATA_L2_N4 );	// L2-N4   Lg1
				case O2:	symbol = "L2-O4";	return linePolyCalc ( z,  DATA_L2_O4 );	// L2-O4   Lg6
				default:	return 0.0;
			};
		case L3:
			switch (lower.index()) {
				case M5:		symbol = "L3-M5";	return 100.0;							// L3-M5    La1
				case N4:	symbol = "L3-N4,5";											// L3-N4,5  Lb2,15
//							Note that the L3_N45 fit is split into two Z ranges
					if ( z < DATA_L3_N45_1[3] ) {
						return linePolyCalc ( z,  DATA_L3_N45_1 );
					} else {
						return linePolyCalc ( z,  DATA_L3_N45_2 );
					};
				case M4:		symbol = "L3-M4";	return linePolyCalc ( z,  DATA_L3_M4 );	// L3-M4   La2
				case O4:	symbol = "L3-O4,5";	return linePolyCalc ( z,  DATA_L3_O45 );	// L3-O4,5 Lb5
				case N1:		symbol = "L3-N1";	return linePolyCalc ( z,  DATA_L3_N1 );	// L3-N1   Lb6
				case M1:		symbol = "L3-M1";											// L3-M1   Ll
//							data is L3-M1 to L3-M4,5 ratio times 100 (Ll/La)*100, in two Z ranges
					float temp;
					if ( z < DATA_Ll_La_1[3] ) {
						temp = linePolyCalc ( z,  DATA_Ll_La_1 );
					} else {
						temp = linePolyCalc ( z,  DATA_Ll_La_2 );
					};
					return temp * ( 100.0 + linePolyCalc ( z,  DATA_L3_M4 ) ) / 100.0;
				default: return 0.0;
			};
//				M lines
		case M2:
			if ( lower.index() == N4) { symbol = "M2-N4"; return .001f; };						// M2-N4
			return 0.0;
		case M3:
			if ( lower.index() == N5) { symbol = "M3-N5"; return .01f; };						// M3-N5  Mg
			return 0.0;
		case M4	:
			if ( lower.index() == N6 ) { symbol = "M4-N6"; return 0.34f;	 };				// M4-N6  Mb
			if ( lower.index() == N2 ) { symbol = "M4,5-N2,3"; return 0.001f; };			// M4,5-N2,3  Mzeta
			return 0.0;
		case M5:
			if ( lower.index() == N6 ) { symbol = "M5-N6,7"; return 0.65f; };			// M5-N6,7 Ma
			return 0.0;
		default: return 0.0;
	};
	return 0.0;

};

float XrayLines::linePolyCalc ( const float z, const float array[] ) {
//	Helper routine to evaluate data in arrays, see intSymbol routine for array definitions

	if ( z < array[2] ) {
//			use linear extrapolation below range of fit
		return ( array[0] + z * array[1] );
	} else {
//			evaluate 4th order polynomial fit
		float sum = 0.0;
		float power = 1.0;
		int i;
		for ( i=0; i<5; i++ ) {
			sum = sum + array[i+4] * power;
			power *= z;
		};
		return sum;
	};
};


void XrayLines::merge_peaks( const float separation, vector <XrayLinesInfo> &lineList_in ) {
    if( lineList_in.size() <= 0 ) return;
    vector <XrayLinesInfo> newLineList;
	vector <bool> lineIncluded( lineList_in.size(), false );
    //		merge lines whose energies are within separation of each other
	bool done = false;
	while( ! done ) {
        //			find strongest line not yet included in merged set
        float largest_int = 0;
        int largest_i = 0;
        int i;
        for( i=0; i<lineList_in.size(); i++ ) {
            if( lineIncluded[i] ) continue;
            if( lineList_in[i].lineIntensity > largest_int ) {
                largest_int = lineList_in[i].lineIntensity;
                largest_i = i;
            };
        };
		float average = lineList_in[largest_i].lineEnergy * lineList_in[largest_i].lineIntensity;
		float sum = lineList_in[largest_i].lineIntensity;
		lineIncluded[largest_i] = true;
		done = true;
        //			look through lines and find any within separation of the strongest line
		for( i=0; i<lineList_in.size(); i++ ) {
			if( ( ! lineIncluded[i] ) && fabs( lineList_in[i].lineEnergy - lineList_in[largest_i].lineEnergy ) < separation ) {
				lineIncluded[i] = true;
				average += lineList_in[i].lineEnergy * lineList_in[i].lineIntensity;
				sum += lineList_in[i].lineIntensity;
			} else {
				if( ! lineIncluded[i] ) done = false;   //  keep looking until all lines in the group are included
			};
		};
		average /= sum;
        XrayLinesInfo temp;
        temp.edgeOccupied = lineList_in[largest_i].edgeOccupied;
        temp.lineIUPAC = lineList_in[largest_i].lineIUPAC;
        temp.lineEnergy = average;
        temp.lineIntensity = sum;
        temp.lineFactor = 1.0;
        newLineList.push_back( temp );
	};

    //		move merged list to input argument
	lineList_in.clear();
	lineList_in.resize( newLineList.size() );
    int i;
	for( i=0; i<newLineList.size(); i++ ) lineList_in[i] = newLineList[i];

	return;
};


string XrayLinesInfo_toString(const XrayLinesInfo &lines)
{
    ostringstream os;
    os << "XrayLinesInfo:" << endl;
    os << "  edgeOccupied:" << lines.edgeOccupied.toString() << endl;

    os << "  lineIUPAC: " << lines.lineIUPAC << endl;
    os << "  lineEnergy: " << lines.lineEnergy << endl;
    os << "  lineFactor: " << lines.lineFactor << endl;

    return os.str();
}

std::string XrayLines::toString() const
{
    ostringstream os;
    os << "XrayLines:" << endl;
    os << "  edgeVacant:" << edgeVacant.toString() << endl;

    os << "  lineCount: " << lineCount << endl;
    os << "  commonFactor_value: " << commonFactor_value << endl;

    os << "  lineList: " << endl;
    int c = 0;
    for(auto it = lineList.begin(); it != lineList.end(); it++)
    {
        os << "lineList[" << c << "]: " << XrayLinesInfo_toString(*it) << endl;
        c++;
    }
    return os.str();
}