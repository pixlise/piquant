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
#include "Element.h"

using namespace std;

//	Elam Ravel Sieber database class

//  Modified Jan. 16, 2017
//      Add check_Z and chack_symbol static member functions

const int Element::MAXZ = 103;

const string Element::SYMBOLS[MAXZ+1] = { " ",
		"H", "He", "Li", "Be", "B",
		"C", "N", "O", "F", "Ne",
		"Na", "Mg", "Al", "Si", "P",
		"S", "Cl", "Ar", "K", "Ca",
		"Sc", "Ti", "V", "Cr", "Mn",
		"Fe", "Co", "Ni", "Cu", "Zn",
		"Ga", "Ge", "As", "Se", "Br",
		"Kr", "Rb", "Sr", "Y", "Zr",
		"Nb", "Mo", "Tc", "Ru", "Rh",
		"Pd", "Ag", "Cd", "In", "Sn",
		"Sb", "Te", "I", "Xe", "Cs",
		"Ba", "La", "Ce", "Pr", "Nd",
		"Pm", "Sm", "Eu", "Gd", "Tb",
		"Dy", "Ho", "Er", "Tm", "Yb",
		"Lu", "Hf", "Ta", "W", "Re",
		"Os", "Ir", "Pt", "Au", "Hg",
		"Tl", "Pb", "Bi", "Po", "At",
		"Rn", "Fr", "Ra", "Ac", "Th",
		"Pa", "U", "Np", "Pu", "Am",
		"Cm", "Bk", "Cf", "Es", "Fm",
		"Md", "No", "Lr"
	};


Element::Element(const int z) {
	if (z>0 && z<=MAXZ) {
		atomicNumberZ = z;
	} else {
		throw string("Atomic number <1 or >103")  ;
	}
};

Element::Element(const string s) {
	int z;
	bool found = false;
	for (z=1; z <= MAXZ; z++ ) {
		if ( s == SYMBOLS[z] ) {
			atomicNumberZ = z;
			found = true;
			break;
		};
	};
	if (! found ) {
		throw string( "Atomic symbol invalid" );
	};
};

Element::Element() { atomicNumberZ = 1; };

Element::Element(const Element& e) { atomicNumberZ = e.atomicNumberZ; };

const string& Element::symbol() const {
	return SYMBOLS[atomicNumberZ];
};

const bool Element::check_symbol ( const string &symbol_in ) {
	int z;
	for (z=1; z <= MAXZ; z++ ) {
		if ( symbol_in == SYMBOLS[z] ) return true;
	};
	return false;
};




float Element::atomicWeight() const {
	static const float ATOMICWEIGHT[MAXZ+1] = { 0.0f,
                 1.0079f,   4.0026f,   6.9410f,   9.0122f,  10.8100f,  12.0110f,  14.0067f,  15.9994f,
                18.9984f,  20.1790f,  22.9898f,  24.3050f,  26.9815f,  28.0855f,  30.9738f,  32.0600f,
                35.4530f,  39.9480f,  39.0983f,  40.0800f,  44.9559f,  47.8800f,  50.9415f,  51.9960f,
                54.9380f,  55.8470f,  58.9332f,  58.6900f,  63.5460f,  65.3800f,  69.7200f,  72.5900f,
                74.9216f,  78.9600f,  79.9040f,  83.8000f,  85.4678f,  87.6200f,  88.9059f,  91.2200f,
                92.9064f,  95.9400f,  97.9070f, 101.0700f, 102.9055f, 106.4200f, 107.8680f, 112.4100f,
               114.8200f, 118.6900f, 121.7500f, 127.6000f, 126.9045f, 131.2900f, 132.9054f, 137.3300f,
               138.9055f, 140.1200f, 140.9077f, 144.2400f, 144.9130f, 150.3600f, 151.9600f, 157.2500f,
               158.9254f, 162.5000f, 164.9304f, 167.2600f, 168.9342f, 173.0400f, 174.9670f, 178.4900f,
               180.9479f, 183.8500f, 186.2070f, 190.2000f, 192.2200f, 195.0800f, 196.9665f, 200.5900f,
               204.3830f, 207.2000f, 208.9804f, 208.9820f, 209.9870f, 222.0180f, 223.0200f, 226.0254f,
               227.0278f, 232.0381f, 231.0359f, 238.0510f, 237.0482f, 239.0520f, 243.0610f, 247.0700f,
               247.0700f, 251.0800f, 252.0830f, 257.0950f,    .0000f,    .0000f,    .0000f,
		};
	return ATOMICWEIGHT[atomicNumberZ] ;
};

float Element::density() const {
	static const float DENSITY[MAXZ+1] = { 0.0f,
                   .071f,     .122f,     .533f,    1.845f,    2.340f,    2.260f,     .810f,    1.140f,
                  1.108f,    1.207f,     .969f,    1.735f,    2.694f,    2.320f,    1.820f,    2.070f,
                  1.560f,    1.400f,     .860f,    1.550f,    2.980f,    4.530f,    6.100f,    7.180f,
                  7.430f,    7.860f,    8.900f,    8.876f,    8.940f,    7.112f,    5.877f,    5.307f,
                  5.720f,    4.780f,    3.110f,    2.600f,    1.529f,    2.540f,    4.456f,    6.494f,
                  8.550f,   10.200f,   11.480f,   12.390f,   12.390f,   12.000f,   10.480f,    8.630f,
                  7.300f,    7.300f,    6.679f,    6.230f,    4.920f,    3.520f,    1.870f,    3.500f,
                  6.127f,    6.637f,    6.761f,    6.994f,    7.200f,    7.510f,    5.228f,    7.877f,
                  8.214f,    8.525f,    8.769f,    9.039f,    9.294f,    6.953f,    9.811f,   13.790f,
                 16.624f,   19.300f,   20.980f,   22.530f,   22.390f,   21.410f,   18.850f,   13.522f,
                 11.830f,   11.330f,    9.730f,    9.300f,     .000f,    4.400f,     .000f,    5.000f,
                 10.050f,   11.700f,   15.340f,   18.920f,   20.210f,   19.800f,   13.640f,   13.490f,
                 14.000f,     .000f,     .000f,     .000f,     .000f,     .000f,     .000f,
		};
	return DENSITY[atomicNumberZ] ;
};

string Element::toString() const
{
    return symbol();
}
