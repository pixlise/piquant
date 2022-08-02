// Copyright (c) 2018-2022 California Institute of Technology (‚ÄúCaltech‚Äù) and
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
#include <vector>
#include <sstream>
#include "interp.h"
#include "Element.h"
#include "XrayEdge.h"

//  Modified June 27, 2012 to update energies from latest Gwyn Williams table and add some comments
//  Checked occupancies vs energy table
//  Modified March 2012 to add level widths from Campbell
//  Modified Oct. 1, 2013 to change K fluor yields to Bambynek values
//      (see note with data below on Hubble erratum)
//` Modified March 31, 2017 to make EdgeIndex and EdgeLevel more consistent
//      Change first EdgeLevel to K and first EdgeIndex to K1

const std::string XrayEdge::EDGE_NAMES[MAXINDEX+1] = {
	"K",  "L1", "L2", "L3", "M1", "M2",
	"M3", "M4", "M5", "N1", "N2", "N3",
	"N4", "N5", "N6", "N7", "O1", "O2",
	"O3", "O4", "O5", "O6", "O7",
	"P1", "P2", "P3", "P4", "Q1"
};

const std::string XrayEdge::EDGE_DESIGNATIONS[MAXINDEX+1] = {
    "1s",    "2s",    "2p1/2", "2p3/2", "3s",    "3p1/2",
    "3p3/2", "3d3/2", "3d5/2", "4s",    "4p1/2", "4p3/2",
    "4d3/2", "4d5/2", "4f5/2", "4f7/2", "5s",    "5p1/2",
    "5p3/2", "5d3/2", "5d5/2", "5f5/2", "5f7/2",
    "6s", "6p1/2", "6p3/2", "6d3/2", "7s"
};


XrayEdge::XrayEdge( const Element e, const EdgeIndex indexIn) {
	if (indexIn>=K1 && indexIn<=Q1) {
		z = e;
		edgeIndex = indexIn;
	} else {
		throw string ("XrayEdge: invalid edge index");
	}
};
/*
XrayEdge::XrayEdge( const Element e, const string s) {
	int i;
	bool found = false;
	for (i=0; i<=MAXINDEX; i++ ) {
		if ( s == EDGE_NAMES[i] ) {
			z=e;
			edgeIndex = i;
			found = true;
			break;
		};
	};
	if ( ! found ) { throw string("XrayEdge: invalid edge name" ); };
};

XrayEdge::XrayEdge( const Element e, const float energy ) {
	int i;
	bool found = false;
	float ee;
	for (i=0; i<=MAXINDEX; i++ ) {
		ee = EDGE_ENERGIES[(e.Z()-1)*(MAXINDEX+1)+i];
		if ( ee == 0.0 ) break;
		if ( energy >= ee ) {
			z=e;
			edgeIndex = i;
			found = true;
			break;
		};
	};
	if ( ! found ) { throw string("XrayEdge: no edges below given energy" ); };
};
*/

XrayEdge::XrayEdge( const XrayEdge& ee ) {
	z = ee.z;
	edgeIndex = ee.index();
};

bool XrayEdge::operator==(const XrayEdge& ee) const {
	if ( z == ee.z && edgeIndex == ee.index() ) {
		return true;
	} else {
		return false;
	}
};

bool XrayEdge::operator<(const XrayEdge& comp) const {
	float eeThis = (*this).energy();
	float eeComp = comp.energy();
//		order by edge energy, then by atomic number if edge energies are identical
	if ( eeThis == eeComp ) return (*this).element().Z() < comp.element().Z();
 	return eeThis < eeComp ;
 };

bool XrayEdge::operator>(const XrayEdge& comp) const {
	float eeThis = (*this).energy();
	float eeComp = comp.energy();
//		order by edge energy, then by atomic number if edge energies are identical
	if ( eeThis == eeComp ) return (*this).element().Z() > comp.element().Z();
 	return eeThis > eeComp ;
 };

const float XrayEdge::degeneracy() const {
	static const float EDGE_DEGENERACIES[MAXINDEX+1] = {
 		2.00, 2.00, 2.00, 4.00, 2.00, 2.00,
		4.00, 4.00, 6.00, 2.00, 2.00, 4.00,
		4.00, 6.00, 6.00, 8.00, 2.00, 2.00,
		4.00, 4.00, 6.00, 6.00, 8.00,
		2.00, 2.00, 4.00, 4.00, 2.00
	};

	return EDGE_DEGENERACIES[edgeIndex];
};

const int XrayEdge::numberOccupied ( vector<EdgeIndex> &edgeList, const Element& el ) {
//	return list of occupied edges for this element
	edgeList.resize(0);
	int indexTemp;
	for (indexTemp=0; indexTemp <= MAXINDEX; indexTemp++ ) {
		int occ = EDGE_OCCUPANCIES[(el.Z()-1)*(MAXINDEX+1)+indexTemp];
		if ( occ > 0 ) {
			edgeList.insert( edgeList.end(), static_cast<EdgeIndex>(indexTemp) );
		};
	};
	return edgeList.size();
};

const int XrayEdge::numberOfEdges ( vector <EdgeIndex> &edgeList, const Element& el, const float excitEnergy ) {
//	return list of occupied edges for this element with energies greater than
//		zero but less than (or equal to) the argument (excitation energy)
	edgeList.resize(0);
	int indexTemp;
	for (indexTemp=0; indexTemp <= MAXINDEX; indexTemp++ ) {
		float ee;
		ee = EDGE_ENERGIES[(el.Z()-1)*(MAXINDEX+1)+indexTemp];
		if ( ee > 0.0 && ( excitEnergy == 0.0 || ee <= excitEnergy ) ) {
			edgeList.insert( edgeList.end(), static_cast<EdgeIndex>(indexTemp) );
		};
	};
	return edgeList.size();
};

int XrayEdge::number () const {
//	return the number of edges for this element with indices greater than
//		this one and energy greater than zero
	int thisMany = 0;
	int indexTemp;
	for (indexTemp=edgeIndex+1; indexTemp <= MAXINDEX; indexTemp++ ) {
		float ee;
		ee = EDGE_ENERGIES[(z.Z()-1)*(MAXINDEX+1)+indexTemp];
		if ( ee > 0.0 ) {
			thisMany += 1;
		};
	};
	return thisMany;
};

const EdgeLevel XrayEdge::level() const {
	static const EdgeLevel levels[MAXINDEX+1] = {
		K, L, L, L, M, M, M, M, M, N, N, N, N, N, N, N, O, O, O, O, O, O, O, P, P, P, P, Q
	};
	return levels[edgeIndex];
};

const EdgeAngularMonmentum XrayEdge::angularMomentum() const {
	static const EdgeAngularMonmentum am[MAXINDEX+1] = {
		s, s, p, p, s, p, p, d, d, s, p, p, d, d, f, f, s, p, p, d, d, f, f, s, p, p, d, s
	};
	return am[edgeIndex];
};

const float XrayEdge::spin() const {
	static const float spins[MAXINDEX+1] = {
		0.0, 0.0, 0.5, 1.5, 0.0, 0.5, 1.5, 1.5, 2.5, 0.0, 0.5, 1.5,
		1.5, 2.5, 2.5, 3.5, 0.0, 0.5, 1.5, 1.5, 2.5, 2.5, 3.5, 0.0, 0.5, 1.5, 1.5, 0.0
	};
	return spins[edgeIndex];
};

const float XrayEdge::cktotal ( const XrayEdge toEdge ) const {
//		total Coster-Kronig transition rate from this edge to given edge
//		including vacancy rattling down through intermediate levels
//			executive decision - only handle up to 10 intermediate levels
	const int N_MAX = 10;
	int zi = z.Z();
//		make sure elements match
	if ( zi != toEdge.element().Z() ) return 0.0;
	int toIndex = toEdge.index();
//		make sure the transition is allowed
	if ( ( toIndex <= edgeIndex ) || ( (*this).level() != toEdge.level() ) ) return 0.0;
//		we need to sum the probabilities of all possible routes for the vacancy to
//			travel from the initial level (this edge) to the final level (index toIndex)
//		the outer loop will run over the number of intermediate levels, from zero to maximum
	int n = toIndex - edgeIndex -1;		//		maximum number of intermediate states
	if ( n > N_MAX ) n = N_MAX;			//		don't do more than size of list (>10 factors=> very small)
	int inter;							// 		actual number of intermediate states
	float sum = 0.0;		//	used to sum probabilities of each route through imtermediates
	float prod;				//	used to accumulate product of individual transition probalities
	int step;				//	used to count transitons in probability product over imtermediate transitions
	int firstIndex;			//	initial state index for each step
	int secondIndex;		//	final state index for each step
	int ii[N_MAX];			//	used to hold indices of intermediate states to generate routes
//		handle the direct transition (zero imtermediate states)
//	cout << endl << "Direct transition indices: " << edgeIndex << "  " << toIndex << endl;
	sum += ck_value ( zi, edgeIndex, toIndex );
//		now loop through the complicated ones
	for ( inter=1; inter<n; inter++ ) {
//				set up list of intermediate states
		for ( step=0; step<inter; step++ ) {
			ii[step] = edgeIndex + step + 1;
		};
//			start of loop to increment indices of intermediate states to cover all
//				possible routes for the vacancy
		bool moreRoutes = true;
		while ( moreRoutes ) {
//			cout << endl << "# intermediates: " << inter << " product indices:" ;
//				accumulate product of transition probabilities through intermediate states
			prod = 1.0;
			for ( step=0; step<=inter; step++ ) {
//					vacancy always starts first at this edge
				if ( step == 0 ) firstIndex = edgeIndex;
				if ( step == inter ) {
					secondIndex = toIndex;	// if last step, end at toEdge
				} else {
					secondIndex = ii[step];	// if not at end, go to next intermediate level
				};
				prod *= ck_value ( zi, firstIndex, secondIndex );
//				cout << "  " << firstIndex << " " << secondIndex << " " << prod;
//					retain index of intermediate level to start next step (vacancy now here)
				firstIndex = secondIndex;
			};
//			cout << endl;
			sum += prod;
//				now increment last index of intermediate states to get next route for vacancy
//					and propagate increment up through list as necessary
			int nextIncr = inter;
			for ( step=0; step<inter; step++ ) {
				ii[inter-1-step] += 1;
				if ( ii[inter-1-step] < toIndex-step ) break;
				nextIncr -= 1;
			};
//				if the first intermediate index in the list wasn't incremented successfully,
//					then we're done
			if ( nextIncr < 1 )  {
				moreRoutes = false;
			} else {
//					re-arrange indices after the last one incremented back to starting order
				for ( step=nextIncr; step<inter; step++ ) {
					ii[step] = ii[step-1] + 1;
				};
			};
		};
	};
//		finally, handle the last one - all intermediate states involved (if any)
	if ( n > 0 ) {
		prod = 1.0;
		for ( step=0; step<=n; step++ ) {
//				must always start first atthis edge and end at toEdge
			if ( step == 0 ) firstIndex = edgeIndex;
			secondIndex = firstIndex + 1;
			prod *= ck_value ( zi, firstIndex, secondIndex );
//			cout << "Final (all states) route, product indices: " << firstIndex << "  " << secondIndex << " " << prod << endl;
//				retain index of intermediate level for next transition
			firstIndex = secondIndex;
		};
		sum += prod;
	};
	return sum;
};

//	large data statements for these functions are at the end of this file
const float XrayEdge::energy() const {
	int i;
	i=z.Z();
	return EDGE_ENERGIES[(i-1)*(MAXINDEX+1)+edgeIndex];
};

const float XrayEdge::occupancy() const {
	int i;
	i=z.Z();
	return EDGE_OCCUPANCIES[(i-1)*(MAXINDEX+1)+edgeIndex];
};

const float XrayEdge::jzero() const {
	int i;
	i=z.Z();
	return EDGE_JZEROS[(i-1)*(MAXINDEX+1)+edgeIndex];
};

//		Level width data from Campbell ATNDT_77_1_2001
const float XrayEdge::width ( ) const {
#include "CampbellWidthData_extrap_XrayEdge_Mar2012.h"
   	int i;
	i=z.Z();
    if( i <= WIDTH_DATA_MAX_Z && edgeIndex <= N7 ) {
        return WIDTH_DATA[ (i-1) * WIDTH_DATA_HORIZ + edgeIndex ];
    } else {
        return 0.00001f;
    }
};

string XrayEdge::toString() const
{
	ostringstream os;
	os << "XrayEdge:" << endl;
    os << "  element: " << element().toString() << endl;
    os << "  edgeIndex: " << name() << endl;
    return os.str();
}

const float XrayEdge::yield_value ( const int zi, const int fromIndex ) const {

//		return fluorescence yield

//	Alternate K-shell fluorescence yields from new fits in J. H. Hubbell et. al.,
//		J. Chem. Phys. Ref. Data, Vol. 23, No. 2, 1994, pp339-364.
//	Fluorescence yields and Coster-Kronig transition rates for K and L shells
//	Krause, J. Phys. Chem. Ref. Data, Vol. 8, No. 2, 1979, pp307-327.
//	values for wK, wL2,and f23 are from Table 1. (values for light atoms in condensed matter)
//	(note that this produces a large step in f23 values at z=30, see discussion in reference
//	section 5.3 L2 Subshell and section 7 last paragraph)

//	Values of wL1 for Z=85-110 and f12 for Z=72-96 from Krause were modified as suggested by
//	W. Jitschin, "Progress in Measurements of L-Subshell Fluorescence, Coster-Kronig,
//	and Auger Values", AIP Conference Proceedings 215, X-ray and Inner-Shell Processes,
//	Knocxville, TN, 1990. T. A. Carlson, M. O. Krause, and S. T. Manson, Eds.
//	(American Institute of Physics, 1990).

float DATA_wK[111] = { 0.0,
0, 0, 9e-05, 3.3e-05, 0.0007, 0.0014, 0.0031, 0.0058, 0.0092, 0.016,
0.021, 0.028, 0.036, 0.048, 0.061, 0.078, 0.097, 0.118, 0.14, 0.163,
0.188, 0.214, 0.243, 0.273, 0.308, 0.34, 0.373, 0.406, 0.44, 0.474,
0.507, 0.535, 0.562, 0.589, 0.618, 0.643, 0.667, 0.69, 0.71, 0.73,
0.747, 0.765, 0.78, 0.794, 0.808, 0.82, 0.831, 0.843, 0.853, 0.862,
0.87, 0.877, 0.884, 0.891, 0.897, 0.902, 0.907, 0.912, 0.917, 0.921,
0.925, 0.929, 0.932, 0.935, 0.938, 0.941, 0.944, 0.947, 0.949, 0.951,
0.953, 0.955, 0.957, 0.958, 0.959, 0.961, 0.962, 0.963, 0.964, 0.965,
0.966, 0.967, 0.968, 0.968, 0.969, 0.969, 0.97, 0.97, 0.971, 0.971,
0.972, 0.972, 0.973, 0.973, 0.974, 0.974, 0.975, 0.975, 0.975, 0.976,
0.976, 0.976, 0.977, 0.977, 0.977, 0.978, 0.978, 0.978, 0.978, 0.979 };


float DATA_wL1[111] = { 0.0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 2e-05, 2.6e-05, 3e-05, 3.9e-05, 7.4e-05, 0.00012, 0.00018, 0.00024, 0.00031,
0.00039, 0.00047, 0.00058, 0.00071, 0.00084, 0.001, 0.0012, 0.0014, 0.0016, 0.0018,
0.0021, 0.0024, 0.0028, 0.0032, 0.0036, 0.0041, 0.0046, 0.0051, 0.0059, 0.0068,
0.0094, 0.01, 0.011, 0.012, 0.013, 0.014, 0.016, 0.018, 0.02, 0.037,
0.039, 0.041, 0.044, 0.046, 0.049, 0.052, 0.055, 0.058, 0.061, 0.064,
0.066, 0.071, 0.075, 0.079, 0.083, 0.089, 0.094, 0.1, 0.106, 0.112,
0.12, 0.128, 0.137, 0.147, 0.144, 0.13, 0.12, 0.114, 0.107, 0.107,
0.107, 0.112, 0.117, 0.122, 0.128, 0.130, 0.130, 0.130, 0.130, 0.130,
0.130, 0.141, 0.150, 0.164, 0.174, 0.182, 0.189, 0.195, 0.202, 0.210,
0.218, 0.224, 0.226, 0.233, 0.240, 0.248, 0.256, 0.265, 0.274, 0.283 };
//	0.107, 0.112, 0.117, 0.122, 0.128, 0.134, 0.139, 0.146, 0.153, 0.161,  	before Jitschin modifications
//	0.162, 0.176, 0.187, 0.205, 0.218, 0.228, 0.236, 0.244, 0.253, 0.263,
//	0.272, 0.28, 0.282, 0.291, 0.3, 0.31, 0.32, 0.331, 0.343, 0.354 };


float DATA_wL2[111] = { 0.0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0.0012, 0.00075, 0.00037, 3.1e-05, 0.00026, 0.00024, 0.00022, 0.00027, 0.00033,
0.00084, 0.0011, 0.0018, 0.0023, 0.0031, 0.0036, 0.0044, 0.0051, 0.0057, 0.0095,
0.012, 0.013, 0.014, 0.016, 0.018, 0.02, 0.022, 0.024, 0.026, 0.028,
0.031, 0.034, 0.037, 0.04, 0.043, 0.047, 0.051, 0.056, 0.061, 0.065,
0.069, 0.074, 0.079, 0.083, 0.09, 0.096, 0.103, 0.11, 0.117, 0.124,
0.132, 0.14, 0.149, 0.158, 0.167, 0.178, 0.189, 0.2, 0.211, 0.222,
0.234, 0.246, 0.258, 0.27, 0.283, 0.295, 0.308, 0.321, 0.334, 0.347,
0.36, 0.373, 0.387, 0.401, 0.415, 0.429, 0.443, 0.456, 0.468, 0.479,
0.472, 0.467, 0.466, 0.464, 0.471, 0.479, 0.485, 0.49, 0.497, 0.506,
0.515, 0.524, 0.533, 0.544, 0.553, 0.562, 0.573, 0.584, 0.59, 0.598 };


float DATA_wL3[111] = { 0.0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 1.2e-13, 0.00075, 0.00038, 3.1e-05, 0.00026, 0.00024, 0.00022, 0.00027, 0.00033,
0.00084, 0.0015, 0.0026, 0.0037, 0.005, 0.0063, 0.0077, 0.0093, 0.011, 0.012,
0.013, 0.015, 0.016, 0.018, 0.02, 0.022, 0.024, 0.026, 0.028, 0.031,
0.034, 0.037, 0.04, 0.043, 0.046, 0.049, 0.052, 0.056, 0.06, 0.064,
0.069, 0.074, 0.079, 0.085, 0.091, 0.097, 0.104, 0.111, 0.118, 0.125,
0.132, 0.139, 0.147, 0.155, 0.164, 0.174, 0.182, 0.192, 0.201, 0.218,
0.22, 0.231, 0.243, 0.255, 0.268, 0.281, 0.294, 0.306, 0.32, 0.333,
0.347, 0.36, 0.373, 0.386, 0.399, 0.411, 0.424, 0.437, 0.45, 0.463,
0.476, 0.489, 0.502, 0.514, 0.526, 0.539, 0.55, 0.56, 0.57, 0.579,
0.588, 0.596, 0.604, 0.611, 0.618, 0.624, 0.63, 0.635, 0.64, 0.644 };

//	Fluorescence yields and Coster-Kronig transition rates for M shells
//	Eugene J. McGuire, "Atomic M-Shell Coster-Kronig, Auger, and Radiative Rates, and Fluorescence
//	Yields for Ca-Th", Physical Review A, Vol. 5, No. 3, March 1972, pp1043-1047.
float DATA_ZM123[29] = { 0.0,
20, 22, 23, 24, 25, 26, 27, 28, 29, 30,
32, 36, 40, 44, 47, 50, 54, 57, 60, 63,
67, 70, 73, 76, 79, 83, 86, 90 };

float DATA_wM1[29] = { 0.0,
8.4e-06, 3.2e-06, 2.9e-06, 2.6e-06, 3.1e-06, 2.8e-06, 2.8e-06, 3.5e-06, 4.1e-06, 4.6e-06,
9.1e-06, 4.9e-05, 7e-05, 0.00012, 0.00017, 0.00025, 0.00047, 0.00084, 0.00081, 0.00087,
0.00108, 0.00115, 0.00145, 0.00165, 0.00213, 0.00289, 0.00395, 0.00453 };

float DATA_wM2[29] = { 0.0,
0, 3.4e-05, 2.3e-05, 1.6e-05, 1.6e-05, 1.6e-05, 1.7e-05, 1.5e-05, 1.6e-05, 2.2e-05,
2.6e-05, 6e-05, 1.4e-05, 0.00026, 0.00039, 0.0007, 0.0009, 0.0011, 0.00132, 0.00147,
0.00185, 0.00197, 0.00264, 0.00325, 0.0423, 0.00652, 0.00975, 0.014 };

float DATA_wM3[29] = { 0.0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 6e-05, 0.00015, 0.00023, 0.00032, 0.00054, 0.00068, 0.00099, 0.00105, 0.00126,
0.00145, 0.00166, 0.00214, 0.0032, 0.0042, 0.00533, 0.0063, 0.0081 };

float DATA_ZM45[19] = { 0.0,
32, 36, 40, 44, 47, 50, 54, 57, 60, 63,
67, 70, 73, 76, 79, 83, 86, 90 };

float DATA_wM4[19] = { 0.0,
0.0027, 0.0027, 0.0027, 0.0029, 0.0027, 0.0027, 0.0027, 0.0027, 0.0026, 0.0041,
0.0067, 0.0086, 0.013, 0.0137, 0.0264, 0.033, 0.0355, 0.0588 };

float DATA_wM5[19] = { 0.0,
0, 0, 0, 0, 0, 0, 0, 0, 0.0032, 0.0059,
0.0106, 0.0149, 0.0205, 0.0232, 0.0256, 0.0325, 0.0362, 0.0497 };


//	Fluorescence yields and Coster-Kronig transition rates for N shells
//	Eugene J. McGuire, "Atomic N-shell Coster-Kronig, Auger, and Radiative Rates and Fluorescence
//	Yields for 38 <= Z <= 103", Physical Review A 9, No. 5, May 1974, pp1840-1851.
//		Values for Z=38 to 50 were adjusted according to instructions on
//		page 1845, at the end of Section IV.a., and the last sentence of the conclusions.

float DATA_ZN123[26] = { 0.0,
38, 40, 42, 44, 47, 50, 54, 57, 58, 60,
63, 65, 67, 70, 73, 74, 77, 79, 83, 86,
90, 92, 96, 100, 103 };

float DATA_wN1[26] = { 0.0,
1.2e-05, 2.02e-06, 7e-07, 3.8e-07, 1.3e-07, 2.5e-08, 1.5e-05, 2.3e-05, 2.2e-05, 2.2e-05,
2.8e-05, 3.3e-05, 3.7e-05, 4.9e-05, 6.1e-05, 7.5e-05, 1e-04, 0.00013, 0.00024, 0.00036,
0.00061, 0.00071, 0.00083, 0.0012, 0.0016 };

float DATA_wN2[26] = { 0.0,
0.013, 2.3e-05, 7.2e-06, 6.2e-06, 8.4e-06, 7e-06, 5.5e-05, 0.00011, 8.1e-05, 7.3e-05,
7.2e-05, 7e-05, 7e-05, 6.5e-05, 7.5e-05, 8e-05, 1e-04, 0.00012, 0.00019, 0.0002,
0.00067, 0.0009, 0.00105, 0.0016, 0.0022 };

float DATA_wN3[26] = { 0.0,
0.013, 2.3e-05, 7.2e-06, 6.2e-06, 8.4e-06, 7e-06, 5.5e-05, 8.5e-05, 5.9e-05, 5.6e-05,
5.6e-05, 5.6e-05, 5.9e-05, 5.2e-05, 5.9e-05, 6.9e-05, 7.7e-05, 9e-05, 0.00014, 0.00012,
0.00045, 0.00059, 0.00069, 0.00103, 0.00132 };

float DATA_ZN4567[21] = { 0.0,
50, 54, 57, 58, 60, 63, 65, 67, 70, 73,
74, 77, 79, 83, 86, 90, 92, 96, 100, 103 };

float DATA_wN4[21] = { 0.0,
4.2e-06, 6.3e-05, 0.00014, 0.00019, 0.00013, 0.00011, 0.00011, 0.00011, 7e-05, 0.00012,
0.00012, 0.00013, 0.00015, 0.00023, 0.00028, 0.00047, 0.0007, 0.00079, 0.0012, 0.0013 };

float DATA_wN5[21] = { 0.0,
4.2e-06, 6.3e-05, 0.00014, 0.00019, 0.00013, 0.00011, 0.00011, 0.00011, 9e-05, 0.00012,
0.00013, 0.00014, 0.00016, 0.00023, 0.00028, 0.00047, 0.00072, 0.00081, 0.0012, 0.0014 };

float DATA_wN67[21] = { 0.0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 3.6e-06,
7.2e-06, 2.1e-05, 3.9e-05, 0.00018, 0.00044, 0.00046, 0.00031, 0.00025, 0.00029, 0.00029 };



// K shell fluorescence yields from Bambynek/Cohen/Burhop/Hubbel
//  As given in "A Review, Bibliography, and Tabulation of K, L, and Higher Atomic Shell X-Ray fluorescence Yields"
//      J.H. Hubbell, P.N. Trehan, Nirmal Singh, and B. Chand, D Mehta, M.L. Garg, R.R. Garg, Surinder Singh, and S. Puri,
//      J. Phys. Chem. Ref. Data, Vol. 23, No. 2, 1994, p. 339-364, Table 8 col. 1 and 2.
//  See especially eratum by Hubble et al. J. Chem. Phys. Ref. Data, Vol. 33, No. 2, 2004, p 621.
float DATA_wK_Bambynek[101] = { 0,
    0, 0, 0.0002928, 0.0006929, 0.001409, 0.002575, 0.004349, 0.006909, 0.01045, 0.01519,
    0.02133, 0.02911, 0.03872, 0.05037, 0.06422, 0.08038, 0.09892, 0.1199, 0.1432, 0.1687,
    0.1962, 0.2256, 0.2564, 0.2885, 0.3213, 0.3546, 0.3880, 0.4212, 0.4538, 0.4857,
    0.5166, 0.5464, 0.5748, 0.6019, 0.6275, 0.6517, 0.6744, 0.6956, 0.7155, 0.7340,
    0.7512, 0.7672, 0.7821, 0.7958, 0.8086, 0.8204, 0.8313, 0.8415, 0.8508, 0.8595,
    0.8676, 0.8750, 0.8819, 0.8883, 0.8942, 0.8997, 0.9049, 0.9096, 0.9140, 0.9181,
    0.9220, 0.9255, 0.9289, 0.9320, 0.9349, 0.9376, 0.9401, 0.9425, 0.9447, 0.9467,
    0.9487, 0.9505, 0.9522, 0.9538, 0.9553, 0.9567, 0.9580, 0.9592, 0.9604, 0.9615,
    0.9625, 0.9634, 0.9643, 0.9652, 0.9659, 0.9667, 0.9674, 0.9680, 0.9686, 0.9691,
    0.9696, 0.9701, 0.9706, 0.9710, 0.9713, 0.9717, 0.9720, 0.9722, 0.9725, 0.9727 };

	float zf = float ( zi ) ;
	switch ( fromIndex ) {
//		case K1:	return DATA_wK[zi];		//		K - Krause
//		case K1: 							//		K - Hubbell et. al.
//			if ( zi < 11 ) return DATA_wK[z.Z()];
//			if ( zi>=11 && zi<=19 ) return 1.434e-1 - 2.5606e-2*zf + 1.3163e-3*zf*zf;
//			if ( zi > 19 )       return -7.6388e-1 + 5.4070e-2*zf - 4.0544e-4*zf*zf
//				- 1.4348e-6*zf*zf*zf + 1.8252e-8*zf*zf*zf*zf;
		case K1: 							//		K - Bambynek
			if ( zi <= 100 ) {
                return DATA_wK_Bambynek[zi];
            } else if ( zi <= 110 ) {
                return DATA_wK[zi];
            } else {
                return 0;
            }
            break;
		case L1:	return DATA_wL1[zi];		//		L1
		case L2:	return DATA_wL2[zi];		//		L2
		case L3:	return DATA_wL3[zi];		//		L3

		case M1:	return interp ( zf, DATA_ZM123, DATA_wM1, 29 ) ;	//		M1
		case M2:	return interp ( zf, DATA_ZM123, DATA_wM2, 29 ) ;	//		M2
		case M3:	return interp ( zf, DATA_ZM123, DATA_wM3, 29 ) ;	//		M3
		case M4:	return interp ( zf, DATA_ZM45, DATA_wM4, 19 ) ;	//		M4
		case M5:	return interp ( zf, DATA_ZM45, DATA_wM5, 19 ) ;	//		M5

		case N1:		return interp ( zf, DATA_ZN123, DATA_wN1, 26 ) ;	//	N1
		case N2:	return interp ( zf, DATA_ZN123, DATA_wN2, 26 ) ;	//	N2
		case N3:	return interp ( zf, DATA_ZN123, DATA_wN3, 26 ) ;	//	N3
		case N4:	return interp ( zf, DATA_ZN4567, DATA_wN4, 21 ) ;	//	N4
		case N5:	return interp ( zf, DATA_ZN4567, DATA_wN5, 21 ) ;	//	N5
		case N6:	return interp ( zf, DATA_ZN4567, DATA_wN67, 21 ) ;//	N6
		case N7:	return interp ( zf, DATA_ZN4567, DATA_wN67, 21 ) ;//	N7

		default:	return 0.0;	break;
	};
};


const float XrayEdge::jump_value ( const int zi, const int fromIndex ) const {
//		return absorption edge jump ratio (also called edge step)

//	all values are calculated by taking x-ray absorption cross section values
//		1 eV above edge divided by value 1 eV below edge
//  edge jump data from cross section database of 4/4/2000 2:24pm
//  Z for first jump value preceeds data declaration, max Z is 98
//  values are averaged where values had little or no Z dependence
//  fits to 1/Z are used where R value of fit is greater than 0.86

//  starts at Z=4.0000
const float J_K_DATA[95] = {
23.530, 21.700, 19.020, 17.440, 15.400, 13.840, 13.610, 11.840, 12.020, 10.950,
10.370, 9.9690, 9.6130, 9.3090, 9.0540, 9.1630, 8.7440, 8.5510, 8.3660, 8.2420,
8.0450, 7.9990, 7.8930, 7.7960, 7.7070, 7.5600, 7.5430, 7.4680, 7.3920, 7.3140,
7.2250, 7.1410, 7.0580, 6.9700, 6.8880, 6.8140, 6.7490, 6.6830, 6.5390, 6.5610,
6.5030, 6.4440, 6.3950, 6.3340, 6.2750, 6.2290, 6.1610, 6.1300, 6.0850, 6.0390,
5.9880, 5.9490, 5.9010, 5.8630, 5.8060, 5.7630, 5.7180, 5.6740, 5.6250, 5.5820,
5.5370, 5.4920, 5.4500, 5.4030, 5.3530, 5.3690, 5.2700, 5.2280, 5.1850, 5.1430,
5.0990, 5.0560, 5.0070, 4.9630, 4.9730, 4.8740, 4.8280, 4.7810, 4.7310, 4.6820,
4.6280, 4.5860, 4.4940, 4.4810, 4.4330, 4.3810, 4.3350, 4.2870, 4.2380, 4.1910,
4.1440, 4.0950, 4.0430, 3.9980, 3.9540 };

// starts at Z=15.000
const float J_L3_DATA[84] = {
6.1900, 5.5460, 4.6460, 4.1830, 10.980, 5.8110, 5.0070, 4.5920, 4.0730, 5.0000,
3.3060, 3.0990, 2.7730, 2.6140, 3.1350, 2.4550, 2.9590, 3.6840, 4.1560, 4.4310,
4.6540, 4.6780, 3.9870, 3.9820, 3.9090, 3.8360, 3.7500, 3.6740, 3.5940, 3.5180,
3.4440, 3.3280, 3.3080, 3.2540, 3.1500, 3.0060, 3.0130, 2.9780, 2.9500, 2.9200,
2.9550, 2.8960, 2.8680, 2.8520, 2.8150, 2.8180, 2.8010, 2.7830, 2.7660, 2.7470,
2.7320, 2.7070, 2.7060, 2.6920, 2.6780, 2.6650, 2.6490, 2.6300, 2.6130, 2.6130,
2.6000, 2.5870, 2.5800, 2.5630, 2.5500, 2.5340, 2.5200, 2.5060, 2.4920, 2.4780,
2.4590, 2.4490, 2.4430, 2.4330, 2.4230, 2.4110, 2.4010, 2.3920, 2.3800, 2.3700,
2.3570, 2.3480, 2.3460, 2.3050 };

// starts at Z=37.000,
const float J_M4_DATA[62] = {
1.9530, 1.1390, 1.1200, 1.1590, 1.1380, 1.0950, 1.1360, 1.0790, 1.1260, 1.0580,
1.1000, 1.0930, 1.1350, 1.0760, 1.0760, 1.1520, 1.0800, 1.1440, 1.0520, 1.1070,
1.0850, 1.0660, 1.1050, 1.2580, 1.2800, 1.2480, 1.2480, 1.2280, 1.2070, 1.1970,
1.1810, 1.1500, 1.1690, 1.1790, 1.1640, 1.1440, 1.1270, 1.1060, 1.0450, 1.1000,
1.0660, 1.0750, 1.0700, 1.0880, 1.1620, 1.2650, 1.3920, 1.5260, 1.6520, 1.6210,
1.5550, 1.4840, 1.4490, 1.4290, 1.4310, 1.4240, 1.4150, 1.3560, 1.4100, 1.3950,
1.3900, 1.3890 };

// starts at Z=37.000,
const float J_M5_DATA[62] = {
1.8430, 1.8080, 1.6310, 1.5650, 1.4210, 1.3020, 1.3210, 1.1930, 1.2160, 1.0840,
1.1420, 1.1340, 1.0790, 1.0750, 1.0610, 1.1300, 1.0630, 1.1490, 1.4510, 1.5690,
1.2710, 1.3340, 1.3110, 1.3500, 1.3480, 1.4390, 1.3570, 1.1420, 1.4700, 1.4220,
1.4420, 1.4750, 1.3950, 1.5220, 1.2620, 1.2550, 1.2330, 1.2030, 1.3020, 1.0500,
1.0870, 1.0940, 1.0920, 1.1380, 1.3670, 1.7570, 2.3250, 2.9690, 3.3890, 3.5210,
3.1350, 2.8030, 2.6020, 2.4810, 2.4800, 2.4330, 2.3870, 2.3950, 2.3380, 2.2170,
2.2350, 2.2400 };

// starts at Z=52.000, ends at Z=76.000
const float J_N3_DATA[25] = {
1.4110, 1.4980, 1.7260, 1.8840, 1.6170, 1.5120, 1.3480, 1.2750, 1.2240, 1.1990,
1.1610, 1.1590, 1.1370, 1.1320, 1.1220, 1.1170, 1.1130, 1.1100, 1.1010, 1.0890,
1.0810, 1.0770, 1.0720, 1.0690, 1.0670 };

// starts at Z=58.000, ends at Z=71.000
const float J_N4_DATA[14] = {
1.5560, 1.3470, 1.2570, 1.1920, 1.1790, 1.0320, 1.1350, 1.1050, 1.0950, 1.0910,
1.0880, 1.0840, 1.0280, 1.0210 };

// starts at Z=58.000
const float J_N5_DATA[41] = {
1.5560, 1.3470, 1.2570, 1.1920, 1.1790, 1.1230, 1.1350, 1.1050, 1.0950, 1.0910,
1.0880, 1.0840, 1.0500, 1.0440, 1.0370, 1.0370, 1.0370, 1.0320, 1.0310, 1.0320,
1.0220, 1.0220, 1.0200, 1.0170, 1.0070, 1.0020, 1.0050, 1.0000, 1.0000, 1.1240,
1.0250, 1.0140, 1.0130, 1.0110, 1.0130, 1.0070, 1.0010, 1.0240, 1.0200, 1.0180,
1.0200 };

// starts at Z=84.000,
const float J_N6_DATA[15] = {
1.2830, 1.3500, 1.3390, 1.2870, 1.1650, 1.1550, 1.0040, 1.0090, 1.0000, 1.0230,
1.0020, 1.0840, 1.0590, 1.0700, 1.0690 };

// starts at Z=81.000,
const float J_N7_DATA[18] = {
1.0470, 1.1680, 1.0000, 1.2830, 1.3500, 1.3390, 1.2870, 1.1650, 1.1550, 1.1130,
1.1100, 1.0960, 1.1110, 1.0910, 1.1330, 1.0960, 1.0900, 1.1150 };

// starts at Z=78.000,
const float J_O1_DATA[21] = {
1.0850, 1.0630, 1.1520, 1.1720, 1.1720, 1.5870, 1.2970, 1.2870, 1.2970, 1.3000,
1.3370, 1.3200, 1.3100, 1.2870, 1.2820, 1.2730, 1.2590, 1.0660, 1.0640, 1.0610,
1.0560 };

	float zf = float ( zi ) ;
	switch ( fromIndex ) {
		case K1: 	if ( zi >= 4 ) return J_K_DATA[zi-4];					//	K - use tabular values
			break;
		case L1:	if ( zi >= 13 ) return 1.1688 - 1.0373/zf;				//	L1	R=0.88
			break;
		case L2:	if ( zi >= 16 ) return 1.40;							//	L2  +/- 0.05
			break;
		case L3:	if ( zi >= 15 && zi <= 98 ) return J_L3_DATA[zi-15];	//	L3 - use tabular values
			break;
		case M1:	if ( zi >= 27 ) return 1.04 ;							//	M1	+/- 0.01
			break;
		case M2:	if ( zi >= 31 ) return 1.058 ;							//	M2	=/- 0.003
			break;
		case M3:	if ( zi >= 31 ) return 1.1882 - 2.4050/zf ;				//	M3	R=0.86 for Z>55
			break;
		case M4:	if ( zi >= 37 && zi <= 98 ) return J_M4_DATA[zi-37];	//	M4 - use tabular values
			break;
		case M5:	if ( zi >= 37 && zi <= 98 ) return J_M5_DATA[zi-37];	//	M5 - use tabular values
			break;
		case N1:	{														//	N1	R=0.91 for Z=48 to 85
						if ( zi >= 48 && zi <= 85 ) return 0.7313 + 24.69/zf;
						if ( zi > 85 ) return 1.017;
						break;
					};
		case N2:	if ( zi >= 52 ) return 1.005 ;							//	N2	except for Z=51, 52, & 61 ??
			break;
		case N3:	{														//	N3 - table Z=52 to 76
						if ( zi >= 52 && zi <= 76 ) return J_N3_DATA[zi-52];
						if ( zi > 76 ) return 1.06 ;
						break;
					};
		case N4:	{														//	N4 - table Z=58 to 71
						if ( zi >= 58 && zi <= 71 ) return J_N4_DATA[zi-58];
						if ( zi > 71 ) return 1.02 ;
						break;
					};
		case N5:	if ( zi >= 58 && zi <= 98 ) return J_N5_DATA[zi-58];	//	N5 - use tabular values
			break;
		case N6:	if ( zi >= 84 && zi <= 98 ) return J_N6_DATA[zi-84];	//	N6 - use tabular values
			break;
		case N7:	if ( zi >= 81 && zi <= 98 ) return J_N7_DATA[zi-81];	//	N7 - use tabular values
			break;
		case O1:	if ( zi >= 78 && zi <= 98 ) return J_O1_DATA[zi-78];	//	O1 - use tabular values
			break;
		case O2:	if ( zi >= 82 ) return 1.16 ;							//	O2	unreliable
			break;
		case O3:	if ( zi >= 84 ) return 1.5  ;							//	O3	unreliable
			break;
		case O4:	if ( zi >= 92 ) return 1.05 ;							//	O4	unreliable
			break;
		case O5:	if ( zi >= 93 ) return 1.10 ;							//	O5	unreliable
			break;

		default:	return 1.00;	break;
	};
	return 1.00;
};


const float XrayEdge::ck_value ( const int zi, const int fromIndex, const int toIndex ) const {

//		Coster-Kronig transition rate from this edge to given edge

//	Fluorescence yields and Coster-Kronig transition rates for K and L shells
//	Krause, J. Phys. Chem. Ref. Data, Vol. 8, No. 2, 1979, pp307-327.
//	values for wK, wL2,and f23 are from Table 1. (values for light atoms in condensed matter)
//	(note that this produces a large step in f23 values at z=30, see discussion in reference
//	section 5.3 L2 Subshell and section 7 last paragraph)

//	Values of wL1 for Z=85-110 and f12 for Z=72-96 from Krause were modified as suggested by
//	W. Jitschin, "Progress in Measurements of L-Subshell Fluorescence, Coster-Kronig,
//	and Auger Values", AIP Conference Proceedings 215, X-ray and Inner-Shell Processes,
//	Knocxville, TN, 1990. T. A. Carlson, M. O. Krause, and S. T. Manson, Eds.
//	(American Institute of Physics, 1990).

float DATA_f12[111] = { 0.0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0.32, 0.32, 0.32, 0.32, 0.32, 0.32, 0.31, 0.31, 0.31,
0.31, 0.31, 0.31, 0.31, 0.3, 0.3, 0.3, 0.3, 0.3, 0.29,
0.29, 0.28, 0.28, 0.28, 0.28, 0.27, 0.27, 0.27, 0.26, 0.26,
0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.17,
0.17, 0.18, 0.18, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19,
0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19,
0.19, 0.18, 0.16, 0.15, 0.13, 0.12, 0.11, 0.09, 0.08, 0.06,
0.05, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04,
0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.03, 0.03, 0.03,
0.02, 0.02, 0.01, 0.01, 0.01, 0, 0, 0, 0, 0 };
//	0.19, 0.18, 0.18, 0.17, 0.16, 0.16, 0.15, 0.14, 0.14, 0.13, 	before Jitschin modifications
//	0.13, 0.12, 0.11, 0.11, 0.1, 0.1, 0.1, 0.09, 0.09, 0.09,
//	0.08, 0.08, 0.07, 0.05, 0.05, 0.04, 0.04, 0.03, 0.03, 0.03,
//	0.02, 0.02, 0.01, 0.01, 0.01, 0, 0, 0, 0, 0 };


float DATA_f13[111] = { 0.0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0.64, 0.64, 0.64, 0.63, 0.62, 0.62, 0.62, 0.62, 0.61,
0.6, 0.59, 0.58, 0.57, 0.58, 0.57, 0.56, 0.55, 0.54, 0.54,
0.53, 0.53, 0.53, 0.52, 0.52, 0.52, 0.52, 0.52, 0.52, 0.52,
0.61, 0.61, 0.61, 0.61, 0.6, 0.6, 0.59, 0.59, 0.59, 0.27,
0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.29, 0.29, 0.29, 0.3,
0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.29, 0.29,
0.28, 0.28, 0.28, 0.28, 0.33, 0.39, 0.45, 0.5, 0.53, 0.56,
0.57, 0.58, 0.58, 0.58, 0.59, 0.58, 0.58, 0.58, 0.58, 0.57,
0.58, 0.57, 0.57, 0.56, 0.55, 0.55, 0.54, 0.54, 0.54, 0.53,
0.53, 0.52, 0.53, 0.52, 0.51, 0.51, 0.5, 0.5, 0.49, 0.48 };

float DATA_f23[111] = { 0.0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0.25, 0.29, 0.35, 0.39, 0.42, 0.43, 0.45, 0.47, 0.24,
0.032, 0.05, 0.063, 0.076, 0.088, 0.1, 0.109, 0.117, 0.126, 0.132,
0.137, 0.141, 0.144, 0.148, 0.15, 0.151, 0.153, 0.155, 0.157, 0.157,
0.156, 0.155, 0.154, 0.154, 0.154, 0.153, 0.153, 0.153, 0.153, 0.152,
0.151, 0.15, 0.149, 0.147, 0.145, 0.143, 0.142, 0.14, 0.139, 0.138,
0.136, 0.135, 0.134, 0.133, 0.13, 0.128, 0.126, 0.124, 0.122, 0.12,
0.118, 0.116, 0.113, 0.111, 0.111, 0.11, 0.109, 0.108, 0.108, 0.108,
0.139, 0.167, 0.192, 0.198, 0.203, 0.2, 0.198, 0.197, 0.196, 0.194,
0.191, 0.189, 0.185, 0.181, 0.178, 0.174, 0.171, 0.165, 0.163, 0.158 };

//	Fluorescence yields and Coster-Kronig transition rates for M shells
//	Eugene J. McGuire, "Atomic M-Shell Coster-Kronig, Auger, and Radiative Rates, and Fluorescence
//	Yields for Ca-Th", Physical Review A, Vol. 5, No. 3, March 1972, pp1043-1047.
float DATA_ZM123[29] = { 0.0,
20, 22, 23, 24, 25, 26, 27, 28, 29, 30,
32, 36, 40, 44, 47, 50, 54, 57, 60, 63,
67, 70, 73, 76, 79, 83, 86, 90 };

float DATA_SM12[29] = { 0.0,
0.328, 0.319, 0.315, 0.319, 0.312, 0.311, 0.308, 0.307, 0.304, 0.283,
0.249, 0.27, 0.278, 0.305, 0.343, 0.315, 0.238, 0.195, 0.236, 0.338,
0.266, 0.272, 0.197, 0.161, 0.148, 0.109, 0.143, 0.072 };

float DATA_SM13[29] = { 0.0,
0.655, 0.639, 0.631, 0.638, 0.623, 0.621, 0.616, 0.614, 0.608, 0.566,
0.52, 0.54, 0.475, 0.457, 0.461, 0.475, 0.505, 0.506, 0.489, 0.485,
0.527, 0.525, 0.561, 0.594, 0.594, 0.65, 0.593, 0.69 };

float DATA_SM14[29] = { 0.0,
0, 0.314, 0.335, 0.397, 0.357, 0.371, 0.376, 0.381, 0.406, 0.374,
0.273, 0.086, 0.108, 0.065, 0.065, 0.067, 0.081, 0.094, 0.092, 0.07,
0.061, 0.056, 0.065, 0.067, 0.067, 0.065, 0.069, 0.063 };

float DATA_SM15[29] = { 0.0,
0, 0.471, 0.503, 0.596, 0.538, 0.556, 0.564, 0.566, 0.61, 0.561,
0.409, 0.127, 0.163, 0.124, 0.097, 0.101, 0.122, 0.14, 0.128, 0.1,
0.09, 0.091, 0.115, 0.109, 0.112, 0.095, 0.1, 0.091 };

float DATA_SM23[29] = { 0.0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0.032, 0.067, 0.073, 0.016, 0.031, 0.034, 0.057, 0.062,
0.106, 0.116, 0.114, 0.107, 0.114, 0.103, 0.128, 0.116 };

float DATA_SM24[29] = { 0.0,
0, 1.057, 1.089, 1.123, 1.108, 1.116, 1.12, 1.122, 1.133, 1.107,
1.085, 0.919, 0.519, 0.55, 0.57, 0.604, 0.612, 0.557, 0.64, 0.514,
0.667, 0.68, 0.674, 0.684, 0.673, 0.662, 0.61, 0.623 };

float DATA_SM25[29] = { 0.0,
0, 0.672, 0.82, 0.834, 0.797, 0.815, 0.817, 0.827, 0.85, 0.811,
0.786, 0.516, 0.309, 0.283, 0.258, 0.252, 0.233, 0.282, 0.172, 0.137,
0.12, 0.105, 0.106, 0.098, 0.095, 0.083, 0.093, 0.088 };

float DATA_SM34[29] = { 0.0,
0, 0.509, 0.558, 0.612, 0.589, 0.6, 0.602, 0.609, 0.623, 0.597,
0.58, 0.395, 0.252, 0.236, 0.223, 0.213, 0.206, 0.198, 0.174, 0.165,
0.145, 0.141, 0.082, 0.106, 0.114, 0.094, 0.072, 0.097 };

float DATA_SM35[29] = { 0.0,
0, 1.22, 1.28, 1.342, 1.317, 1.329, 1.355, 1.341, 1.36, 1.32,
1.292, 1.039, 0.677, 0.672, 0.689, 0.678, 0.688, 0.678, 0.712, 0.72,
0.751, 0.761, 0.81, 0.764, 0.782, 0.75, 0.768, 0.725 };

float DATA_ZM45[19] = { 0.0,
32, 36, 40, 44, 47, 50, 54, 57, 60, 63,
67, 70, 73, 76, 79, 83, 86, 90 };

float DATA_fM45[19] = { 0.0,
0, 0, 0, 0, 0, 0, 0, 0, 0.267, 0.369,
0.408, 0.479, 0.411, 0.418, 0.046, 0.035, 0.065, 0.066 };

//	Fluorescence yields and Coster-Kronig transition rates for N shells
//	Eugene J. McGuire, "Atomic N-shell Coster-Kronig, Auger, and Radiative Rates and Fluorescence
//	Yields for 38 <= Z <= 103", Physical Review A 9, No. 5, May 1974, pp1840-1851.
//		Values for Z=38 to 50 were adjusted according to instructions on
//		page 1845, at the end of Section IV.a., and the last sentence of the conclusions.

float DATA_ZN123[26] = { 0.0,
38, 40, 42, 44, 47, 50, 54, 57, 58, 60,
63, 65, 67, 70, 73, 74, 77, 79, 83, 86,
90, 92, 96, 100, 103 };

float DATA_SN12[26] = { 0.0,
0.33, 0.3, 0.18, 0.11, 0.05, 0.01, 0.2, 0.22, 0.23, 0.24,
0.25, 0.25, 0.25, 0.25, 0.2, 0.19, 0.19, 0.27, 0.26, 0.22,
0.31, 0.36, 0.42, 0.53, 0.46 };

float DATA_SN13[26] = { 0.0,
0.66, 0.61, 0.36, 0.23, 0.11, 0.03, 0.4, 0.43, 0.47, 0.48,
0.5, 0.49, 0.5, 0.5, 0.56, 0.56, 0.57, 0.5, 0.46, 0.47,
0.33, 0.33, 0.25, 0.18, 0.2 };

float DATA_SN14[26] = { 0.0,
0, 0.35, 0.22, 0.14, 0.07, 0.02, 0.23, 0.2, 0.15, 0.14,
0.11, 0.11, 0.1, 0.08, 0.073, 0.062, 0.059, 0.049, 0.055, 0.053,
0.059, 0.049, 0.051, 0.048, 0.048 };

float DATA_SN15[26] = { 0.0,
0, 0.51, 0.32, 0.21, 0.1, 0.02, 0.35, 0.3, 0.23, 0.21,
0.16, 0.16, 0.15, 0.12, 0.11, 0.092, 0.091, 0.074, 0.083, 0.08,
0.088, 0.074, 0.077, 0.072, 0.073 };

float DATA_SN167[26] = { 0.0,
0, 0, 0, 0, 0, 0, 0, 0, 0.17, 0.38,
0.48, 0.59, 0.64, 0.65, 0.63, 0.59, 0.56, 0.45, 0.38, 0.41,
0.27, 0.25, 0.23, 0.19, 0.24 };

float DATA_SN23[26] = { 0.0,
0, 0, 0, 0, 0, 0, 0, 0.018, 0.061, 0.065,
0.134, 0.127, 0.142, 0.122, 0.118, 0.114, 0.119, 0.09, 0.144, 0.047,
0.092, 0.111, 0.113, 0.138, 0.174 };

float DATA_SN24[26] = { 0.0,
0, 1.08, 1.12, 1.12, 1.12, 1.11, 0.76, 0.73, 0.72, 0.72,
0.68, 0.68, 0.67, 0.69, 0.68, 0.67, 0.68, 0.67, 0.58, 0.72,
0.53, 0.49, 0.5, 0.39, 0.4 };

float DATA_SN25[26] = { 0.0,
0, 0.69, 0.85, 0.86, 0.87, 0.85, 0.21, 0.21, 0.16, 0.15,
0.11, 0.11, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.11,
0.11, 0.11, 0.1, 0.1, 0.1 };

float DATA_SN267[26] = { 0.0,
0, 0, 0, 0, 0, 0, 0, 0, 0.44, 0.59,
0.7, 0.75, 0.81, 0.84, 0.84, 0.79, 0.77, 0.68, 0.53, 0.63,
0.3, 0.26, 0.22, 0.25, 0.25 };

float DATA_SN34[26] = { 0.0,
0, 0.52, 0.62, 0.63, 0.63, 0.62, 0.2, 0.21, 0.18, 0.17,
0.16, 0.16, 0.15, 0.14, 0.15, 0.14, 0.15, 0.14, 0.14, 0.14,
0.13, 0.13, 0.13, 0.12, 0.12 };

float DATA_SN35[26] = { 0.0,
0, 1.25, 1.35, 1.36, 1.36, 1.34, 0.77, 0.745, 0.76, 0.75,
0.76, 0.75, 0.76, 0.75, 0.74, 0.73, 0.74, 0.71, 0.65, 0.72,
0.57, 0.55, 0.55, 0.49, 0.48 };

float DATA_SN367[26] = { 0.0,
0, 0, 0, 0, 0, 0, 0, 0, 0.41, 0.57,
0.71, 0.77, 0.8, 0.86, 0.87, 0.82, 0.8, 0.74, 0.62, 0.66,
0.34, 0.29, 0.25, 0.31, 0.32 };

float DATA_ZN4567[21] = { 0.0,
50, 54, 57, 58, 60, 63, 65, 67, 70, 73,
74, 77, 79, 83, 86, 90, 92, 96, 100, 103 };

float DATA_SN45[21] = { 0.0,
0, 0, 0, 0, 0, 0, 0, 0, 0.45, 0.015,
0.021, 0.034, 0.049, 0.003, 0.004, 0.002, 0.016, 0.021, 0.032, 0.03 };

float DATA_SN467[21] = { 0.0,
0, 0, 0, 1.02, 1.4, 1.62, 1.72, 1.75, 1.41, 1.74,
1.7, 1.61, 1.51, 1.63, 1.07, 1.02, 0.871, 0.849, 0.813, 0.832 };

float DATA_SN567[21] = { 0.0,
0, 0, 0, 1.02, 1.4, 1.67, 1.72, 1.75, 1.82, 1.77,
1.73, 1.8, 1.59, 1.63, 1.07, 1.02, 0.886, 0.868, 0.841, 0.857 };

//		select the transition rate based on the indices of the two levels
	float zf = float ( zi );
	switch ( fromIndex ) {
		case K1:	return 0.0;				// no transitions for K edge
		case L1:							//	L1 - transitions allowed to L2 and L3
			switch ( toIndex ) {
				case L2:		return DATA_f12[zi];		// L1-L2
				case L3:		return DATA_f13[zi];		// L1-L3
				default:	return 0.0;
			};
		case L2:							//	L2 - transition allowed to L3
			switch ( toIndex ) {
				case L3:		return DATA_f23[zi];		// L2-L3
				default:	return 0.0;
			};
		case L3:							// no transitions for L3 edge
			return 0.0;
		case M1:							//	M1 - transitions allowed to all other M edges
			switch ( toIndex ) {
				case M2:		return interp ( zf, DATA_ZM123, DATA_SM12, 29 );	// M1-M2
				case M3:		return interp ( zf, DATA_ZM123, DATA_SM13, 29 );	// M1-M3
				case M4:		return interp ( zf, DATA_ZM123, DATA_SM14, 29 );	// M1-M4
				case M5:		return interp ( zf, DATA_ZM123, DATA_SM15, 29 );	// M1-M5
				default:	return 0.0;
			};
		case M2:							//	M2 - transitions allowed to M3, M4, and M5 edges
			switch ( toIndex ) {
				case M3:		return interp ( zf, DATA_ZM123, DATA_SM23, 29 );	// M2-M3
				case M4:		return interp ( zf, DATA_ZM123, DATA_SM24, 29 );	// M2-M4
				case M5:		return interp ( zf, DATA_ZM123, DATA_SM25, 29 );	// M2-M5
				default:	return 0.0;
			};
		case M3:							//	M3 - transitions allowed to M4 and M5 edges
			switch ( toIndex ) {
				case M4:		return interp ( zf, DATA_ZM123, DATA_SM34, 29 );	// M3-M4
				case M5:		return interp ( zf, DATA_ZM123, DATA_SM35, 29 );	// M3-M5
				default:	return 0.0;
			};
		case M4:							//	M4 - transitions allowed to M5 edge
			switch ( toIndex ) {
				case M5:		return interp ( zf, DATA_ZM45, DATA_fM45, 19 );		// M4-M5
				default:	return 0.0;
			};
		case M5:		return 0.0;			//	M5 - no transitions allowed

		case N1:							//	N1 - transitions allowed to all other N edges
			switch ( toIndex ) {
				case N2:		return interp ( zf, DATA_ZN123, DATA_SN12, 26 );	// N1-N2
				case N3:		return interp ( zf, DATA_ZN123, DATA_SN13, 26 );	// N1-N3
				case N4:		return interp ( zf, DATA_ZN123, DATA_SN14, 26 );	// N1-N4
				case N5:		return interp ( zf, DATA_ZN123, DATA_SN15, 26 );	// N1-N5
				case N6:	return 0.5 * interp ( zf, DATA_ZN123, DATA_SN167, 26 );	// N1-N6
				case N7:	return 0.5 * interp ( zf, DATA_ZN123, DATA_SN167, 26 );	// N1-N7
				default:	return 0.0;
			};
		case N2:							//	N2 - transitions allowed to N3 thru N7 edges
			switch ( toIndex ) {
				case N3:		return interp ( zf, DATA_ZN123, DATA_SN23, 26 );	// N2-N3
				case N4:		return interp ( zf, DATA_ZN123, DATA_SN24, 26 );	// N2-N4
				case N5:		return interp ( zf, DATA_ZN123, DATA_SN25, 26 );	// N2-N5
				case N6:	return 0.5 * interp ( zf, DATA_ZN123, DATA_SN267, 26 );	// N2-N6
				case N7:	return 0.5 * interp ( zf, DATA_ZN123, DATA_SN267, 26 );	// N2-N7
				default:	return 0.0;
			};
		case N3:							//	N3 - transitions allowed to N4 thru N7 edges
			switch ( toIndex ) {
				case N4:		return interp ( zf, DATA_ZN123, DATA_SN34, 26 );	// N3-N4
				case N5:		return interp ( zf, DATA_ZN123, DATA_SN35, 26 );	// N3-N5
				case N6:	return 0.5 * interp ( zf, DATA_ZN123, DATA_SN367, 26 );	// N3-N6
				case N7:	return 0.5 * interp ( zf, DATA_ZN123, DATA_SN367, 26 );	// N3-N7
				default:	return 0.0;
			};
		case N4:							//	N4 - transitions allowed to N5, N6, and N7 edges
			switch ( toIndex ) {
				case N5:		return interp ( zf, DATA_ZN4567, DATA_SN45, 21 );	// N4-N5
				case N6:	return 0.5 * interp ( zf, DATA_ZN4567, DATA_SN467, 21 );	// N4-N6
				case N7:	return 0.5 * interp ( zf, DATA_ZN4567, DATA_SN467, 21 );	// N4-N7
				default:	return 0.0;
			};
		case N5:							//	N5 - transitions allowed to N6 and N7 edges
			switch ( toIndex ) {
				case N6:	return 0.5 * interp ( zf, DATA_ZN4567, DATA_SN567, 21 );	// N5-N6
				case N7:	return 0.5 * interp ( zf, DATA_ZN4567, DATA_SN567, 21 );	// N5-N7
				default:	return 0.0;
			};
		default: return 0.0;
	};
	return 0.0;
};



//https://userweb.jlab.org/~gwyn/ebindene.html
//downloaded June 28, 2012
//Values compiled by Gwyn Williams
//Last updated June 21, 2000
//[Additional values for Z=93-103 from Browne & Firestone "Table of Radioactive Isotopes" John Wiley 1986.]
//[Some low levels added to make emission lines agree with Bearden and Burr, Rev. Mod. Phys. 39 (1967) 78-124,
// as quoted in CRC Handbook 51st Ed. (1970-71) E126-E164.]

//    K       L1       L2       L3       M1       M2       M3       M4       M5       N1       N2       N3       N4       N5
//    N6       N7       O1       O2       O3       O4       O5       O6       O7       P1       P2       P3       P4       Q1
const float XrayEdge::EDGE_ENERGIES[MAXZ*(MAXINDEX+1)] = {
//    K,      L-I,    L-II,   L-III,  M-I,    M-II,   M-III,  M-IV,   M-V,    N-I,    N-II,   N-III,  N-IV,   N-V,    N-VI,   N-VII,  O-I,    O-II,   O-III,  O-IV,   O-V,    O6,     O7,     P-I,    P-II,   P-III,  P4,     Q1,     ,       ,       Element
//    1s,     2s,     2p1/2,  2p3/2,  3s,     3p1/2,  3p3/2,  3d3/2,  3d5/2,  4s,     4p1/2,  4p3/2,  4d3/2,  4d5/2,  4f5/2,  4f7/2,  5s,     5p1/2,  5p3/2,  5d3/2,  5d5/2,  ,       ,       6s,     6p1/2,  6p3/2,  ,       ,       ,       ,
    13.6,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     1,      H
    24.6,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     2,      He
    54.7,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     3,      Li
    111.5,  0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     4,      Be
    188,    0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     5,      B
    284.2,  0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     6,      C
    409.9,  37.3,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     7,      N
    543.1,  41.6,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     8,      O
    696.7,  0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     9,      F
    870.2,  48.5,   21.7,   21.6,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     10,     Ne
    1070.8, 63.5,   30.65,  30.81,  0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     11,     Na
    1303,   88.7,   49.78,  49.5,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     12,     Mg
    1559.6, 117.8,  72.95,  72.55,  0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     13,     Al
    1839,   149.7,  99.82,  99.42,  0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     14,     Si
    2145.5, 189,    136,    135,    0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     15,     P
    2472,   230.9,  163.6,  162.5,  0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     16,     S
    2822.4, 270,    202,    200,    0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     17,     Cl
    3205.9, 326.3,  250.6,  248.4,  29.3,   15.9,   15.7,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     18,     Ar
    3608.4, 378.6,  297.3,  294.6,  34.8,   18.3,   18.3,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     19,     K
    4038.5, 438.4,  349.7,  346.2,  44.3,   25.4,   25.4,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     20,     Ca
    4492,   498,    403.6,  398.7,  51.1,   28.3,   28.3,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     21,     Sc
    4966,   560.9,  460.2,  453.8,  58.7,   32.6,   32.6,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     22,     Ti
    5465,   626.7,  519.8,  512.1,  66.3,   37.2,   37.2,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     23,     V
    5989,   696,    583.8,  574.1,  74.1,   42.2,   42.2,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     24,     Cr
    6539,   769.1,  649.9,  638.7,  82.3,   47.2,   47.2,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     25,     Mn
    7112,   844.6,  719.9,  706.8,  91.3,   52.7,   52.7,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     26,     Fe
    7709,   925.1,  793.2,  778.1,  101,    58.9,   59.9,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     27,     Co
    8333,   1008.6, 870,    852.7,  110.8,  68,     66.2,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     28,     Ni
    8979,   1096.7, 952.3,  932.7,  122.5,  77.3,   75.1,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     29,     Cu
    9659,   1196.2, 1044.9, 1021.8, 139.8,  91.4,   88.6,   10.2,   10.1,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     30,     Zn
    10367,  1299,   1143.2, 1116.4, 159.5,  103.5,  100,    18.7,   18.7,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     31,     Ga
    11103,  1414.6, 1248.1, 1217,   180.1,  124.9,  120.8,  29.8,   29.2,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     32,     Ge
    11867,  1527,   1359.1, 1323.6, 204.7,  146.2,  141.2,  41.7,   41.7,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     33,     As
    12658,  1652,   1474.3, 1433.9, 229.6,  166.5,  160.7,  55.5,   54.6,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     34,     Se
    13474,  1782,   1596,   1550,   257,    189,    182,    70,     69,     0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     35,     Br
    14326,  1921,   1730.9, 1678.4, 292.8,  222.2,  214.4,  95,     93.8,   27.5,   14.1,   14.1,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     36,     Kr
    15200,  2065,   1864,   1804,   326.7,  248.7,  239.1,  113,    112,    30.5,   16.3,   15.3,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     37,     Rb
    16105,  2216,   2007,   1940,   358.7,  280.3,  270,    136,    134.2,  38.9,   21.3,   20.1,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     38,     Sr
    17038,  2373,   2156,   2080,   392,    310.6,  298.8,  157.7,  155.8,  43.8,   24.4,   23.1,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     39,     Y
    17998,  2532,   2307,   2223,   430.3,  343.5,  329.8,  181.1,  178.8,  50.6,   28.5,   27.1,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     40,     Zr
    18986,  2698,   2465,   2371,   466.6,  376.1,  360.6,  205,    202.3,  56.4,   32.6,   30.8,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     41,     Nb
    20000,  2866,   2625,   2520,   506.3,  411.6,  394,    231.1,  227.9,  63.2,   37.6,   35.5,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     42,     Mo
    21044,  3043,   2793,   2677,   544,    447.6,  417.7,  257.6,  253.9,  69.5,   42.3,   39.9,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     43,     Tc
    22117,  3224,   2967,   2838,   586.1,  483.5,  461.4,  284.2,  280,    75,     46.3,   43.2,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     44,     Ru
    23220,  3412,   3146,   3004,   628.1,  521.3,  496.5,  311.9,  307.2,  81.4,   50.5,   47.3,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     45,     Rh
    24350,  3604,   3330,   3173,   671.6,  559.9,  532.3,  340.5,  335.2,  87.1,   55.7,   50.9,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     46,     Pd
    25514,  3806,   3524,   3351,   719,    603.8,  573,    374,    368.3,  97,     63.7,   58.3,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     47,     Ag
    26711,  4018,   3727,   3538,   772,    652.6,  618.4,  411.9,  405.2,  109.8,  63.9,   63.9,   11.7,   10.7,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     48,     Cd
    27940,  4238,   3938,   3730,   827.2,  703.2,  665.3,  451.4,  443.9,  122.9,  73.5,   73.5,   17.7,   16.9,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     49,     In
    29200,  4465,   4156,   3929,   884.7,  756.5,  714.6,  493.2,  484.9,  137.1,  83.6,   83.6,   24.9,   23.9,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     50,     Sn
    30491,  4698,   4380,   4132,   946,    812.7,  766.4,  537.5,  528.2,  153.2,  95.6,   95.6,   33.3,   32.1,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     51,     Sb
    31814,  4939,   4612,   4341,   1006,   870.8,  820,    583.4,  573,    169.4,  103.3,  103.3,  41.9,   40.4,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     52,     Te
    33169,  5188,   4852,   4557,   1072,   931,    875,    630.8,  619.3,  186,    123,    123,    50.6,   48.9,   0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     53,     I
    34561,  5453,   5107,   4786,   1148.7, 1002.1, 940.6,  689,    676.4,  213.2,  146.7,  145.5,  69.5,   67.5,   0,      0,      23.3,   13.4,   12.1,   0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     54,     Xe
    35985,  5714,   5359,   5012,   1211,   1071,   1003,   740.5,  726.6,  232.3,  172.4,  161.3,  79.8,   77.5,   0,      0,      22.7,   14.2,   12.1,   0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     55,     Cs
    37441,  5989,   5624,   5247,   1293,   1137,   1063,   795.7,  780.5,  253.5,  192,    178.6,  92.6,   89.9,   0,      0,      30.3,   17,     14.8,   0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     56,     Ba
    38925,  6266,   5891,   5483,   1362,   1209,   1128,   853,    836,    274.7,  205.8,  196,    105.3,  102.5,  0,      0,      34.3,   19.3,   16.8,   0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     57,     La
    40443,  6549,   6164,   5723,   1436,   1274,   1187,   902.4,  883.8,  291,    223.2,  206.5,  109,    109,    0.1,    0.1,    37.8,   19.8,   17,     0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     58,     Ce
    41991,  6835,   6440,   5964,   1511,   1337,   1242,   948.3,  928.8,  304.5,  236.3,  217.6,  115.1,  115.1,  2,      2,      37.4,   22.3,   22.3,   0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     59,     Pr
    43569,  7126,   6722,   6208,   1575,   1403,   1297,   1003.3, 980.4,  319.2,  243.3,  224.6,  120.5,  120.5,  1.5,    1.5,    37.5,   21.1,   21.1,   0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     60,     Nd
    45184,  7428,   7013,   6459,   1655,   1471,   1357,   1052,   1027,   332.7,  242,    242,    120,    120,    1.5,    1.5,    37.5,   21.1,   21.1,   0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     61,     Pm
    46834,  7737,   7312,   6716,   1723,   1541,   1420,   1110.9, 1083.4, 347.2,  265.6,  247.4,  129,    129,    5.2,    5.2,    37.4,   21.3,   21.3,   0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     62,     Sm
    48519,  8052,   7617,   6977,   1800,   1614,   1481,   1158.6, 1127.5, 360,    284,    257,    133,    127.7,  0,      0,      32,     22,     22,     0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     63,     Eu
    50239,  8376,   7930,   7243,   1881,   1688,   1544,   1221.9, 1189.6, 378.6,  286,    271,    142.6,  142.6,  8.6,    8.6,    36,     28,     21,     0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     64,     Gd
    51996,  8708,   8252,   7514,   1968,   1768,   1611,   1276.9, 1241.1, 396,    322.4,  284.1,  150.5,  150.5,  7.7,    2.4,    45.6,   28.7,   22.6,   0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     65,     Tb
    53789,  9046,   8581,   7790,   2047,   1842,   1676,   1333,   1292.6, 414.2,  333.5,  293.2,  153.6,  153.6,  8,      4.3,    49.9,   26.3,   26.3,   0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     66,     Dy
    55618,  9394,   8918,   8071,   2128,   1923,   1741,   1392,   1351,   432.4,  343.5,  308.2,  160,    160,    8.6,    5.2,    49.3,   30.8,   24.1,   0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     67,     Ho
    57486,  9751,   9264,   8358,   2207,   2006,   1812,   1453,   1409,   449.8,  366.2,  320.2,  167.6,  167.6,  8,      4.7,    50.6,   31.4,   24.7,   0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     68,     Er
    59390,  10116,  9617,   8648,   2307,   2090,   1885,   1515,   1468,   470.9,  385.9,  332.6,  175.5,  175.5,  8,      4.6,    54.7,   31.8,   25,     0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     69,     Tm
    61332,  10486,  9978,   8944,   2398,   2173,   1950,   1576,   1528,   480.5,  388.7,  339.7,  191.2,  182.4,  2.5,    1.3,    52,     30.3,   24.1,   0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     70,     Yb
    63314,  10870,  10349,  9244,   2491,   2264,   2024,   1639,   1589,   506.8,  412.4,  359.2,  206.1,  196.3,  8.9,    7.5,    57.3,   33.6,   26.7,   0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     71,     Lu
    65351,  11271,  10739,  9561,   2601,   2365,   2108,   1716,   1662,   538,    438.2,  380.7,  220,    211.5,  15.9,   14.2,   64.2,   38,     29.9,   0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     72,     Hf
    67416,  11682,  11136,  9881,   2708,   2469,   2194,   1793,   1735,   563.4,  463.4,  400.9,  237.9,  226.4,  23.5,   21.6,   69.7,   42.2,   32.7,   0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     73,     Ta
    69525,  12100,  11544,  10207,  2820,   2575,   2281,   1872,   1809,   594.1,  490.4,  423.6,  255.9,  243.5,  33.6,   31.4,   75.6,   45.3,   36.8,   0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     74,     W
    71676,  12527,  11959,  10535,  2932,   2682,   2367,   1949,   1883,   625.4,  518.7,  446.8,  273.9,  260.5,  42.9,   40.5,   83,     45.6,   34.6,   0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     75,     Re
    73871,  12968,  12385,  10871,  3049,   2792,   2457,   2031,   1960,   658.2,  549.1,  470.7,  293.1,  278.5,  53.4,   50.7,   84,     58,     44.5,   0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     76,     Os
    76111,  13419,  12824,  11215,  3174,   2909,   2551,   2116,   2040,   691.1,  577.8,  495.8,  311.9,  296.3,  63.8,   60.8,   95.2,   63,     48,     0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     77,     Ir
    78395,  13880,  13273,  11564,  3296,   3027,   2645,   2202,   2122,   725.4,  609.1,  519.4,  331.6,  314.6,  74.5,   71.2,   101.7,  65.3,   51.7,   0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     78,     Pt
    80725,  14353,  13734,  11919,  3425,   3148,   2743,   2291,   2206,   762.1,  642.7,  546.3,  353.2,  335.1,  87.6,   84,     107.2,  74.2,   57.2,   0,      0,      0,      0,      0,      0,      0,      0,      0,      //,     79,     Au
    83102,  14839,  14209,  12284,  3562,   3279,   2847,   2385,   2295,   802.2,  680.2,  576.6,  378.2,  358.8,  104,    99.9,   127,    83.1,   64.5,   9.6,    7.8,    0,      0,      0,      0,      0,      0,      0,      //,     80,     Hg
    85530,  15347,  14698,  12658,  3704,   3416,   2957,   2485,   2389,   846.2,  720.5,  609.5,  405.7,  385,    122.2,  117.8,  136,    94.6,   73.5,   14.7,   12.5,   0,      0,      0,      0,      0,      0,      0,      //,     81,     Tl
    88005,  15861,  15200,  13035,  3851,   3554,   3066,   2586,   2484,   891.8,  761.9,  643.5,  434.3,  412.2,  141.7,  136.9,  147,    106.4,  83.3,   20.7,   18.1,   0,      0,      0,      0,      0,      0,      0,      //,     82,     Pb
    90526,  16388,  15711,  13419,  3999,   3696,   3177,   2688,   2580,   939,    805.2,  678.8,  464,    440.1,  162.3,  157,    159.3,  119,    92.6,   26.9,   23.8,   0,      0,      0,      0,      0,      0,      0,      //,     83,     Bi
    93105,  16939,  16244,  13814,  4149,   3854,   3302,   2798,   2683,   995,    851,    705,    500,    473,    184,    184,    177,    132,    104,    31,     31,     0,      0,      0,      0,      0,      0,      0,      //,     84,     Po
    95730,  17493,  16785,  14214,  4317,   4008,   3426,   2909,   2787,   1042,   886,    740,    533,    507,    210,    210,    195,    148,    115,    40,     40,     0,      0,      0,      0,      0,      0,      0,      //,     85,     At
    98404,  18049,  17337,  14619,  4482,   4159,   3538,   3022,   2892,   1097,   929,    768,    567,    541,    238,    238,    214,    164,    127,    48,     48,     0,      0,      26,     0,      0,      0,      0,      //,     86,     Rn
    101137, 18639,  17907,  15031,  4652,   4327,   3663,   3136,   3000,   1153,   980,    810,    603,    577,    268,    268,    234,    182,    140,    58,     58,     0,      0,      34,     15,     15,     0,      0,      //,     87,     Fr
    103922, 19237,  18484,  15444,  4822,   4490,   3792,   3248,   3105,   1208,   1058,   879,    636,    603,    299,    299,    254,    200,    153,    68,     68,     0,      0,      44,     19,     19,     0,      0,      //,     88,     Ra
    106755, 19840,  19083,  15871,  5002,   4656,   3909,   3370,   3219,   1269,   1080,   890,    675,    639,    319,    319,    272,    215,    167,    80,     80,     0,      0,      0,      0,      0,      0,      0,      //,     89,     Ac
    109651, 20472,  19693,  16300,  5182,   4830,   4046,   3491,   3332,   1330,   1168,   966.4,  712.1,  675.2,  342.4,  333.1,  290,    229,    182,    92.5,   85.4,   0,      0,      41.4,   24.5,   16.6,   0,      0,      //,     90,     Th
    112601, 21105,  20314,  16733,  5367,   5001,   4174,   3611,   3442,   1387,   1224,   1007,   743,    708,    371,    360,    310,    232,    232,    94,     94,     0,      0,      0,      0,      0,      0,      0,      //,     91,     Pa
    115606, 21757,  20948,  17166,  5548,   5182,   4303,   3728,   3552,   1439,   1271,   1043,   778.3,  736.2,  388.2,  377.4,  321,    257,    192,    102.8,  94.2,   0,      0,      43.9,   26.8,   16.8,   0,      0,      //,     92,     U
    118669.0, 22427.0, 21600.0, 17610.0,  5739.0,  5366.0,  4435.0,  3849.0,  3664.0,  1501.0,  1328.0,  1085.0,   816.0,   771.0, 	//	Np 93
    414.0,   403.0,   338.0,   274.0,   206.0,   109.0,   101.0,     0.0,     0.0,    47.0,    29.0,    18.0,     0.0,     0.0,
    121791.0, 23104.0, 22266.0, 18057.0,  5933.0,  5547.0,  4563.0,  3970.0,  3775.0,  1559.0,  1380.0,  1123.0,   846.0,   798.0, 	//	Pu 94
    436.0,   424.0,   350.0,   283.0,   213.0,   113.0,   102.0,     0.0,     0.0,    46.0,    29.0,    16.0,     0.0,     0.0,
    124982.0, 23808.0, 22952.0, 18510.0,  6133.0,  5739.0,  4698.0,  4096.0,  3890.0,  1620.0,  1438.0,  1165.0,   880.0,   829.0, 	//	Am 95
    461.0,   446.0,   365.0,   298.0,   219.0,   116.0,   106.0,     0.0,     0.0,    48.0,    29.0,    16.0,     0.0,     0.0,
    128241.0, 24526.0, 23651.0, 18970.0,  6337.0,  5937.0,  4838.0,  4224.0,  4009.0,  1684.0,  1498.0,  1207.0,   916.0,   862.0, 	//	Cm 96
    484.0,   470.0,   383.0,   313.0,   229.0,   124.0,   110.0,     0.0,     0.0,    50.0,    30.0,    16.0,     0.0,     0.0,
    131556.0, 25256.0, 24371.0, 19435.0,  6545.0,  6138.0,  4976.0,  4353.0,  4127.0,  1748.0,  1558.0,  1249.0,   955.0,   898.0, 	//	Bk 97
    511.0,   495.0,   399.0,   326.0,   237.0,   130.0,   117.0,     0.0,     0.0,    52.0,    32.0,    16.0,     0.0,     0.0,
    134939.0, 26010.0, 25108.0, 19907.0,  6761.0,  6345.0,  5116.0,  4484.0,  4247.0,  1813.0,  1620.0,  1292.0,   991.0,   930.0, 	//	Cf 98
    538.0,   520.0,   416.0,   341.0,   245.0,   137.0,   122.0,     0.0,     0.0,    54.0,    33.0,    17.0,     0.0,     0.0,
    138396.0, 26782.0, 25865.0, 20384.0,  6981.0,  6558.0,  5259.0,  4617.0,  4368.0,  1883.0,  1630.0,  1336.0,  1029.0,   965.0, 	//	Es 99
    564.0,   546.0,   434.0,   357.0,   255.0,   142.0,   127.0,     0.0,     0.0,    57.0,    35.0,    17.0,     0.0,     0.0,
    141926.0, 27574.0, 26641.0, 20868.0,  7208.0,  6776.0,  5405.0,  4752.0,  4491.0,  1952.0,  1749.0,  1379.0,  1067.0,  1000.0, 	//	Fm 100
    591.0,   572.0,   452.0,   373.0,   262.0,   149.0,   133.0,     0.0,     0.0,    59.0,    36.0,    17.0,     0.0,     0.0,
};

/*
//  Replaced with updated values    Jun 28, 2012    see above

//  Gwyn P. Williams, "Electron Binding Energies of the Elements",
//   CRC Handbook of Chemistry and Physics, R. C. Weast, Ed. 66th Edition F170-3 (1985).
//[Additional values for Z=93-103 from Browne & Firestone "Table of Radioactive Isotopes" John Wiley 1986.]
//[Some low levels added to make emission lines agree with Bearden and Burr, Rev. Mod. Phys. 39 (1967) 78-124,
// as quoted in CRC Handbook 51st Ed. (1970-71) E126-E164.]

//    K       L1       L2       L3       M1       M2       M3       M4       M5       N1       N2       N3       N4       N5
//   N6       N7       O1       O2       O3       O4       O5       O6       O7       P1       P2       P3       P4       Q1
	const float XrayEdge::EDGE_ENERGIES[MAXZ*(MAXINDEX+1)] = {
   13.6,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, 	//	H 1
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
   24.6,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, 	//	He 2
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
   54.7,     5.3,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, 	//	Li 3
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
  111.5,     8.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, 	//	Be 4
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
  188.0,    12.6,     4.7,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, 	//	B 5
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
  284.2,    18.0,     7.2,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, 	//	C 6
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
  409.9,    37.3,    17.5,    17.5,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, 	//	N 7
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
  543.1,    41.6,    18.2,    18.2,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, 	//	O 8
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
  696.7,    45.0,    19.9,    19.9,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, 	//	F 9
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
  870.2,    48.5,    21.7,    21.6,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, 	//	Ne 10
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
 1070.8,    63.5,    30.4,    30.5,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, 	//	Na 11
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
 1303.0,    88.6,    49.6,    49.2,     2.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, 	//	Mg 12
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
 1559.0,   117.8,    72.9,    72.5,     4.0,     2.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, 	//	Al 13
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
 1839.0,   149.7,    99.8,    99.2,     8.0,     2.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, 	//	Si 14
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
 2145.5,   189.0,   136.0,   135.0,    12.0,     7.0,     6.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, 	//	P 15
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
 2472.0,   230.9,   163.6,   162.5,    14.0,     8.0,     7.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, 	//	S 16
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
 2822.0,   270.0,   202.0,   200.0,    18.0,    10.0,    10.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, 	//	Cl 17
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
 3205.9,   326.3,   250.6,   248.4,    29.3,    15.9,    15.7,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, 	//	Ar 18
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
 3608.4,   378.6,   297.3,   294.6,    34.8,    18.3,    18.3,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, 	//	K 19
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
 4038.5,   438.4,   349.7,   346.2,    44.3,    25.4,    25.4,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, 	//	Ca 20
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
 4492.0,   498.0,   403.6,   398.7,    51.1,    28.3,    28.3,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, 	//	Sc 21
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
 4966.0,   560.9,   460.2,   453.8,    58.7,    32.6,    32.6,     0.0,     2.0,     0.0,     0.0,     0.0,     0.0,     0.0, 	//	Ti 22
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
 5465.0,   626.7,   519.8,   512.1,    66.3,    37.2,    37.2,     0.0,     2.0,     0.0,     0.0,     0.0,     0.0,     0.0, 	//	V 23
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
 5989.0,   696.0,   583.8,   574.1,    74.1,    42.2,    42.2,     2.0,     2.0,     0.0,     0.0,     0.0,     0.0,     0.0, 	//	Cr 24
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
 6539.0,   769.1,   649.9,   638.7,    82.3,    47.2,    47.2,     2.0,     2.0,     0.0,     0.0,     0.0,     0.0,     0.0, 	//	Mn 25
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
 7112.0,   844.6,   719.9,   706.8,    91.3,    52.7,    52.7,     2.0,     2.0,     0.0,     0.0,     0.0,     0.0,     0.0, 	//	Fe 26
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
 7709.0,   925.1,   793.2,   778.1,   101.0,    58.9,    59.9,     3.0,     3.0,     0.0,     0.0,     0.0,     0.0,     0.0, 	//	Co 27
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
 8333.0,  1008.6,   870.0,   852.7,   110.8,    68.0,    66.2,     4.0,     4.0,     0.0,     0.0,     0.0,     0.0,     0.0, 	//	Ni 28
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
 8979.0,  1096.7,   952.3,   932.7,   122.5,    77.3,    75.1,     5.0,     5.0,     0.0,     0.0,     0.0,     0.0,     0.0, 	//	Cu 29
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
 9659.0,  1196.2,  1044.9,  1021.8,   139.8,    91.4,    88.6,    10.2,    10.1,     0.0,     0.0,     0.0,     0.0,     0.0, 	//	Zn 30
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
10367.0,  1299.0,  1143.2,  1116.4,   159.5,   103.5,   100.0,    18.7,    18.7,     1.0,     0.0,     2.0,     0.0,     0.0, 	//	Ga 31
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
11103.0,  1414.6,  1248.1,  1217.0,   180.1,   124.9,   120.8,    29.8,    29.2,     5.0,     3.0,     0.0,     0.0,     0.0, 	//	Ge 32
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
11867.0,  1527.0,  1359.1,  1323.6,   204.7,   146.2,   141.2,    41.7,    41.7,     8.0,     3.0,     3.0,     0.0,     0.0, 	//	As 33
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
12658.0,  1652.0,  1474.3,  1433.9,   229.6,   166.5,   160.7,    55.5,    54.6,    12.0,     3.0,     3.0,     0.0,     0.0, 	//	Se 34
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
13474.0,  1782.0,  1596.0,  1550.0,   257.0,   189.0,   182.0,    70.0,    69.0,    27.0,     3.0,     3.0,     0.0,     0.0, 	//	Br 35
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
14326.0,  1921.0,  1730.9,  1678.4,   292.8,   222.2,   214.4,    95.0,    93.8,    27.5,    14.1,    14.1,     0.0,     0.0, 	//	Kr 36
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
15200.0,  2065.0,  1864.0,  1804.0,   326.7,   248.7,   239.1,   113.0,   112.0,    30.5,    16.3,    15.3,     0.0,     0.0, 	//	Rb 37
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
16105.0,  2216.0,  2007.0,  1940.0,   358.7,   280.3,   270.0,   136.0,   134.2,    38.9,    21.6,    20.1,     0.0,     0.0, 	//	Sr 38
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
17038.0,  2373.0,  2156.0,  2080.0,   392.0,   310.6,   298.8,   157.7,   155.8,    43.8,    24.4,    23.1,     0.0,     0.0, 	//	Y 39
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
17998.0,  2532.0,  2307.0,  2223.0,   430.3,   343.5,   329.8,   181.1,   178.8,    50.6,    28.5,    27.1,     0.0,     0.0, 	//	Zr 40
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
18986.0,  2698.0,  2465.0,  2371.0,   466.6,   376.1,   360.6,   205.0,   202.3,    56.4,    32.6,    30.8,     0.0,     0.0, 	//	Nb 41
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
20000.0,  2866.0,  2625.0,  2520.0,   506.3,   411.6,   394.0,   231.1,   227.9,    63.2,    37.6,    35.5,     0.0,     0.0, 	//	Mo 42
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
21044.0,  3043.0,  2793.0,  2677.0,   544.0,   447.6,   417.7,   257.6,   253.9,    69.5,    42.3,    39.9,     0.0,     0.0, 	//	Tc 43
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
22117.0,  3224.0,  2967.0,  2838.0,   586.1,   483.3,   461.5,   284.2,   280.0,    75.0,    46.3,    43.2,     0.0,     0.0, 	//	Ru 44
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
23220.0,  3412.0,  3146.0,  3004.0,   628.1,   521.3,   496.5,   311.9,   307.2,    81.4,    50.5,    47.3,     2.0,     2.0, 	//	Rh 45
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
24350.0,  3604.0,  3330.0,  3173.0,   671.6,   559.9,   532.3,   340.5,   335.2,    87.1,    55.7,    50.9,     2.0,     2.0, 	//	Pd 46
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
25514.0,  3806.0,  3524.0,  3351.0,   719.0,   603.8,   573.0,   374.0,   368.3,    97.0,    63.7,    58.3,     4.0,     4.0, 	//	Ag 47
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
26711.0,  4018.0,  3727.0,  3538.0,   772.0,   652.6,   618.4,   411.9,   405.2,   109.8,    63.9,    63.9,    11.7,    10.7, 	//	Cd 48
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
27940.0,  4238.0,  3938.0,  3730.0,   827.2,   703.2,   665.3,   451.4,   443.9,   122.9,    73.5,    73.5,    17.7,    16.9, 	//	In 49
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
29200.0,  4465.0,  4156.0,  3929.0,   884.7,   756.5,   714.6,   493.2,   484.9,   137.1,    83.6,    83.6,    24.9,    23.9, 	//	Sn 50
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
30491.0,  4698.0,  4380.0,  4132.0,   940.0,   812.7,   766.4,   537.5,   528.2,   153.2,    95.6,    95.6,    33.3,    32.1, 	//	Sb 51
    0.0,     0.0,     7.0,     2.0,     2.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
31814.0,  4939.0,  4612.0,  4341.0,  1006.0,   870.8,   820.8,   583.4,   573.0,   169.4,   103.3,   103.3,    41.9,    40.4, 	//	Te 52
    0.0,     0.0,    12.0,     2.0,     2.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
33169.0,  5188.0,  4852.0,  4557.0,  1072.0,   931.0,   875.0,   630.8,   619.3,   186.0,   123.0,   123.0,    50.6,    48.9, 	//	I 53
    0.0,     0.0,    14.0,     3.0,     3.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
34561.0,  5453.0,  5107.0,  4786.0,  1148.7,  1002.1,   940.6,   689.0,   676.4,   213.2,   146.7,   145.5,    69.5,    67.5, 	//	Xe 54
    0.0,     0.0,    23.3,    13.4,    12.1,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
35985.0,  5714.0,  5359.0,  5012.0,  1211.0,  1071.0,  1003.0,   740.5,   726.6,   232.3,   172.4,   161.3,    79.8,    77.5, 	//	Cs 55
    0.0,     0.0,    22.7,    14.2,    12.1,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
37441.0,  5989.0,  5624.0,  5247.0,  1293.0,  1137.0,  1063.0,   795.7,   780.5,   253.5,   192.0,   178.6,    92.6,    89.9, 	//	Ba 56
    0.0,     0.0,    30.3,    17.0,    14.8,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
38925.0,  6266.0,  5891.0,  5483.0,  1362.0,  1209.0,  1128.0,   853.0,   836.0,   274.7,   205.8,   196.0,   105.3,   102.5, 	//	La 57
    0.0,     0.0,    34.3,    19.3,    16.8,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
40443.0,  6548.0,  6164.0,  5723.0,  1436.0,  1274.0,  1187.0,   902.4,   883.8,   291.0,   223.2,   206.5,   109.0,   109.0, 	//	Ce 58
    0.1,     0.0,    37.8,    19.8,    17.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
41991.0,  6835.0,  6440.0,  5964.0,  1511.0,  1337.0,  1242.0,   948.3,   928.8,   304.5,   236.3,   217.6,   115.1,   115.1, 	//	Pr 59
    2.0,     0.0,    37.4,    22.3,    22.3,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
43569.0,  7126.0,  6722.0,  6208.0,  1575.0,  1403.0,  1297.0,  1003.3,   980.4,   319.2,   243.3,   224.6,   120.5,   120.5, 	//	Nd 60
    1.5,     0.0,    37.5,    21.1,    21.1,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
45184.0,  7428.0,  7013.0,  6459.0,  1650.0,  1471.4,  1357.0,  1052.0,  1027.0,   331.0,   242.0,   242.0,   120.0,   120.0, 	//	Pm 61
    4.0,     4.0,    38.0,    22.0,    22.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
46834.0,  7737.0,  7312.0,  6716.0,  1723.0,  1541.0,  1419.8,  1110.9,  1083.4,   347.2,   265.6,   247.4,   129.0,   129.0, 	//	Sm 62
    5.2,     5.2,    37.4,    21.3,    21.3,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
48519.0,  8052.0,  7617.0,  6977.0,  1800.0,  1614.0,  1481.0,  1158.6,  1127.5,   360.0,   284.0,   257.0,   133.0,   127.7, 	//	Eu 63
    6.0,     6.0,    32.0,    22.0,    22.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
50239.0,  8376.0,  7930.0,  7243.0,  1881.0,  1688.0,  1544.0,  1221.9,  1189.6,   378.6,   286.0,   271.0,   142.6,   142.6, 	//	Gd 64
    8.6,     8.6,    36.0,    20.0,    20.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
51996.0,  8708.0,  8252.0,  7514.0,  1968.0,  1768.0,  1611.0,  1276.9,  1241.1,   396.0,   322.4,   284.1,   150.5,   150.5, 	//	Tb 65
    7.7,     2.4,    45.6,    28.7,    22.6,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
53789.0,  9046.0,  8581.0,  7790.0,  2047.0,  1842.0,  1676.0,  1333.0,  1292.0,   414.2,   333.5,   293.2,   153.6,   153.6, 	//	Dy 66
    8.0,     4.3,    49.9,    26.3,    26.3,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
55618.0,  9394.0,  8918.0,  8071.0,  2128.0,  1923.0,  1741.0,  1392.0,  1351.0,   432.4,   343.5,   308.2,   160.0,   160.0, 	//	Ho 67
    8.6,     5.2,    49.3,    30.8,    24.1,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
57486.0,  9751.0,  9264.0,  8358.0,  2206.0,  2006.0,  1812.0,  1453.0,  1409.0,   449.8,   366.2,   320.2,   167.6,   167.6, 	//	Er 68
    4.7,     4.7,    50.6,    31.4,    24.7,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
59390.0, 10116.0,  9617.0,  8648.0,  2307.0,  2090.0,  1885.0,  1515.0,  1468.0,   470.9,   385.9,   332.6,   175.5,   175.5, 	//	Tm 69
    4.6,     4.6,    54.7,    31.8,    25.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
61332.0, 10486.0,  9978.0,  8944.0,  2398.0,  2173.0,  1950.0,  1576.0,  1528.0,   480.5,   388.7,   339.7,   191.2,   182.4, 	//	Yb 70
    2.5,     1.3,    52.0,    30.3,    24.1,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
63314.0, 10870.0, 10349.0,  9244.0,  2491.0,  2264.0,  2024.0,  1639.0,  1589.0,   506.8,   412.4,   359.2,   206.1,   196.3, 	//	Lu 71
    8.9,     7.5,    57.3,    33.6,    26.7,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
65351.0, 11271.0, 10739.0,  9561.0,  2601.0,  2365.0,  2107.0,  1716.0,  1662.0,   538.0,   438.2,   380.7,   220.0,   211.5, 	//	Hf 72
   15.9,    14.2,    64.2,    38.0,    29.9,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
67416.0, 11682.0, 11136.0,  9881.0,  2708.0,  2469.0,  2194.0,  1793.0,  1735.0,   563.4,   463.4,   400.9,   237.9,   226.4, 	//	Ta 73
   23.5,    21.6,    69.7,    42.2,    32.7,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
69525.0, 12100.0, 11544.0, 10207.0,  2820.0,  2575.0,  2281.0,  1872.0,  1809.0,   594.1,   490.4,   423.6,   255.9,   243.5, 	//	W 74
   33.6,    31.4,    75.6,    45.3,    36.8,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
71676.0, 12527.0, 11959.0, 10535.0,  2932.0,  2682.0,  2367.0,  1949.0,  1883.0,   625.4,   518.7,   446.8,   273.9,   260.5, 	//	Re 75
   42.9,    40.5,    83.0,    45.6,    34.6,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
73871.0, 12968.0, 12385.0, 10871.0,  3049.0,  2792.0,  2457.0,  2031.0,  1960.0,   658.2,   549.1,   470.7,   293.1,   278.5, 	//	Os 76
   53.4,    50.7,    84.0,    58.0,    44.5,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
76111.0, 13419.0, 12824.0, 11215.0,  3174.0,  2909.0,  2551.0,  2116.0,  2040.0,   691.1,   577.8,   495.8,   311.9,   296.3, 	//	Ir 77
   63.8,    60.8,    95.2,    63.0,    48.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
78395.0, 13880.0, 13273.0, 11564.0,  3296.0,  3027.0,  2645.0,  2202.0,  2122.0,   725.4,   609.1,   519.4,   331.6,   314.6, 	//	Pt 78
   74.5,    71.2,   101.7,    65.3,    51.7,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
80725.0, 14353.0, 13734.0, 11919.0,  3425.0,  3148.0,  2743.0,  2291.0,  2206.0,   762.1,   642.7,   546.3,   353.2,   335.1, 	//	Au 79
   87.6,    83.9,   107.2,    74.2,    57.2,     5.0,     5.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
83102.0, 14839.0, 14209.0, 12284.0,  3562.0,  3279.0,  2847.0,  2385.0,  2295.0,   802.2,   680.2,   576.6,   378.2,   358.8, 	//	Hg 80
  104.0,    99.9,   127.0,    83.1,    64.5,     9.6,     7.8,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
85530.0, 15347.0, 14698.0, 12658.0,  3704.0,  3416.0,  2957.0,  2485.0,  2389.0,   846.2,   720.5,   609.5,   405.7,   385.0, 	//	Tl 81
  122.2,   117.8,   136.0,    94.6,    73.5,    14.7,    12.5,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
88005.0, 15861.0, 15200.0, 13035.0,  3851.0,  3554.0,  3066.0,  2586.0,  2484.0,   891.8,   761.9,   643.5,   434.3,   412.2, 	//	Pb 82
  141.7,   136.9,   147.0,   106.4,    83.3,    20.7,    18.1,     0.0,     0.0,     3.0,     1.0,     0.0,     0.0,     0.0,
90526.0, 16388.0, 15711.0, 13419.0,  3999.0,  3696.0,  3177.0,  2688.0,  2580.0,   939.0,   805.2,   678.8,   464.0,   440.1, 	//	Bi 83
  162.3,   157.0,   159.3,   119.0,    92.6,    26.9,    23.8,     0.0,     0.0,     8.0,     3.0,     3.0,     0.0,     0.0,
93105.0, 16939.0, 16244.0, 13814.0,  4149.0,  3854.0,  3302.0,  2798.0,  2683.0,   995.0,   851.0,   705.0,   500.0,   473.0, 	//	Po 84
  184.0,   184.0,   177.0,   132.0,   104.0,    31.0,    31.0,     0.0,     0.0,     9.0,     4.0,     1.0,     0.0,     0.0,
95730.0, 17493.0, 16785.0, 14214.0,  4317.0,  4008.0,  3426.0,  2909.0,  2787.0,  1042.0,   886.0,   740.0,   533.0,   507.0, 	//	At 85
  210.0,   210.0,   195.0,   148.0,   115.0,    40.0,    40.0,     0.0,     0.0,    13.0,     6.0,     1.0,     0.0,     0.0,
98404.0, 18049.0, 17337.0, 14619.0,  4482.0,  4159.0,  3538.0,  3022.0,  2892.0,  1097.0,   929.0,   768.0,   567.0,   541.0, 	//	Rn 86
  238.0,   238.0,   214.0,   164.0,   127.0,    48.0,    48.0,     0.0,     0.0,    16.0,     8.0,     2.0,     0.0,     0.0,
101137.0, 18639.0, 17907.0, 15031.0,  4652.0,  4327.0,  3663.0,  3136.0,  3000.0,  1153.0,   980.0,   810.0,   603.0,   577.0, 	//	Fr 87
  268.0,   268.0,   234.0,   182.0,   140.0,    58.0,    58.0,     0.0,     0.0,    24.0,    14.0,     7.0,     0.0,     0.0,
103922.0, 19237.0, 18484.0, 15444.0,  4822.0,  4490.0,  3792.0,  3248.0,  3105.0,  1208.0,  1058.0,   879.0,   636.0,   603.0, 	//	Ra 88
  299.0,   299.0,   254.0,   200.0,   153.0,    68.0,    68.0,     0.0,     0.0,    31.0,    20.0,    12.0,     0.0,     0.0,
106755.0, 19840.0, 19083.0, 15871.0,  5002.0,  4656.0,  3909.0,  3370.0,  3219.0,  1269.0,  1080.0,   890.0,   675.0,   639.0, 	//	Ac 89
  319.0,   319.0,   272.0,   215.0,   167.0,    80.0,    80.0,     0.0,     0.0,    37.0,    24.0,    15.0,     0.0,     0.0,
109651.0, 20472.0, 19693.0, 16300.0,  5182.0,  4830.0,  4046.0,  3491.0,  3332.0,  1330.0,  1168.0,   966.4,   712.1,   675.2, 	//	Th 90
  342.4,   333.1,   290.0,   229.0,   182.0,    92.5,    85.4,     0.0,     0.0,    41.4,    24.5,    16.6,     0.0,     0.0,
112601.0, 21105.0, 20314.0, 16733.0,  5367.0,  5001.0,  4174.0,  3611.0,  3442.0,  1387.0,  1224.0,  1007.0,   743.0,   708.0, 	//	Pa 91
  371.0,   360.0,   310.0,   232.0,   187.0,    94.0,    94.0,     0.0,     0.0,    43.0,    27.0,    17.0,     0.0,     0.0,
115606.0, 21757.0, 20948.0, 17166.0,  5548.0,  5182.0,  4303.0,  3728.0,  3552.0,  1439.0,  1271.0,  1043.0,   778.3,   736.2, 	//	U 92
  388.2,   377.4,   321.0,   257.0,   192.0,   102.8,    94.2,     0.0,     0.0,    43.9,    26.8,    16.8,     0.0,     0.0,
118669.0, 22427.0, 21600.0, 17610.0,  5739.0,  5366.0,  4435.0,  3849.0,  3664.0,  1501.0,  1328.0,  1085.0,   816.0,   771.0, 	//	Np 93
  414.0,   403.0,   338.0,   274.0,   206.0,   109.0,   101.0,     0.0,     0.0,    47.0,    29.0,    18.0,     0.0,     0.0,
121791.0, 23104.0, 22266.0, 18057.0,  5933.0,  5547.0,  4563.0,  3970.0,  3775.0,  1559.0,  1380.0,  1123.0,   846.0,   798.0, 	//	Pu 94
  436.0,   424.0,   350.0,   283.0,   213.0,   113.0,   102.0,     0.0,     0.0,    46.0,    29.0,    16.0,     0.0,     0.0,
124982.0, 23808.0, 22952.0, 18510.0,  6133.0,  5739.0,  4698.0,  4096.0,  3890.0,  1620.0,  1438.0,  1165.0,   880.0,   829.0, 	//	Am 95
  461.0,   446.0,   365.0,   298.0,   219.0,   116.0,   106.0,     0.0,     0.0,    48.0,    29.0,    16.0,     0.0,     0.0,
128241.0, 24526.0, 23651.0, 18970.0,  6337.0,  5937.0,  4838.0,  4224.0,  4009.0,  1684.0,  1498.0,  1207.0,   916.0,   862.0, 	//	Cm 96
  484.0,   470.0,   383.0,   313.0,   229.0,   124.0,   110.0,     0.0,     0.0,    50.0,    30.0,    16.0,     0.0,     0.0,
131556.0, 25256.0, 24371.0, 19435.0,  6545.0,  6138.0,  4976.0,  4353.0,  4127.0,  1748.0,  1558.0,  1249.0,   955.0,   898.0, 	//	Bk 97
  511.0,   495.0,   399.0,   326.0,   237.0,   130.0,   117.0,     0.0,     0.0,    52.0,    32.0,    16.0,     0.0,     0.0,
134939.0, 26010.0, 25108.0, 19907.0,  6761.0,  6345.0,  5116.0,  4484.0,  4247.0,  1813.0,  1620.0,  1292.0,   991.0,   930.0, 	//	Cf 98
  538.0,   520.0,   416.0,   341.0,   245.0,   137.0,   122.0,     0.0,     0.0,    54.0,    33.0,    17.0,     0.0,     0.0,
138396.0, 26782.0, 25865.0, 20384.0,  6981.0,  6558.0,  5259.0,  4617.0,  4368.0,  1883.0,  1630.0,  1336.0,  1029.0,   965.0, 	//	Es 99
  564.0,   546.0,   434.0,   357.0,   255.0,   142.0,   127.0,     0.0,     0.0,    57.0,    35.0,    17.0,     0.0,     0.0,
141926.0, 27574.0, 26641.0, 20868.0,  7208.0,  6776.0,  5405.0,  4752.0,  4491.0,  1952.0,  1749.0,  1379.0,  1067.0,  1000.0, 	//	Fm 100
  591.0,   572.0,   452.0,   373.0,   262.0,   149.0,   133.0,     0.0,     0.0,    59.0,    36.0,    17.0,     0.0,     0.0,
		};
*/

//		ocupancies from Jzero reference, see below
//K  L1  L2  L3  M1  M2  M3  M4  M5  N1  N2  N3  N4  N5  N6  N7  O1  O2 O3  O4  O5  O6  O7  P1  P2  P3  P4  Q1
const int XrayEdge::EDGE_OCCUPANCIES[MAXZ*(MAXINDEX+1)] = {
 1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	H 1
 2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	He 2
 2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Li 3
 2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Be 4
 2,  2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	B 5
 2,  2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	C 6
 2,  2,  2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	N 7
 2,  2,  2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	O 8
 2,  2,  2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	F 9
 2,  2,  2,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Ne 10
 2,  2,  2,  4,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Na 11
 2,  2,  2,  4,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Mg 12
 2,  2,  2,  4,  2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Al 13
 2,  2,  2,  4,  2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Si 14
 2,  2,  2,  4,  2,  2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	P 15
 2,  2,  2,  4,  2,  2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	S 16
 2,  2,  2,  4,  2,  2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Cl 17
 2,  2,  2,  4,  2,  2,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Ar 18
 2,  2,  2,  4,  2,  2,  4,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	K 19
 2,  2,  2,  4,  2,  2,  4,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Ca 20
 2,  2,  2,  4,  2,  2,  4,  1,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Sc 21
 2,  2,  2,  4,  2,  2,  4,  0,  2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Ti 22
 2,  2,  2,  4,  2,  2,  4,  0,  3,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	V 23
 2,  2,  2,  4,  2,  2,  4,  4,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Cr 24
 2,  2,  2,  4,  2,  2,  4,  4,  1,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Mn 25
 2,  2,  2,  4,  2,  2,  4,  4,  2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Fe 26
 2,  2,  2,  4,  2,  2,  4,  4,  3,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Co 27
 2,  2,  2,  4,  2,  2,  4,  4,  4,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Ni 28
 2,  2,  2,  4,  2,  2,  4,  4,  6,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Cu 29
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Zn 30
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Ga 31
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Ge 32
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	As 33
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Se 34
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Br 35
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Kr 36
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Rb 37
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Sr 38
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Y 39
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  2,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Zr 40
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Nb 41
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Mo 42
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  1,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Tc 43
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  3,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Ru 44
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  4,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Rh 45
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Pd 46
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Ag 47
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Cd 48
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  0,  0,  2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	In 49
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  0,  0,  2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Sn 50
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  0,  0,  2,  2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Sb 51
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  0,  0,  2,  2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Te 52
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  0,  0,  2,  2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	I 53
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  0,  0,  2,  2,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0, 	//	Xe 54
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  0,  0,  2,  2,  4,  0,  0,  0,  0,  1,  0,  0,  0,  0, 	//	Cs 55
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  0,  0,  2,  2,  4,  0,  0,  0,  0,  2,  0,  0,  0,  0, 	//	Ba 56
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  0,  0,  2,  2,  4,  1,  0,  0,  0,  2,  0,  0,  0,  0, 	//	La 57
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  1,  0,  2,  2,  4,  1,  0,  0,  0,  2,  0,  0,  0,  0, 	//	Ce 58
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  3,  0,  2,  2,  4,  0,  0,  0,  0,  2,  0,  0,  0,  0, 	//	Pr 59
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  4,  0,  2,  2,  4,  0,  0,  0,  0,  2,  0,  0,  0,  0, 	//	Nd 60
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  4,  1,  2,  2,  4,  0,  0,  0,  0,  2,  0,  0,  0,  0, 	//	Pm 61
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  4,  2,  2,  2,  4,  0,  0,  0,  0,  2,  0,  0,  0,  0, 	//	Sm 62
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  1,  2,  2,  4,  0,  0,  0,  0,  2,  0,  0,  0,  0, 	//	Eu 63
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  1,  2,  2,  4,  1,  0,  0,  0,  2,  0,  0,  0,  0, 	//	Gd 64
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  3,  2,  2,  4,  0,  0,  0,  0,  2,  0,  0,  0,  0, 	//	Tb 65
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  4,  2,  2,  4,  0,  0,  0,  0,  2,  0,  0,  0,  0, 	//	Dy 66
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  5,  2,  2,  4,  0,  0,  0,  0,  2,  0,  0,  0,  0, 	//	Ho 67
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  6,  2,  2,  4,  0,  0,  0,  0,  2,  0,  0,  0,  0, 	//	Er 68
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  7,  2,  2,  4,  0,  0,  0,  0,  2,  0,  0,  0,  0, 	//	Tm 69
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  0,  0,  0,  0,  2,  0,  0,  0,  0, 	//	Yb 70
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  1,  0,  0,  0,  2,  0,  0,  0,  0, 	//	Lu 71
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  2,  0,  0,  0,  2,  0,  0,  0,  0, 	//	Hf 72
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  3,  0,  0,  0,  2,  0,  0,  0,  0, 	//	Ta 73
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  0,  0,  0,  2,  0,  0,  0,  0, 	//	W 74
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  1,  0,  0,  2,  0,  0,  0,  0, 	//	Re 75
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  2,  0,  0,  2,  0,  0,  0,  0, 	//	Os 76
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  3,  0,  0,  2,  0,  0,  0,  0, 	//	Ir 77
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  5,  0,  0,  1,  0,  0,  0,  0, 	//	Pt 78
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  0,  0,  1,  0,  0,  0,  0, 	//	Au 79
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  0,  0,  2,  0,  0,  0,  0, 	//	Hg 80
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  0,  0,  2,  1,  0,  0,  0, 	//	Tl 81
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  0,  0,  2,  2,  0,  0,  0, 	//	Pb 82
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  0,  0,  2,  2,  1,  0,  0, 	//	Bi 83
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  0,  0,  2,  2,  2,  0,  0, 	//	Po 84
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  0,  0,  2,  2,  3,  0,  0, 	//	At 85
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  0,  0,  2,  2,  4,  0,  0, 	//	Rn 86
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  0,  0,  2,  2,  4,  0,  1, 	//	Fr 87
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  0,  0,  2,  2,  4,  0,  2, 	//	Ra 88
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  0,  0,  2,  2,  4,  1,  2, 	//	Ac 89
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  0,  0,  2,  2,  4,  2,  2, 	//	Th 90
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  2,  0,  2,  2,  4,  1,  2, 	//	Pa 91
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  3,  0,  2,  2,  4,  1,  2, 	//	U 92
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  4,  0,  2,  2,  4,  1,  2, 	//	Np 93
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  6,  0,  2,  2,  4,  0,  2, 	//	Pu 94
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  6,  1,  2,  2,  4,  0,  2, 	//	Am 95
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  6,  1,  2,  2,  4,  1,  2, 	//	Cm 96
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  6,  3,  2,  2,  4,  0,  2, 	//	Bk 97
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  6,  4,  2,  2,  4,  0,  2, 	//	Cf 98
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  6,  5,  2,  2,  4,  0,  2, 	//	Es 99
 2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  6,  6,  2,  2,  4,  0,  2, 	//	Fm 100
		};

//		Taken from F. Biggs, L. Mendelsohn, and J. Mann, At. Data and Nucl. Data Tables Vol. 16, pp201-309 (1975).
	const float XrayEdge::EDGE_JZEROS[MAXZ*(MAXINDEX+1)] = {
0.84900, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 	//	H 1
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.53500, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 	//	He 2
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.32900, 1.94000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 	//	Li 3
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.23700, 1.34000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 	//	Be 4
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.18600, 1.00000, 0.61500, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 	//	B 5
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.15300, 0.80400, 0.48800, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 	//	C 6
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.13000, 0.67200, 0.40700, 0.40700, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 	//	N 7
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.11300, 0.57900, 0.35000, 0.35000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 	//	O 8
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.10000, 0.50800, 0.30700, 0.30700, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 	//	F 9
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.09000, 0.45300, 0.27400, 0.27400, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 	//	Ne 10
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.08150, 0.39000, 0.22500, 0.22500, 2.07000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 	//	Na 11
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.07450, 0.34200, 0.19200, 0.19200, 1.59000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 	//	Mg 12
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.06860, 0.30500, 0.16800, 0.16800, 1.24000, 0.91900, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 	//	Al 13
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.06350, 0.27500, 0.14900, 0.14900, 1.04000, 0.74400, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 	//	Si 14
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.05920, 0.25100, 0.13400, 0.13400, 0.89700, 0.63100, 0.63100, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 	//	P 15
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.05530, 0.23000, 0.12200, 0.12200, 0.79400, 0.55100, 0.55100, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 	//	S 16
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.05200, 0.21300, 0.11200, 0.11200, 0.71300, 0.49000, 0.49000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 	//	Cl 17
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.04900, 0.19800, 0.10400, 0.10400, 0.64900, 0.44200, 0.44200, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 	//	Ar 18
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.04640, 0.18500, 0.09650, 0.09650, 0.57000, 0.37600, 0.37600, 0.00000, 0.00000, 2.46000, 0.00000, 0.00000, 0.00000, 0.00000, 	//	K 19
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.04400, 0.17400, 0.09020, 0.09020, 0.50800, 0.33000, 0.33000, 0.00000, 0.00000, 1.95000, 0.00000, 0.00000, 0.00000, 0.00000, 	//	Ca 20
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.04180, 0.16400, 0.08470, 0.08470, 0.47100, 0.30400, 0.30400, 0.31000, 0.00000, 1.84000, 0.00000, 0.00000, 0.00000, 0.00000, 	//	Sc 21
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.03990, 0.15500, 0.07980, 0.07980, 0.44100, 0.28200, 0.28200, 0.00000, 0.27500, 1.76000, 0.00000, 0.00000, 0.00000, 0.00000, 	//	Ti 22
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.03810, 0.14700, 0.07550, 0.07550, 0.41500, 0.26400, 0.26400, 0.00000, 0.82500, 1.69000, 0.00000, 0.00000, 0.00000, 0.00000, 	//	V 23
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.03650, 0.14000, 0.07160, 0.07160, 0.39800, 0.25300, 0.25300, 0.25100, 0.25100, 1.85000, 0.00000, 0.00000, 0.00000, 0.00000, 	//	Cr 24
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.03500, 0.13300, 0.06810, 0.06810, 0.37300, 0.23500, 0.23500, 0.21500, 0.21500, 1.59000, 0.00000, 0.00000, 0.00000, 0.00000, 	//	Mn 25
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.03600, 0.12700, 0.06500, 0.06500, 0.35500, 0.22300, 0.22300, 0.21000, 0.21000, 1.55000, 0.00000, 0.00000, 0.00000, 0.00000, 	//	Fe 26
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.03230, 0.12200, 0.06210, 0.06210, 0.33900, 0.21200, 0.21200, 0.19000, 0.19000, 1.51000, 0.00000, 0.00000, 0.00000, 0.00000, 	//	Co 27
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.03110, 0.11700, 0.05950, 0.05950, 0.32500, 0.20200, 0.20200, 0.17900, 0.17900, 1.47000, 0.00000, 0.00000, 0.00000, 0.00000, 	//	Ni 28
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.03010, 0.11200, 0.05710, 0.05710, 0.31500, 0.19600, 0.19600, 0.18500, 0.18500, 1.64000, 0.00000, 0.00000, 0.00000, 0.00000, 	//	Cu 29
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.02900, 0.10800, 0.05490, 0.05490, 0.29900, 0.18600, 0.18600, 0.16200, 0.16200, 1.41000, 0.00000, 0.00000, 0.00000, 0.00000, 	//	Zn 30
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.02810, 0.10400, 0.05280, 0.05280, 0.28500, 0.17600, 0.17600, 0.14600, 0.14600, 1.18000, 0.00000, 0.91500, 0.00000, 0.00000, 	//	Ga 31
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.02720, 0.10100, 0.05090, 0.05090, 0.27100, 0.16700, 0.16700, 0.13300, 0.13300, 1.03000, 0.76900, 0.00000, 0.00000, 0.00000, 	//	Ge 32
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.02630, 0.09710, 0.04910, 0.04910, 0.25800, 0.15800, 0.15800, 0.12300, 0.12300, 0.92200, 0.67400, 0.67400, 0.00000, 0.00000, 	//	As 33
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.02550, 0.09390, 0.04750, 0.04750, 0.24700, 0.15100, 0.15100, 0.11400, 0.11400, 0.83900, 0.60400, 0.60400, 0.00000, 0.00000, 	//	Se 34
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.02480, 0.09090, 0.04590, 0.04590, 0.23600, 0.14400, 0.14400, 0.10700, 0.10700, 0.77200, 0.54900, 0.54900, 0.00000, 0.00000, 	//	Br 35
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.02340, 0.08620, 0.04330, 0.04430, 0.22200, 0.13400, 0.13700, 0.10100, 0.10100, 0.70500, 0.49600, 0.50800, 0.00000, 0.00000, 	//	Kr 36
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.02270, 0.08350, 0.04190, 0.04290, 0.21300, 0.12900, 0.13100, 0.09500, 0.09570, 0.63100, 0.43300, 0.44300, 0.00000, 0.00000, 	//	Rb 37
0.00000, 0.00000, 2.56000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.02210, 0.08090, 0.04060, 0.04160, 0.20400, 0.12300, 0.12500, 0.09010, 0.09070, 0.57300, 0.38900, 0.39700, 0.00000, 0.00000, 	//	Sr 38
0.00000, 0.00000, 2.06000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.02150, 0.07850, 0.03940, 0.04040, 0.19600, 0.11800, 0.12000, 0.08560, 0.08630, 0.53300, 0.36100, 0.36600, 0.45400, 0.00000, 	//	Y 39
0.00000, 0.00000, 1.90000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.02090, 0.07620, 0.03820, 0.03920, 0.18900, 0.11400, 0.11600, 0.08160, 0.08230, 0.50100, 0.33900, 0.34200, 0.34200, 0.00000, 	//	Zr 40
0.00000, 0.00000, 1.80000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.02030, 0.07410, 0.03710, 0.03810, 0.18200, 0.10900, 0.11200, 0.07790, 0.07870, 0.48000, 0.32400, 0.32600, 0.38500, 0.00000, 	//	Nb 41
0.00000, 0.00000, 1.88000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01980, 0.07200, 0.03600, 0.03710, 0.17600, 0.10500, 0.10800, 0.07470, 0.07540, 0.45600, 0.30600, 0.30900, 0.34700, 0.35700, 	//	Mo 42
0.00000, 0.00000, 1.82000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01930, 0.07000, 0.03500, 0.03610, 0.17000, 0.10200, 0.10400, 0.07170, 0.07240, 0.43000, 0.28700, 0.29100, 0.29900, 0.30500, 	//	Tc 43
0.00000, 0.00000, 1.62000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01880, 0.06820, 0.03410, 0.03520, 0.16400, 0.09820, 0.10100, 0.06900, 0.06970, 0.41600, 0.27600, 0.28100, 0.29600, 0.30300, 	//	Ru 44
0.00000, 0.00000, 1.74000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01840, 0.06640, 0.03320, 0.03430, 0.15900, 0.09490, 0.09730, 0.06650, 0.06710, 0.39900, 0.26300, 0.27000, 0.27700, 0.28300, 	//	Rh 45
0.00000, 0.00000, 1.70000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01790, 0.06470, 0.03230, 0.03350, 0.15400, 0.09190, 0.09430, 0.06420, 0.06470, 0.38600, 0.25400, 0.26100, 0.27800, 0.28500, 	//	Pd 46
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01750, 0.06310, 0.03150, 0.03270, 0.14900, 0.08910, 0.09150, 0.06200, 0.06250, 0.36900, 0.24200, 0.24900, 0.24700, 0.25100, 	//	Ag 47
0.00000, 0.00000, 1.64000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01710, 0.06150, 0.03070, 0.03190, 0.14500, 0.08640, 0.08880, 0.06000, 0.06050, 0.35200, 0.23000, 0.23700, 0.22400, 0.22800, 	//	Cd 48
0.00000, 0.00000, 1.45000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01670, 0.06000, 0.02990, 0.03110, 0.14100, 0.08380, 0.08630, 0.05810, 0.05860, 0.33600, 0.21900, 0.22600, 0.20600, 0.20800, 	//	In 49
0.00000, 0.00000, 1.24000, 0.96200, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01630, 0.05860, 0.02920, 0.03040, 0.13700, 0.08140, 0.08390, 0.05630, 0.05690, 0.32200, 0.20900, 0.21600, 0.19200, 0.19300, 	//	Sn 50
0.00000, 0.00000, 1.10000, 0.83600, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01590, 0.05720, 0.02850, 0.02980, 0.13300, 0.07910, 0.08160, 0.05460, 0.05520, 0.30900, 0.20000, 0.20700, 0.17900, 0.18100, 	//	Sb 51
0.00000, 0.00000, 1.00000, 0.73700, 0.76400, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01560, 0.05590, 0.02780, 0.02910, 0.13000, 0.07690, 0.07940, 0.05300, 0.05360, 0.29600, 0.19100, 0.19800, 0.16800, 0.17000, 	//	Te 52
0.00000, 0.00000, 0.92000, 0.66500, 0.69300, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01530, 0.05470, 0.02720, 0.02850, 0.12600, 0.07480, 0.07740, 0.05150, 0.05210, 0.28500, 0.18400, 0.19000, 0.15900, 0.16100, 	//	I 53
0.00000, 0.00000, 0.85300, 0.60800, 0.63700, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01490, 0.05350, 0.02650, 0.02790, 0.12300, 0.07290, 0.07540, 0.05010, 0.05070, 0.27400, 0.17600, 0.18200, 0.15100, 0.15300, 	//	Xe 54
0.00000, 0.00000, 0.79700, 0.56200, 0.59200, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01460, 0.05230, 0.02590, 0.02730, 0.12000, 0.07100, 0.07350, 0.04880, 0.04940, 0.26400, 0.17000, 0.17600, 0.14400, 0.14500, 	//	Cs 55
0.00000, 0.00000, 0.72000, 0.49800, 0.52200, 0.00000, 0.00000, 0.00000, 0.00000, 2.74000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01430, 0.05120, 0.02540, 0.02680, 0.11700, 0.06920, 0.07180, 0.04750, 0.04810, 0.25500, 0.16300, 0.16900, 0.13700, 0.13900, 	//	Ba 56
0.00000, 0.00000, 0.65700, 0.45200, 0.47300, 0.00000, 0.00000, 0.00000, 0.00000, 2.23000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01400, 0.05010, 0.02480, 0.02630, 0.11400, 0.06750, 0.07010, 0.04630, 0.04690, 0.24600, 0.15800, 0.16300, 0.13100, 0.13300, 	//	La 57
0.00000, 0.00000, 0.61600, 0.42400, 0.44000, 0.50900, 0.00000, 0.00000, 0.00000, 2.07000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01370, 0.04900, 0.02430, 0.02570, 0.11100, 0.06580, 0.06840, 0.04520, 0.04580, 0.23900, 0.15300, 0.15900, 0.12700, 0.12900, 	//	Ce 58
0.14200, 0.00000, 0.60200, 0.41400, 0.43100, 0.49800, 0.00000, 0.00000, 0.00000, 2.05000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01350, 0.04800, 0.02380, 0.02520, 0.10900, 0.06430, 0.06690, 0.04410, 0.04470, 0.23400, 0.15000, 0.15600, 0.12500, 0.12600, 	//	Pr 59
0.14900, 0.00000, 0.60900, 0.41700, 0.44200, 0.00000, 0.00000, 0.00000, 0.00000, 2.15000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01320, 0.04700, 0.02330, 0.02480, 0.10600, 0.06280, 0.06540, 0.04310, 0.04360, 0.22800, 0.14600, 0.15200, 0.12100, 0.12300, 	//	Nd 60
0.14100, 0.00000, 0.59700, 0.40700, 0.43400, 0.00000, 0.00000, 0.00000, 0.00000, 2.13000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01290, 0.04610, 0.02280, 0.02430, 0.10400, 0.06130, 0.06400, 0.04210, 0.04260, 0.22300, 0.14200, 0.14800, 0.11800, 0.11900, 	//	Pm 61
0.13400, 0.13400, 0.58600, 0.39900, 0.42700, 0.00000, 0.00000, 0.00000, 0.00000, 2.11000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01270, 0.04520, 0.02230, 0.02390, 0.10100, 0.05990, 0.06260, 0.04120, 0.04170, 0.21700, 0.13900, 0.14500, 0.11500, 0.11600, 	//	Sm 62
0.12900, 0.12900, 0.57500, 0.39100, 0.42000, 0.00000, 0.00000, 0.00000, 0.00000, 2.09000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01240, 0.04430, 0.02190, 0.02340, 0.09920, 0.05860, 0.06130, 0.04030, 0.04080, 0.21200, 0.13500, 0.14200, 0.11200, 0.11400, 	//	Eu 63
0.12400, 0.13000, 0.56400, 0.38400, 0.41300, 0.00000, 0.00000, 0.00000, 0.00000, 2.07000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01220, 0.04340, 0.02140, 0.02300, 0.09710, 0.05730, 0.06010, 0.03940, 0.03990, 0.20700, 0.13200, 0.13800, 0.10900, 0.11000, 	//	Gd 64
0.11300, 0.11600, 0.53800, 0.36600, 0.39100, 0.46800, 0.00000, 0.00000, 0.00000, 0.19300, 0.00000, 0.00000, 0.00000, 0.00000,
0.01200, 0.04260, 0.02100, 0.02260, 0.09500, 0.05610, 0.05890, 0.03860, 0.03910, 0.20300, 0.12900, 0.13600, 0.10700, 0.10800, 	//	Tb 65
0.11600, 0.12100, 0.54400, 0.37000, 0.40000, 0.00000, 0.00000, 0.00000, 0.00000, 2.03000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01170, 0.04180, 0.02060, 0.02220, 0.09310, 0.05490, 0.05770, 0.03780, 0.03840, 0.19900, 0.12600, 0.13300, 0.10400, 0.10600, 	//	Dy 66
0.11300, 0.11700, 0.53500, 0.36400, 0.39400, 0.00000, 0.00000, 0.00000, 0.00000, 2.01000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01150, 0.04100, 0.02020, 0.02190, 0.09120, 0.05370, 0.05660, 0.03700, 0.03760, 0.19400, 0.12300, 0.13000, 0.10200, 0.10400, 	//	Ho 67
0.11000, 0.11400, 0.52600, 0.35800, 0.38800, 0.00000, 0.00000, 0.00000, 0.00000, 1.99000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01130, 0.04020, 0.01980, 0.02150, 0.08930, 0.05260, 0.05560, 0.03620, 0.03690, 0.19000, 0.12100, 0.12800, 0.09960, 0.10100, 	//	Er 68
0.10700, 0.11000, 0.51700, 0.35300, 0.38300, 0.00000, 0.00000, 0.00000, 0.00000, 1.97000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01110, 0.03950, 0.01940, 0.02110, 0.08760, 0.05160, 0.05450, 0.03550, 0.03620, 0.18600, 0.11800, 0.12500, 0.09750, 0.09940, 	//	Tm 69
0.10400, 0.10700, 0.50900, 0.34700, 0.37800, 0.00000, 0.00000, 0.00000, 0.00000, 1.95000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01090, 0.03880, 0.01900, 0.02080, 0.08590, 0.05050, 0.05350, 0.03480, 0.03550, 0.18300, 0.11600, 0.12300, 0.09540, 0.09740, 	//	Yb 70
0.10200, 0.10500, 0.50100, 0.34200, 0.37300, 0.00000, 0.00000, 0.00000, 0.00000, 1.94000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01070, 0.03810, 0.01870, 0.02050, 0.08420, 0.04950, 0.05260, 0.03420, 0.03490, 0.17900, 0.11300, 0.12000, 0.09310, 0.09510, 	//	Lu 71
0.09530, 0.09670, 0.48000, 0.32700, 0.35400, 0.48900, 0.00000, 0.00000, 0.00000, 1.77000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01050, 0.03740, 0.01830, 0.02010, 0.08260, 0.04860, 0.05170, 0.03360, 0.03430, 0.17500, 0.11000, 0.11800, 0.09080, 0.09270, 	//	Hf 72
0.08980, 0.09070, 0.46100, 0.31400, 0.33800, 0.42300, 0.00000, 0.00000, 0.00000, 1.68000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01030, 0.03670, 0.01800, 0.01980, 0.08110, 0.04760, 0.05080, 0.03300, 0.03370, 0.17100, 0.10800, 0.11500, 0.08850, 0.09050, 	//	Ta 73
0.08530, 0.08580, 0.44300, 0.30200, 0.32400, 0.38200, 0.00000, 0.00000, 0.00000, 1.61000, 0.00000, 0.00000, 0.00000, 0.00000,
0.01010, 0.03610, 0.01770, 0.01950, 0.07960, 0.04670, 0.04990, 0.03240, 0.03310, 0.16700, 0.10500, 0.11300, 0.08630, 0.08830, 	//	W 74
0.08130, 0.08170, 0.42700, 0.29100, 0.31100, 0.35200, 0.00000, 0.00000, 0.00000, 1.55000, 0.00000, 0.00000, 0.00000, 0.00000,
0.00997, 0.03540, 0.01730, 0.01920, 0.07810, 0.04580, 0.04910, 0.03180, 0.03250, 0.16300, 0.10300, 0.11000, 0.08420, 0.08620, 	//	Re 75
0.07770, 0.07820, 0.41200, 0.27900, 0.30000, 0.32700, 0.34600, 0.00000, 0.00000, 1.49000, 0.00000, 0.00000, 0.00000, 0.00000,
0.00980, 0.03480, 0.01700, 0.01890, 0.07670, 0.04500, 0.04820, 0.03120, 0.03200, 0.16000, 0.10100, 0.10800, 0.08220, 0.08420, 	//	Os 76
0.07440, 0.07510, 0.39700, 0.26800, 0.29000, 0.30700, 0.34200, 0.00000, 0.00000, 1.44000, 0.00000, 0.00000, 0.00000, 0.00000,
0.00962, 0.03420, 0.01670, 0.01870, 0.07530, 0.04420, 0.04750, 0.03070, 0.03140, 0.15600, 0.09860, 0.10600, 0.08030, 0.08220, 	//	Ir 77
0.07160, 0.07230, 0.38400, 0.25800, 0.28100, 0.29100, 0.30600, 0.00000, 0.00000, 1.40000, 0.00000, 0.00000, 0.00000, 0.00000,
0.00946, 0.03360, 0.01640, 0.01840, 0.07400, 0.04340, 0.04670, 0.03020, 0.03090, 0.15300, 0.09640, 0.10400, 0.07840, 0.08030, 	//	Pt 78
0.06900, 0.06980, 0.37500, 0.25000, 0.27400, 0.28700, 0.30400, 0.00000, 0.00000, 1.45000, 0.00000, 0.00000, 0.00000, 0.00000,
0.00929, 0.03310, 0.01610, 0.01810, 0.07270, 0.04260, 0.04590, 0.02970, 0.03040, 0.15000, 0.09430, 0.10200, 0.07670, 0.07850, 	//	Au 79
0.06660, 0.06750, 0.36300, 0.24200, 0.26600, 0.27300, 0.28800, 0.00000, 0.00000, 1.41000, 0.00000, 0.00000, 0.00000, 0.00000,
0.00913, 0.03250, 0.01580, 0.01790, 0.07140, 0.04180, 0.04520, 0.02920, 0.03000, 0.14700, 0.09230, 0.09950, 0.07500, 0.07680, 	//	Hg 80
0.06450, 0.06530, 0.35000, 0.23300, 0.25600, 0.25200, 0.26500, 0.00000, 0.00000, 1.29000, 0.00000, 0.00000, 0.00000, 0.00000,
0.00898, 0.03200, 0.01550, 0.01760, 0.07020, 0.04110, 0.04450, 0.02870, 0.02950, 0.14400, 0.09030, 0.09760, 0.07330, 0.07510, 	//	Tl 81
0.06250, 0.06330, 0.33700, 0.22400, 0.24600, 0.23500, 0.24400, 0.00000, 0.00000, 1.15000, 0.91300, 0.00000, 0.00000, 0.00000,
0.00882, 0.03140, 0.01530, 0.01740, 0.06890, 0.04030, 0.04380, 0.02830, 0.02900, 0.14100, 0.08840, 0.09570, 0.07180, 0.07350, 	//	Pb 82
0.06060, 0.06140, 0.32500, 0.21500, 0.23700, 0.22100, 0.22700, 0.00000, 0.00000, 1.04000, 0.80800, 0.00000, 0.00000, 0.00000,
0.00867, 0.03090, 0.01500, 0.01710, 0.06780, 0.03960, 0.04320, 0.02780, 0.02860, 0.13800, 0.08660, 0.09390, 0.07030, 0.07200, 	//	Bi 83
0.05890, 0.05970, 0.31400, 0.20700, 0.22800, 0.20800, 0.21400, 0.00000, 0.00000, 0.95500, 0.71300, 0.83200, 0.00000, 0.00000,
0.00852, 0.03040, 0.01470, 0.01690, 0.06660, 0.03900, 0.04250, 0.02740, 0.02820, 0.13500, 0.08480, 0.09210, 0.06880, 0.07060, 	//	Po 84
0.05730, 0.05800, 0.30300, 0.20000, 0.22000, 0.19700, 0.20300, 0.00000, 0.00000, 0.88700, 0.64500, 0.75300, 0.00000, 0.00000,
0.00838, 0.02990, 0.01450, 0.01670, 0.06550, 0.03830, 0.04190, 0.02700, 0.02780, 0.13200, 0.08310, 0.09040, 0.06740, 0.06910, 	//	At 85
0.05580, 0.05650, 0.29300, 0.19300, 0.21200, 0.18800, 0.19400, 0.00000, 0.00000, 0.83000, 0.59300, 0.69200, 0.00000, 0.00000,
0.00823, 0.02940, 0.01420, 0.01650, 0.06440, 0.03760, 0.04130, 0.02660, 0.02740, 0.13000, 0.08140, 0.08870, 0.06610, 0.06780, 	//	Rn 86
0.05430, 0.05500, 0.28400, 0.18700, 0.20500, 0.17900, 0.18500, 0.00000, 0.00000, 0.78200, 0.55100, 0.64400, 0.00000, 0.00000,
0.00809, 0.02890, 0.01400, 0.01620, 0.06340, 0.03700, 0.04070, 0.02620, 0.02700, 0.12700, 0.07980, 0.08710, 0.06480, 0.06650, 	//	Fr 87
0.05300, 0.05370, 0.27500, 0.18100, 0.19900, 0.17200, 0.17700, 0.00000, 0.00000, 0.71400, 0.49600, 0.57100, 0.00000, 2.65000,
0.00796, 0.02850, 0.01370, 0.01600, 0.06230, 0.03640, 0.04010, 0.02580, 0.02660, 0.12500, 0.07820, 0.08560, 0.06350, 0.06520, 	//	Ra 88
0.05170, 0.05240, 0.26600, 0.17500, 0.19200, 0.16500, 0.17000, 0.00000, 0.00000, 0.65700, 0.45600, 0.52000, 0.00000, 0.26600,
0.00782, 0.02800, 0.01350, 0.01580, 0.06130, 0.03580, 0.03950, 0.02540, 0.02620, 0.12200, 0.07670, 0.08410, 0.06230, 0.06400, 	//	Ac 89
0.05050, 0.05120, 0.25800, 0.16900, 0.18700, 0.15800, 0.16300, 0.00000, 0.00000, 0.61700, 0.42900, 0.48400, 0.60800, 1.99000,
0.00769, 0.02750, 0.01330, 0.01560, 0.06030, 0.03520, 0.03900, 0.02500, 0.02590, 0.12000, 0.07520, 0.08260, 0.06110, 0.06280, 	//	Th 90
0.04930, 0.05000, 0.25100, 0.16400, 0.18100, 0.15200, 0.15700, 0.00000, 0.00000, 0.58500, 0.40800, 0.45600, 0.53500, 1.86000,
0.00756, 0.02710, 0.01300, 0.01540, 0.05930, 0.03460, 0.03840, 0.02470, 0.02550, 0.11800, 0.07380, 0.08120, 0.06000, 0.06170, 	//	Pa 91
0.04820, 0.04890, 0.24500, 0.16100, 0.17700, 0.14900, 0.15300, 0.20700, 0.00000, 0.58000, 0.40100, 0.45900, 0.55800, 1.93000,
0.00743, 0.02670, 0.01280, 0.01520, 0.05840, 0.03400, 0.03790, 0.02430, 0.02520, 0.11500, 0.07240, 0.07980, 0.05890, 0.06060, 	//	U 92
0.04720, 0.04790, 0.23800, 0.15700, 0.17300, 0.14500, 0.14900, 0.19300, 0.00000, 0.56500, 0.39000, 0.45000, 0.54300, 1.91000,
0.00730, 0.02620, 0.01260, 0.01500, 0.05750, 0.03340, 0.03740, 0.02400, 0.02480, 0.11300, 0.07100, 0.07850, 0.05790, 0.05960, 	//	Np 93
0.04620, 0.04690, 0.23300, 0.15300, 0.16900, 0.14200, 0.14500, 0.18200, 0.00000, 0.55200, 0.38000, 0.44200, 0.53200, 1.89000,
0.00718, 0.02580, 0.01240, 0.01490, 0.05650, 0.03290, 0.03690, 0.02370, 0.02450, 0.11100, 0.06970, 0.07720, 0.05690, 0.05850, 	//	Pu 94
0.04520, 0.04600, 0.22800, 0.15000, 0.16600, 0.13900, 0.14200, 0.18300, 0.00000, 0.55000, 0.37600, 0.44900, 0.00000, 1.99000,
0.00705, 0.02540, 0.01210, 0.01470, 0.05570, 0.03240, 0.03640, 0.02330, 0.02420, 0.10900, 0.06840, 0.07600, 0.05590, 0.05760, 	//	Am 95
0.04430, 0.04500, 0.22200, 0.14600, 0.16200, 0.13600, 0.13900, 0.17400, 0.18700, 0.53800, 0.36700, 0.44100, 0.00000, 1.97000,
0.00693, 0.02500, 0.01190, 0.01450, 0.05480, 0.03180, 0.03590, 0.02300, 0.02390, 0.10700, 0.06720, 0.07480, 0.05490, 0.05660, 	//	Cm 96
0.04350, 0.04420, 0.21700, 0.14200, 0.15900, 0.13200, 0.13600, 0.15900, 0.16600, 0.51600, 0.35300, 0.42100, 0.51600, 1.83000,
0.00681, 0.02460, 0.01170, 0.01430, 0.05390, 0.03130, 0.03550, 0.02270, 0.02360, 0.10500, 0.06600, 0.07360, 0.05400, 0.05570, 	//	Bk 97
0.04270, 0.04340, 0.21300, 0.13900, 0.15600, 0.13000, 0.13300, 0.16000, 0.16900, 0.51500, 0.35100, 0.42700, 0.00000, 1.93000,
0.00670, 0.02420, 0.01150, 0.01410, 0.05310, 0.03080, 0.03500, 0.02240, 0.02330, 0.10400, 0.06480, 0.07250, 0.05310, 0.05480, 	//	Cf 98
0.04190, 0.04260, 0.20800, 0.13600, 0.15300, 0.12700, 0.13100, 0.15400, 0.16200, 0.50500, 0.34400, 0.42100, 0.00000, 1.91000,
0.00658, 0.02380, 0.01130, 0.01400, 0.05230, 0.03030, 0.03460, 0.02210, 0.02300, 0.10200, 0.06360, 0.07140, 0.05230, 0.05400, 	//	Es 99
0.04110, 0.04180, 0.20400, 0.13300, 0.15000, 0.12400, 0.12800, 0.14900, 0.15600, 0.49500, 0.33700, 0.41500, 0.00000, 1.89000,
0.00646, 0.02340, 0.01110, 0.01380, 0.05150, 0.02980, 0.03410, 0.02180, 0.02270, 0.09990, 0.06250, 0.07030, 0.05150, 0.05320, 	//	Fm 100
0.04040, 0.04110, 0.19900, 0.13000, 0.14800, 0.12100, 0.12600, 0.14400, 0.15100, 0.48500, 0.33100, 0.41000, 0.00000, 1.87000,
		};


