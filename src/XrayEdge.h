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

#ifndef XrayEdge_h
#define XrayEdge_h

#include <vector>
#include <string>
#include "Element.h"

//	Modified Oct. 1, 2015 to make maxz and maxIndex functions static
//  Modified Feb. 11, 2017 to include an EdgeLevel that matches nothing (for spectrum components)
//` Modified March 31, 2017 to make EdgeIndex and EdgeLevel more consistent
//      Change first EdgeLevel to K and first EdgeIndex to K1

//	enum for edge indices
 enum EdgeIndex {
	K1=0, L1, L2, L3, M1, M2,   //  Use K1 for K edge to avoid conflict with K level below
	M3, M4, M5, N1, N2, N3,
	N4, N5, N6, N7, O1, O2,
	O3, O4, O5, O6, O7, P1, P2, P3, P4, Q1 };
//	enums for principal quantum number (energy level) and angular momentum
enum EdgeLevel { NO_EDGE = -1, K, L, M, N, O, P, Q};   //  Modified to use for spectrum components
enum EdgeAngularMonmentum {s, p, d, f};
class XrayEdge {
friend class XrayLines;
public:
	XrayEdge ( const Element e, const EdgeIndex indexIn = K1 );
//	XrayEdge ( const Element e, const std::string s );
//	XrayEdge ( const Element e, const float energy = 1.0e99 );
	XrayEdge ( const XrayEdge& ee );
//		must have default constructor to declare arrays
	XrayEdge() { edgeIndex=K1; }
	bool operator==(const XrayEdge&) const;
	bool operator<(const XrayEdge& comp) const;
	bool operator>(const XrayEdge& comp) const;
	const float energy() const;
	const std::string &name() const { return EDGE_NAMES[edgeIndex]; };
	const std::string &symbol() const { return EDGE_NAMES[edgeIndex]; };
	const std::string &designation() const { return EDGE_DESIGNATIONS[edgeIndex]; };
	const float degeneracy() const;
	const float occupancy() const;
	const EdgeLevel level() const;
	const EdgeAngularMonmentum angularMomentum() const;
	const float spin() const;
//		fluorescence yield
	const float yield() const { return yield_value ( z.Z(), edgeIndex ); };
//		absorption edge jump ratio (also called edge step)
	const float jump() const { return jump_value ( z.Z(), edgeIndex ); };
	const float step() const { return jump_value ( z.Z(), edgeIndex ); };
//		Coster-Kronig transition rate from this edge to given edge
	const float ck ( const XrayEdge toEdge) const { return z == toEdge.element() ?
		ck_value ( z.Z(), edgeIndex, toEdge.index() ) : 0.0; } ;
//		total Coster-Kronig transition rate from this edge to given edge including intermediate levels
	const float cktotal ( const XrayEdge toEdge ) const;
//		Compton profile parameter Jzero (J of p=0)
	const float jzero ( ) const;
//		Level width data from Campbell ATNDT_77_1_2001
	const float width ( ) const;
//	these are static functions so they can be used before any objects are created
//		returns number of occupied edges
	static const int numberOccupied ( std::vector <EdgeIndex> &edgeList,  const Element&  el ) ;
//		returns number of occupied edges with non-zero energy < excitation energy (0.0 => all)
	static const int numberOfEdges ( std::vector <EdgeIndex> &edgeList,  const Element&  el, const float excitEnergy = 0.0 ) ;
//		private data access functions
	int number() const;
	const Element&  element() const { return z; };
	const EdgeIndex index() const { return edgeIndex; };
	static const int maxz() { return MAXZ; };
	static const int maxIndex() { return MAXINDEX; };

    std::string toString() const;

private:
	Element z;
//		many routines here and in XrayLines assume the indices are in a particular order
	EdgeIndex edgeIndex;
	static const int MAXINDEX = 27;
	static const int MAXZ = 100;
	static const std::string EDGE_NAMES[];
	static const std::string EDGE_DESIGNATIONS[];
//		calculate for index because compiler doesn't like static const with 2 indices
	static const float EDGE_ENERGIES[];
	static const int EDGE_OCCUPANCIES[];
	static const float EDGE_JZEROS[];
	const float yield_value ( const int zi, const int fromIndex ) const;
	const float jump_value ( const int zi, const int fromIndex ) const;
	const float ck_value ( const int zi, const int fromIndex, const int toIndex ) const;
};
#endif
