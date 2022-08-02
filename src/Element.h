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

#ifndef Element_h
#define Element_h

#include <string>

//	Elam Ravel Sieber database class

class Element {

public:
	Element(const int z);
	Element();
	Element(const std::string s);
	Element(const Element& e);
	float atomicWeight() const;
	float density() const;
	bool operator==(const Element& e) const { return atomicNumberZ == e.atomicNumber(); };
	bool operator<(const Element& e) const { return atomicNumberZ < e.atomicNumber(); };
	bool operator>(const Element& e) const { return atomicNumberZ > e.atomicNumber(); };
	int Z() const { return atomicNumberZ; };
	int atomicNumber() const { return atomicNumberZ; };
	const std::string& symbol() const;
	static const int maxZ() { return MAXZ; };
	static const bool check_Z ( const int z_in ) { return ( 1 <= z_in && z_in <= MAXZ ); };
	static const bool check_symbol ( const std::string &symbol_in );

	std::string toString() const;

private:
	int atomicNumberZ;
	static const int MAXZ;
	static const std::string SYMBOLS[];
};

#endif
