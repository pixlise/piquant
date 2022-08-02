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

#ifndef ScatterXsectTable_h
#define ScatterXsectTable_h

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <math.h>
#include "Element.h"

using namespace std;

//	Angle-dependent form factors for coherent and incoherent scattering - database class

class ScatterXsectTable {
//	contains and evaluates a table of angle-dependent form factors for scattering cross sections 
//		for a particular element
//		member functions return cross section values at a single energy and angle
//		class contains entire table for all energies
//		data read from database file (currently hard-wired in implementation file)
//	As of Jan. 30, 2003 data sources are:
//		Incoherent scattering function S(x): J. H. Hubbell, Wm. J. Vieglie, E. A. Briggs, R. T. Brown, 
//			D. T. Cromer, and R. J. Howerton, J. Chem. Phys. Ref. Data, Vol. 4, No. 3, 1975, pp471-538.  
//			Erratum Vol. 6, No. 2, 1977, pp615-616.
//		Coherent atomic structure factors (form factors) F(x):  RTAB: The Rayleigh Scattering Database, 
//			Lynn Kissel, LLNL, Version 2.1, 2000-09-29.  See L. Kissell, Rad. Phys. Chem. Vol. 59, 
//			No. 30, 2000, pp185-200.
//	Calculations of the angle-dependent differential cross sections follow R. Tertian and F. Claisse, 
//			"Principles of Quantitative X-ray Fluorescence Analysis" (Heyden and Sons, 
//			London), 1982, LC 545.836, Dewey QC482.S6, ISBN 0-85501-709-0 (pp 26-31).
//	Calculations of doubly-differential Compton cross sections (vs scatter angle and scatered photon energy)
//			are calcualted via the relativistic impulse approximation (RIA) of Ribberfors and Berggren.
//			Compton profiles are the analytic profiles of Brusa et. al., with Jzero valures from Biggs, Mandlesohn,and Mann.
//		Gudrun Alm Carlsson, Carl. A. Carlsson, karl-Frederik Berggren, and Roland Ribberfors, 
//			Med. Phys. 9(6), Nov./Dec. 1982, 868-879. (This paper provides a good, brief summary of the RIA.)
//		D. Brusa, G. Stutz, J. A. Riveros, J. M. Fernandez-Varea, and F. Salvat, Nuclear Instruments and 
//			Methods in Physics Research A 379 (1996) 167-175.  (Analytic approximation to Compton profiles 
//			plus the full set of equations for the RIA with correct units for calculating the cross section 
//			in barns per atom.)
//		F.Biggs, L.Mendelsohn, and J.Mann, At. Data and Nucl. Data Tables Vol. 16, pp201-309 (1975).  (The 
//			values for calculated Compton profiles at zero momentum, used as a parameter in the analytic profiles
//			of Brusa et. al.  These are stored in the XrayEdge class, file XrayEdge.cpp.)




public:
	ScatterXsectTable(const Element& el);
	ScatterXsectTable();
//		equality operators
	bool operator == ( const ScatterXsectTable& test ) const {
		return ( thisElement == test.element() ); };
	const Element& element() const { return thisElement; };
	float coherent(const float energy_in, const float angle_in) const { return cohCalc(energy_in, angle_in); };
	float Rayleigh(const float energy_in, const float angle_in) const { return cohCalc(energy_in, angle_in); };
//			single-diferential Compton cross sections (vs incident energy and angle)
	float incoherent(const float energy, const float angle_in) const { return incohCalc(energy, angle_in); };
	float Compton(const float energy_in, const float angle_in) const { return incohCalc(energy_in, angle_in); };
//			doubly-differential Compton cross sections (vs incident energy, angle, and scattered energy)
	float incoherent(const float energy, const float angle_in, const float ePrime_in ) const { 
		return incohCalc(energy, angle_in, ePrime_in); };
	float Compton(const float energy_in, const float angle_in, const float ePrime_in ) const { 
		return incohCalc(energy_in, angle_in, ePrime_in); };
	static float eCompton ( const float energy_in, const float angle_in ) {
			return eC ( energy_in, angle_in ); };
	static float eComptonUp ( const float energy_in, const float angle_in );
	float shift ( const float energy_in, const float angle_in ) {
		return eC ( energy_in, angle_in ) - energy_in; };
	int cohCount() const { return numberCoherent; };
	int incohCount() const { return numberIncoherent; };
	float coherentX(const int index) const { return xCoh.at(index); };
	float coherentValue(const int index) const { return fofx.at(index); };
	float incoherentX(const int index) const { return xIncoh.at(index); };
	float incoherentValue(const int index) const { return sofx.at(index); };

private:
	float cohCalc(const float energy_in, const float angle_in) const;
	static float eC (const float energy_in, const float angle_in);
	float incohCalc(const float energy_in, const float angle_in) const;	
	float incohCalc(const float energy_in, const float angle_in, const float ePrime_in) const;	

//individual data
	Element thisElement;
	int numberCoherent;
	vector <float> xCoh;
	vector <float> fofx;
	int numberIncoherent;
	vector <float> xIncoh;
	vector <float> sofx;
};
#endif
