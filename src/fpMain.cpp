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

#include <sstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include "Element.h"
#include "XrayEdge.h"
#include "XrayLines.h"
#include "XraySource.h"
#include "XrayXsectTable.h"
#include "ScatterXsectTable.h"
#include "XRFconstants.h"
#include "sampleEdgeList.h"
#include "fpExcitation.h"
#include "fpEdgeAbsorption.h"
#include "fpPrimary.h"
#include "fpSecondary.h"
#include "fpCK.h"
#include "fpMain.h"
#include "fpBeams.h"
#include "XRFcontrols.h"

// NY edit 9-9-2013
#include "XrayOptic.h"
#include "toStringHelpers.h"

//	added X-ray tube current     Nov. 30, 2011
//	changed contEn[0] < 0 to contEn[0] <= 0 in fpContScat   Dec. 13, 2011
//	also check for energy <= 0 IN fpLineScat

//  Modified, Sep. 30, 2013    W. T. Elam and Nick Yang
//  To include effects of incident beam optic on lines and continuum

//  Modified Dec. 18, 2015     W. T. Elam      fpCalc
//  Fix several bugs with secondary fluorescence calculation
//  One was serious: factor in XrayLines was changed for element producing secondary excitation, but element emitting it
//  Others were minor, skipping secondary fluorescence calculation in some cases when it should not have been skipped
//  Change AmpTekDet to XrayDetector class   May 11, 2016
//  Modified Dec. 22, 2016 to include dust on optic (or anywhere in incident beam path) and dust on specimen (treat same as window)
//  Re-written Feb. 2, 2017
//      Use XrayMaterial class for specimen composition, thickness, and X-ray parameters
//      Use new conditions structure and setup for fp calculations
//      Use XrayEnergyCal class (in XraySpectrum files) for energy to channel conversions
//      Get number of channels from output vector (it must be properly sized before call)
//  Modified Feb. 10, 2017
//      Put FPstorage in this file at file scope
//  Modified July 25, 2018
//      Write out some useful information if calculated intensity is zero or nan
//  Modified May 25, 2019
//      Temporary fix for dependence of calculated intensity on element list
//          Include all elements in list of edge energies and intensities for excitation integrals
//      Fix wrong index in secEdgeIndex loop in fpCalc, was incorrectly changed on Dec. 18, 2015
//  Modified Nov. 24, 2020
//      Add pure element X-rayLines in fpPrep (so that matrix effect factor can be calculated)
//      Add matrix effect factor to sample XrayLines in fpCalc
//  Modified Jan. 7, 2021
//      Implement SEC_FLUOR_THRESHOLD from XRFcontrols.h in fpCalc (and re-arrange sec fluor criteria)


using namespace std;

// performs fundamental parameters calculation of x-ray fluorescence line intensities
//		using formulas developed by Sherman, Gillam and Heal, and Shiraiwa and Fujino.
//		For details and references see R. Tertian and F. Claisse, "Principles of Quantitative
//			X-ray Fluorescence Analysis" (Heyden and Sons, London), 1982, LC 545.836,
//			Dewey QC482.S6, ISBN 0-85501-709-0.
//     Copyright 2006  W. T. Elam


string FPstorage_toString(const FPstorage &storage)
{
	ostringstream os;
	os << "FPstorage:" << endl;
	os << "  sampleElements:" << endl;
	int c = 0;
	for(auto it = storage.sampleElements.begin(); it != storage.sampleElements.end(); it++)
	{
		os << "  [" << c << "]: " << it->toString() << endl;
		c++;
	}

	c = 0;
	for(auto it = storage.sampleEdges.begin(); it != storage.sampleEdges.end(); it++)
	{
		os << "  [" << c << "]: " << it->toString() << endl;
		c++;
	}

	c = 0;
	for(auto it = storage.elementIndices.begin(); it != storage.elementIndices.end(); it++)
	{
		os << "  [" << c << "]: " << *it << endl;
		c++;
	}

	os << "  excitEnergies: " << floatVecToString(storage.excitEnergies) << endl;
	os << "  excitIntensities: " << floatVecToString(storage.excitIntensities) << endl;

	os << "  sinExcit: " << storage.sinExcit << endl;
	os << "  sinEmerg: " << storage.sinEmerg << endl;
	os << "  geometry: " << storage.geometry << endl;

	c = 0;
	for(auto it = storage.pureLines.begin(); it != storage.pureLines.end(); it++)
	{
		os << "  [" << c << "]: " << it->toString() << endl;
		c++;
	}

	return os.str();
}

void fpIDlist ( const float detRes, const float eMin, const float eMax,
			   vector <float> &energies, vector <Element> &elements ) {

	energies.resize(0);
	elements.resize(0);
	if ( detRes <= 0 ) return;
	if ( eMax <= 0 ) return;
	if ( eMax <= eMin ) return;
	float fwhm2 = detRes / 2;
//		loop over all elements
//	int maximumZ = Element::maxZ();
	int maximumZ = 94;	//	end of cross-section tables
	int z;
	for ( z=1; z<=maximumZ; z++ ) {
//			eliminate some elements that are rare and cause problems
		if ( z == 49 ) continue;
		Element el(z);
//			generate list of absorption edges between eMin and eMax
		vector <EdgeIndex> edgeIndexList;
		XrayEdge::numberOfEdges ( edgeIndexList, el, eMax );
		int edgeIndex;
		for ( edgeIndex=0; edgeIndex<edgeIndexList.size(); edgeIndex++ ) {
			XrayEdge edge ( el, edgeIndexList[edgeIndex] );
			if ( eMin > 0 && edge.energy() < eMin ) continue;
//				examine all lines emitted from a vacancy in this edge
			XrayLines lines ( edge );
			if ( lines.numberOfLines() <= 0 ) continue;
//				combine emission lines within FWHM/2 of most intense line and iterate
//				find the most intense line
			float maxInt = 0;
			int maxIndex = 0;
			int j;
			for ( j=0; j<lines.numberOfLines(); j++ ) {
				if ( lines.relative(j) > maxInt ) {
					maxInt = lines.relative(j);
					maxIndex = j;
				};
			};
			float maxEnergy = lines.energy( maxIndex );
			float avgEnergy = 0;
			float averageInt = 0;
//				find and combine all lines within FWHM/2 of most intense line
			for ( j=0; j<lines.numberOfLines(); j++ ) {
				float en = lines.energy(j);
				if ( fabs ( en - maxEnergy ) > fwhm2 ) continue;
				float inten = lines.relative(j);
				avgEnergy += en * inten;
				averageInt += inten;
			};
			avgEnergy /= averageInt;
			energies.push_back( avgEnergy );
			elements.push_back( el );
		};
	};
};



void fpPrep (FPstorage &storage, const XrayMaterial sample, const XRFconditions &conditions_in,
            vector <XrayLines> &pureLines ) {

//		reset storage to match this specimen and conditions

	int i;
	float temp;

	pureLines.resize ( 0 );
//		save sample element list
	storage.sampleElements.resize( sample.element_list().size() );
	for ( i=0; i<sample.element_list().size(); i++ ) storage.sampleElements[i] = sample.element_list()[i];

//		create list of x-ray absorption edges excited by this xray source
	sampleEdgeList ( storage.sampleElements, conditions_in.source, storage.sampleEdges, conditions_in.eMin );

//		sort by decreasing edge energy (sort in standard template library is increasing)
	sort( storage.sampleEdges.begin(), storage.sampleEdges.end() );
	reverse( storage.sampleEdges.begin(), storage.sampleEdges.end() );

//		generate vector of lines emitted by vacancy at each edge (if any)
	for ( i=0; i<storage.sampleEdges.size(); i++ ) {
		XrayLines thisLine ( storage.sampleEdges[i] );
		if ( thisLine.numberOfLines() <= 0 ) continue;
//			set intensity factor to zero for each line
		int lineIndex;
		for ( lineIndex=0; lineIndex<thisLine.numberOfLines(); lineIndex++ ) {
			thisLine.factor ( lineIndex, 0.0 );
		};
		pureLines.push_back ( thisLine );
	};
//		keep track of index in element list associated with each edge so that
//			corresponding fractions and absorption tables can be found easily
	storage.elementIndices.resize( pureLines.size() );
	for ( i=0; i<pureLines.size(); i++ ) {
//			generate list of element indices
		int elementIndex;
		for ( elementIndex=0; elementIndex<storage.sampleElements.size(); elementIndex++ ) {
			if ( storage.sampleElements[elementIndex] == pureLines[i].edge().element() ) {
				storage.elementIndices[i] = elementIndex;
			};
		};
	};

//		calculate some quantities which don't depend on individual lines
	storage.sinExcit = sin ( conditions_in.excitAngle * RADDEG );
	if ( storage.sinExcit < 1.0e-6f ) storage.sinExcit = 1.0e-6f;
	storage.sinEmerg = sin ( conditions_in.emergAngle * RADDEG );
	if ( storage.sinEmerg < 1.0e-6f ) storage.sinEmerg = 1.0e-6f;
	storage.geometry =  ( storage.sinExcit / storage.sinEmerg ) / ( 4.0 * PI );

//		include optic center energy in excitation energies list
	float opticCenterEnergy[2] = { -1, -1 };
	if( !conditions_in.optic.DefaultCheck() ) {
		opticCenterEnergy[0] = conditions_in.optic.GetCenterEnergy() - conditions_in.optic.GetBandwidth()/2;
		opticCenterEnergy[1] = conditions_in.optic.GetCenterEnergy() + conditions_in.optic.GetBandwidth()/2;
	}

//		generate continuum and characteristic line excitation energies and intensities
//		(intensities will be pre-multiplied by energy intervals and integration coefficients)
//      Temporary fix for dependence of calculated intensity on element list
//          Include all elements in list of edge energies for excitation integrals
	int maximumZ = 94;	//	end of cross-section tables
	vector <Element> all_elements;
	int z;
	for ( z=1; z<=maximumZ; z++ ) {
        Element new_el( z );
        all_elements.push_back( new_el );
    }
    vector <XrayEdge> all_edges;
	sampleEdgeList ( all_elements, conditions_in.source, all_edges, conditions_in.eMin );
	fpExcitation( all_edges, conditions_in.source, opticCenterEnergy, conditions_in.eMin, storage.excitEnergies, storage.excitIntensities );
//	cout << "excit n all elements" << storage.excitEnergies.size() << endl;
//	fpExcitation( storage.sampleEdges, conditions_in.source, opticCenterEnergy, conditions_in.eMin, storage.excitEnergies, storage.excitIntensities );
//	cout << "excit n std elements" << storage.excitEnergies.size() << endl;

//		apply incident beam corrections
	fpIncidentBeam( conditions_in, storage.excitEnergies, storage.excitIntensities );

//		calculate x-ray fluorescence intensity for each pure element emission line

	int edgeIndex;
	for ( edgeIndex=0; edgeIndex<pureLines.size() ;edgeIndex++ ) {
//			get index of corresponding info in element and absorption table vectors
		int ePri = storage.elementIndices[edgeIndex];
//			calculate subshell absorption for this edge at excitation energies
		vector <float> edgeAbs;
		fpEdgeAbsorption ( pureLines[edgeIndex].edge(),
            sample.cross_section_table ( storage.sampleElements[ePri] ),
            storage.excitEnergies, edgeAbs );
//			calculate pure element absorption at incient energies (for pure element emission)
		vector <float> pureIncAbs(storage.excitEnergies.size());
		for ( i=0; i<storage.excitEnergies.size(); i++ ) {
			pureIncAbs[i] = sample.cross_section( storage.sampleElements[ePri], storage.excitEnergies[i] );
		};
		int lineIndex;
		for ( lineIndex=0; lineIndex<pureLines[edgeIndex].numberOfLines() ;lineIndex++ ) {
//	primary fluorescence
//				fluorescence of line in pure element under same measurement conditions as sample
			float muSp = sample.cross_section( storage.sampleElements[ePri], pureLines[edgeIndex].energy(lineIndex) );
			float pure = fpPrimary ( pureLines[edgeIndex], edgeAbs,
			 	1.0, storage.excitEnergies, storage.excitIntensities, muSp,
			 	pureIncAbs, storage.sinExcit, storage.sinEmerg, storage.geometry );
//				add primary fluorescence into line intensity factor for pure element lines
			 temp = pureLines[edgeIndex].factor( lineIndex );
			 pureLines[edgeIndex].factor( lineIndex, pure + temp );
//	no secondary fluorescence calculated for pure elements
		};
//	Coster_Kronig transitions from primary to secondary edges, for sample lines and pure lines
//						only need to check edges which are higher energy => lower indices in list
//						make sure elements and energy levels match
		int secEdgeIndex;
		for ( secEdgeIndex=edgeIndex+1; secEdgeIndex<pureLines.size(); secEdgeIndex++ ) {
			float cktemp = pureLines[edgeIndex].edge().cktotal( pureLines[secEdgeIndex].edge() );
			if ( cktemp <= 0.0 ) continue;
			int secLineIndex;
			for ( secLineIndex=0; secLineIndex<pureLines[secEdgeIndex].numberOfLines(); secLineIndex++ ) {
//					Conter-Kronig transitions for pure element
				float muSp = sample.cross_section( storage.sampleElements[ePri], pureLines[secEdgeIndex].energy(secLineIndex) );
				float cksum = fpCK ( pureLines[secEdgeIndex], edgeAbs, pureLines[edgeIndex].edge(),
					1.0, storage.excitEnergies, storage.excitIntensities, muSp,
					 pureIncAbs, storage.sinExcit, storage.sinEmerg, storage.geometry );
//						add Coster_Kronig transitions into line intensity factor for pure lines
				temp = pureLines[secEdgeIndex].factor(secLineIndex);
				pureLines[secEdgeIndex].factor(secLineIndex, cksum + temp);
			};  //	end of loop over Conter-Kronig secondary emission lines
		};  //	end of loop over Conter-Kronig secondary absorption edges

	};  //	end of loop over primary emission lines

//		apply incident beam corrections
//		apply detector response correction
	for ( edgeIndex=0; edgeIndex<pureLines.size() ;edgeIndex++ ) {
		int lineIndex;
		for ( lineIndex=0; lineIndex<pureLines[edgeIndex].numberOfLines() ;lineIndex++ ) {
			float lineEnergy = pureLines[edgeIndex].energy( lineIndex );
			float emergCorr = fpEmergentBeam( lineEnergy, conditions_in );
			float detResp = conditions_in.detector.response( lineEnergy );
			temp = pureLines[edgeIndex].factor( lineIndex );
			pureLines[edgeIndex].factor( lineIndex, temp * emergCorr * detResp );
		};
	};
	//  Add vector of pure element lines to fpStorage (for matrix efect factor calculation in fpCalc)
    storage.pureLines = pureLines;
};


void fpCalc(const FPstorage &storage, const XrayMaterial sample, const XRFconditions &conditions_in,
            std::vector <XrayLines> &sampleLines ) {

	int i;
	float temp;

	sampleLines.resize ( 0 );

	if ( sample.number_of_elements() != storage.sampleElements.size() ) {
		cout << "fpCalc  Bad Element list   " << sample.number_of_elements() << "  " << storage.sampleElements.size() << endl;
		return;
	};

	for ( i=0; i<storage.sampleEdges.size(); i++ ) {
//			generate vector of lines emitted by vacancy at each edge
		XrayLines thisLine ( storage.sampleEdges[i] );
		if ( thisLine.numberOfLines() <= 0 ) continue;
//			set intensity factor to zero for each line
		int lineIndex;
		for ( lineIndex=0; lineIndex<thisLine.numberOfLines(); lineIndex++ ) {
			thisLine.factor ( lineIndex, 0 );
		};
		sampleLines.push_back ( thisLine );
	};

//		load vector with sample absorption at each excitation energy
	vector <float> sampleIncAbs( storage.excitEnergies.size() );
	for ( i=0; i<sampleIncAbs.size(); i++ ) {
		sampleIncAbs[i] = sample.cross_section( storage.excitEnergies[i] );
	};


//		calculate x-ray fluorescence intensity for each sample emission line

	int edgeIndex;
	for ( edgeIndex=0; edgeIndex<sampleLines.size() ;edgeIndex++ ) {
//			get index of corresponding info in element and absorption table vectors
		int ePri = storage.elementIndices[edgeIndex];
		float fPri = sample.fraction( storage.sampleElements[ePri] );
		if ( fPri <= 0 ) continue;
//			calculate subshell absorption for this edge at excitation energies
		vector <float> edgeAbs;
		fpEdgeAbsorption ( sampleLines[edgeIndex].edge(),
            sample.cross_section_table ( storage.sampleElements[ePri] ),
            storage.excitEnergies, edgeAbs );
		int lineIndex;
		for ( lineIndex=0; lineIndex<sampleLines[edgeIndex].numberOfLines() ;lineIndex++ ) {
//				calculate sample absorption at emission line energy
			float muSpri = sample.cross_section( sampleLines[edgeIndex].energy(lineIndex) );
//	primary fluorescence
			float pri = fpPrimary ( sampleLines[edgeIndex], edgeAbs, fPri,
                storage.excitEnergies, storage.excitIntensities, muSpri,
			 	sampleIncAbs, storage.sinExcit, storage.sinEmerg, storage.geometry, sample.mass_thickness() );
//				add primary fluorescence into line intensity factor
			temp = sampleLines[edgeIndex].factor(lineIndex) + pri;
			if( temp  <= 0 || isnan( temp ) ) {
                cout << "Warning - emission line calculated intensity is zero for ";
                cout << sampleLines[edgeIndex].edge().element().symbol();
                cout << "   " << sampleLines[edgeIndex].symbolSiegbahn( lineIndex );
                cout << "   " << temp << "  fraction " << fPri<< "  energy " << ePri;
                cout << endl;
            }
			sampleLines[edgeIndex].factor(lineIndex, temp);
//            cout << "pri fluor  " << sampleLines[edgeIndex].edge().element().symbol() << "   " << sampleLines[edgeIndex].symbolIUPAC(lineIndex) << "   " << pri << endl;

//	secondary fluorescence induced by this line on lower-energy lines
//					(since edges are ordered by energy, only need to check
//					those which are lower in the list)
            if( fPri < SEC_FLUOR_THRESHOLD ) continue;
			int secEdgeIndex;
			for ( secEdgeIndex=edgeIndex+1; secEdgeIndex<sampleLines.size(); secEdgeIndex++ ) {
//					skip if this line is below the minimum energy
				if ( sampleLines[secEdgeIndex].edge().energy() < conditions_in.eMin ) continue;
//					skip if primary line can't excite this edge
				if ( sampleLines[edgeIndex].energy(lineIndex) < sampleLines[secEdgeIndex].edge().energy() ) continue;
//					get element list index for exciter
				int eSec = storage.elementIndices[secEdgeIndex];
                float fSec = sample.fraction( storage.sampleElements[eSec] );
//				    skip secondary fluorescence calculations if this is a minor element
//                if( fSec < 0.001f ) continue;
//					skip on some other conditions to speed things up
//				if( fSec * 1000 < fPri ) continue;
				int secLineIndex;
//						calculate subshell absorption for secondary edge at primary line energy
				vector <float> lineEnergy(1);
				vector <float> secAbs(1);
				lineEnergy[0] = sampleLines[edgeIndex].energy(lineIndex);
				fpEdgeAbsorption ( sampleLines[secEdgeIndex].edge(),
                        sample.cross_section_table ( storage.sampleElements[eSec] ),
                        lineEnergy, secAbs );
//					cout << sampleLines[edgeIndex].symbolIUPAC(lineIndex) << "  " << secEdgeIndex << "  " << sampleLines[secEdgeIndex].numberOfLines() << endl;
                float sec_total = 0;
				for ( secLineIndex=0; secLineIndex<sampleLines[secEdgeIndex].numberOfLines(); secLineIndex++ ) {
					float muSsec = sample.cross_section( sampleLines[secEdgeIndex].energy(secLineIndex) );
					float sec = fpSecondary ( sampleLines[secEdgeIndex], secAbs[0],
			 			fSec, sampleLines[edgeIndex], lineIndex, edgeAbs,
			 			fPri, storage.excitEnergies, storage.excitIntensities, muSsec, muSpri,
			 			sampleIncAbs, storage.sinExcit, storage.sinEmerg, storage.geometry, sample.mass_thickness() );
//					cout << edgeIndex << "  " << lineIndex << "  " << secEdgeIndex << "  " << secLineIndex << "  " << sec << endl;
//							add secondary fluorescence into line intensity factor for secondary line
					temp = sampleLines[secEdgeIndex].factor(secLineIndex);
					sampleLines[secEdgeIndex].factor(secLineIndex, sec + temp);
					sec_total += sec;
//                    cout << "sec fluor  from " << sampleLines[edgeIndex].edge().element().symbol() << "   " << sampleLines[edgeIndex].symbolIUPAC(lineIndex);
//                    cout << "   to " << sampleLines[secEdgeIndex].edge().element().symbol() << "  " << sampleLines[secEdgeIndex].symbolIUPAC(secLineIndex) << "   s " << sec_total << endl;
				};
//					end of loop over secondary emission lines
			 };
//				end of loop over secondary absorption edges
		};
//			end of loop over primary emission lines

//	Coster_Kronig transitions from primary to secondary edges, for sample lines and pure lines
//						only need to check edges which are higher energy => lower indices in list
//						make sure elements and energy levels match
		int secEdgeIndex;
		for ( secEdgeIndex=edgeIndex+1; secEdgeIndex<sampleLines.size(); secEdgeIndex++ ) {
			float cktemp = sampleLines[edgeIndex].edge().cktotal( sampleLines[secEdgeIndex].edge() );
			if ( cktemp <= 0.0 ) continue;
			int secLineIndex;
			for ( secLineIndex=0; secLineIndex<sampleLines[secEdgeIndex].numberOfLines(); secLineIndex++ ) {
				float muSsec = sample.cross_section( sampleLines[secEdgeIndex].energy(secLineIndex) );
				float cksum = fpCK ( sampleLines[secEdgeIndex], edgeAbs, sampleLines[edgeIndex].edge(),
					fPri, storage.excitEnergies, storage.excitIntensities, muSsec,
					 sampleIncAbs, storage.sinExcit, storage.sinEmerg, storage.geometry, sample.mass_thickness() );
//					cout << sampleLines[edgeIndex].symbolIUPAC(lineIndex) << "  " << sampleLines[secEdgeIndex].edge().symbol();
//					cout << "   " << edgeIndex << "  " << lineIndex << "  " << secEdgeIndex << "  " << cktemp << "   " << cksum << endl;
//						add Coster_Kronig transitions into line intensity factor for sample lines
				temp = sampleLines[secEdgeIndex].factor(secLineIndex);
				sampleLines[secEdgeIndex].factor(secLineIndex, cksum + temp);
			};
//				end of loop over Conter-Kronig secondary emission lines
		};
//			end of loop over Conter-Kronig secondary absorption edges
	};
//		end of loop over primary absorption edges

//		apply emergent beam corrections
//		apply detector response correction
	for ( edgeIndex=0; edgeIndex<sampleLines.size(); edgeIndex++ ) {
		int lineIndex;
		for ( lineIndex=0; lineIndex<sampleLines[edgeIndex].numberOfLines(); lineIndex++ ) {
			float lineEnergy = sampleLines[edgeIndex].energy( lineIndex );
			float emergCorr = fpEmergentBeam( lineEnergy, conditions_in );
			float detResp = conditions_in.detector.response( lineEnergy );
			temp = sampleLines[edgeIndex].factor(lineIndex);
			sampleLines[edgeIndex].factor( lineIndex, temp * emergCorr * detResp );
/*			if( temp * emergCorr * detResp <= 0 ) {
                cout << sampleLines[edgeIndex].edge().element().symbol();
                cout << "   " << sampleLines[edgeIndex].symbolSiegbahn( lineIndex );
                cout << "   " << temp << "  " << emergCorr << "  " << detResp;
                cout << endl;
            }                       */
            //  Add matrix effect factor
            unsigned int ePri = storage.elementIndices[edgeIndex];
            float fPri = sample.fraction( storage.sampleElements[ePri] );
            float pureInt = storage.pureLines[edgeIndex].intensity( lineIndex );
            float mf = sampleLines[edgeIndex].intensity( lineIndex ) / ( fPri * pureInt );
            sampleLines[edgeIndex].matrix( lineIndex, mf );
//            cout << "fpCalc  M " << sampleLines[edgeIndex].edge().element().symbol();
//            cout << " " << sampleLines[edgeIndex].symbolIUPAC( lineIndex );
//            cout << "   i " << sampleLines[edgeIndex].intensity( lineIndex );
//            cout << "   p " << pureInt << "   f " << fPri << "   mf " << mf;
//            cout << "   M " << sampleLines[edgeIndex].matrix( lineIndex );
//            cout << endl;
		};
	};

};


void fpRayleigh(const FPstorage &storage, const XrayMaterial &sample, const XRFconditions &conditions_in,
            std::vector <XrayLines> &scatterLines ) {

//	calculates Rayleigh scatter of tube characteristic lines & sets factor in XrayLInes objects

	if ( sample.number_of_elements() != storage.sampleElements.size() ) {
		cout << "Bad Element list   " << sample.number_of_elements() << "  " << storage.sampleElements.size() << endl;
		return;
	};
//		some things that don't depend on energy
	float theta = conditions_in.excitAngle * RADDEG + conditions_in.emergAngle * RADDEG;
//		get list with intensities of tube characteristic lines
	conditions_in.source.lines( scatterLines, conditions_in.eMin );
	int edgeIndex;
	for ( edgeIndex=0; edgeIndex<scatterLines.size(); edgeIndex++ ) {
		int lineIndex;
		for ( lineIndex=0; lineIndex<scatterLines[edgeIndex].numberOfLines(); lineIndex++ ) {
//				calculate Rayleigh scatter for each line
			float lineEn = scatterLines[edgeIndex].energy( lineIndex );
			float lineInt = scatterLines[edgeIndex].factor( lineIndex );
//			    apply incident beam corrections
			lineInt *= fpIncidentBeam ( lineEn, conditions_in );
//				calculate sample absorption at desired energy
			float muSamp = sample.cross_section( lineEn );
//				calculate Rayleigh cross section at given energy and angle
			float sigmaCoh = sample.coherent( lineEn, theta );
			float denominator = conditions_in.excitCosecant * muSamp + conditions_in.emergCosecant * muSamp;
//				include factor for finite thickness here, since it does not depend on spectrum energy (only on line energy)
			if( sample.mass_thickness() > 0 ) {
				float expArg = denominator * sample.mass_thickness();
				if( expArg < EXP_FLOAT_TEST ) denominator /= ( 1 - exp( - expArg ) );
			};
			float cohInt = lineInt * conditions_in.emergCosecant * sigmaCoh / denominator;
//		    apply emerging beam corrections
            float emergCorr = fpEmergentBeam( lineEn, conditions_in );
//		    apply detector response correction
			float detResp = conditions_in.detector.response( lineEn );
			cohInt *= emergCorr * detResp;
//          put the calculated Rayleigh intensity into the factor for this emission line
			scatterLines[edgeIndex].factor( lineIndex, cohInt );
        };
	};
};


void fpContScat(const FPstorage &storage, const XrayEnergyCal &cal_in, const XrayMaterial &sample,
				const XRFconditions &conditions_in, vector <float> &continuumSpec ) {

//	calculates background from Compton and Rayleigh scatter of tube continuum

	if ( sample.number_of_elements() != storage.sampleElements.size() ) {
		cout << "Bad Element list   " << sample.number_of_elements() << "  " << storage.sampleElements.size() << endl;
		return;
	};
	int nChan = continuumSpec.size();
	int iChan;
	for( iChan=0; iChan<nChan; iChan++ ) {
		vector <float> contEn(1);
		vector <float> contInt(1);
		contEn[0] = cal_in.energy( iChan );
		if( contEn[0] <= 0 ) {
			continuumSpec[iChan] = 0;
			continue;
		};
//			find continuum intensity at desired energy
		contInt[0] = conditions_in.source.continuum( contEn[0] );
//			apply incident beam corrections
		fpIncidentBeam ( conditions_in, contEn, contInt );
//			calculate sample absorption at desired energy
		float muSamp = sample.cross_section( contEn[0] );
//			calculate Compton and Rayleigh cross section at given energy and angle
		float theta = conditions_in.excitAngle * RADDEG + conditions_in.emergAngle * RADDEG;
		float sigmaIncoh = sample.incoherent( contEn[0], theta );
		float sigmaCoh = sample.coherent( contEn[0], theta );
//   ***** should probably add window scatter here, someday   ******     and dust scatter    ******************
//			ignore Compton shift and use same energy for incident and scattered, Compton and Rayleigh
		float denominator = conditions_in.excitCosecant * muSamp + conditions_in.emergCosecant * muSamp;
		float bkgEst = contInt[0] * conditions_in.emergCosecant * ( sigmaCoh + sigmaIncoh ) / denominator;
		if( sample.mass_thickness() > 0 ) {
			float expArg = denominator * sample.mass_thickness();
			if( expArg < EXP_FLOAT_TEST ) bkgEst *= ( 1 - exp( - expArg ) );
		};
//		apply incident beam corrections
//		apply detector response correction
        int i;
		for ( i=0; i<contEn.size(); i++ ) {
            float emergCorr = fpEmergentBeam( contEn[i], conditions_in );
            float detResp = conditions_in.detector.response( contEn[i] );
			bkgEst *= emergCorr * detResp;
		};
//			result is per keV, so multiply by channel width in keV to get counts in each channel
		bkgEst *= cal_in.energyPerChannel() / 1000;
		continuumSpec[iChan] = bkgEst;
	};
};


void fpCompton(const FPstorage &storage, const XrayEnergyCal &cal_in, const XrayMaterial &sample,
				const XRFconditions &conditions_in, SpectrumComponent &component_out ) {

//	calculates peaks from Compton and Rayleigh scatter of tube characteristic lines

	if ( sample.number_of_elements() != storage.sampleElements.size() ) {
		cout << "Bad Element list   " << sample.number_of_elements() << "  " << storage.sampleElements.size() << endl;
		return;
	};
//		need vectors for the functions that calculate absorption corrections
	vector <float> lineEn(1);
	vector <float> lineInt(1);
//		some things that don't depend on energy
	float theta = conditions_in.excitAngle * RADDEG + conditions_in.emergAngle * RADDEG;
	const int nChan = component_out.spectrum.size();
//		get tube characteristic lines
	vector <XrayLines> sourceLines;
	conditions_in.source.lines( sourceLines, conditions_in.eMin );
	int edgeIndex;
	for ( edgeIndex=0; edgeIndex<sourceLines.size(); edgeIndex++ ) {
		int lineIndex;
		for ( lineIndex=0; lineIndex<sourceLines[edgeIndex].numberOfLines(); lineIndex++ ) {
            //  Check to see if this emission line should be included in the component
            if( ! checkComponent( component_out, sourceLines[edgeIndex], lineIndex ) ) continue;
//      Remove all lines except L alpha (L3-M4,5) (or L beta 1 l2-M4) line
/*            if( sourceLines[edgeIndex].edge().index() == L3 ) {
//                if( ! ( sourceLines[edgeIndex].edgeSource(lineIndex).index() == M4 || sourceLines[edgeIndex].edgeSource(lineIndex).index() == M5 ) ) {
                    continue;
//                }
            }
            if( sourceLines[edgeIndex].edge().index() == L2 ) {
                if( ! ( sourceLines[edgeIndex].edgeSource(lineIndex).index() == M4 ) ) {
                    continue;
                }
            }
            if( sourceLines[edgeIndex].edge().index() == L1 ) {
                continue;
            }   */
//				calculate Rayleigh scatter for each line
			lineEn[0] = sourceLines[edgeIndex].energy( lineIndex );
			lineInt[0] = sourceLines[edgeIndex].intensity( lineIndex );
//			    apply incident beam corrections
			fpIncidentBeam ( conditions_in, lineEn, lineInt );
//				calculate sample absorption at desired energy
			float muSamp = sample.cross_section( lineEn[0] );
//				calculate Compton scatter for this line with correct Compton profile
			float enC = ScatterXsectTable::eCompton( lineEn[0], theta );
//				calculate sample absorption at the Compton shifted line energy - use only this energy to save compute time
			float muSampC = sample.cross_section( enC );
			float denominator = conditions_in.excitCosecant * muSamp + conditions_in.emergCosecant * muSampC;
//				include factor for finite thickness here, since it does not depend on spectrum energy (only on line energy)
			if( sample.mass_thickness() > 0 ) {
				float expArg = denominator * sample.mass_thickness();
				if( expArg < EXP_FLOAT_TEST ) denominator /= ( 1 - exp( - expArg ) );
			};
//		    apply incident beam corrections
//		    apply detector response correction
            float emergCorr = fpEmergentBeam( enC, conditions_in );
            float detResp = conditions_in.detector.response( enC );
			int iChanMin = cal_in.channel( lineEn[0] - 3 * ( lineEn[0] - enC ) ) - 1;
			if ( iChanMin < 0 ) iChanMin = 0;
			int iChanMax = cal_in.channel( lineEn[0] ) + 2;
			if ( iChanMax > nChan ) iChanMax = nChan;
			//int iChanMid = ( iChanMin + iChanMax ) / 2;
			//float e = cal_in.energy(iChanMid );
			int iChan;
			for ( iChan=iChanMin; iChan<iChanMax; iChan++ ) {
//   ***** should probably add window scatter here, someday   ******            and dust *****************
				float en = cal_in.energy(iChan );
				if( en <= 0 ) continue;
//					calculate Compton cross section at given energy and angle
				float sigmaIncoh = sample.incoherent( lineEn[0], theta, en );
//					doubly-differential Compton cross section is per eV,
//							so multiply by channel width in eV to get counts in each channel
				component_out.spectrum[iChan] += lineInt[0] * conditions_in.emergCosecant * sigmaIncoh
                    / denominator * emergCorr * detResp * cal_in.energyPerChannel( iChan );
			};
		};
	};
};
