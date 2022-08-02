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

#include <math.h>
#include "XrayDetector.h"
#include "XrayEdge.h"
#include "XrayLines.h"
#include "XRFconstants.h"
#include <sstream>


using namespace std;

//	this class does resolution and response calculations for an AmpTek XR-100CR detector
//		Modified to calculate resolution from Fano factor, separately for Si, Ge, & CZT     April 12, 2011
//  Modified to add escape calculations    Oct. 29, 2013   W. T. Elam   APL/UW
//  Modified May 4, 2016
//  	Changed name of class from AmpTekDet to XrayDetector
//  	Added ability to retrieve and change Fano factor, energy per pair, and resolution reference energy
//  	Moved Mn Ka (usual resolution reference energy) to define above
//  Modified Aug. 10, 2016
//      Prevent resolution from becoming not-a-number at zero energy for bad values of reference resolution or fano factor
//      Change actually made in file XrayDetector.h
//  Modified Oct. 5, 2016
//      Fixed bug in escape peak fraction calculation, wrong MAC used in equation
//  Modified Feb. 1, 2017
//      Include default values for each detector type here (instead of where conditions are set up)
//      Use XrayMaterial class for window, dead layer, and active layer
//      Note: response was calculated incorrectly for CdTe before this change, element fractions were left out
//  Modified Jan. 24, 2018
//      Add Si K alpha 3 and K alpha 4 satellite lines to escape peaks for Si detectors
//  Modified Apr. 11, 2018
//      Calculate resolution as electronic noise plus Fano instead of Mn Ka resolution adjusted using Fano
//          This makes it easier and more reasonable to avoid square root of a negative number or zero resolution
//      Moved setResolution into this file as part of this change
//  Modified Aug. 26, 2020
//      Include Shcole and Procop model for incomplete charge collection tail on peaks
//      F. Scholze and M. Procop, X-Ray Spectrometry 2009, 39, 312-321.
//  Modified Sept. 15, 2020
//      Add Compton escape calculations of shelf at low energies
//  Modified Sept. 30, 2020
//      Test and debug Compton escape calculation with X-ray tube data
//  Modified Dec. 7, 2020
//      Added detector shelf calculation from photoelectron and Auger electron escape (active volume and front contact)
//      Fix bug in escape satellite lines, wrong place in loop making them included multiple times
//      Include K edge jump ratio in L lines when photon energy is above K edge
//  Modified Mar. 8, 2021   Change tail Czero to be 0.75 below Si K edge and 0.92 above
//  Modified Mar. 8, 2021   Change tail Czero to be 0.75 always and Zzero to follow linear formula below Si K edge
//  Modified Mar. 15, 2021  Front contact shelf disabled to improve Na results
//  Modified May 14, 2021   Move shelf factor and slope to XrayDetector and control via -T option (header change only)
//  Modified June 17, 2021  Peter Nemere found and fixed bug in detector shelf calculation - uninitialized variable new_Auger.binding_energy
//  Modified June 19, 2021  Fixing above bug had repercussions for calibration, lots of flailing around to check shelf factor, slope, & tail parameters
//                          Final choice was shelf factor 1 & slope 0, tail unchanged, and front contact shelf enabled with thickness 150 um
//  Modified July 10, 2021  Add pulse resolving time for simple pulse pileup calculation (header change only for now)


// Principal Auger electron energies for selected elements added for shelf calculations
//      From R. N. Yasko, and R. D. Whitmoyer, Journal of Vacuum Science and Technology 8, 733 (1971); https://doi.org/10.1116/1.1315385
const int number_energies_Auger_KLL = 19;
//                                                                  H   He  Li  Be  B   C   N   O   F   Ne   Na   Mg   Al   Si   P    S    Cl   Ar
//                                                          Z   0   1   2   3   4   5   6   7   8   9   10   11   12   13   14   15   16   17   18
const float energies_Auger_KLL[number_energies_Auger_KLL] = {   0,  0,  0,  0,  0,  176,268,383,516,659,818, 1039,1180,1478,1730,1850,2105,2375,2660};

XrayDetector::XrayDetector( float detR_in, float detW_in, float detD_in, float detA_in, DetectorType detType_in ) {
	detType = detType_in;
	Element be(4);	//	beryllium
	XrayMaterial be_win( be );
	window = be_win;
	res_fwhm_energy = RESOLUTION_REFERENCE_ENERGY;
	//  Set up detector materials and set default values (different for each detector type)
//			*****  need metal contact ?   **************
	switch ( detType ) {
		default:
		case SI_PIN:    {
//				instantiate AmpTek Detector XR-100 PIN diode
//				resolution 250 eV, Be window 0.5 mil (12.5 micron), dead layer 1 micron, active layer .5 mm
//				resolution changed to 180 eV to better match AmpTek XR-100CR detector # N5270
//			    XrayDetector det_PIN( 180.0f, 12.5f, 0.1f, 0.5f, Si );
            float detW = ( detW_in==0 ? 12.5f : detW_in ) * CM_MICRON;
            window.thickness( detW );
			Element si(14);	//	Silicon
			XrayMaterial si_det( si );
			deadLayer = si_det;
			float detD = ( detD_in==0 ? 0.1f : detD_in ) * CM_MICRON;
			deadLayer.thickness( detD );
            activeLayer = si_det;
			float detA = ( detA_in==0 ? 0.5f : detA_in ) * CM_MM;
			activeLayer.thickness( detA );
			pair_energy = 3.86f;	//	Handbook of X-ray Spectrometry, 2002, p216.
			fano_factor = 0.12f;
            setResolution( detR_in==0 ? 180.0f : detR_in );
            //  Set up the front contact
			Element al(13);	//	Aluminum
			Element ox(8);	//	Oxygen
			vector <Element> front_contact_elements;
			front_contact_elements.push_back( al );
			front_contact_elements.push_back( si );
			front_contact_elements.push_back( ox );
			vector <float> front_contact_fractions = { 0.31f, 0.32f, 0.37f }; //  Roughly AlSiO2, Al contact plus SiO2 passivation layer (not active)
			XrayMaterial al_contact( front_contact_elements, front_contact_fractions );
			frontContact = al_contact;
			frontContact.thickness( CONTACT_SDD_UM ); //  defined in XrayDetector.h, for SDD from Scholze and Procop, need better value for individual detectors
			break;
			};
		case SI_SDD:    {
//				instantiate AmpTek Detector SDD
//				resolution 150 eV, Be window 0.5 mil (12.5 micron), dead layer 1 micron, active layer .5 mm
//			    XrayDetector det_SDD( 155.0f, 12.5f, 0.1f, 0.5f, Si );
            float detW = ( detW_in==0 ? 12.5f : detW_in ) * CM_MICRON;
            window.thickness( detW );
			Element si(14);	//	silicon
			XrayMaterial si_det( si );
			deadLayer = si_det;
			float detD = ( detD_in==0 ? 0.1f : detD_in ) * CM_MICRON;
			deadLayer.thickness( detD );
            activeLayer = si_det;
			float detA = ( detA_in==0 ? 0.5f : detA_in ) * CM_MM;
			activeLayer.thickness( detA );
			pair_energy = 3.86f;	//	Handbook of X-ray Spectrometry, 2002, p216.
			fano_factor = 0.12f;
            setResolution( detR_in==0 ? 155.0f : detR_in );
            //  Set up the front contact
			Element al(13);	//	Aluminum
			Element ox(8);	//	Oxygen
			vector <Element> front_contact_elements;
			front_contact_elements.push_back( al );
			front_contact_elements.push_back( si );
			front_contact_elements.push_back( ox );
			//  Roughly AlSiO2, Al contact plus SiO2 passivation layer (which is not part of active volume)
			vector <float> front_contact_fractions = { 0.31f, 0.32f, 0.37f };
			XrayMaterial al_contact( front_contact_elements, front_contact_fractions );
			frontContact = al_contact;
			frontContact.thickness( CONTACT_SDD_UM ); //  defined in XrayDetector.h, for SDD from Scholze and Procop, need better value for individual detectors
			break;
			};
		case CD_TE: {
//				CdTe detector, resolution 290 eV, Be window 250 micron, dead layer 1 micron, active layer 1 mm
//				See http://www.amptek.com/xrcdtaps.html   April 23, 2009
//			    XrayDetector det_CdTe( 500.0f, 12.5f, 0.1f, 1.0f, CdTe );
            float detW = ( detW_in==0 ? 12.5f : detW_in ) * CM_MICRON;
            window.thickness( detW );
            const int ne_CdTe = 2;
            const int z_CdTe[ne_CdTe] = { 48, 52 };
            float f_CdTe[ne_CdTe];
//				use atomic weights to convert from formula (1:1 ratio) to weight fractions
			Element cd(48);	//	cadmium
			Element te(52);	//	tellurium
			f_CdTe[0] = cd.atomicWeight() / ( cd.atomicWeight() + te.atomicWeight() );
			f_CdTe[1] = te.atomicWeight() / ( cd.atomicWeight() + te.atomicWeight() );
			XrayMaterial cd_te( ne_CdTe, z_CdTe, f_CdTe );
			cd_te.density( 5.85f );//				5.85 is CdTe density in gm/cm3
			deadLayer = cd_te;
			float detD = ( detD_in==0 ? 0.1f : detD_in ) * CM_MICRON;
			deadLayer.thickness( detD );
			activeLayer = cd_te;
			float detA = ( detA_in==0 ? 1.0f : detA_in ) * CM_MM;
			activeLayer.thickness( detA );
			pair_energy = 5.0f;	//	Redus et al., MRS Bulletin,
			fano_factor = 0.089f;
            setResolution( detR_in==0 ? 500.0f : detR_in );
            //  Set up the front contact
			Element al(13);	//	Aluminum
			XrayMaterial al_contact( al );  //  Just Al for now, may need to modify this if one of these detectors is better characterized
			frontContact = al_contact;
			frontContact.thickness( CONTACT_SDD_UM ); //  value for SDD, need data for CdTe detector
			break;
			}
		case HP_GE: {
//				HPGe detector, resolution 295 eV, Be window 75 micron, dead layer 1 micron, active layer 10 mm
//				Canberra Model GL0110S Serial 09024876       Dec. 7, 2010
//				actual resolution in measurements was closer to 650 eV related to Mn Ka
//  			XrayDetector det_HPGe( 650.0f, 75.0f, 0.1f, 10.0f, HPGe );
            float detW = ( detW_in==0 ? 75.0f : detW_in ) * CM_MICRON;
            window.thickness( detW );
			Element ge(32);	//	germanium
			XrayMaterial ge_det( ge );
			deadLayer = ge_det;
			float detD = ( detD_in==0 ? 0.1f : detD_in ) * CM_MICRON;
			deadLayer.thickness( detD );
			activeLayer = ge_det;
			float detA = ( detA_in==0 ? 10.0f : detA_in ) * CM_MM;
			activeLayer.thickness( detA );
			pair_energy = 2.96f;	//	Handbook of X-ray Spectrometry, 2002, p216.
//			fano_factor = 0.08f;
			fano_factor = 0.15f;   //  from APS beamline 1ID spectrum, GHSR1002301.001  DEC 04, 2010 09:38:20.174
            setResolution( detR_in==0 ? 649.95f : detR_in );
			break;
			};
	};
	initialize_shelf();
};

void XrayDetector::setResolution ( const float res_in, const float ref_energy ) {
    if( ref_energy >= 0 ) res_fwhm_energy = ref_energy; //  Resolution is given at this energy
    if ( res_in <= 0 ) return;
    //  Convert resolution at reference energy to zero energy (electronic noise contribution)
	float factor = EIGHT_LN_2 * fano_factor * pair_energy;   //  8 ln2, defined in XRFconstants.h
    float en2 = res_in*res_in - factor * res_fwhm_energy;
    if( en2 > 0 ) electronic_noise = sqrt( en2 );
    else electronic_noise = 0;
    return;
};



float XrayDetector::resolution ( const float energy ) const {
//		calculate detector resolution as a function of energy
//  ratio of fwhm to sigma in Gaussian sqrt( 8 ln2 ) ~= 2.35482
//  8 ln2 = 5.545177
//  Put check in member functions to set fano and energy_per_pair
//      to be sure argument to sqrt does not become negative at zero energy
    if( energy < 0 ) return electronic_noise;
	float factor = EIGHT_LN_2 * fano_factor * pair_energy;   //  8 ln2, defined in XRFconstants.h
	float res = sqrt ( electronic_noise*electronic_noise + factor * energy );
	return res;
};


float XrayDetector::response ( const float energy ) const {
//		calculate detector response at a given energy
	float resp = 1;
//			*****  need metal contact ?   **************
	resp = window.transmission( energy )
            * frontContact.transmission( energy )
            * deadLayer.transmission( energy )
            * activeLayer.absorption( energy );
	return resp;
};

float XrayDetector::escape ( const float energy, vector<EscapeLines> &escape_line_vector ) const {
    //  This calculation follows S. J. B. Reed, "Electron Microprobe Analysis" 2nd ed. (Cambridge, 1997) ISBN 0-521-41956-5 p115-117.
    //  It covers front escape at normal incidence only.
    //  He gives reference to (Dalton 1974) for back escape and to (Statham 1976a) for non-normal detector incidence.
    //  This routine also only checks emission lines from the K and L3 edges.
    //  Satellite lines associated with the Si K alpha lines were added Jan. 24, 3018 (only K alpha 3 and K alpha 4)
    //      Intensities and energies are from Table I (page 401) in:
    //      J.L. Campbell, G. Cauchon, M.-C. LÈpy, L. McDonald, J. Plagnard, P. Stemmler, W.J. Teesdale, G. White,
    //      A quantitative explanation of low-energy tailing features of Si(Li) and Ge X-ray detectors, using synchrotron radiation,
    //      Nuclear Instruments and Methods in Physics Research Section A: Accelerators, Spectrometers, Detectors and Associated Equipment,
    //      Volume 418, Issues 2ñ3, 1998, Pages 394-404, ISSN 0168-9002,
    //      https://doi.org/10.1016/S0168-9002(98)00889-4. (http://www.sciencedirect.com/science/article/pii/S0168900298008894)
    escape_line_vector.clear();
    int ie;
    vector <Element> detAE = activeLayer.element_list();
    //  loop over elements in detector active area
    for( ie=0; ie<detAE.size(); ie++ ) {
        float totAbs_incEnergy = activeLayer.cross_section( energy );
        //  check K lines
        XrayEdge edge_K( detAE[ie], K1 );
        float rk = edge_K.jump();
        //  see if photon has enough energy to produce escape
        if ( energy < edge_K.energy() ) continue;
        XrayLines lines_K ( edge_K );
        float satellite_fraction = 0;   //  Used to normalize total escape fraction when satellites lies are included
        if( detAE[ie].Z() == 14 ) {
            //  Add K alpha 3 satellite line from silicon
            float energy_K_alpha3 = 1751.0f;
            float relative_K_alpha3 = 0.057f;
            //  calculate escape energy using line energy of escaped photon
            float escEnergy = energy - energy_K_alpha3;
            EscapeLines new_escape;
            new_escape.energy = escEnergy;
            //  calculate fraction for this escape peak
            float totAbs_escEnergy = activeLayer.cross_section( energy_K_alpha3 );
            float absRatio = totAbs_escEnergy / totAbs_incEnergy;
            float x = 1.0f - absRatio * log ( 1.0f + (1.0f/absRatio) );
            new_escape.fraction = activeLayer.fraction( detAE[ie] ) * relative_K_alpha3 * ( ( rk - 1 ) / rk ) * edge_K.yield() * x / 2.0f;
            escape_line_vector.push_back( new_escape );
//  if( 3680 < energy && energy < 3700 ) cout << "      Ca esc Ka3 satellite " << new_escape.energy << "   f=" << new_escape.fraction * 100 << " %" << endl;
            satellite_fraction += new_escape.fraction;
            //  Add K alpha 4 satellite line from silicon
            float energy_K_alpha4 = 1753.0f;
            float relative_K_alpha4 = 0.030f;
            //  calculate escape energy using line energy of escaped photon
            escEnergy = energy - energy_K_alpha4;
            new_escape.energy = escEnergy;
            //  calculate fraction for this escape peak
            totAbs_escEnergy = activeLayer.cross_section( energy_K_alpha4 );
            absRatio = totAbs_escEnergy / totAbs_incEnergy;
            x = 1.0f - absRatio * log ( 1.0f + (1.0f/absRatio) );
            new_escape.fraction = activeLayer.fraction( detAE[ie] ) * relative_K_alpha4 * ( ( rk - 1 ) / rk ) * edge_K.yield() * x / 2.0f;
            escape_line_vector.push_back( new_escape );
//  if( 3680 < energy && energy < 3700 ) cout << "      Ca esc Ka4 satellite   en=" << new_escape.energy << "   f=" << new_escape.fraction * 100 << " %" << endl;
           satellite_fraction += new_escape.fraction;
        }
        int il;
        for( il=0; il<lines_K.numberOfLines(); il++ ) {
            //  calculate escape energy using line energy of escaped photon
            float escEnergy = energy - lines_K.energy( il );
            //  instantiate new structure to hold info on this escape peak
            EscapeLines new_escape;
            new_escape.energy = escEnergy;
            //  calculate fraction for this escape peak
            float totAbs_escEnergy = activeLayer.cross_section( lines_K.energy( il ) );
            float absRatio = totAbs_escEnergy / totAbs_incEnergy;
            float x = 1.0f - absRatio * log ( 1.0f + (1.0f/absRatio) );
            new_escape.fraction = activeLayer.fraction( detAE[ie] ) * lines_K.relative( il ) * ( ( rk - 1 ) / rk ) * edge_K.yield() * x / 2.0f;
            //  Subtract satellite fraction proportionally from each K line
            new_escape.fraction -= lines_K.relative( il ) * satellite_fraction;
            escape_line_vector.push_back( new_escape );
// if( 3680 < energy && energy < 3700 ) cout << "      Ca esc  " << lines_K.symbolIUPAC(il) << "  en=" << new_escape.energy << "   f=" << new_escape.fraction * 100 << " %" << endl;
        }
        //  check L3 lines
        XrayEdge edge_L3( detAE[ie], L3 );
        //  see if photon has enough energy to produce escape
        if ( energy < edge_L3.energy() ) continue;
        XrayLines lines_L3 ( edge_L3 );
        for( il=0; il<lines_L3.numberOfLines(); il++ ) {
            //  calculate escape energy using line energy of escaped photon
            float escEnergy = energy - lines_L3.energy( il );
            //  instantiate new structure to hold info on this escape peak
            EscapeLines new_escape;
            new_escape.energy = escEnergy;
            //  calculate fraction for this escape peak
            float totAbs_escEnergy = activeLayer.cross_section( lines_L3.energy( il ) );
            float absRatio = totAbs_escEnergy / totAbs_incEnergy;
            float x = 1.0f - absRatio * log ( 1.0f + (1.0f/absRatio) );
            float rk = edge_L3.jump();
            new_escape.fraction = activeLayer.fraction( detAE[ie] ) * lines_L3.relative( il ) * ( ( rk - 1 ) / rk ) * edge_L3.yield() * x / 2.0f;
            //  Correct for K jump ratio when photon energy is above K edge
            if( energy > edge_K.energy() ) new_escape.fraction /= edge_K.jump();
            escape_line_vector.push_back( new_escape );
//       if( 3680 < energy && energy < 3700 ) cout << "      Ca esc  " << lines_L3.symbolIUPAC(il) << "  en=" << new_escape.energy << "   f=" << new_escape.fraction * 100 << " %" << endl;
        }
    }
    float total_escape_fraction = 0;
    for( ie=0; ie<escape_line_vector.size(); ie++ ) {
        total_escape_fraction += escape_line_vector[ie].fraction;
    }
// if( 3680 < energy && energy < 3700 ) cout << "      Ca esc  total   " << total_escape_fraction * 100 << " %" << endl;
    return 1 - total_escape_fraction;
};

//  Tail calculations added August 26, 2020
//      See F. Scholze and M. Procop, X-Ray Spectrometry 2009, 39, 312-321.

const float XrayDetector::tail_C0( const float energy ) const {
    const Element Si( 14 );
    const XrayEdge Si_K_edge( Si, K1 );
    float c0 = 0.75f;
//    if( activeLayer.fraction( Si ) > 0 && energy > Si_K_edge.energy() ) c0 = 0.92f;
//    if( activeLayer.fraction( Si ) > 0 && energy > Si_K_edge.energy() ) c0 = 0.65f;
    return c0;
};

float XrayDetector::z_of_C( const float photon_energy, const float channel_energy ) const {
    //  Calculates fraction of charge collected for a given spectrum channel energy
    //      and the incident photon energy
    //  The returns the depth corresponding to that fraction
    if( photon_energy <= 0 ) return 0;
    float charge_ratio = channel_energy / photon_energy;
//    return z_of_C( charge_ratio );
//}

//float XrayDetector::z_of_C( const float charge_ratio ) const {
    //  Calculates depth from front surface corresponding to a given fraction of charge collection
    //      in incomplete charge collection region near surface
    float tail_z0 = 35e-7f; //  35 nanometers (expressed in cm)
    const Element Si( 14 );
    const XrayEdge Si_K_edge( Si, K1 );
    if( activeLayer.fraction( Si ) > 0 && photon_energy > 1250 && photon_energy < Si_K_edge.energy() ) {
        //  Formula from matching tail for Mg (in MgCO3), Si (in SiO2), and Al (in Al2O3) from PIXL FM elemental calibration
        tail_z0 = 35e-7f + ( photon_energy - 1250 ) * 210e-7f / 500;   //  240 nm at 1740 eV, ~120 nm at 1480 eV, and 35 nm at 1250 eV
    }
    if( charge_ratio >= 1.0f ) return tail_z0;
    float lost_fraction = charge_ratio - tail_C0( photon_energy );
    if( lost_fraction <= 0 || tail_C0( photon_energy ) >= 1 ) return 0;
    float depth_fraction = pow( lost_fraction / ( 1 - tail_C0( photon_energy ) ), 1/tail_a );
    return tail_z0 * depth_fraction;
}

float XrayDetector::exp_term_of_C( const float photon_energy, const float channel_energy ) const {
    //  Calculates one term in the fraction of incident intensity in a single channel of the peak tail
    //  Using this form saves double calculation of this quantity (which is in an innermost loop)
    float z_new = z_of_C( photon_energy, channel_energy );
    float exponent = activeLayer.density() * activeLayer.photo( photon_energy ) * z_new;
    float exp_term = exp( - exponent );
    return exp_term;
}

float XrayDetector::tail_fraction( const float photon_energy, const float channel_E1, const float channel_E2 ) const {
    //  Returns fraction of incident intensity for a given photon energy that appears between
    //      channel_E1 energy and channel_E2 energy in the spectrum
    float z1 = z_of_C( photon_energy, channel_E1 );
    float z2 = z_of_C( photon_energy, channel_E2 );
    float center_z = ( z1 + z2 ) /2;
    float delta_z = z2 - z1;
    if( center_z < 0 || delta_z == 0 ) return 0;
    float alpha = activeLayer.density() * activeLayer.photo( photon_energy );
    float beta =  alpha * center_z;
    if( beta > EXP_FLOAT_TEST ) return 0;
    if( delta_z < 0 ) delta_z = -delta_z;
    float fraction = delta_z * alpha * exp( - beta );
    return fraction;
};

float XrayDetector::ce_minimum( const float channel_energy ) const {
    //  Calculate the minimum photon energy that can yield a Compton shift of the input energy
    //  This is the minimum energy that can produce Compton escape intensity at the input energy
    //  It occurs for 180 degree scattering angle (cos(theta)=-1)
    if( channel_energy <= 0 ) return 0;
    float min_energy = 0.5 * ( channel_energy + sqrt( channel_energy * ( channel_energy + 2 * ME ) ) );
    return min_energy;
};

float XrayDetector::ce_cos_angle( const float photon_energy, const float channel_energy ) const {
    //  Compton scattering angle that produces a Compton electron of channel_energy for an
    //      incident photon of at photon_energy
    if( photon_energy <= channel_energy ) return 0;
    float ec = photon_energy - channel_energy;
    float cos_theta = 1 - ( channel_energy * ME ) / ( ec * photon_energy );
    if( cos_theta < -1 || cos_theta > 1 ) cos_theta = 0;
    return cos_theta;
};

float XrayDetector::ce_fraction( const float photon_energy, const float channel_energy ) const {
    //  Calculates the fraction of incident photons that pass through the detector and undergo
    //      Compton scatter while passing through, such that the Compton electron deposits energy
    //      channel_energy into the detector active volume
    //  Such an event appears at channel_energy in the spectrum (well below the incident photon energy)
    float cos_th = ce_cos_angle( photon_energy, channel_energy );
    if( cos_th == 0 ) return 0;
    //  Calculate some things for either case below
    float th = acos( cos_th );  //  Scatter angle for this incident energy and Compton electron energy
    float ec = ScatterXsectTable::eCompton(  photon_energy, th );   //  Compton-shifted photon energy
    float xsect = activeLayer.incoherent( photon_energy, th );   //  Compton integrated cross-section at incident energy
    float mu_inc = activeLayer.photo( photon_energy );   //  Absorption cross-section at incident energy
    float mu_c = activeLayer.photo( ec );   //  Absorption cross-section at Compton-shifted energy
    float mu_star = mu_inc - mu_c / cos_th;  //  Note that if cos_th is < 0 then this is always > 0
    float prefactor = xsect / mu_star;  //  Density cancels out
    float thickness = activeLayer.mass_thickness(); //  Use mass thickness to go with cross section in cm2/gm
    float exp_term = 0;
    if( cos_th > 0 ) {
        //  Forward scatter, photon emerges from the back of the detector active volume
        float term1_exponent = thickness * mu_c / cos_th;
        if( term1_exponent > EXP_FLOAT_TEST ) return 0;
        //  Use this form for better numeric performance when mu_star < 0
        exp_term = exp( - term1_exponent ) * ( 1 - exp( - thickness * mu_star ) );
    } else {
        //  Back scatter, photon emerges from the front of the detector active volume
        exp_term = 1 - exp( - thickness * mu_star );
    }
    return prefactor * exp_term;
};

//  Electron-escape shelf calculations added Dec. 7, 2020
//      These are also from the Scholze and Procop paper listed above

void XrayDetector::initialize_shelf() {
    //  Sets up calculations of detector shelf from photoelectron and Auger electron escape,
    //      from photons absorbed in active volume and in front contact
    //  One-time computations that don't depend on photon energy
    shelf_constants.clear();
    int ie;
    //  Shelf from photons absorbed in the active volume
    vector <Element> detAE = activeLayer.element_list();
    //  loop over elements in detector active area
    for( ie=0; ie<detAE.size(); ie++ ) {
        //	Create list of x-ray absorption edges for this element
		vector <EdgeIndex> edgeIndexList;
		XrayEdge::numberOfEdges ( edgeIndexList, detAE[ie], 1e6 );  //  Make voltage high to get all edges,= (any not excited will be skipped)
		int edgeIndex;
		for ( edgeIndex=0; edgeIndex<edgeIndexList.size(); edgeIndex++ ) {
			XrayEdge thisEdge ( detAE[ie], edgeIndexList[edgeIndex] );
			if( thisEdge.level() != K ) continue;   //  Disable everything but K edges
            //  Add the photoelectric contributions to the shelf
            ShelfConstants new_photo;
            new_photo.type = PHOTO_ACTIVE_VOLUME;
            new_photo.element = detAE[ie];
            new_photo.binding_energy = thisEdge.energy();
            new_photo.energy = thisEdge.energy();
            //  Compute some things that do not depend on photon energy and store in prefactor
            float beta = 1;
            if( thisEdge.angularMomentum() == s ) beta = 2;
            float rk = thisEdge.jump();
            float abs_ratio = ( rk - 1 ) / rk;
            new_photo.prefactor = abs_ratio * activeLayer.density() * ( 1 - ( beta / 8 ) ) / 4;
            shelf_constants.push_back( new_photo );
            //  Include only principal Augeer lines for selected light elements, and only for K shell
            if( thisEdge.level() != K || detAE[ie].Z() > number_energies_Auger_KLL-1 ) continue;
            //  Add the Auger electron contributions to the shelf
            ShelfConstants new_Auger;
            new_Auger.type = AUGER_ACTIVE_VOLUME;
            new_Auger.element = detAE[ie];
            new_Auger.binding_energy = thisEdge.energy();
            new_Auger.energy = energies_Auger_KLL[ detAE[ie].Z() ];
            //  Compute some things that do not depend on photon energy and store in prefactor
            new_Auger.prefactor = abs_ratio * activeLayer.density() * ( 1 - thisEdge.yield() ) / 4;
            shelf_constants.push_back( new_Auger );
		}
    }
    //  Shelf from photons absorbed in the front contact
    detAE = frontContact.element_list();
    //  loop over elements in detector active area
    for( ie=0; ie<detAE.size(); ie++ ) {
        //	Create list of x-ray absorption edges for this element
		vector <EdgeIndex> edgeIndexList;
		XrayEdge::numberOfEdges ( edgeIndexList, detAE[ie], 1e6 );  //  Make voltage high to get all edges,= (any not excited will be skipped)
		int edgeIndex;
		for ( edgeIndex=0; edgeIndex<edgeIndexList.size(); edgeIndex++ ) {
			XrayEdge thisEdge ( detAE[ie], edgeIndexList[edgeIndex] );
            //  Add the photoelectric contributions to the shelf
            ShelfConstants new_photo;
            new_photo.type = PHOTO_FRONT_CONTACT;
            new_photo.element = detAE[ie];
            new_photo.binding_energy = thisEdge.energy();
            new_photo.energy = thisEdge.energy();
            //  Compute some things that do not depend on photon energy and store in prefactor and terms
            float beta = 1;
            if( thisEdge.angularMomentum() == s ) beta = 2;
            float rk = thisEdge.jump();
            float abs_ratio = ( rk - 1 ) / rk;
            float big_D = frontContact.thickness();
            new_photo.prefactor = abs_ratio * frontContact.density() / 4;
            new_photo.term1 = 2 * big_D;
            new_photo.term2 = - ( 1 + (beta/4) ) * big_D*big_D;
            new_photo.term3 = ( beta / 8 ) * pow(big_D,4);
            shelf_constants.push_back( new_photo );
            //  Include only principal Auger lines for selected light elements, and only for K shell
            if( thisEdge.level() != K || detAE[ie].Z() > number_energies_Auger_KLL-1 ) continue;
            //  Add the Auger electron contributions to the shelf
            ShelfConstants new_Auger;
            new_Auger.type = AUGER_FRONT_CONTACT;
            new_Auger.element = detAE[ie];
            new_Auger.binding_energy = thisEdge.energy();
            new_Auger.energy = energies_Auger_KLL[ detAE[ie].Z() ];
            //  Compute some things that do not depend on photon energy and store in prefactor and terms
            new_Auger.prefactor = abs_ratio * frontContact.density() * ( 1 - thisEdge.yield() ) / 4;
            new_Auger.term1 = 2 * big_D;
            new_Auger.term2 = - big_D*big_D;
            shelf_constants.push_back( new_Auger );
//            cout.setf( ios::scientific, ios::floatfield );
//            cout.precision(2);
//            cout << "Shelf const  pre " << new_Auger.prefactor << "  ratio " << abs_ratio << "   den " << frontContact.density() << "   1-w " << ( 1 - thisEdge.yield() ) << "   D (nm) " << 1e7*big_D << endl;
//            cout.setf( ios::fixed, ios::floatfield );
//            cout.precision(4);
		};
    }
    return;
};

int XrayDetector::electron_shelf( const float photon_energy, std::vector <ShelfStruct> &shelf_contributions ) const {
    //  Calculates detector shelf from photoelectron and Auger electron escape,
    //      from photons absorbed in active volume and in front contact
    shelf_contributions.clear();
    unsigned int i_sh;
    for( i_sh=0; i_sh<shelf_constants.size(); i_sh++ ) {
        //  Ignore if not enough energy to excite this mechanism
        if( photon_energy < shelf_constants[i_sh].binding_energy ) continue;
        //  Finish the calculation for this contribution and add it to the output argument
        if( shelf_constants[i_sh].type == PHOTO_ACTIVE_VOLUME ) {
            float electron_energy = photon_energy - shelf_constants[i_sh].energy;
            float big_R = electron_range( electron_energy, activeLayer.density() );
            float abs_active = activeLayer.photo( shelf_constants[i_sh].element, photon_energy );
            ShelfStruct new_shelf;
            new_shelf.type = shelf_constants[i_sh].type;
            new_shelf.energy_start = photon_energy - electron_energy;
            new_shelf.energy_end = photon_energy;
            new_shelf.probability = shelf_constants[i_sh].prefactor * abs_active * big_R;
            shelf_contributions.push_back( new_shelf );
        } else if( shelf_constants[i_sh].type == AUGER_ACTIVE_VOLUME ) {
            float electron_energy = shelf_constants[i_sh].energy;
            float big_R = electron_range( electron_energy, activeLayer.density() );
            float abs_active = activeLayer.photo( shelf_constants[i_sh].element, photon_energy );
            ShelfStruct new_shelf;
            new_shelf.type = shelf_constants[i_sh].type;
            new_shelf.energy_start = photon_energy - electron_energy;
            new_shelf.energy_end = photon_energy;
            new_shelf.probability = shelf_constants[i_sh].prefactor * abs_active * big_R;
            shelf_contributions.push_back( new_shelf );
        } else if( shelf_constants[i_sh].type == PHOTO_FRONT_CONTACT ) {
            float big_D = frontContact.thickness();
            float electron_energy = photon_energy - shelf_constants[i_sh].energy;
            float big_R = electron_range( electron_energy, frontContact.density() );
            float abs_contact = frontContact.photo( shelf_constants[i_sh].element, photon_energy );
            ShelfStruct new_shelf;
            new_shelf.type = shelf_constants[i_sh].type;
            new_shelf.energy_start = 0;
            new_shelf.energy_end = electron_energy;
            float new_prefactor = shelf_constants[i_sh].prefactor * abs_contact;
            if( big_R < big_D ) {
                new_shelf.probability = shelf_constants[i_sh].prefactor * abs_contact * big_R;
            } else {
                float t2 = shelf_constants[i_sh].term2 / big_R;
                float t3 = shelf_constants[i_sh].term3 / big_R*big_R*big_R;
                new_shelf.probability = new_prefactor * ( shelf_constants[i_sh].term1 + t2 + t3 );
            }
            shelf_contributions.push_back( new_shelf );
        } else if( shelf_constants[i_sh].type == AUGER_FRONT_CONTACT ) {
            float big_D = frontContact.thickness();
            float electron_energy = shelf_constants[i_sh].energy;
            float big_R = electron_range( electron_energy, frontContact.density() );
            float abs_contact = frontContact.photo( shelf_constants[i_sh].element, photon_energy );
            ShelfStruct new_shelf;
            new_shelf.type = shelf_constants[i_sh].type;
            new_shelf.energy_start = 0;
            new_shelf.energy_end = electron_energy;
            if( big_R < big_D ) {
                new_shelf.probability = shelf_constants[i_sh].prefactor * abs_contact * big_R;
            } else {
                float new_prefactor = shelf_constants[i_sh].prefactor * abs_contact;
                float t2 = shelf_constants[i_sh].term2 / big_R;
                new_shelf.probability = new_prefactor * ( shelf_constants[i_sh].term1 + t2 );
            }
            shelf_contributions.push_back( new_shelf );
//            cout.precision(0);
//            cout << "AUGER_FRONT_CONTACT   ee " << electron_energy << "   st " << new_shelf.energy_start << "   end " << new_shelf.energy_end;
//            cout.setf( ios::scientific, ios::floatfield );
//            cout.precision(2);
//            cout << "  pre " << shelf_constants[i_sh].prefactor << "  R (nm) " << 1e7*big_R << "   abs " << abs_contact << "   p " << new_shelf.probability << endl;
//            cout.setf( ios::fixed, ios::floatfield );
//            cout.precision(4);
        }
    }
    return shelf_contributions.size();
};

float XrayDetector::electron_range( const float electron_energy, const float density ) const {
    //  Computes the range of an electron using the high-energy expression of Fittings
    //  Original formula gives range in nanometers with density in gm/cm3 and energy in keV
    return 90.0f * pow( density, -0.8 ) * pow( electron_energy/1000, 1.7 )   * 1e-7;    //  Convert to cm
};



string XrayDetector::toString() const
{
    std::ostringstream os;
    os << "XrayDetector:" << endl;

    os << "  detType: ";
    switch(detType)
    {
    case NO_DETECTOR:
        os << "NO_DETECTOR";
        break;
    case SI_PIN:
        os << "SI_PIN";
        break;
    case SI_SDD:
        os << "SI_SDD";
        break;
    case CD_TE:
        os << "CD_TE";
        break;
    case HP_GE:
        os << "HP_GE";
        break;
    case BAD_DETECTOR:
        os << "BAD_DETECTOR";
        break;
    }
    os << endl;

    os << "  electronic_noise: " << electronic_noise << endl;
    os << "  window: " << window.toString() << endl;
    os << "  deadLayer: " << deadLayer.toString() << endl;
    os << "  activeLayer: " << activeLayer.toString() << endl;


    os << "  fano_factor: " << fano_factor << endl;
    os << "  pair_energy: " << pair_energy << endl;
    os << "  res_fwhm_energy: " << res_fwhm_energy << endl;
    os << "  tail_a: " << tail_a << endl;
//    os << "  tail_C0: " << tail_C0 << endl;
//    os << "  tail_z0: " << tail_z0 << endl;

    return os.str();
}
