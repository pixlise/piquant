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

#ifndef XrayDetector_h
#define XrayDetector_h

#include <vector>
#include "XRFconstants.h"
#include "XrayMaterial.h"
#include "XrayEdge.h"

//  Default values, actual values controlled via -T option (or maybe in configuration file in the future)
//  Changed when bug fixed, unitialized variable in electron escape shelf from active volume for Auger electrons
//  Fixing bug changed shelf and allowed factor to be unity and slope disabled, had to recalibrate with these values
//  Front contact sheld enabled with default thickness of 150 microns, that value also included here
#define SHELF_FACTOR 1         //  Multiplicative factor to increase (or decrease) overall shelf size
#define SHELF_SLOPE 0           //  Slope proportional to loss energy (loss energy is negative)
#define SHELF_SLOPE_START 0     //   Fraction of electron energy
//  Value is 150 nanometers (expressed in cm)  from Scholze and Procop Table 1
#define CONTACT_SDD_UM 150e-7f  //  May need better values for individual detectors, put in config someday


enum DetectorType { NO_DETECTOR = 0, SI_PIN, SI_SDD, CD_TE, HP_GE, BAD_DETECTOR };
//  Additions here must also be made in the functions that read the configuration files

struct EscapeLines {
    float energy;
    float fraction;
};

enum Shelf_Type { PHOTO_ACTIVE_VOLUME = 1, AUGER_ACTIVE_VOLUME, PHOTO_FRONT_CONTACT, AUGER_FRONT_CONTACT };

struct ShelfConstants {
    Shelf_Type type;
    Element element;
    float binding_energy;
    float energy;
    float prefactor;
    float term1;
    float term2;
    float term3;
};

struct ShelfStruct {
    Shelf_Type type;
    float energy_start;
    float energy_end;
    float probability;
};

class XrayDetector {

//	resolution and response calculations for an energy-dispersive detector
//		Modified to calculate resolution from Fano factor, separately for Si, Ge, & CZT     April 12, 2011
//          Also remove Using namespace std from header file
//  Changed class name from AmpTekDet to XrayDetector   May 11, 2016
//      Removed resolution function to set resolution to avoid name conflict with resolution function vs energy
//  Modified Aug. 10, 2016
//      Prevent resolution from becoming not-a-number at zero energy for bad values of reference resolution or fano factor
//  Modified Dec. 6, 2016
//      Change enum DetectorType to match with interpretation of EMSA format files (kept compatibility with xsp files)
//  Modified Feb. 1, 2017
//      Use XrayMaterial class for window, dead layer, and active layer
//  Modified Feb. 20, 2017
//      Set al defaults here to zero so defaults in constructor will be used
//  Modified Apr. 11, 2018
//      Calculate resolution as electronic noise plus Fano instead of Mn Ka resolution adjusted using Fano
//          This makes it easier and more reasonable to avoid square root of a negative number or zero resolution
//      Moved setResolution out of this file as part of this change
//  Modified July 10, 2021  Add pulse resolving time for simple pulse pileup calculation

public:
//		must have default constructor to declare arrays
//		Main constructor has all arguments defaulted, so it serves as the default constructor

//		resolution in eV at Mn Ka (5984 eV)
//			or another energy as set using fwhn_energy member function
//		Be window thickness in microns
//		dead layer thickness in microns
//		active layer thickness in mm
	XrayDetector( float detR_in = 0, float detW_in = 0, float detD_in = 0, float detA_in = 0, DetectorType detType_in = SI_PIN );
//		Member functions
	DetectorType type () const { return detType; };
	float resolution ( const float energy = RESOLUTION_REFERENCE_ENERGY ) const;
	float response ( const float energy ) const;
	float escape ( const float energy, std::vector<EscapeLines> &escape_line_vector ) const;
//		set and retrieve private data
	void setResolution ( const float res_in, const float ref_energy = -1 );
	float fano ( ) const { return fano_factor; };
	float energy_per_pair ( ) const { return pair_energy; };
	float fwhm_energy ( ) const { return res_fwhm_energy; };
	void fano ( const float fano_in ) {
        if ( fano_in > 0 && fano_in < 1.0f ) fano_factor = fano_in;
        return; };
	void energy_per_pair ( const float energy_per_pair_in ) {
        if ( energy_per_pair_in > 0 ) pair_energy = energy_per_pair_in;
        return; };
	void fwhm_energy ( const float res_fwhm_energy_in )
		{ if ( res_fwhm_energy_in >= 0 ) res_fwhm_energy = res_fwhm_energy_in; return; };
    const XrayMaterial &window_material() const { return window; };
    const XrayMaterial &dead_layer() const { return deadLayer; };
    const XrayMaterial &active_layer() const { return activeLayer; };
	float z_of_C( const float photon_energy, const float channel_energy ) const;
//	float z_of_C( const float charge_ratio ) const;
	float energy_for_C0( const float photon_energy ) const { return photon_energy * tail_C0( photon_energy ); };
	float exp_term_of_C( const float photon_energy, const float channel_energy ) const;
	float tail_fraction( const float photon_energy, const float channel_E1, const float channel_E2 ) const;
	float ce_minimum( const float channel_energy ) const;
	float ce_cos_angle( const float photon_energy, const float channel_energy ) const;
	float ce_fraction( const float photon_energy, const float channel_energy ) const;
	int electron_shelf( const float photon_energy, std::vector <ShelfStruct> &shelf_contributions ) const;
	float get_shelf_factor() const { return shelf_factor; };
	float get_shelf_slope() const { return shelf_slope; }
	float get_shelf_slope_start() const { return shelf_slope_start; }
	void set_shelf_factor( const float shelf_factor_in ) { shelf_factor = shelf_factor_in; };
	void set_shelf_slope( const float shelf_slope_in ) { shelf_slope = shelf_slope_in; };
	void set_shelf_slope_start( const float shelf_slope_start_in ) { shelf_slope_start = shelf_slope_start_in; };
	float pileup_time() const { return pulse_resolving_time; };
	void pileup_time( const float time_in ) { pulse_resolving_time = time_in; };

    std::string toString() const;

private:

	DetectorType detType;
//		resolution at Mn Ka line (5984 eV) in eV
//      changed April 11, 2018 to be resolution extrapolated to zero energy
//          which is electronic noise contribution to resolution
	float electronic_noise;
//		detector window is beryllium
	XrayMaterial window;
//		detector dead layer
    XrayMaterial deadLayer;
//		detector active layer
    XrayMaterial activeLayer;
    XrayMaterial frontContact;  //  Added Dec. 7, 2020 for electron loss shelf calculations
	float fano_factor;
	float pair_energy;
	float res_fwhm_energy;
	//  Tail calculations added August 26, 2020       [Scholze and Procop values for SDD, Table 1:  a=0.5, C0=0.9, z0=50 nm]
	float tail_a = 0.4f;        //  From ATLO Fe55 and ElCal SiO2 (Fe55 may be different peaking time?)
//	float tail_C0 = 0.8f;       //  Changed in tail code, different values below and above Si K edge
//	float tail_z0 = 120e-7f;    //  from ElCal BHVO, COQ-1, and GYP-B   120 nanometers (expressed in cm) [was 220 from Fe55 & SiO2]
//	float tail_z0 = 240e-7f;    //  from ElCal SiO2
//	float tail_z0 = 35e-7f;    //  from ElCal ZnS (S peak)
//	float tail_z0 = 300e-7f;    //  from ElCal
	vector <ShelfConstants> shelf_constants;
	float shelf_factor = SHELF_FACTOR;
	float shelf_slope = SHELF_SLOPE;
	float shelf_slope_start = SHELF_SLOPE_START;
	//  Pulse resolving time for simple pulse pileup calculation added July 10, 2021
	float pulse_resolving_time = 0.1e-6;    //  0.1 microsecond (integration time for fast channel used for pileup rejection)

	void initialize_shelf();
	float electron_range( const float electron_energy, const float density ) const;
	const float tail_C0( const float energy ) const;

};

#endif
