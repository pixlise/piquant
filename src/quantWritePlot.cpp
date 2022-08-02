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

//#include <iostream>
//#include <iomanip>
#include "quantWritePlot.h"
#include "XRFconstants.h"
#include "XRFcontrols.h"

//  Written July 25, 2018
//      Consolidate code from main program to generate plot file
//  Modified Sep. 19, 2018
//      Make plot with all-zero measured spectrum work OK
//  Modified May 14, 2019
//`     Write version and file name for calculated spectra (configuration file name put into spectrum for this purpose)
//  Modified Feb. 26, 2021  Add DETECTOR_SHELF and DETECTOR_CE spectrum component types

using namespace std;

int quantWritePlot( const XraySpectrum &singleSpectrum, const std::string plotPathName,
                const PIQUANT_SUBCOMMAND cmd, const int detector_select,
                const std::vector <XraySpectrum> &spectrum_vec, const std::string version ) {

    //      Write a CSV file with plotting information

        ofstream plotFile ( plotPathName.c_str() );
        if( !plotFile ) {
            cout << "Error opening plot file " << plotPathName << endl;
            return -1;
        }
        if( cmd == PRIMARY ) plotFile << "Calculated Primary Spectrum";
        else if( cmd == CALCULATE ) plotFile << "Calculated Full Spectrum";
		plotFile << "   PIQUANT " << version << "  " << singleSpectrum.file_name() << endl;
		bool channels = false;
        if( singleSpectrum.calibration().good() ) {
            channels = false;
            plotFile << "Energy (keV)";
        } else {
            channels = true;
            plotFile << "Channel";
		}
        //  See if there is a measured spectrum to include in plot
		bool measured_present = true;
        int is;
        float sum = 0;
        for( is=0; is<singleSpectrum.numberOfChannels(); is++ )
                sum += singleSpectrum.meas()[is];
        if( sum <= 0 ) measured_present = false;
		if( measured_present && cmd != BULK_SUM_MAX ) plotFile << ", meas";
		if( measured_present && cmd == BULK_SUM_MAX ) plotFile << ", sum";
		bool calc_present = singleSpectrum.calc().size() >= singleSpectrum.numberOfChannels();
		//  Make plots for spectrum of all zeros possible
		if( cmd == PLOT && !calc_present && sum == 0 && singleSpectrum.meas().size() == singleSpectrum.numberOfChannels() ) measured_present = true;
		bool bkg_present = singleSpectrum.bkg().size() >= singleSpectrum.numberOfChannels();
		//  Find number of background components (will need to plot them if more than one)

		bool max_value_present = singleSpectrum.max_value().size() >= singleSpectrum.numberOfChannels();
        if( calc_present ) plotFile << ", calc";
        if( max_value_present ) plotFile << ", max_value";
		if( bkg_present ) plotFile << ", bkg";
		if( measured_present ) plotFile << ", sigma";
		bool residual_present = measured_present && calc_present && singleSpectrum.residual().size() >= singleSpectrum.numberOfChannels();
		if( residual_present ) plotFile << ", residual";
		if( cmd == PLOT ) {
            if( spectrum_vec.size() > 1 && detector_select < 0 ) {
                int isv;
                for( isv=0; isv<spectrum_vec.size(); isv++ ) {
                    if( spectrum_vec[isv].numberOfChannels() < singleSpectrum.numberOfChannels() ) continue;
                    plotFile << ", " << "Det_" << isv;
                }
            }
		} else {
            int ic;
            for( ic=0; ic<singleSpectrum.numberOfComponents(); ic++ ) {
                if( ! singleSpectrum.component( ic ).plot ) continue;
                if( singleSpectrum.component( ic ).spectrum.size() < singleSpectrum.numberOfChannels() ) continue;
                plotFile << ", " << componentDescription( singleSpectrum.component( ic ) );
            }
		}
		plotFile << endl;
        for ( is=0; is<singleSpectrum.numberOfChannels(); is++ ) {
            float e = singleSpectrum.energy( is );
            if( ! channels ) e /= 1000; //  convert to keV if not channels
            plotFile << e;
            if( measured_present ) plotFile << ", " << singleSpectrum.meas()[is];
            //  Plot calculation + background for better comparison to measured spectrum
            //  Now background is fit so already included in calculation
            if( calc_present ) plotFile << ", " << singleSpectrum.calc()[is];
            if( max_value_present ) plotFile << ", " << singleSpectrum.max_value()[is];
            if( bkg_present ) plotFile << ", " << singleSpectrum.bkg()[is];
            if( measured_present ) plotFile << ", " << singleSpectrum.sigma()[is];
            if( residual_present ) plotFile << ", " << singleSpectrum.residual()[is];
            if( cmd == PLOT ) {
                if( spectrum_vec.size() > 1 && detector_select < 0 ) {
                    int isv;
                    for( isv=0; isv<spectrum_vec.size(); isv++ ) {
                        if( spectrum_vec[isv].numberOfChannels() < singleSpectrum.numberOfChannels() ) continue;
                        plotFile << ", "  << spectrum_vec[isv].meas()[is];
                    }
                }
            } else {
                int ic;
                for( ic=0; ic<singleSpectrum.numberOfComponents(); ic++ ) {
                    if( singleSpectrum.component( ic ).spectrum.size() < singleSpectrum.numberOfChannels() ) continue;
                    if( ! singleSpectrum.component( ic ).plot ) continue;
                    if( bkg_present && !singleSpectrum.component( ic ).bkg ) {
                        //      Plot with background added so it sits better on spectrum visually in plot
                        plotFile << ", " << singleSpectrum.component( ic ).coefficient * singleSpectrum.component( ic ).spectrum[is] + singleSpectrum.bkg()[is];
                    } else {
                        //  Plot without background since it is part of the background (or bkg is missing)
                        plotFile << ", " << singleSpectrum.component( ic ).coefficient * singleSpectrum.component( ic ).spectrum[is];
                    }
                }
            }
            plotFile << endl;
        };
        plotFile.close();

	return 0;
};
